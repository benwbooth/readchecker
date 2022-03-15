use std::env;
use sorted_vec::SortedSet;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufReader;
use std::io::Read;
use sha1::Sha1;
use sha1::digest::Digest;
use bio::io::fastq;
use bio::io::fastq::FastqRead;
use anyhow::Result;
use anyhow::bail;
use flate2::bufread::MultiGzDecoder;
use byte_unit::Byte;
use rayon::prelude::*;
use dashmap::DashMap;
use atomic_counter::AtomicCounter;
use atomic_counter::RelaxedCounter;
use generic_array::GenericArray;
use generic_array::typenum::U20;
use std::sync::Arc;
use std::sync::Mutex;
use rand::seq::SliceRandom;

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 6 {
        eprintln!("Usage: {} to_index_fastq_files.tsv query_fastq_files.tsv found_reads_out.tsv uncategorized_reads_out.tsv missing_reads_out.fastq", args[0]);
        std::process::exit(1);
    }
    let found_reads = File::create(&args[3])?;
    let fr = Arc::new(Mutex::new(csv::WriterBuilder::new().
        delimiter(b'\t').
        has_headers(false).
        from_writer(found_reads)));

    let uncategorized_reads = File::create(&args[4])?;
    let ur = Arc::new(Mutex::new(csv::WriterBuilder::new().
        delimiter(b'\t').
        has_headers(false).
        from_writer(uncategorized_reads)));

    let missing_reads = Arc::new(Mutex::new(fastq::Writer::to_file(&args[5])?));

    let threads = match std::env::var("THREADS") {
        Ok(t) => t.parse::<usize>()?,
        _ => num_cpus::get(),
    };
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global()?;

    let bytes_per_page = procfs::page_size()?;

    let mut to_index_fastq_reader = csv::ReaderBuilder::new().
        delimiter(b'\t').
        has_headers(false).
        from_path(&args[1])?;
    let mut to_index_fastq = Vec::<String>::new();
    let file2category = DashMap::<usize,String>::new();
    for (i, record) in to_index_fastq_reader.records().enumerate() {
        let record = record?;
        let file = &record[0];
        let category= if record.len() > 1 { &record[1] } else { "" };
        to_index_fastq.push(file.to_string());
        file2category.insert(i, category.to_string());
    }

    let seqsha2fileidxset = DashMap::<GenericArray<u8, U20>, SortedSet<u32>>::new();
    let counter = RelaxedCounter::new(0usize);
    let result: Result<(), anyhow::Error> = to_index_fastq.par_iter().with_max_len(1).enumerate().try_for_each(|(i, to_index)| {

        let br: Box<dyn Read> = if to_index.ends_with(".gz") {
            Box::new(MultiGzDecoder::new(BufReader::new(File::open(to_index)?)))
        } else {
            Box::new(File::open(to_index)?)
        };
        let mut reader = fastq::Reader::new(br);

        let mut record = fastq::Record::new();
        let mut record_num = 0usize;
        record_num += 1;
        while let Err(e) = reader.read(&mut record) { 
            eprintln!("Could not read file {to_index} at record {record_num}, line {record_line}: {e:?}", record_line=record_num*4);
            record_num += 1;
        }
        while !record.is_empty() {
            let mut seed = Sha1::new();
            seed.update(record.seq());
            let bytes  = seed.finalize();
            match seqsha2fileidxset.get_mut(&bytes) {
                Some(mut fileidxset) => {
                    (*fileidxset).insert(i.try_into()?);
                    ()
                },
                _ => {
                    let mut fileidxset = SortedSet::new();
                    fileidxset.insert(i.try_into()?);
                    seqsha2fileidxset.insert(bytes, fileidxset);
                }
            }
            record_num += 1;
            while let Err(e) = reader.read(&mut record) { 
                eprintln!("Could not read file {to_index} at record {record_num}, line {record_line}: {e:?}", record_line=record_num*4);
                record_num += 1;
            }
        }
        let mem_usage_bytes = procfs::process::Process::myself()?.stat()?.rss * bytes_per_page;
        let c = counter.inc();
        eprintln!("Indexing file {to_index} ({c}/{to_index_length}), total reads={reads}, memory usage={memusage}", 
            c=c+1, 
            to_index_length=to_index_fastq.len(), 
            reads=seqsha2fileidxset.len(),
            memusage=Byte::from_bytes(mem_usage_bytes.try_into()?).get_appropriate_unit(false).to_string());
        Ok(())
    });
    result?;

    let mut query_fastq = Vec::<String>::new();
    let mut query_fastq_reader = csv::ReaderBuilder::new().
        delimiter(b'\t').
        has_headers(false).
        from_path(&args[2])?;
    let expected_category = DashMap::<usize,SortedSet<String>>::new();
    for (i, record) in query_fastq_reader.records().enumerate() {
        let record = record?;
        let file = &record[0];
        let categories = &record[1];

        let mut rdr = csv::ReaderBuilder::new().from_reader(categories.as_bytes());
        let cat_records = rdr.records().collect::<std::result::Result<Vec<csv::StringRecord>,csv::Error>>()?;
        for cat_record in cat_records {
            for cat in cat_record.iter() {
                expected_category.entry(i).or_insert_with(|| SortedSet::new()).insert(cat.to_string());
            }
        }
        query_fastq.push(file.to_string());
    }
    let counter = RelaxedCounter::new(0usize);
    let result: Result<(), anyhow::Error> = query_fastq.par_iter().enumerate().with_max_len(1).
        try_for_each_with((fr, ur, missing_reads), 
        |(fr, ur, missing_reads), (i, file)| 
    {
        let c = counter.inc();
        let mem_usage_bytes = procfs::process::Process::myself()?.stat()?.rss * bytes_per_page;
        eprintln!("Processing file {file} ({c}/{query_fastq_length}), memory usage={memusage}", 
            c=c+1,
            query_fastq_length=query_fastq.len(), 
            memusage=Byte::from_bytes(mem_usage_bytes.try_into()?).get_appropriate_unit(false).to_string());

        let br: Box<dyn Read> = if file.ends_with(".gz") {
            Box::new(MultiGzDecoder::new(BufReader::new(File::open(file)?)))
        } else {
            Box::new(File::open(file)?)
        };
        let mut reader = fastq::Reader::new(br);

        let mut category2count = BTreeMap::<String,usize>::new();
        let mut record = fastq::Record::new();
        let mut record_num = 0usize;
        record_num += 1;
        let mut read_count = 0usize;
        while let Err(e) = reader.read(&mut record) { 
            eprintln!("Could not parse fastq file {file} at record {record_num}, line {record_line}: {e:?}", record_line=record_num*4);
            record_num += 1;
        }
        while !record.is_empty() {
            let mut seed = Sha1::new();
            seed.update(record.seq());
            let bytes  = seed.finalize();
            match seqsha2fileidxset.get(&bytes) {
                Some(fileidxset) => {
                    let mut expected = Vec::<String>::new();
                    let mut unexpected = Vec::<String>::new();
                    for j in fileidxset.iter() {
                        if let Some(cat) = file2category.get(&(*j as usize)) {
                            if let Some(exp) = expected_category.get(&i) {
                                if exp.contains(&*cat) {
                                    expected.push(cat.to_string());
                                }
                                else {
                                    unexpected.push(cat.to_string());
                                }
                            }
                            else {
                                unexpected.push(cat.to_string());
                            }
                        }
                    }
                    if !expected.is_empty() {
                        let selected = expected.choose_multiple(&mut rand::thread_rng(), 1).collect::<Vec<&String>>()[0];
                        match category2count.get_mut(selected.as_str()) {
                            Some(c2c) => *c2c += 1,
                            None => { category2count.insert(selected.to_string(), 1); }
                        }
                    }
                    else if !unexpected.is_empty() {
                        let selected = unexpected.choose_multiple(&mut rand::thread_rng(), 1).collect::<Vec<&String>>()[0];
                        match category2count.get_mut(selected.as_str()) {
                            Some(c2c) => *c2c += 1,
                            None => { category2count.insert(selected.to_string(), 1); }
                        }
                    }
                    else {
                        let mut file_wtr = csv::WriterBuilder::new().from_writer(vec![]);
                        file_wtr.write_record(fileidxset.iter().map(|idx| to_index_fastq[*idx as usize].to_string()).collect::<Vec<_>>())?;
                        let files_str= String::from_utf8(file_wtr.into_inner()?)?;
                        match category2count.get_mut("") {
                            Some(c2c) => *c2c += 1,
                            None => { category2count.insert("".to_string(), 1); },
                        }
                        eprintln!("Read {read_id} matches uncategorized file(s): {files_str}", read_id=record.id());

                        let mut ur = ur.lock().or_else(|_| bail!("Could not write to uncategorized_reads file {f}!", f=args[4]))?;
                        (*ur).write_record([file, record.id(), files_str.trim_end_matches(['\r','\n'])])?;
                    }
                },
                None => {
                    let id = record.id();
                    eprintln!("{file}: read {id} not found in index");
                    let mut mr = missing_reads.lock().or_else(|_| bail!("Could not write to missing_reads file {f}!", f=args[4]))?;
                    (*mr).write_record(&record)?;
                }
            }
            record_num += 1;
            read_count += 1;
            while let Err(e) = reader.read(&mut record) { 
                eprintln!("Could not parse fastq file {file} at record {record_num}, line {record_line}: {e:?}", record_line=record_num*4);
                record_num += 1;
            }
        }
        let mut summary = vec![];
        for (category, count) in category2count.iter() {
            let category = if category == "" { "other" } else { category };
            summary.push(format!("{category} ({count}/{read_count}, {pct:0.2}%)", pct=(*count as f64 / read_count as f64)*100.0));
        }
        let mut file_wtr = csv::WriterBuilder::new().from_writer(vec![]);
        file_wtr.write_record(summary)?;
        let files_str= String::from_utf8(file_wtr.into_inner()?)?;
        let mut fr = fr.lock().or_else(|_| bail!("Could not write to found_reads file {f}!", f=args[3]))?;
        (*fr).write_record([file, files_str.trim_end_matches(['\r','\n'])])?;
        Ok(())
    });
    result?;
    Ok(())
}
