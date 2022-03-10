use std::env;
use std::collections::BTreeSet;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;
use std::io::Read;
//use sha2::Sha256;
//use sha2::digest::Digest;
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
use csv::WriterBuilder;
use std::sync::Arc;
use std::sync::Mutex;

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        eprintln!("Usage: {} to_index_fastq_files.list query_fastq_files.list found_reads_out.tsv missing_reads_out.fastq", args[0]);
        std::process::exit(1);
    }
    let found_reads = File::create(&args[3])?;
    let fr = Arc::new(Mutex::new(WriterBuilder::new().delimiter(b'\t').from_writer(found_reads)));

    let missing_reads = Arc::new(Mutex::new(fastq::Writer::to_file(&args[4])?));

    let threads = match std::env::var("THREADS") {
        Ok(t) => t.parse::<usize>()?,
        _ => 15usize,
    };
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global()?;

    let bytes_per_page = procfs::page_size()?;

    let to_index_fastq_files = File::open(&args[1])?;
    let mut to_index_fastq_files_br = BufReader::new(to_index_fastq_files);
    let mut line = String::new();
    let mut to_index_fastq = Vec::<String>::new();
    while to_index_fastq_files_br.read_line(&mut line)? > 0 {
        let file = line.trim_end_matches(['\r','\n']).to_string();
        to_index_fastq.push(file);
        line.clear();
    }

    let seqsha2fileidxset = DashMap::<GenericArray<u8, U20>, BTreeSet<u32>>::new();
    let counter = RelaxedCounter::new(0usize);
    let result: Result<(), anyhow::Error> = to_index_fastq.par_iter().enumerate().try_for_each(|(i, to_index)| {

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
            //eprintln!("seq={seq}, bytes={bytes:?}", seq=std::str::from_utf8(record.seq())?);
            match seqsha2fileidxset.get_mut(&bytes) {
                Some(mut fileidxset) => {
                    (*fileidxset).insert(i as u32);
                    ()
                },
                _ => {
                    let mut fileidxset = BTreeSet::new();
                    fileidxset.insert(i as u32);
                    seqsha2fileidxset.insert(bytes, fileidxset);
                }
            }
            record_num += 1;
            while let Err(e) = reader.read(&mut record) { 
                eprintln!("Could not read file {to_index} at record {record_num}, line {record_line}: {e:?}", record_line=record_num*4);
                record_num += 1;
            }
        }
        let mem_usage_bytes = procfs::process::Process::myself()?.stat()?.rss as u64 * bytes_per_page as u64;
        let c = counter.inc();
        eprintln!("Indexing file {to_index} ({c}/{to_index_length}), total reads={reads}, memory usage={memusage}", 
            c=c+1, 
            to_index_length=to_index_fastq.len(), 
            reads=seqsha2fileidxset.len(),
            memusage=Byte::from_bytes(mem_usage_bytes as u128).get_appropriate_unit(false).to_string());
        Ok(())
    });
    result?;

    let query_fastq_files = File::open(&args[2])?;
    let mut query_fastq_files_br = BufReader::new(query_fastq_files);
    let mut query_fastq = Vec::<String>::new();
    line.clear();
    while query_fastq_files_br.read_line(&mut line)? > 0 {
        let file = line.trim_end_matches(['\r','\n']).to_string();
        query_fastq.push(file);
        line.clear();
    }
    let counter = RelaxedCounter::new(0usize);
    let result: Result<(), anyhow::Error> = query_fastq.par_iter().try_for_each_with((fr, missing_reads), |(fr, missing_reads), file| {
        let c = counter.inc();
        let mem_usage_bytes = procfs::process::Process::myself()?.stat()?.rss as u64 * bytes_per_page as u64;
        eprintln!("Processing file {file} ({c}/{query_fastq_length}), memory usage={memusage}", 
            c=c+1,
            query_fastq_length=query_fastq.len(), 
            memusage=Byte::from_bytes(mem_usage_bytes as u128).get_appropriate_unit(false).to_string());

        let br: Box<dyn Read> = if file.ends_with(".gz") {
            Box::new(MultiGzDecoder::new(BufReader::new(File::open(file)?)))
        } else {
            Box::new(File::open(file)?)
        };
        let mut reader = fastq::Reader::new(br);

        let mut fileidxset2count = BTreeMap::<BTreeSet<u32>,usize>::new();
        let mut record = fastq::Record::new();
        let mut record_num = 0usize;
        record_num += 1;
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
                    match fileidxset2count.get_mut(&*fileidxset) {
                        Some(count) => *count += 1,
                        _ => {
                            fileidxset2count.insert((*fileidxset).clone(), 1usize);
                            ()
                        },
                    };
                },
                None => {
                    let id = record.id();
                    eprintln!("{file}: read {id} not found in index");
                    let mut mr = missing_reads.lock().or_else(|_| bail!("Could not write to missing_reads file {f}!", f=args[4]))?;
                    (*mr).write_record(&record)?;
                }
            }
            record_num += 1;
            while let Err(e) = reader.read(&mut record) { 
                eprintln!("Could not parse fastq file {file} at record {record_num}, line {record_line}: {e:?}", record_line=record_num*4);
                record_num += 1;
            }
        }
        for (fileidxset, count) in fileidxset2count.iter() {
            let mut files = Vec::new();
            for i in fileidxset.iter() {
                files.push(to_index_fastq[*i as usize].clone());
            }
            let mut file_wtr = WriterBuilder::new().from_writer(vec![]);
            file_wtr.write_record(files)?;
            let files_str= String::from_utf8(file_wtr.into_inner()?)?;
            let mut fr = fr.lock().or_else(|_| bail!("Could not write to found_reads file {f}!", f=args[3]))?;
            (*fr).write_record([file, files_str.trim_end_matches(['\r','\n']), &count.to_string()])?
        }
        Ok(())
    });
    result?;
    Ok(())
}
