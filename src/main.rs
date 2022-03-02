use std::env;
use std::collections::BTreeSet;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;
use std::io::Read;
use std::io::Write;
use sha2::Sha256;
use sha2::digest::Digest;
use bio::io::fastq;
use bio::io::fastq::FastqRead;
use anyhow::Result;
use flate2::read::GzDecoder;
use byte_unit::Byte;
use rayon::prelude::*;
use dashmap::DashMap;

fn main() -> Result<()> {
    let threads = match std::env::var("THREADS") {
        Ok(t) => t.parse::<usize>()?,
        _ => 15usize,
    };
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global()?;

    let bytes_per_page = procfs::page_size()?;

    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} to_index_fastq_files.list query_fastq_files.list >found_reads.list 2> missing_reads.list", args[0]);
        std::process::exit(1);
    }

    let to_index_fastq_files = File::open(&args[1])?;
    let mut to_index_fastq_files_br = BufReader::new(to_index_fastq_files);
    let mut line = String::new();
    let mut to_index_fastq = Vec::<String>::new();
    while to_index_fastq_files_br.read_line(&mut line)? > 0 {
        let file = line.trim_end_matches(['\r','\n']).to_string();
        to_index_fastq.push(file);
        line.clear();
    }

    let seqsha2fileidxset = DashMap::<[u8; 32], BTreeSet<u32>>::new();
    let result: Result<(), anyhow::Error> = to_index_fastq.par_iter().enumerate().try_for_each(|(i, to_index)| {
        let mut seed = Sha256::new();
        let mut bytes = [0u8; 32];
        let mem_usage_bytes = procfs::process::Process::myself()?.stat()?.rss as u64 * bytes_per_page as u64;
        eprintln!("Indexing file {to_index} ({ip1}/{to_index_length}), memory usage={memusage}", 
            ip1=i+1, 
            to_index_length=to_index_fastq.len(), 
            memusage=Byte::from_bytes(mem_usage_bytes as u128).get_appropriate_unit(false).to_string());
        let mut reader: fastq::Reader<BufReader<Box<dyn Read>>> = if to_index.ends_with(".gz") {
            fastq::Reader::new(Box::new(GzDecoder::new(File::open(to_index)?)))
        } else {
            fastq::Reader::new(Box::new(File::open(to_index)?))
        };
        let mut record = fastq::Record::new();
        while reader.read(&mut record).is_err() { }
        while !record.is_empty() {
            seed.reset();
            seed.update(record.seq());
            seed.write(&mut bytes)?;
            match seqsha2fileidxset.get_mut(&bytes) {
                Some(mut fileidxset) => {
                    (*fileidxset).insert(i as u32);
                    ()
                },
                _ => {
                    let mut fileidxset = BTreeSet::new();
                    fileidxset.insert(i as u32);
                    seqsha2fileidxset.insert(bytes.clone(), fileidxset);
                }
            }
            while reader.read(&mut record).is_err() { }
        }
        Ok(())
    });
    result?;

    let query_fastq_files = File::open(&args[1])?;
    let mut query_fastq_files_br = BufReader::new(query_fastq_files);
    let mut query_fastq = Vec::<String>::new();
    line.clear();
    while query_fastq_files_br.read_line(&mut line)? > 0 {
        let file = line.trim_end_matches(['\r','\n']).to_string();
        query_fastq.push(file);
        line.clear();
    }
    let result: Result<(), anyhow::Error> = query_fastq.par_iter().enumerate().try_for_each(|(i, file)| {
        let mem_usage_bytes = procfs::process::Process::myself()?.stat()?.rss as u64 * bytes_per_page as u64;
        eprintln!("Processing file {file} ({ip1}/{query_fastq_length}), memory usage={memusage}", 
            ip1=i+1, 
            query_fastq_length=query_fastq.len(), 
            memusage=Byte::from_bytes(mem_usage_bytes as u128).get_appropriate_unit(false).to_string());
        let mut seed = Sha256::new();
        let mut bytes = [0u8; 32];
        let mut reader: fastq::Reader<BufReader<Box<dyn Read>>> = if file.ends_with(".gz") {
            fastq::Reader::new(Box::new(GzDecoder::new(File::open(file)?)))
        } else {
            fastq::Reader::new(Box::new(File::open(file)?))
        };
        let mut fileidxset2count = BTreeMap::<BTreeSet<u32>,usize>::new();
        let mut record = fastq::Record::new();
        while reader.read(&mut record).is_err() { }
        while !record.is_empty() {
            seed.reset();
            seed.update(record.seq());
            seed.write(&mut bytes)?;
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
                }
            }
            while reader.read(&mut record).is_err() { }
        }
        for (fileidxset, count) in fileidxset2count.iter() {
            let mut files = Vec::new();
            for i in fileidxset.iter() {
                files.push(to_index_fastq[*i as usize].clone());
            }
            println!("{file}\t{files_str}\t{count}", files_str=files.join(","));
        }
        Ok(())
    });
    result?;
    Ok(())
}
