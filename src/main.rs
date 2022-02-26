use std::env;
use std::collections::HashMap;
use std::collections::HashSet;
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

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} to_index_fastq_files.list query_fastq_files.list > missing_reads.fastq", args[0]);
        std::process::exit(1);
    }

    let to_index_fastq_files = File::open(&args[1])?;
    let mut to_index_fastq_files_br = BufReader::new(to_index_fastq_files);
    let mut line = String::new();
    let mut to_index_fastq = Vec::<String>::new();
    while to_index_fastq_files_br.read_line(&mut line)? > 0 {
        to_index_fastq.push(line.trim_end_matches(['\r','\n']).to_string());
    }

    let mut seed = Sha256::new();
    let mut bytes = [0u8; 32];
    let mut hashmap = HashMap::<[u8; 32], Vec<usize>>::new();
    for (i, to_index) in to_index_fastq.iter().enumerate() {
        eprintln!("Indexing file {to_index}");
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
            match hashmap.get_mut(&bytes) {
                Some(vec) => vec.push(i),
                None => {
                    let mut vec = Vec::new();
                    vec.push(i);
                    hashmap.insert(bytes.clone(), vec);
                }
            }
            while reader.read(&mut record).is_err() { }
        }
    }

    let mut line = String::new();
    let query_fastq_files = File::open(&args[1])?;
    let mut query_fastq_files_br = BufReader::new(query_fastq_files);
    while query_fastq_files_br.read_line(&mut line)? > 0 {
        let file = line.trim_end_matches(['\r','\n']);
        let mut reader: fastq::Reader<BufReader<Box<dyn Read>>> = if file.ends_with(".gz") {
            fastq::Reader::new(Box::new(GzDecoder::new(File::open(file)?)))
        } else {
            fastq::Reader::new(Box::new(File::open(file)?))
        };
        let mut fileset = HashSet::<usize>::new();
        let mut record = fastq::Record::new();
        while reader.read(&mut record).is_err() { }
        while !record.is_empty() {
            seed.reset();
            seed.update(record.seq());
            seed.write(&mut bytes)?;
            match hashmap.get(&bytes) {
                Some(vec) => {
                    fileset.extend(vec);
                },
                None => {
                    let id = record.id();
                    eprintln!("{file}: read {id} not found in index");
                }
            }
            while reader.read(&mut record).is_err() { }
        }
        let mut found_in = Vec::new();
        for i in fileset.iter() {
            found_in.push(to_index_fastq[*i].clone());
        }
        let found_in_str = found_in.join(",");
        println!("{file}\t{found_in_str}");
    }
    Ok(())
}
