use std::env;
use std::io;
use std::collections::HashSet;
use sha2::Sha256;
use bio::io::fastq;
use bio::io::fastq::FastqRead;
use anyhow::Result;
use sha2::Digest;

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} query.fastq ... < to_index.fastq > missing_reads.fastq", args[0]);
        std::process::exit(1);
    }
    let mut hashset = HashSet::new();
    let mut reader = fastq::Reader::new(io::stdin());
    let mut record = fastq::Record::new();
    reader.read(&mut record)?;
    while !record.is_empty() {
        let hash = Sha256::new().chain_update(record.seq()).finalize();
        hashset.insert(hash);
        reader.read(&mut record)?;
    }

    let mut total = 0u64;
    let mut notfound = 0u64;
    let mut writer = fastq::Writer::new(io::stdout());
    for arg in &args[1..args.len()] {
        let mut record = fastq::Record::new();
        let mut reader2 = fastq::Reader::from_file(&arg)?;
        reader2.read(&mut record)?;
        while !record.is_empty() {
            total += 1;
            let hash = Sha256::new().chain_update(record.seq()).finalize();
            if !hashset.contains(&hash) {
                writer.write_record(&record)?;
                notfound += 1;
            }
            reader2.read(&mut record)?;
        }
    }
    eprintln!("Missing {} of {} reads", notfound, total);
    Ok(())
}
