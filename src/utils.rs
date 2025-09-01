use anyhow::Result;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::str;

pub type ChrMap = FxHashMap<String, u8>;

#[derive(Debug, Clone)]
pub struct Pair {
    pub chr1: u8,
    pub pos1: u32,
    pub chr2: u8,
    pub pos2: u32,
}

pub fn create_chr_map(chrom_size_file: Option<&str>) -> ChrMap {
    if let Some(filename) = chrom_size_file {
        create_chr_map_from_file(filename).unwrap_or_else(|_| {
            eprintln!(
                "Warning: Could not load {}, using default human chromosome map",
                filename
            );
            create_default_chr_map()
        })
    } else {
        create_default_chr_map()
    }
}

pub fn create_chr_map_from_file(filename: &str) -> Result<ChrMap> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut map = ChrMap::default();
    let mut chr_index = 1u8;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 {
            let chr_name = parts[0].to_string();
            map.insert(chr_name, chr_index);
            chr_index = chr_index.saturating_add(1);
        }
    }

    Ok(map)
}

pub fn create_default_chr_map() -> ChrMap {
    let mut map = ChrMap::default();

    // Human chromosomes 1-22, X, Y
    for i in 1..=22 {
        map.insert(i.to_string(), i as u8);
    }
    map.insert("X".to_string(), 23);
    map.insert("Y".to_string(), 24);
    map.insert("23".to_string(), 23);
    map.insert("24".to_string(), 24);

    // Add common prefixes
    for i in 1..=22 {
        map.insert(format!("chr{}", i), i as u8);
    }
    map.insert("chrX".to_string(), 23);
    map.insert("chrY".to_string(), 24);
    map.insert("chr23".to_string(), 23);
    map.insert("chr24".to_string(), 24);

    map
}

pub fn get_genome_lengths(chrom_size_file: Option<&str>) -> Vec<u32> {
    if let Some(filename) = chrom_size_file {
        get_genome_lengths_from_file(filename).unwrap_or_else(|_| {
            eprintln!(
                "Warning: Could not load {}, using default hg19 lengths",
                filename
            );
            get_default_genome_lengths()
        })
    } else {
        get_default_genome_lengths()
    }
}

pub fn get_genome_lengths_from_file(filename: &str) -> Result<Vec<u32>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut lengths = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 {
            if let Ok(length) = parts[1].parse::<u32>() {
                lengths.push(length);
            }
        }
    }

    Ok(lengths)
}

pub fn get_default_genome_lengths() -> Vec<u32> {
    // hg19 chromosome lengths (from UCSC)
    vec![
        249250621, // chr1
        242193529, // chr2
        198295559, // chr3
        191154276, // chr4
        180915260, // chr5
        171115067, // chr6
        159138663, // chr7
        146364022, // chr8
        141213431, // chr9
        135534747, // chr10
        135006516, // chr11
        133851895, // chr12
        115169878, // chr13
        107349540, // chr14
        102531392, // chr15
        90354753,  // chr16
        81195210,  // chr17
        78077248,  // chr18
        59128983,  // chr19
        63025520,  // chr20
        48129895,  // chr21
        51304566,  // chr22
        155270560, // chrX
        59373566,  // chrY
    ]
}

#[inline]
pub fn parse_u32_fast(s: &[u8]) -> Option<u32> {
    if s.is_empty() {
        return None;
    }

    let mut result = 0u32;
    for &byte in s {
        if byte >= b'0' && byte <= b'9' {
            result = result * 10 + (byte - b'0') as u32;
        } else {
            return None;
        }
    }
    Some(result)
}

#[inline]
pub fn parse_chr(s: &[u8], chr_map: &ChrMap) -> Option<u8> {
    let s_str = str::from_utf8(s).ok()?;
    chr_map.get(s_str).copied()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn repo_file(name: &str) -> String {
        let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        p.push(name);
        p.to_string_lossy().into_owned()
    }

    #[test]
    fn reads_chrom_size_lengths_and_map() {
        let path = repo_file("chrom.size");

        // Validate lengths parsing
        let lengths = get_genome_lengths_from_file(&path).expect("should read chrom.size lengths");
        assert!(!lengths.is_empty(), "chrom.size produced no lengths");
        // Spot-check first contig length from repo file
        assert_eq!(lengths[0], 138_735_004, "first contig length mismatch");

        // Validate chromosome map parsing
        let map = create_chr_map_from_file(&path).expect("should read chrom.size map");
        assert!(map.contains_key("ptg000001l"), "missing first contig key");
        assert!(map.contains_key("ptg000040l"), "missing expected contig key");
    }
}
