use crate::utils::{Pair, ChrMap, parse_u32_fast, parse_chr};
use memchr::memchr;
use std::io::{BufRead, BufReader};
use std::io::Read;
use flate2::read::MultiGzDecoder;
use anyhow::Result;

pub struct PairIterator<R: BufRead> {
    reader: R,
    chr_map: ChrMap,
    buffer: String,
}

impl<R: BufRead> PairIterator<R> {
    pub fn new(reader: R, chr_map: ChrMap) -> Self {
        Self {
            reader,
            chr_map,
            buffer: String::with_capacity(1024),
        }
    }
}

impl<R: BufRead> Iterator for PairIterator<R> {
    type Item = Result<Pair>;

    fn next(&mut self) -> Option<Self::Item> {
        static mut LINE_COUNT: u64 = 0;
        static mut PARSED_COUNT: u64 = 0;
        static mut DEBUG_SHOWN: bool = false;
        
        loop {
            self.buffer.clear();
            match self.reader.read_line(&mut self.buffer) {
                Ok(0) => {
                    unsafe {
                        if LINE_COUNT > 0 {
                            eprintln!("Debug: Processed {} lines, parsed {} pairs", LINE_COUNT, PARSED_COUNT);
                        }
                    }
                    return None; // EOF
                }
                Ok(_) => {
                    unsafe {
                        LINE_COUNT += 1;
                        if !DEBUG_SHOWN && LINE_COUNT <= 3 {
                            eprintln!("Debug line {}: {}", LINE_COUNT, self.buffer.trim());
                        }
                        if LINE_COUNT == 3 {
                            DEBUG_SHOWN = true;
                        }
                    }
                    
                    if let Some(pair) = parse_line(&self.buffer, &self.chr_map) {
                        unsafe {
                            PARSED_COUNT += 1;
                            if PARSED_COUNT <= 3 {
                                eprintln!("Debug: Parsed pair {}: chr{}:{} - chr{}:{}", 
                                         PARSED_COUNT, pair.chr1, pair.pos1, pair.chr2, pair.pos2);
                            }
                        }
                        return Some(Ok(pair));
                    }
                    // Invalid line, continue to next
                }
                Err(e) => return Some(Err(e.into())),
            }
        }
    }
}

fn parse_line(line: &str, chr_map: &ChrMap) -> Option<Pair> {
    let line = line.trim_end();
    
    // Split by whitespace (spaces, not tabs in this format)
    let fields: Vec<&str> = line.split_whitespace().collect();
    
    unsafe {
        static mut PARSE_ATTEMPT_COUNT: u64 = 0;
        PARSE_ATTEMPT_COUNT += 1;
        if PARSE_ATTEMPT_COUNT <= 3 {
            eprintln!("Debug parse_line {}: {} fields, line='{}'", PARSE_ATTEMPT_COUNT, fields.len(), line);
        }
    }
    
    if fields.len() < 9 {
        unsafe {
            static mut FIELD_ERROR_COUNT: u64 = 0;
            FIELD_ERROR_COUNT += 1;
            if FIELD_ERROR_COUNT <= 3 {
                eprintln!("Debug: Line rejected - only {} fields (need 9+)", fields.len());
            }
        }
        return None;
    }
    
    // Actual field mapping from demo.txt analysis:
    // Field 1: chr1, Field 2: pos1, Field 3: frag1, Field 4: str1
    // Field 5: chr2, Field 6: pos2, Field 7: frag2, Field 8: str2  
    // Field 9: mapq1, Field 12: mapq2
    // Original script: chr1=$2, pos1=$3, frag1=$4, chr2=$6, pos2=$7, frag2=$8, mapq1=$9, mapq2=$12
    
    let chr1_str = fields[1];   // Field 2 in awk = index 1
    let pos1_str = fields[2];   // Field 3 in awk = index 2
    let frag1_str = fields[3];  // Field 4 in awk = index 3
    let chr2_str = fields[5];   // Field 6 in awk = index 5
    let pos2_str = fields[6];   // Field 7 in awk = index 6
    let frag2_str = fields[7];  // Field 8 in awk = index 7
    let mapq1_str = fields[8];  // Field 9 in awk = index 8
    let mapq2_str = if fields.len() > 11 { fields[11] } else { "0" }; // Field 12 in awk = index 11
    
    // Parse values with detailed error reporting
    let chr1 = match chr_map.get(chr1_str).copied() {
        Some(c) => c,
        None => {
            unsafe {
                static mut CHR_ERROR_COUNT: u64 = 0;
                CHR_ERROR_COUNT += 1;
                if CHR_ERROR_COUNT <= 3 {
                    eprintln!("Debug: Failed to parse chr1='{}' (not in map)", chr1_str);
                }
            }
            return None;
        }
    };
    
    let pos1 = match pos1_str.parse() {
        Ok(p) => p,
        Err(_) => {
            unsafe {
                static mut POS1_ERROR_COUNT: u64 = 0;
                POS1_ERROR_COUNT += 1;
                if POS1_ERROR_COUNT <= 3 {
                    eprintln!("Debug: Failed to parse pos1='{}'", pos1_str);
                }
            }
            return None;
        }
    };
    
    let frag1 = frag1_str.parse::<u32>().ok()?;
    let chr2 = chr_map.get(chr2_str).copied()?;
    let pos2 = pos2_str.parse().ok()?;
    let frag2 = frag2_str.parse::<u32>().ok()?;
    let mapq1 = mapq1_str.parse::<u32>().ok()?;
    let mapq2 = mapq2_str.parse::<u32>().ok()?;
    
    // Debug: Show what we parsed for first few lines
    unsafe {
        static mut DEBUG_PARSE_COUNT: u64 = 0;
        DEBUG_PARSE_COUNT += 1;
        if DEBUG_PARSE_COUNT <= 3 {
            eprintln!("Debug parse {}: chr1={}, pos1={}, frag1={}, chr2={}, pos2={}, frag2={}, mapq1={}, mapq2={}", 
                     DEBUG_PARSE_COUNT, chr1_str, pos1_str, frag1_str, chr2_str, pos2_str, frag2_str, mapq1_str, mapq2_str);
            eprintln!("  Filter check: mapq1={}>0? mapq2={}>0? frag1={}!=frag2={}?", mapq1, mapq2, frag1, frag2);
        }
    }
    
    // Apply filters from original script: $9>0 && $12>0 && $4!=$8
    if mapq1 > 0 && mapq2 > 0 && frag1 != frag2 {
        unsafe {
            static mut ACCEPTED_COUNT: u64 = 0;
            ACCEPTED_COUNT += 1;
            if ACCEPTED_COUNT <= 3 {
                eprintln!("Debug: Accepted pair {}: chr{}:{} - chr{}:{} (mapq1={}, mapq2={}, frag1={}, frag2={})", 
                         ACCEPTED_COUNT, chr1, pos1, chr2, pos2, mapq1, mapq2, frag1, frag2);
            }
        }
        Some(Pair { chr1, pos1, chr2, pos2 })
    } else {
        unsafe {
            static mut FILTERED_COUNT: u64 = 0;
            FILTERED_COUNT += 1;
            if FILTERED_COUNT <= 10 {
                eprintln!("Debug: Filtered pair {} - mapq1={}>0? mapq2={}>0? frag1={}!=frag2={}?", 
                         FILTERED_COUNT, mapq1, mapq2, frag1, frag2);
            }
        }
        None
    }
}

pub fn open_file<R: Read>(reader: R) -> Result<PairIterator<BufReader<MultiGzDecoder<R>>>> {
    let decoder = MultiGzDecoder::new(reader);
    let buf_reader = BufReader::with_capacity(64 * 1024, decoder);
    let chr_map = crate::utils::create_chr_map();
    Ok(PairIterator::new(buf_reader, chr_map))
}

pub fn open_file_uncompressed<R: Read>(reader: R) -> Result<PairIterator<BufReader<R>>> {
    let buf_reader = BufReader::with_capacity(64 * 1024, reader);
    let chr_map = crate::utils::create_chr_map();
    Ok(PairIterator::new(buf_reader, chr_map))
}