use crate::utils::{ChrMap, Pair};
use anyhow::Result;
use flate2::read::MultiGzDecoder;
use std::io::Read;
use std::io::{BufRead, BufReader};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

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
        static LINE_COUNT: AtomicU64 = AtomicU64::new(0);
        static PARSED_COUNT: AtomicU64 = AtomicU64::new(0);
        static DEBUG_SHOWN: AtomicBool = AtomicBool::new(false);

        loop {
            self.buffer.clear();
            match self.reader.read_line(&mut self.buffer) {
                Ok(0) => {
                    let line_count = LINE_COUNT.load(Ordering::Relaxed);
                    let parsed_count = PARSED_COUNT.load(Ordering::Relaxed);
                    if line_count > 0 {
                        eprintln!(
                            "Debug: Processed {} lines, parsed {} pairs",
                            line_count, parsed_count
                        );
                    }
                    return None; // EOF
                }
                Ok(_) => {
                    let line_count = LINE_COUNT.fetch_add(1, Ordering::Relaxed) + 1;
                    if !DEBUG_SHOWN.load(Ordering::Relaxed) && line_count <= 3 {
                        eprintln!("Debug line {}: {}", line_count, self.buffer.trim());
                    }
                    if line_count == 3 {
                        DEBUG_SHOWN.store(true, Ordering::Relaxed);
                    }

                    if let Some(pair) = parse_line(&self.buffer, &self.chr_map) {
                        let parsed_count = PARSED_COUNT.fetch_add(1, Ordering::Relaxed) + 1;
                        if parsed_count <= 3 {
                            eprintln!(
                                "Debug: Parsed pair {}: chr{}:{} - chr{}:{}",
                                parsed_count, pair.chr1, pair.pos1, pair.chr2, pair.pos2
                            );
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

    static PARSE_ATTEMPT_COUNT: AtomicU64 = AtomicU64::new(0);
    let parse_count = PARSE_ATTEMPT_COUNT.fetch_add(1, Ordering::Relaxed) + 1;
    if parse_count <= 3 {
        eprintln!(
            "Debug parse_line {}: {} fields, line='{}'",
            parse_count,
            fields.len(),
            line
        );
    }

    if fields.len() < 9 {
        static FIELD_ERROR_COUNT: AtomicU64 = AtomicU64::new(0);
        let error_count = FIELD_ERROR_COUNT.fetch_add(1, Ordering::Relaxed) + 1;
        if error_count <= 3 {
            eprintln!(
                "Debug: Line rejected - only {} fields (need 9+)",
                fields.len()
            );
        }
        return None;
    }

    // Actual field mapping from demo.txt analysis:
    // Field 1: chr1, Field 2: pos1, Field 3: frag1, Field 4: str1
    // Field 5: chr2, Field 6: pos2, Field 7: frag2, Field 8: str2
    // Field 9: mapq1, Field 12: mapq2
    // Original script: chr1=$2, pos1=$3, frag1=$4, chr2=$6, pos2=$7, frag2=$8, mapq1=$9, mapq2=$12

    let chr1_str = fields[1]; // Field 2 in awk = index 1
    let pos1_str = fields[2]; // Field 3 in awk = index 2
    let frag1_str = fields[3]; // Field 4 in awk = index 3
    let chr2_str = fields[5]; // Field 6 in awk = index 5
    let pos2_str = fields[6]; // Field 7 in awk = index 6
    let frag2_str = fields[7]; // Field 8 in awk = index 7
    let mapq1_str = fields[8]; // Field 9 in awk = index 8
    let mapq2_str = if fields.len() > 11 { fields[11] } else { "0" }; // Field 12 in awk = index 11

    // Parse values with detailed error reporting
    let chr1 = match chr_map.get(chr1_str).copied() {
        Some(c) => c,
        None => {
            static CHR_ERROR_COUNT: AtomicU64 = AtomicU64::new(0);
            let error_count = CHR_ERROR_COUNT.fetch_add(1, Ordering::Relaxed) + 1;
            if error_count <= 3 {
                eprintln!("Debug: Failed to parse chr1='{}' (not in map)", chr1_str);
            }
            return None;
        }
    };

    let pos1 = match pos1_str.parse() {
        Ok(p) => p,
        Err(_) => {
            static POS1_ERROR_COUNT: AtomicU64 = AtomicU64::new(0);
            let error_count = POS1_ERROR_COUNT.fetch_add(1, Ordering::Relaxed) + 1;
            if error_count <= 3 {
                eprintln!("Debug: Failed to parse pos1='{}'", pos1_str);
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
    static DEBUG_PARSE_COUNT: AtomicU64 = AtomicU64::new(0);
    let debug_count = DEBUG_PARSE_COUNT.fetch_add(1, Ordering::Relaxed) + 1;
    if debug_count <= 3 {
        eprintln!("Debug parse {}: chr1={}, pos1={}, frag1={}, chr2={}, pos2={}, frag2={}, mapq1={}, mapq2={}", 
                 debug_count, chr1_str, pos1_str, frag1_str, chr2_str, pos2_str, frag2_str, mapq1_str, mapq2_str);
        eprintln!(
            "  Filter check: mapq1={}>0? mapq2={}>0? frag1={}!=frag2={}?",
            mapq1, mapq2, frag1, frag2
        );
    }

    // Apply filters from original script: $9>0 && $12>0 && $4!=$8
    if mapq1 > 0 && mapq2 > 0 && frag1 != frag2 {
        static ACCEPTED_COUNT: AtomicU64 = AtomicU64::new(0);
        let accepted_count = ACCEPTED_COUNT.fetch_add(1, Ordering::Relaxed) + 1;
        if accepted_count <= 3 {
            eprintln!("Debug: Accepted pair {}: chr{}:{} - chr{}:{} (mapq1={}, mapq2={}, frag1={}, frag2={})", 
                     accepted_count, chr1, pos1, chr2, pos2, mapq1, mapq2, frag1, frag2);
        }
        Some(Pair {
            chr1,
            pos1,
            chr2,
            pos2,
        })
    } else {
        static FILTERED_COUNT: AtomicU64 = AtomicU64::new(0);
        let filtered_count = FILTERED_COUNT.fetch_add(1, Ordering::Relaxed) + 1;
        if filtered_count <= 10 {
            eprintln!(
                "Debug: Filtered pair {} - mapq1={}>0? mapq2={}>0? frag1={}!=frag2={}?",
                filtered_count, mapq1, mapq2, frag1, frag2
            );
        }
        None
    }
}

pub fn open_file<R: Read>(
    reader: R,
    chrom_size_file: Option<&str>,
) -> Result<PairIterator<BufReader<MultiGzDecoder<R>>>> {
    let decoder = MultiGzDecoder::new(reader);
    let buf_reader = BufReader::with_capacity(64 * 1024, decoder);
    let chr_map = crate::utils::create_chr_map(chrom_size_file);
    Ok(PairIterator::new(buf_reader, chr_map))
}

pub fn open_file_uncompressed<R: Read>(
    reader: R,
    chrom_size_file: Option<&str>,
) -> Result<PairIterator<BufReader<R>>> {
    let buf_reader = BufReader::with_capacity(64 * 1024, reader);
    let chr_map = crate::utils::create_chr_map(chrom_size_file);
    Ok(PairIterator::new(buf_reader, chr_map))
}
