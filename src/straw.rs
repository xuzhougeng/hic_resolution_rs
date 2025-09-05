use anyhow::{anyhow, Context, Result};
use flate2::read::ZlibDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::collections::BTreeMap;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

// Magic string for slice files (no NUL terminator)
const HICSLICE_MAGIC: &[u8] = b"HICSLICE";

// Minimal structures
#[derive(Clone, Debug)]
struct IndexEntry { size: i64, position: i64 }

#[derive(Clone, Debug)]
struct Chromosome { name: String, index: i32, length: i64 }

#[allow(dead_code)]
#[derive(Debug)]
struct HicFile {
    file: BufReader<File>,
    version: i32,
    master: i64,
    genome_id: String,
    nvi_pos: i64,
    nvi_len: i64,
    chromosomes: Vec<Chromosome>,
    resolutions: Vec<i32>,
    path: PathBuf,
}

impl HicFile {
    fn open(path: &Path) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("Open {:?}", path))?;
        let mut reader = BufReader::new(file);
        if !read_magic(&mut reader)? { return Err(anyhow!("Not a .hic file: missing HIC magic")); }
        let version = read_i32(&mut reader)?;
        if version < 6 { return Err(anyhow!("Unsupported .hic version {} (<6)", version)); }
        let master = read_i64(&mut reader)?;
        let genome_id = read_cstring(&mut reader)?;
        let (nvi_pos, nvi_len) = if version > 8 { (read_i64(&mut reader)?, read_i64(&mut reader)?) } else { (0, 0) };
        let nattr = read_i32(&mut reader)?;
        for _ in 0..nattr { let _ = read_cstring(&mut reader)?; let _ = read_cstring(&mut reader)?; }
        let num_chromosomes = read_i32(&mut reader)? as usize;
        let mut chromosomes = Vec::with_capacity(num_chromosomes);
        for i in 0..num_chromosomes {
            let name = read_cstring(&mut reader)?;
            let length = if version > 8 { read_i64(&mut reader)? } else { read_i32(&mut reader)? as i64 };
            chromosomes.push(Chromosome { name, index: i as i32, length });
        }
        let nbp = read_i32(&mut reader)? as usize;
        let mut resolutions = Vec::with_capacity(nbp);
        for _ in 0..nbp { resolutions.push(read_i32(&mut reader)?); }
        Ok(HicFile { file: reader, version, master, genome_id, nvi_pos, nvi_len, chromosomes, resolutions, path: path.to_path_buf() })
    }

    fn get_matrix_zoom_data(&mut self, chr1_idx: i32, chr2_idx: i32, unit: &str, resolution: i32) -> Result<Option<MatrixZoomData>> {
        let (c1, c2) = if chr1_idx <= chr2_idx { (chr1_idx, chr2_idx) } else { (chr2_idx, chr1_idx) };
        self.file.seek(SeekFrom::Start(self.master as u64))?;
        if self.version > 8 { let _ = read_i64(&mut self.file)?; } else { let _ = read_i32(&mut self.file)?; }
        let key = format!("{}_{}", c1, c2);
        let nentries = read_i32(&mut self.file)?;
        let mut my_file_pos: Option<i64> = None;
        for _ in 0..nentries {
            let k = read_cstring(&mut self.file)?;
            let fpos = read_i64(&mut self.file)?;
            let _size = read_i32(&mut self.file)?;
            if k == key { my_file_pos = Some(fpos); }
        }
        let my_file_pos = match my_file_pos { Some(p) => p, None => return Ok(None) };
        let (block_map, sum_counts, block_bin_count, block_col_count) = read_matrix(&mut self.file, my_file_pos, unit, resolution)?;
        Ok(Some(MatrixZoomData {
            version: self.version,
            resolution,
            is_intra: c1 == c2,
            num_bins1: (self.chromosomes[c1 as usize].length / resolution as i64) as i32,
            num_bins2: (self.chromosomes[c2 as usize].length / resolution as i64) as i32,
            block_map,
            sum_counts,
            block_bin_count,
            block_col_count,
            c1,
            c2,
        }))
    }
}

#[allow(dead_code)]
#[derive(Debug)]
struct MatrixZoomData {
    version: i32,
    resolution: i32,
    is_intra: bool,
    num_bins1: i32,
    num_bins2: i32,
    sum_counts: f32,
    block_bin_count: i32,
    block_col_count: i32,
    block_map: BTreeMap<i32, IndexEntry>,
    c1: i32,
    c2: i32,
}

fn read_matrix<R: Read + Seek>(r: &mut R, my_file_pos: i64, unit: &str, resolution: i32) -> Result<(BTreeMap<i32, IndexEntry>, f32, i32, i32)> {
    r.seek(SeekFrom::Start(my_file_pos as u64))?;
    let _c1 = read_i32(r)?;
    let _c2 = read_i32(r)?;
    let nres = read_i32(r)?;
    let mut found = false;
    let mut block_map = BTreeMap::new();
    let mut sum_counts = 0f32;
    let mut block_bin_count = 0;
    let mut block_col_count = 0;
    for _ in 0..nres {
        let (bm, sum, bbc, bcc, is_match) = read_matrix_zoom_data(r, unit, resolution)?;
        if is_match {
            block_map = bm; sum_counts = sum; block_bin_count = bbc; block_col_count = bcc; found = true; break;
        }
    }
    if !found { return Err(anyhow!("Resolution {} at unit {} not found in matrix", resolution, unit)); }
    Ok((block_map, sum_counts, block_bin_count, block_col_count))
}

fn read_matrix_zoom_data<R: Read + Seek>(r: &mut R, my_unit: &str, my_binsize: i32) -> Result<(BTreeMap<i32, IndexEntry>, f32, i32, i32, bool)> {
    let unit = read_cstring(r)?;
    let _old_zoom = read_i32(r)?;
    let sum_counts = read_f32(r)?;
    let _occupied = read_f32(r)?;
    let _stddev = read_f32(r)?;
    let _p95 = read_f32(r)?;
    let bin_size = read_i32(r)?;
    let block_bin_count = read_i32(r)?;
    let block_col_count = read_i32(r)?;
    let is_match = unit == my_unit && bin_size == my_binsize;
    let nblocks = read_i32(r)?;
    let mut block_map = BTreeMap::new();
    if is_match {
        for _ in 0..nblocks {
            let block_number = read_i32(r)?;
            let file_position = read_i64(r)?;
            let size_in_bytes = read_i32(r)? as i64;
            block_map.insert(block_number, IndexEntry { size: size_in_bytes, position: file_position });
        }
    } else {
        let skip = nblocks as i64 * (4 + 8 + 4) as i64;
        r.seek(SeekFrom::Current(skip))?;
    }
    Ok((block_map, sum_counts, block_bin_count, block_col_count, is_match))
}

#[derive(Clone, Debug)]
struct ContactRecord { bin_x: i32, bin_y: i32, counts: f32 }

fn read_block(path: &Path, idx: &IndexEntry, version: i32) -> Result<Vec<ContactRecord>> {
    if idx.size <= 0 { return Ok(Vec::new()); }
    let mut f = File::open(path).with_context(|| format!("Open {:?}", path))?;
    let mut comp = vec![0u8; idx.size as usize];
    f.seek(SeekFrom::Start(idx.position as u64))?;
    f.read_exact(&mut comp)?;
    let mut dec = ZlibDecoder::new(&comp[..]);
    let mut buf = Vec::new();
    dec.read_to_end(&mut buf)?;
    let mut cur = std::io::Cursor::new(buf);

    let n_records = read_i32(&mut cur)? as usize;
    let mut out = Vec::with_capacity(n_records);
    if version < 7 {
        for _ in 0..n_records {
            let bin_x = read_i32(&mut cur)?;
            let bin_y = read_i32(&mut cur)?;
            let counts = read_f32(&mut cur)?;
            out.push(ContactRecord { bin_x, bin_y, counts });
        }
        return Ok(out);
    }

    let bin_x_offset = read_i32(&mut cur)?;
    let bin_y_offset = read_i32(&mut cur)?;
    let use_short = read_u8(&mut cur)? == 0;
    let mut use_short_bin_x = true;
    let mut use_short_bin_y = true;
    if version > 8 {
        use_short_bin_x = read_u8(&mut cur)? == 0;
        use_short_bin_y = read_u8(&mut cur)? == 0;
    }
    let typ = read_u8(&mut cur)?;
    match typ {
        1 => {
            if use_short_bin_x && use_short_bin_y {
                let row_count = read_i16(&mut cur)? as i32;
                for _ in 0..row_count {
                    let bin_y = bin_y_offset + read_i16(&mut cur)? as i32;
                    let col_count = read_i16(&mut cur)? as i32;
                    for _ in 0..col_count {
                        let bin_x = bin_x_offset + read_i16(&mut cur)? as i32;
                        let counts = if use_short { read_i16(&mut cur)? as f32 } else { read_f32(&mut cur)? };
                        out.push(ContactRecord { bin_x, bin_y, counts });
                    }
                }
            } else if use_short_bin_x && !use_short_bin_y {
                let row_count = read_i32(&mut cur)?;
                for _ in 0..row_count {
                    let bin_y = bin_y_offset + read_i32(&mut cur)?;
                    let col_count = read_i16(&mut cur)? as i32;
                    for _ in 0..col_count {
                        let bin_x = bin_x_offset + read_i16(&mut cur)? as i32;
                        let counts = if use_short { read_i16(&mut cur)? as f32 } else { read_f32(&mut cur)? };
                        out.push(ContactRecord { bin_x, bin_y, counts });
                    }
                }
            } else if !use_short_bin_x && use_short_bin_y {
                let row_count = read_i16(&mut cur)? as i32;
                for _ in 0..row_count {
                    let bin_y = bin_y_offset + read_i16(&mut cur)? as i32;
                    let col_count = read_i32(&mut cur)?;
                    for _ in 0..col_count {
                        let bin_x = bin_x_offset + read_i32(&mut cur)?;
                        let counts = if use_short { read_i16(&mut cur)? as f32 } else { read_f32(&mut cur)? };
                        out.push(ContactRecord { bin_x, bin_y, counts });
                    }
                }
            } else {
                let row_count = read_i32(&mut cur)?;
                for _ in 0..row_count {
                    let bin_y = bin_y_offset + read_i32(&mut cur)?;
                    let col_count = read_i32(&mut cur)?;
                    for _ in 0..col_count {
                        let bin_x = bin_x_offset + read_i32(&mut cur)?;
                        let counts = if use_short { read_i16(&mut cur)? as f32 } else { read_f32(&mut cur)? };
                        out.push(ContactRecord { bin_x, bin_y, counts });
                    }
                }
            }
        }
        2 => {
            let n_pts = read_i32(&mut cur)?;
            let w = read_i16(&mut cur)? as i32;
            for i in 0..n_pts {
                let row = i / w;
                let col = i - row * w;
                let bin_x = bin_x_offset + col;
                let bin_y = bin_y_offset + row;
                if use_short {
                    let c = read_i16(&mut cur)?;
                    if c != -32768 { out.push(ContactRecord { bin_x, bin_y, counts: c as f32 }); }
                } else {
                    let counts = read_f32(&mut cur)?;
                    if !counts.is_nan() { out.push(ContactRecord { bin_x, bin_y, counts }); }
                }
            }
        }
        _ => {}
    }
    Ok(out)
}

pub fn dump_hic_genome_wide(input: &Path, binsize: i32, output: &Path) -> Result<()> {
    let mut hic = HicFile::open(input)?;
    // Build chromosome keys (skip index <= 0 per C++ code)
    let mut chr_keys: BTreeMap<String, i16> = BTreeMap::new();
    let mut key_counter: i16 = 0;
    for chr in &hic.chromosomes {
        if chr.index > 0 { chr_keys.insert(chr.name.clone(), key_counter); key_counter += 1; }
    }

    // Open output .slc.gz
    let out = File::create(output).with_context(|| format!("Create {:?}", output))?;
    let mut enc = GzEncoder::new(BufWriter::new(out), Compression::default());

    // Write header
    enc.write_all(HICSLICE_MAGIC)?;
    enc.write_all(&(binsize as i32).to_le_bytes())?;
    enc.write_all(&(chr_keys.len() as i32).to_le_bytes())?;
    for (name, key) in &chr_keys {
        let nb = name.as_bytes();
        enc.write_all(&(nb.len() as i32).to_le_bytes())?;
        enc.write_all(nb)?;
        enc.write_all(&(*key as i16).to_le_bytes())?;
    }

    // Iterate chromosome pairs
    let n = hic.chromosomes.len();
    for i in 0..n {
        let c1_idx = hic.chromosomes[i].index;
        if c1_idx <= 0 { continue; }
        for j in i..n {
            let c2_idx = hic.chromosomes[j].index;
            if c2_idx <= 0 { continue; }
            if let Some(mzd) = hic.get_matrix_zoom_data(c1_idx, c2_idx, "BP", binsize)? {
                let key1 = *chr_keys.get(&hic.chromosomes[mzd.c1 as usize].name).unwrap();
                let key2 = *chr_keys.get(&hic.chromosomes[mzd.c2 as usize].name).unwrap();
                for (_, idx) in mzd.block_map.iter() {
                    let records = read_block(&hic.path, idx, mzd.version)?;
                    for rec in records {
                        if rec.counts > 0.0 && rec.counts.is_finite() {
                            enc.write_all(&key1.to_le_bytes())?;
                            enc.write_all(&rec.bin_x.to_le_bytes())?;
                            enc.write_all(&key2.to_le_bytes())?;
                            enc.write_all(&rec.bin_y.to_le_bytes())?;
                            enc.write_all(&rec.counts.to_le_bytes())?;
                        }
                    }
                }
            }
        }
    }

    enc.finish()?.flush()?;
    Ok(())
}

// ----------------- low-level readers -----------------
fn read_magic<R: Read>(r: &mut R) -> Result<bool> { let s = read_cstring(r)?; Ok(s.starts_with("HIC")) }
fn read_u8<R: Read>(r: &mut R) -> Result<u8> { let mut b=[0u8;1]; r.read_exact(&mut b)?; Ok(b[0]) }
fn read_i16<R: Read>(r: &mut R) -> Result<i16> { let mut b=[0u8;2]; r.read_exact(&mut b)?; Ok(i16::from_le_bytes(b)) }
fn read_i32<R: Read>(r: &mut R) -> Result<i32> { let mut b=[0u8;4]; r.read_exact(&mut b)?; Ok(i32::from_le_bytes(b)) }
fn read_i64<R: Read>(r: &mut R) -> Result<i64> { let mut b=[0u8;8]; r.read_exact(&mut b)?; Ok(i64::from_le_bytes(b)) }
fn read_f32<R: Read>(r: &mut R) -> Result<f32> { let mut b=[0u8;4]; r.read_exact(&mut b)?; Ok(f32::from_le_bytes(b)) }
fn _read_f64<R: Read>(r: &mut R) -> Result<f64> { let mut b=[0u8;8]; r.read_exact(&mut b)?; Ok(f64::from_le_bytes(b)) }
fn read_cstring<R: Read>(r: &mut R) -> Result<String> {
    let mut buf = Vec::new();
    let mut byte = [0u8;1];
    loop {
        r.read_exact(&mut byte)?;
        if byte[0] == 0 { break; }
        buf.push(byte[0]);
    }
    Ok(String::from_utf8(buf).unwrap_or_default())
}

pub fn list_hic_chromosomes(input: &Path) -> Result<()> {
    let hic = HicFile::open(input)?;
    // Print available BP resolutions
    let mut res = hic.resolutions.clone();
    res.sort_unstable();
    println!("# Resolutions (BP): {}", res.iter().map(|r| r.to_string()).collect::<Vec<_>>().join(", "));
    // Print chromosomes table
    println!("# Chromosomes (name\tlength)");
    for chr in hic.chromosomes.iter() {
        if chr.index > 0 {
            println!("{}\t{}", chr.name, chr.length);
        }
    }
    Ok(())
}

pub fn effres_hic(input: &Path, chrom_req: Option<&str>, thr: i32, pct: f64) -> Result<()> {
    let mut hic = HicFile::open(input)?;
    // If no chromosome provided, compute min/mean/max coverage across chromosomes per resolution
    if chrom_req.is_none() {
        println!("# File: {}", input.display());
        println!("# Mode: all chromosomes coverage summary");
        println!("# Filters: length >= 2,500,000 bp; exclude no-signal contigs per resolution");
        println!("# Threshold per bin: {} contacts", thr);
        println!("resolution_bp\tmin_cov\tmean_cov\tmax_cov");

        let mut resolutions = hic.resolutions.clone();
        resolutions.sort_unstable();

        // Collect usable chromosome indices: index>0 and length >= 2,500,000 bp
        let chr_idxs: Vec<i32> = hic
            .chromosomes
            .iter()
            .filter(|c| c.index > 0 && c.length >= 2_500_000)
            .map(|c| c.index)
            .collect();

        for res in resolutions {
            let mut covs: Vec<f64> = Vec::with_capacity(chr_idxs.len());
            for &ci in &chr_idxs {
                let cov_opt = match hic.get_matrix_zoom_data(ci, ci, "BP", res)? {
                    None => None,
                    Some(mzd) => {
                        let mut counts: HashMap<i32, f64> = HashMap::new();
                        for (_, idx) in mzd.block_map.iter() {
                            let records = read_block(&hic.path, idx, mzd.version)?;
                            for rec in records {
                                *counts.entry(rec.bin_x).or_insert(0.0) += rec.counts as f64;
                                *counts.entry(rec.bin_y).or_insert(0.0) += rec.counts as f64;
                            }
                        }
                        if counts.is_empty() {
                            None // exclude no-signal contig for this resolution
                        } else {
                            let covered = counts.values().filter(|&&v| v >= thr as f64).count();
                            Some(covered as f64 / counts.len() as f64)
                        }
                    }
                };
                if let Some(cov) = cov_opt { covs.push(cov); }
            }
            if covs.is_empty() {
                println!("{}\t{:.3}\t{:.3}\t{:.3}", res, 0.0, 0.0, 0.0);
            } else {
                let min = covs
                    .iter()
                    .copied()
                    .fold(f64::INFINITY, f64::min);
                let max = covs
                    .iter()
                    .copied()
                    .fold(f64::NEG_INFINITY, f64::max);
                let mean = covs.iter().sum::<f64>() / (covs.len() as f64);
                println!("{}\t{:.3}\t{:.3}\t{:.3}", res, min, mean, max);
            }
        }
        return Ok(());
    }

    // Single chromosome path (original: resolution vs coverage and effective resolution thresholding)
    // Collect available chromosomes with owned names to avoid borrowing conflicts
    let avail: Vec<(String, i32)> = hic
        .chromosomes
        .iter()
        .filter(|c| c.index > 0)
        .map(|c| (c.name.clone(), c.index))
        .collect();
    // Flexible name matching: case-insensitive, optional "chr" prefix
    let req_s = chrom_req.unwrap().to_lowercase();
    let req_trim = req_s.trim_start_matches("chr").to_string();
    let mut c_idx_opt: Option<i32> = None;
    for (name, idx) in &avail {
        let nm = name.to_lowercase();
        if nm == req_s || nm.trim_start_matches("chr") == &req_trim {
            c_idx_opt = Some(*idx);
            break;
        }
    }
    let c_idx = match c_idx_opt {
        Some(i) => i,
        None => {
            eprintln!(
                "[ERROR] 未找到染色体 '{}', 可选值: {}",
                chrom_req.unwrap(),
                avail
                    .iter()
                    .map(|(n, _)| n.as_str())
                    .collect::<Vec<_>>()
                    .join(", ")
            );
            return Ok(());
        }
    };

    println!("# File: {}", input.display());
    let cname = hic.chromosomes[c_idx as usize].name.clone();
    println!("# Chromosome: {}", cname);
    println!("# Threshold per bin: {} contacts", thr);
    println!("# Required coverage: {:.1}% bins\n", pct * 100.0);
    println!("resolution_bp\tcoverage");

    let mut resolutions = hic.resolutions.clone();
    resolutions.sort_unstable();
    let mut eff_res: Option<i32> = None;
    for res in resolutions {
        match hic.get_matrix_zoom_data(c_idx, c_idx, "BP", res)? {
            None => {
                println!("{}\t{:.3}", res, 0.0);
            }
            Some(mzd) => {
                // Accumulate per-bin counts using a sparse map to mirror the Python reference
                let mut counts: HashMap<i32, f64> = HashMap::new();
                for (_, idx) in mzd.block_map.iter() {
                    let records = read_block(&hic.path, idx, mzd.version)?;
                    for rec in records {
                        *counts.entry(rec.bin_x).or_insert(0.0) += rec.counts as f64;
                        *counts.entry(rec.bin_y).or_insert(0.0) += rec.counts as f64;
                    }
                }
                let mut cov = 0.0f64;
                if !counts.is_empty() {
                    let covered = counts.values().filter(|&&v| v >= thr as f64).count();
                    cov = covered as f64 / counts.len() as f64;
                }
                println!("{}\t{:.3}", res, cov);
                if eff_res.is_none() && cov >= pct {
                    eff_res = Some(res);
                }
            }
        }
    }

    if let Some(r) = eff_res {
        println!(
            "\nEffective resolution on {}: {} bp (≥{:.0}% bins ≥ {} contacts)",
            cname, r, pct * 100.0, thr
        );
    } else {
        println!(
            "\nNo resolution met the {:.0}% / {} contacts criterion.",
            pct * 100.0, thr
        );
    }
    Ok(())
}
