use crate::utils::{get_genome_lengths, Pair};
use rayon::prelude::*;

pub struct Coverage {
    pub bins: Vec<Vec<u32>>,
    pub bin_width: u32,
    pub chr_lengths: Vec<u32>,
}

impl Coverage {
    pub fn new(bin_width: u32, chrom_size_file: Option<&str>) -> Self {
        let chr_lengths = get_genome_lengths(chrom_size_file);
        let bins: Vec<Vec<u32>> = chr_lengths
            .iter()
            .map(|&len| {
                let num_bins = (len / bin_width) + 1;
                vec![0u32; num_bins as usize]
            })
            .collect();

        Self {
            bins,
            bin_width,
            chr_lengths,
        }
    }

    pub fn from_lengths(bin_width: u32, chr_lengths: Vec<u32>) -> Self {
        let bins: Vec<Vec<u32>> = chr_lengths
            .iter()
            .map(|&len| {
                let num_bins = (len / bin_width) + 1;
                vec![0u32; num_bins as usize]
            })
            .collect();

        Self {
            bins,
            bin_width,
            chr_lengths,
        }
    }

    pub fn increment(&mut self, chr: u8, pos: u32) {
        let chr_idx = (chr as usize).saturating_sub(1);
        if chr_idx >= self.bins.len() {
            return;
        }

        if pos >= self.chr_lengths[chr_idx] {
            return;
        }

        let bin_idx = (pos / self.bin_width) as usize;
        if bin_idx < self.bins[chr_idx].len() {
            let v = &mut self.bins[chr_idx][bin_idx];
            *v = v.saturating_add(1);
        }
    }

    pub fn add_pair(&mut self, pair: &Pair) {
        self.increment(pair.chr1, pair.pos1);
        self.increment(pair.chr2, pair.pos2);
    }

    pub fn get_counts(&self, bin_size: u32) -> Vec<Vec<u32>> {
        let bins_per_chunk = bin_size / self.bin_width;

        self.bins
            .par_iter()
            .map(|chr_bins| {
                let mut result = Vec::new();
                let mut i = 0;

                while i < chr_bins.len() {
                    let end = (i + bins_per_chunk as usize).min(chr_bins.len());
                    let sum: u32 = chr_bins[i..end].iter().copied().sum();
                    result.push(sum);
                    i += bins_per_chunk as usize;
                }

                result
            })
            .collect()
    }

    pub fn count_good_bins(&self, bin_size: u32, threshold: u32) -> u64 {
        let bins_per_chunk = bin_size / self.bin_width;

        // For very large bin sizes, use optimized approach
        if bins_per_chunk >= 1000 {
            return self.count_good_bins_large(bin_size, threshold);
        }

        self.bins
            .par_iter()
            .map(|chr_bins| {
                let mut count = 0u64;
                let chunk_size = bins_per_chunk as usize;

                // Process in chunks
                for chunk in chr_bins.chunks(chunk_size) {
                    let sum: u32 = chunk.iter().copied().sum();

                    if sum >= threshold {
                        count += 1;
                    }
                }

                count
            })
            .sum()
    }

    // Optimized version for large bin sizes
    fn count_good_bins_large(&self, bin_size: u32, threshold: u32) -> u64 {
        let bins_per_chunk = bin_size / self.bin_width;

        self.bins
            .par_iter()
            .map(|chr_bins| {
                let mut count = 0u64;
                let chunk_size = bins_per_chunk as usize;
                let num_chunks = (chr_bins.len() + chunk_size - 1) / chunk_size;

                // Use batched processing for better cache performance
                for chunk_idx in 0..num_chunks {
                    let start = chunk_idx * chunk_size;
                    let end = (start + chunk_size).min(chr_bins.len());

                    let sum: u32 = chr_bins[start..end].iter().copied().sum();

                    if sum >= threshold {
                        count += 1;
                    }
                }

                count
            })
            .sum()
    }

    pub fn total_genome_size(&self) -> u64 {
        self.chr_lengths.iter().map(|&x| x as u64).sum()
    }

    pub fn get_total_contacts(&self) -> u64 {
        self.bins
            .par_iter()
            .map(|chr_bins| {
                chr_bins.iter().map(|&x| x as u64).sum::<u64>()
            })
            .sum()
    }

    pub fn get_non_zero_bins(&self) -> u64 {
        self.bins
            .par_iter()
            .map(|chr_bins| {
                chr_bins.iter().filter(|&&x| x > 0).count() as u64
            })
            .sum()
    }
}
