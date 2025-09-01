use crate::coverage::Coverage;

pub fn find_resolution(
    coverage: &Coverage,
    prop: f64,
    count_threshold: u32,
    step_size: u32,
) -> u32 {
    let genome_size = coverage.total_genome_size();

    let mut low = coverage.bin_width;
    let mut high = coverage.bin_width;

    println!("Starting resolution search...");
    println!("Genome size: {} bp", genome_size);

    // Analyze data sparsity to set reasonable bounds
    let total_contacts = coverage.get_total_contacts();
    let non_zero_bins = coverage.get_non_zero_bins();
    let total_base_bins = genome_size / coverage.bin_width as u64;

    println!("Data analysis:");
    println!("  Total contacts: {}", total_contacts);
    println!(
        "  Non-zero 50bp bins: {} / {} ({:.2}%)",
        non_zero_bins,
        total_base_bins,
        non_zero_bins as f64 * 100.0 / total_base_bins as f64
    );

    // If data is very sparse, adjust search strategy
    let sparsity = non_zero_bins as f64 / total_base_bins as f64;
    let adjusted_step_size = if sparsity < 0.01 {
        println!(
            "  Detected sparse data ({:.4}% coverage), using larger step size",
            sparsity * 100.0
        );
        step_size * 10
    } else {
        step_size
    };

    // Find reasonable upper bound with large steps, but limit maximum
    let max_reasonable_bin = 10_000_000; // 10 Mb maximum
    let limit = max_reasonable_bin
        .min((genome_size.min(u64::from(u32::MAX))) as u32);
    let mut iteration = 0;
    let mut found_upper = false;

    loop {
        iteration += 1;
        if iteration % 10 == 0 {
            println!(
                "  Coarse search iteration {}: testing bin size {}",
                iteration, high
            );
        }

        let good_bins = coverage.count_good_bins(high, count_threshold);
        let total_bins = genome_size / high as u64;
        let required_bins = (prop * total_bins as f64) as u64;

        if iteration <= 5 {
            println!(
                "  Bin size: {}, Good bins: {}, Total bins: {}, Required: {}",
                high, good_bins, total_bins, required_bins
            );
        }

        if good_bins >= required_bins {
            println!(
                "Found upper bound: {} bp (good bins: {}/{})",
                high, good_bins, total_bins
            );
            found_upper = true;
            break;
        }

        // If we've gone too far without finding an upper bound, stop
        if high >= limit {
            println!(
                "Warning: Reached search limit ({} bp) without meeting requirement.",
                limit
            );
            // Ensure 'high' is within limit and aligned to bin multiple
            high = round_to_bin_multiple(limit, coverage.bin_width);
            break;
        }

        low = high;
        // Increase and align to multiple of base bin width
        let mut next = high.saturating_add(adjusted_step_size);
        next = round_to_bin_multiple(next, coverage.bin_width);
        if next == high { // avoid stalling if step < bin width
            next = next.saturating_add(coverage.bin_width);
        }
        high = next;
    }

    if !found_upper {
        println!(
            "Error: No bin size up to {} bp satisfies >= {:.1}% bins with >= {} contacts.",
            high,
            prop * 100.0,
            count_threshold
        );
        println!(
            "Returning upper limit ({} bp). Result does not satisfy the target proportion.",
            high
        );
        return high;
    }

    println!("Binary search range: {} - {} bp", low, high);

    // Binary search for exact resolution
    let mut binary_iteration = 0;
    while high > low + coverage.bin_width {
        binary_iteration += 1;
        let mid = round_to_bin_multiple(low + (high - low) / 2, coverage.bin_width);

        if binary_iteration % 5 == 0 || binary_iteration <= 3 {
            println!(
                "  Binary search iteration {}: testing {}",
                binary_iteration, mid
            );
        }

        let good_bins = coverage.count_good_bins(mid, count_threshold);
        let total_bins = genome_size / mid as u64;
        let required_bins = (prop * total_bins as f64) as u64;

        if good_bins >= required_bins {
            high = mid;
            if binary_iteration <= 3 {
                println!(
                    "    Success: {} good bins >= {} required",
                    good_bins, required_bins
                );
            }
        } else {
            low = mid;
            if binary_iteration <= 3 {
                println!(
                    "    Failed: {} good bins < {} required",
                    good_bins, required_bins
                );
            }
        }

        // Safety check to prevent infinite loop
        if binary_iteration > 100 {
            println!(
                "Warning: Binary search taking too long, stopping at iteration {}",
                binary_iteration
            );
            break;
        }
    }

    println!("Final resolution: {} bp", high);
    high
}

fn round_to_bin_multiple(value: u32, bin_width: u32) -> u32 {
    ((value + bin_width - 1) / bin_width) * bin_width
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_round_to_bin_multiple() {
        assert_eq!(round_to_bin_multiple(75, 50), 100);
        assert_eq!(round_to_bin_multiple(100, 50), 100);
        assert_eq!(round_to_bin_multiple(125, 50), 150);
        assert_eq!(round_to_bin_multiple(1, 50), 50);
    }
}
