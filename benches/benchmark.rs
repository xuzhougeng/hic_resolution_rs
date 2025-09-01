use criterion::{black_box, criterion_group, criterion_main, Criterion};
use hic_resolution_rs::coverage::Coverage;
use hic_resolution_rs::utils::Pair;

fn benchmark_coverage_build(c: &mut Criterion) {
    c.bench_function("coverage_build_1M_pairs", |b| {
        b.iter(|| {
            let coverage = Coverage::new(50);

            // Simulate 1M pairs
            for i in 0..1_000_000 {
                let pair = Pair {
                    chr1: ((i % 22) + 1) as u8,
                    pos1: (i * 1000) % 100_000_000,
                    chr2: ((i % 22) + 1) as u8,
                    pos2: ((i * 1000) + 500) % 100_000_000,
                };
                coverage.add_pair(&pair);
            }

            black_box(coverage)
        })
    });
}

fn benchmark_resolution_search(c: &mut Criterion) {
    // Pre-build coverage with some data
    let coverage = Coverage::new(50);
    for i in 0..100_000 {
        let pair = Pair {
            chr1: ((i % 22) + 1) as u8,
            pos1: (i * 1000) % 100_000_000,
            chr2: ((i % 22) + 1) as u8,
            pos2: ((i * 1000) + 500) % 100_000_000,
        };
        coverage.add_pair(&pair);
    }

    c.bench_function("resolution_search", |b| {
        b.iter(|| {
            hic_resolution_rs::resolution::find_resolution(
                black_box(&coverage),
                black_box(0.8),
                black_box(1000),
                black_box(1000),
            )
        })
    });
}

criterion_group!(
    benches,
    benchmark_coverage_build,
    benchmark_resolution_search
);
criterion_main!(benches);
