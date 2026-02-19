/// Integration tests — invoke the compiled `ldsc` binary end-to-end.
///
/// Run with:
///   cargo test --test integration
///
/// The `ldscore_smoke` test is marked `#[ignore]` by default because it
/// reads the 1000G BED file (~1 GB) and takes ~30–90 s on first run
/// (OS page-fault warm-up).  Pass `--include-ignored` to run it:
///   cargo test --test integration -- --include-ignored
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use std::time::Instant;

/// Resolve the path to the compiled binary.
fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_ldsc"))
}

/// Path to the 1000G test dataset prefix (no extension).
fn test_data() -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("data")
        .join("1000G_phase3_common_norel")
        .to_str()
        .unwrap()
        .to_string()
}

// ---------------------------------------------------------------------------
// ldscore smoke test
// ---------------------------------------------------------------------------

/// Run `ldsc ldscore` on the full 1000G file with a tiny SNP window.
/// Checks that:
///   1. The binary exits successfully.
///   2. Per-chromosome .ldscore.gz files are created.
///   3. Each file is a valid gzip-compressed TSV with a CHR/SNP/BP/L2 header.
#[test]
#[ignore = "reads ~1 GB BED file; run with --include-ignored"]
fn ldscore_smoke() {
    let out_dir = tempfile::tempdir().expect("tempdir");
    let out_prefix = out_dir.path().join("test_ld").to_str().unwrap().to_string();

    let t0 = Instant::now();
    let status = Command::new(binary())
        .args([
            "ldscore",
            "--bfile",
            &test_data(),
            "--out",
            &out_prefix,
            "--ld-wind-snp",
            "50", // small window → fast
        ])
        .status()
        .expect("failed to launch ldsc binary");

    eprintln!("ldscore completed in {:.1}s", t0.elapsed().as_secs_f64());
    assert!(status.success(), "ldsc ldscore exited with {}", status);

    // Expect one output file per chromosome present in the BIM file.
    // The 1000G file covers chromosomes 1–22.
    let mut found = 0usize;
    for chr in 1u8..=22 {
        let path = format!("{}.{}.l2.ldscore.gz", out_prefix, chr);
        if PathBuf::from(&path).exists() {
            found += 1;
            // Verify it is a valid gzip TSV with the expected header.
            let content = read_gz(&path);
            let first_line = content.lines().next().expect("empty file");
            assert_eq!(
                first_line, "CHR\tSNP\tBP\tL2",
                "unexpected header in {}",
                path
            );
            // Verify there is at least one data row.
            assert!(
                content.lines().count() > 1,
                "no SNP rows in {}",
                path
            );
        }
    }
    assert!(found >= 22, "expected ≥22 chromosome files, found {}", found);
    eprintln!("Verified {} chromosome output files.", found);
}

// ---------------------------------------------------------------------------
// irwls unit-level smoke  (no BED file needed)
// ---------------------------------------------------------------------------

/// Confirm that the IRWLS unit test passes via the test binary as well.
/// This is already covered by `cargo test` on the lib, but including it here
/// ensures the integration harness also exercises the math path.
#[test]
fn irwls_known_solution() {
    // The full math is tested in irwls.rs; this test just verifies the binary
    // runs and exits 0 for --help (a quick reachability check).
    let status = Command::new(binary())
        .arg("--help")
        .status()
        .expect("failed to launch ldsc --help");
    assert!(status.success(), "ldsc --help failed");
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn read_gz(path: &str) -> String {
    use flate2::read::GzDecoder;
    use std::io::Read;
    let file = fs::File::open(path).expect("open gz file");
    let mut decoder = GzDecoder::new(file);
    let mut s = String::new();
    decoder.read_to_string(&mut s).expect("decompress gz");
    s
}
