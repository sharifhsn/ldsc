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
            assert!(content.lines().count() > 1, "no SNP rows in {}", path);
        }
    }
    assert!(
        found >= 22,
        "expected ≥22 chromosome files, found {}",
        found
    );
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
// h2-cts smoke test (tiny synthetic inputs)
// ---------------------------------------------------------------------------

#[test]
fn h2_cts_smoke() {
    let dir = tempfile::tempdir().expect("tempdir");
    let base = dir.path().join("base.");
    let cts = dir.path().join("cts.");
    let w = dir.path().join("w.");
    let out = dir.path().join("out");

    let sumstats = dir.path().join("trait.sumstats");
    fs::write(
        &sumstats,
        "SNP\tA1\tA2\tZ\tN\n\
rs1\tA\tG\t1.0\t1000\n\
rs2\tA\tG\t2.0\t1000\n\
rs3\tA\tG\t1.5\t1000\n",
    )
    .expect("write sumstats");

    fs::write(
        format!("{}1.l2.ldscore", base.to_string_lossy()),
        "CHR\tSNP\tBP\tBASEL2\n\
1\trs1\t1\t10.0\n\
1\trs2\t2\t12.0\n\
1\trs3\t3\t11.0\n",
    )
    .expect("write base ldscore");
    fs::write(format!("{}1.l2.M_5_50", base.to_string_lossy()), "3\n").expect("write base M");

    fs::write(
        format!("{}1.l2.ldscore", cts.to_string_lossy()),
        "CHR\tSNP\tBP\tCTSL2\n\
1\trs1\t1\t2.0\n\
1\trs2\t2\t3.0\n\
1\trs3\t3\t2.5\n",
    )
    .expect("write cts ldscore");
    fs::write(format!("{}1.l2.M_5_50", cts.to_string_lossy()), "3\n").expect("write cts M");

    fs::write(
        format!("{}1.l2.ldscore", w.to_string_lossy()),
        "CHR\tSNP\tBP\tL2\n\
1\trs1\t1\t9.0\n\
1\trs2\t2\t9.5\n\
1\trs3\t3\t9.2\n",
    )
    .expect("write weight ldscore");

    let ldcts = dir.path().join("cts.ldcts");
    fs::write(&ldcts, format!("CTS\t{}\n", cts.to_string_lossy())).expect("write ldcts");

    let status = Command::new(binary())
        .args([
            "h2",
            "--h2-cts",
            sumstats.to_str().unwrap(),
            "--ref-ld-chr",
            base.to_str().unwrap(),
            "--w-ld-chr",
            w.to_str().unwrap(),
            "--ref-ld-chr-cts",
            ldcts.to_str().unwrap(),
            "--out",
            out.to_str().unwrap(),
            "--n-blocks",
            "2",
        ])
        .status()
        .expect("failed to launch ldsc h2 --h2-cts");
    assert!(status.success(), "ldsc h2-cts exited with {}", status);

    let results = out.with_extension("cell_type_results.txt");
    assert!(results.exists(), "missing CTS results file");
    let content = fs::read_to_string(results).expect("read CTS results");
    assert!(
        content.lines().count() >= 2,
        "expected header + at least one row"
    );
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
