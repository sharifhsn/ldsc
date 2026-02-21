/// Continuous-annotation binning into annot files (Python --cts-bin).
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use crate::cli::CtsAnnotArgs;
use crate::ldscore::parse_bim;
use crate::parse::resolve_text_path;

struct BinSpec {
    labels: Vec<String>,
    lower_bounds: Vec<f64>,
    cut_breaks: Vec<f64>,
}

pub fn run(args: CtsAnnotArgs) -> Result<()> {
    let snps =
        parse_bim(&args.bimfile).with_context(|| format!("reading BIM '{}'", args.bimfile))?;
    println!("Read {} SNPs from '{}'", snps.len(), args.bimfile);
    let snp_ids: Vec<String> = snps.iter().map(|s| s.snp.clone()).collect();

    let cts_files = split_paths(&args.cts_bin);
    anyhow::ensure!(
        !cts_files.is_empty(),
        "--cts-bin must list at least one file"
    );
    let breaks = parse_breaks(&args.cts_breaks, cts_files.len())?;
    let cts_names = parse_cts_names(args.cts_names.as_deref(), cts_files.len())?;

    let mut bin_specs: Vec<BinSpec> = Vec::new();
    let mut bin_indices: Vec<Vec<usize>> = Vec::new();

    for (i, path) in cts_files.iter().enumerate() {
        let values = read_cts_values(path, &snp_ids)?;
        let spec = compute_bins(&values, &breaks[i])?;
        let indices = assign_bins(&values, &spec.cut_breaks)?;
        bin_specs.push(spec);
        bin_indices.push(indices);
    }

    let sizes: Vec<usize> = bin_specs.iter().map(|b| b.labels.len()).collect();
    let mut combos = cartesian_indices(&sizes);
    combos.sort_by(|a, b| compare_combo(a, b, &bin_specs));

    let col_names = build_col_names(&combos, &bin_specs, &cts_names);
    let mut combo_map: HashMap<Vec<usize>, usize> = HashMap::new();
    for (idx, combo) in combos.iter().enumerate() {
        combo_map.insert(combo.clone(), idx);
    }

    let mut row_buf = vec![0u8; col_names.len()];
    let writer = open_writer(&args.annot_file)?;
    let mut w = BufWriter::new(writer);
    writeln!(w, "CHR\tBP\tSNP\tCM\t{}", col_names.join("\t"))?;

    for (row_idx, snp) in snps.iter().enumerate() {
        row_buf.fill(0);
        let combo: Vec<usize> = bin_indices.iter().map(|v| v[row_idx]).collect();
        let col_idx = combo_map.get(&combo).context("missing bin combination")?;
        row_buf[*col_idx] = 1;
        write!(w, "{}\t{}\t{}\t{}", snp.chr, snp.bp, snp.snp, snp.cm)?;
        for v in &row_buf {
            write!(w, "\t{}", v)?;
        }
        writeln!(w)?;
    }

    println!(
        "Wrote CTS annot matrix ({} columns) to '{}'",
        col_names.len(),
        args.annot_file
    );
    Ok(())
}

fn split_paths(raw: &str) -> Vec<String> {
    raw.split(',')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .collect()
}

fn parse_breaks(raw: &str, expected: usize) -> Result<Vec<Vec<f64>>> {
    let fixed = raw.replace('N', "-");
    let parts: Vec<&str> = fixed.split('x').collect();
    anyhow::ensure!(
        parts.len() == expected,
        "Need one set of breaks for each --cts-bin file"
    );
    let mut out = Vec::new();
    for part in parts {
        let vals: Vec<f64> = part
            .split(',')
            .filter(|s| !s.trim().is_empty())
            .map(|s| {
                s.parse::<f64>()
                    .with_context(|| format!("parsing break '{}'", s))
            })
            .collect::<Result<_>>()?;
        anyhow::ensure!(!vals.is_empty(), "Empty break list in --cts-breaks");
        out.push(vals);
    }
    Ok(out)
}

fn parse_cts_names(raw: Option<&str>, expected: usize) -> Result<Vec<String>> {
    if let Some(raw) = raw {
        let names: Vec<String> = raw.split(',').map(|s| s.trim().to_string()).collect();
        anyhow::ensure!(
            names.len() == expected,
            "Must provide one --cts-names entry per --cts-bin file"
        );
        Ok(names)
    } else {
        Ok((0..expected).map(|i| format!("ANNOT{}", i)).collect())
    }
}

fn read_cts_values(path: &str, snp_ids: &[String]) -> Result<Vec<f64>> {
    let resolved = resolve_text_path(path)?;
    let file = File::open(&resolved).with_context(|| format!("opening CTS file '{}'", path))?;
    let reader = BufReader::new(file);
    let mut values = Vec::with_capacity(snp_ids.len());
    for (i, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("reading CTS line {}", i + 1))?;
        let cols: Vec<&str> = line.split_whitespace().collect();
        if cols.is_empty() {
            continue;
        }
        anyhow::ensure!(
            cols.len() >= 2,
            "CTS line {} must have at least 2 columns",
            i + 1
        );
        let expected = snp_ids
            .get(values.len())
            .context("CTS file has extra rows")?;
        anyhow::ensure!(
            cols[0] == expected,
            "CTS SNP mismatch at row {}: expected {}, found {}",
            values.len() + 1,
            expected,
            cols[0]
        );
        let v: f64 = cols[1]
            .parse()
            .with_context(|| format!("parsing CTS value '{}' on line {}", cols[1], i + 1))?;
        values.push(v);
    }
    anyhow::ensure!(
        values.len() == snp_ids.len(),
        "CTS file '{}' has {} rows; expected {}",
        path,
        values.len(),
        snp_ids.len()
    );
    Ok(values)
}

fn compute_bins(values: &[f64], breaks: &[f64]) -> Result<BinSpec> {
    let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min = values.iter().cloned().fold(f64::INFINITY, f64::min);
    let mut cut_breaks = breaks.to_vec();
    let mut name_breaks = cut_breaks.clone();

    let all_ge_max = cut_breaks.iter().all(|&b| b >= max);
    let all_le_min = cut_breaks.iter().all(|&b| b <= min);
    if all_ge_max || all_le_min {
        anyhow::bail!("All breaks lie outside the range of the cts variable");
    }

    if cut_breaks.iter().all(|&b| b <= max) {
        name_breaks.push(max);
        cut_breaks.push(max + 1.0);
    }
    if cut_breaks.iter().all(|&b| b >= min) {
        name_breaks.push(min);
        cut_breaks.push(min - 1.0);
    }

    name_breaks.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    cut_breaks.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    if let Some(first) = name_breaks.first_mut() {
        *first = f64::NEG_INFINITY;
    }
    if let Some(last) = name_breaks.last_mut() {
        *last = f64::INFINITY;
    }

    let mut labels = Vec::new();
    let mut lower_bounds = Vec::new();
    for i in 0..cut_breaks.len() - 1 {
        let lb = name_breaks[i];
        let ub = name_breaks[i + 1];
        let lb_str = if lb.is_infinite() && lb.is_sign_negative() {
            "min".to_string()
        } else {
            lb.to_string()
        };
        let ub_str = if ub.is_infinite() && ub.is_sign_positive() {
            "max".to_string()
        } else {
            ub.to_string()
        };
        labels.push(format!("{}_{}", lb_str, ub_str));
        lower_bounds.push(lb);
    }

    Ok(BinSpec {
        labels,
        lower_bounds,
        cut_breaks,
    })
}

fn assign_bins(values: &[f64], cut_breaks: &[f64]) -> Result<Vec<usize>> {
    let mut out = Vec::with_capacity(values.len());
    for &v in values {
        let mut found = None;
        for i in 0..cut_breaks.len() - 1 {
            if v > cut_breaks[i] && v <= cut_breaks[i + 1] {
                found = Some(i);
                break;
            }
        }
        let idx = found.context("Some SNPs have no annotation in --cts-bin")?;
        out.push(idx);
    }
    Ok(out)
}

fn cartesian_indices(sizes: &[usize]) -> Vec<Vec<usize>> {
    let mut combos: Vec<Vec<usize>> = vec![Vec::new()];
    for &size in sizes {
        let mut next = Vec::new();
        for combo in &combos {
            for i in 0..size {
                let mut c = combo.clone();
                c.push(i);
                next.push(c);
            }
        }
        combos = next;
    }
    combos
}

fn compare_combo(a: &Vec<usize>, b: &Vec<usize>, specs: &[BinSpec]) -> std::cmp::Ordering {
    for (i, (ai, bi)) in a.iter().zip(b.iter()).enumerate() {
        let ka = specs[i].lower_bounds[*ai];
        let kb = specs[i].lower_bounds[*bi];
        match ka.partial_cmp(&kb).unwrap_or(std::cmp::Ordering::Equal) {
            std::cmp::Ordering::Equal => {}
            other => return other,
        }
    }
    std::cmp::Ordering::Equal
}

fn build_col_names(combos: &[Vec<usize>], specs: &[BinSpec], names: &[String]) -> Vec<String> {
    let mut out = Vec::with_capacity(combos.len());
    for combo in combos {
        if combo.len() == 1 {
            let label = &specs[0].labels[combo[0]];
            out.push(format!("{}_{}", names[0], label));
        } else {
            let mut parts = Vec::with_capacity(combo.len());
            for (i, idx) in combo.iter().enumerate() {
                parts.push(format!("{}_{}", names[i], specs[i].labels[*idx]));
            }
            out.push(parts.join("_"));
        }
    }
    out
}

fn open_writer(path: &str) -> Result<Box<dyn Write>> {
    if path.ends_with(".gz") {
        use flate2::Compression;
        use flate2::write::GzEncoder;
        let file = File::create(path).with_context(|| format!("creating '{}'", path))?;
        Ok(Box::new(GzEncoder::new(file, Compression::default())))
    } else {
        let file = File::create(path).with_context(|| format!("creating '{}'", path))?;
        Ok(Box::new(file))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_assign_bins_basic() {
        let values = vec![0.1, 0.2, 0.3];
        let spec = compute_bins(&values, &[0.15]).unwrap();
        let bins = assign_bins(&values, &spec.cut_breaks).unwrap();
        assert_eq!(bins, vec![0, 1, 1]);
        assert_eq!(
            spec.labels,
            vec!["min_0.15".to_string(), "0.15_max".to_string()]
        );
    }

    #[test]
    fn test_combo_sort_and_names() {
        let specs = vec![
            BinSpec {
                labels: vec!["min_0.0".into(), "0.0_max".into()],
                lower_bounds: vec![f64::NEG_INFINITY, 0.0],
                cut_breaks: vec![],
            },
            BinSpec {
                labels: vec!["min_1.0".into(), "1.0_max".into()],
                lower_bounds: vec![f64::NEG_INFINITY, 1.0],
                cut_breaks: vec![],
            },
        ];
        let mut combos = cartesian_indices(&[2, 2]);
        combos.sort_by(|a, b| compare_combo(a, b, &specs));
        let names = vec!["A".to_string(), "B".to_string()];
        let cols = build_col_names(&combos, &specs, &names);
        assert_eq!(
            cols,
            vec![
                "A_min_0.0_B_min_1.0",
                "A_min_0.0_B_1.0_max",
                "A_0.0_max_B_min_1.0",
                "A_0.0_max_B_1.0_max"
            ]
        );
    }
}
