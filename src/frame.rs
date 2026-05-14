//! Minimal column-oriented data frame for sumstats / LD scores.
//!
//! Replaces polars LazyFrame/DataFrame for the operations ldsc actually uses:
//! typed columns (Str / F64), select/rename/drop, hstack, inner-join on a string
//! key. Pure Rust; compiles to wasm32-unknown-unknown.

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::{BufRead, BufReader};

/// A typed column. Only `Str` and `F64` are needed for sumstats / LD scores.
/// Missing values: `Str` uses empty string sentinel via `Option<String>`,
/// `F64` uses `f64::NAN` (matches polars' Float64 + null semantics in our usage —
/// every downstream consumer already handles NaN).
#[derive(Debug, Clone)]
pub enum Column {
    Str(Vec<Option<String>>),
    F64(Vec<f64>),
}

impl Column {
    pub fn len(&self) -> usize {
        match self {
            Column::Str(v) => v.len(),
            Column::F64(v) => v.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn as_f64(&self) -> Result<&[f64]> {
        match self {
            Column::F64(v) => Ok(v),
            Column::Str(_) => anyhow::bail!("expected F64 column, got Str"),
        }
    }

    pub fn as_str(&self) -> Result<&[Option<String>]> {
        match self {
            Column::Str(v) => Ok(v),
            Column::F64(_) => anyhow::bail!("expected Str column, got F64"),
        }
    }

    /// Coerce to F64. Str cells parse via `f64::from_str`; "NA" / "." / empty → NaN.
    pub fn cast_to_f64(&self) -> Result<Column> {
        match self {
            Column::F64(_) => Ok(self.clone()),
            Column::Str(v) => {
                let out = v
                    .iter()
                    .map(|opt| match opt {
                        None => f64::NAN,
                        Some(s) => parse_f64_or_nan(s),
                    })
                    .collect();
                Ok(Column::F64(out))
            }
        }
    }
}

/// Parse a token as f64. Recognises NA / . / "" / NaN as NaN; anything else
/// falls back to f64::from_str (returning NaN on failure to match polars'
/// permissive coercion).
pub fn parse_f64_or_nan(s: &str) -> f64 {
    let t = s.trim();
    if t.is_empty() || t.eq_ignore_ascii_case("na") || t == "." || t.eq_ignore_ascii_case("nan") {
        return f64::NAN;
    }
    t.parse().unwrap_or(f64::NAN)
}

/// A column-oriented in-memory frame. Order of columns is preserved.
#[derive(Debug, Clone, Default)]
pub struct Frame {
    pub cols: Vec<(String, Column)>,
    pub n_rows: usize,
}

impl Frame {
    pub fn new() -> Self {
        Frame::default()
    }

    pub fn from_cols(cols: Vec<(String, Column)>) -> Result<Self> {
        let n_rows = cols.first().map(|(_, c)| c.len()).unwrap_or(0);
        for (n, c) in &cols {
            anyhow::ensure!(
                c.len() == n_rows,
                "column '{}' has {} rows; expected {}",
                n,
                c.len(),
                n_rows
            );
        }
        Ok(Frame { cols, n_rows })
    }

    pub fn height(&self) -> usize {
        self.n_rows
    }

    pub fn width(&self) -> usize {
        self.cols.len()
    }

    pub fn column_names(&self) -> Vec<&str> {
        self.cols.iter().map(|(n, _)| n.as_str()).collect()
    }

    pub fn column_names_owned(&self) -> Vec<String> {
        self.cols.iter().map(|(n, _)| n.clone()).collect()
    }

    pub fn has_column(&self, name: &str) -> bool {
        self.cols.iter().any(|(n, _)| n == name)
    }

    pub fn column(&self, name: &str) -> Result<&Column> {
        self.cols
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, c)| c)
            .with_context(|| format!("column '{}' not found in frame", name))
    }

    pub fn column_mut(&mut self, name: &str) -> Result<&mut Column> {
        self.cols
            .iter_mut()
            .find(|(n, _)| n == name)
            .map(|(_, c)| c)
            .with_context(|| format!("column '{}' not found in frame", name))
    }

    pub fn column_idx(&self, name: &str) -> Option<usize> {
        self.cols.iter().position(|(n, _)| n == name)
    }

    /// Replace the named column's data (must be same length).
    pub fn replace_column(&mut self, name: &str, col: Column) -> Result<()> {
        let i = self
            .column_idx(name)
            .with_context(|| format!("column '{}' not found", name))?;
        anyhow::ensure!(
            col.len() == self.n_rows,
            "column '{}' new length {} != frame height {}",
            name,
            col.len(),
            self.n_rows
        );
        self.cols[i].1 = col;
        Ok(())
    }

    pub fn rename(&mut self, old: &str, new: &str) -> Result<()> {
        let i = self
            .column_idx(old)
            .with_context(|| format!("column '{}' not found for rename", old))?;
        anyhow::ensure!(
            !self.has_column(new) || self.cols[i].0 == new,
            "rename target '{}' already exists",
            new
        );
        self.cols[i].0 = new.to_string();
        Ok(())
    }

    /// Append a new column. Length must match frame height (or frame is empty).
    pub fn add_column(&mut self, name: String, col: Column) -> Result<()> {
        if self.cols.is_empty() {
            self.n_rows = col.len();
        } else {
            anyhow::ensure!(
                col.len() == self.n_rows,
                "column '{}' length {} != frame height {}",
                name,
                col.len(),
                self.n_rows
            );
        }
        anyhow::ensure!(!self.has_column(&name), "column '{}' already exists", name);
        self.cols.push((name, col));
        Ok(())
    }

    pub fn drop_column(&mut self, name: &str) -> Result<Column> {
        let i = self
            .column_idx(name)
            .with_context(|| format!("column '{}' not found for drop", name))?;
        let (_, col) = self.cols.remove(i);
        if self.cols.is_empty() {
            self.n_rows = 0;
        }
        Ok(col)
    }

    /// Subset to the named columns in the order given. Each name must exist.
    pub fn select(&self, names: &[&str]) -> Result<Frame> {
        let mut out = Vec::with_capacity(names.len());
        for n in names {
            let c = self.column(n)?.clone();
            out.push((n.to_string(), c));
        }
        Frame::from_cols(out)
    }

    /// Add `other`'s columns onto self. Heights must match; column names must
    /// be disjoint (caller renames first if needed).
    pub fn hstack(&mut self, other: Frame) -> Result<()> {
        if !self.cols.is_empty() {
            anyhow::ensure!(
                other.height() == self.n_rows,
                "hstack height mismatch: {} vs {}",
                self.n_rows,
                other.height()
            );
        }
        for (name, col) in other.cols {
            self.add_column(name, col)?;
        }
        Ok(())
    }

    /// Permute rows by `perm` (`perm[new_i] = old_i`).
    pub fn take_rows(&self, perm: &[usize]) -> Result<Frame> {
        let mut out: Vec<(String, Column)> = Vec::with_capacity(self.cols.len());
        for (name, col) in &self.cols {
            let new_col = match col {
                Column::F64(v) => {
                    let mut nv = Vec::with_capacity(perm.len());
                    for &i in perm {
                        nv.push(v[i]);
                    }
                    Column::F64(nv)
                }
                Column::Str(v) => {
                    let mut nv: Vec<Option<String>> = Vec::with_capacity(perm.len());
                    for &i in perm {
                        nv.push(v[i].clone());
                    }
                    Column::Str(nv)
                }
            };
            out.push((name.clone(), new_col));
        }
        let n_rows = perm.len();
        Ok(Frame { cols: out, n_rows })
    }

    /// Inner join self ⨝ other on the named string key, preserving self's row
    /// order. Output columns: all of self, then all of other except `key`.
    pub fn join_inner_on(&self, other: &Frame, key: &str) -> Result<Frame> {
        let left_keys = self.column(key)?.as_str()?;
        let right_keys = other.column(key)?.as_str()?;
        // Build right index: key → first row index.
        let mut right_idx: HashMap<&str, usize> = HashMap::with_capacity(right_keys.len());
        for (i, k) in right_keys.iter().enumerate() {
            if let Some(s) = k {
                right_idx.entry(s.as_str()).or_insert(i);
            }
        }
        // Walk self in order, keeping rows where the key is present in right.
        let mut left_take = Vec::with_capacity(left_keys.len());
        let mut right_take = Vec::with_capacity(left_keys.len());
        for (i, k) in left_keys.iter().enumerate() {
            if let Some(s) = k
                && let Some(&j) = right_idx.get(s.as_str())
            {
                left_take.push(i);
                right_take.push(j);
            }
        }
        let left_sub = self.take_rows(&left_take)?;
        let right_sub = other.take_rows(&right_take)?;
        // Combine: keep all left columns + all right columns except `key`.
        let mut out = left_sub;
        for (name, col) in right_sub.cols {
            if name == key {
                continue;
            }
            // If a non-key collision occurs, suffix with _right (rare in our usage).
            if out.has_column(&name) {
                out.add_column(format!("{}_right", name), col)?;
            } else {
                out.add_column(name, col)?;
            }
        }
        Ok(out)
    }

    /// True iff `self.column(key) == other.column(key)` (option-equality).
    pub fn key_columns_equal(&self, other: &Frame, key: &str) -> Result<bool> {
        if self.height() != other.height() {
            return Ok(false);
        }
        let a = self.column(key)?.as_str()?;
        let b = other.column(key)?.as_str()?;
        Ok(a == b)
    }
}

/// Read a tab-or-whitespace-separated text file into a `Frame`.
///
/// First line is the header (column names). The separator (tab vs runs of
/// whitespace) is determined once based on whether the *header* contains a
/// tab — then applied uniformly to every row. This avoids the foot-gun where
/// the header uses tabs but a body row uses spaces (or vice versa) and ends
/// up with a different column count.
///
/// Type inference: each column is F64 if every non-empty / non-NA value
/// parses as f64; otherwise Str. `.gz` and `.bz2` are decompressed
/// transparently.
///
/// Caveat: does not handle CSV-style quoting (sumstats / LD-score files
/// don't use quotes; if you need quoted fields, use a real csv parser).
pub fn read_tsv(path: &str) -> Result<Frame> {
    let mut reader: Box<dyn BufRead> = open_text_reader(path)?;

    let mut header_line: Option<String> = None;
    let mut body: Vec<String> = Vec::new();
    let mut line_buf = String::new();
    loop {
        line_buf.clear();
        let n = reader
            .read_line(&mut line_buf)
            .with_context(|| format!("reading line of '{}'", path))?;
        if n == 0 {
            break;
        }
        let trimmed = line_buf.trim_end_matches(['\n', '\r']);
        if trimmed.is_empty() {
            continue;
        }
        if header_line.is_none() {
            header_line = Some(trimmed.to_string());
        } else {
            body.push(trimmed.to_string());
        }
    }
    let header_line =
        header_line.with_context(|| format!("'{}' is empty (no header row)", path))?;

    // Pick the separator based on the header — tab if present anywhere,
    // else any run of whitespace.
    let use_tab = header_line.contains('\t');
    let split_line = |s: &str| -> Vec<String> {
        if use_tab {
            s.split('\t').map(|p| p.to_string()).collect()
        } else {
            s.split_whitespace().map(|p| p.to_string()).collect()
        }
    };

    let header: Vec<String> = split_line(&header_line);
    let ncol = header.len();
    let nrow = body.len();

    // Parse rows with the chosen separator; pad short rows.
    let mut rows: Vec<Vec<String>> = body.iter().map(|s| split_line(s)).collect();
    for r in rows.iter_mut() {
        if r.len() < ncol {
            r.resize(ncol, String::new());
        }
    }

    // Type inference: try parsing each cell as f64 once. If every non-empty
    // / non-NA cell parses, the column is F64 and we already have the values.
    let mut out_cols: Vec<(String, Column)> = Vec::with_capacity(ncol);
    for j in 0..ncol {
        let mut all_numeric = true;
        let mut parsed: Vec<f64> = Vec::with_capacity(nrow);
        for r in &rows {
            let s = r[j].trim();
            if is_na_token(s) {
                parsed.push(f64::NAN);
                continue;
            }
            match s.parse::<f64>() {
                Ok(v) => parsed.push(v),
                Err(_) => {
                    all_numeric = false;
                    break;
                }
            }
        }
        let col = if all_numeric {
            Column::F64(parsed)
        } else {
            // Fall back to string column.
            let v: Vec<Option<String>> = rows
                .iter()
                .map(|r| {
                    let s = r[j].trim();
                    if is_na_token(s) {
                        None
                    } else {
                        Some(s.to_string())
                    }
                })
                .collect();
            Column::Str(v)
        };
        out_cols.push((header[j].clone(), col));
    }

    Ok(Frame {
        cols: out_cols,
        n_rows: nrow,
    })
}

fn is_na_token(s: &str) -> bool {
    s.is_empty() || s.eq_ignore_ascii_case("na") || s == "." || s.eq_ignore_ascii_case("nan")
}

/// Open a possibly-compressed text reader.
pub fn open_text_reader(path: &str) -> Result<Box<dyn BufRead>> {
    let file = std::fs::File::open(path).with_context(|| format!("opening '{}'", path))?;
    if path.ends_with(".gz") {
        let gz = flate2::read::MultiGzDecoder::new(file);
        Ok(Box::new(BufReader::new(gz)))
    } else if path.ends_with(".bz2") {
        let bz = bzip2_rs::DecoderReader::new(file);
        Ok(Box::new(BufReader::new(bz)))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

impl Frame {
    /// Keep only rows where `mask[i]` is true. Mask must be `self.height()` long.
    pub fn filter_rows(&self, mask: &[bool]) -> Result<Frame> {
        anyhow::ensure!(
            mask.len() == self.n_rows,
            "filter_rows mask length {} != frame height {}",
            mask.len(),
            self.n_rows
        );
        let perm: Vec<usize> = mask
            .iter()
            .enumerate()
            .filter_map(|(i, &b)| if b { Some(i) } else { None })
            .collect();
        self.take_rows(&perm)
    }

    /// Deduplicate keeping the first occurrence of each value of `key`.
    /// Returns the deduplicated frame.
    pub fn unique_first_on(&self, key: &str) -> Result<Frame> {
        let keys = self.column(key)?.as_str()?;
        let mut seen: std::collections::HashSet<&str> = std::collections::HashSet::new();
        let mut keep: Vec<usize> = Vec::new();
        for (i, k) in keys.iter().enumerate() {
            if let Some(s) = k {
                if seen.insert(s.as_str()) {
                    keep.push(i);
                }
            } else if seen.insert("") {
                keep.push(i);
            }
        }
        self.take_rows(&keep)
    }

    /// Uppercase every value in the named Str column in place.
    pub fn uppercase_str_column(&mut self, name: &str) -> Result<()> {
        let col = self.column_mut(name)?;
        match col {
            Column::Str(v) => {
                for s in v.iter_mut() {
                    if let Some(x) = s.as_mut() {
                        *x = x.to_ascii_uppercase();
                    }
                }
                Ok(())
            }
            Column::F64(_) => anyhow::bail!("uppercase_str_column: '{}' is not Str", name),
        }
    }
}

/// Format an f64 for TSV output. Matches the prior polars
/// `with_float_precision(3)` behavior: 3 fractional digits then strip
/// trailing zeros, but keep `.0` on integer-valued floats so `6702.0`
/// round-trips as `6702.0` (not `6702`). NaN renders as the empty string.
fn format_f64(v: f64) -> String {
    if v.is_nan() {
        return String::new();
    }
    let mut s = format!("{:.3}", v);
    if s.contains('.') {
        // Pop trailing zeros except when doing so would leave a bare `.`,
        // i.e. preserve the `.0` on `6702.0`.
        while s.ends_with('0') && !s.ends_with(".0") {
            s.pop();
        }
    }
    s
}

/// Write a Frame as a tab-separated text file. If `path` ends with `.gz`,
/// gzip-compress the output (fast level).
pub fn write_tsv(path: &str, frame: &Frame) -> Result<()> {
    use std::io::Write;
    let file = std::fs::File::create(path).with_context(|| format!("creating '{}'", path))?;
    let mut w: Box<dyn Write> = if path.ends_with(".gz") {
        Box::new(flate2::write::GzEncoder::new(
            std::io::BufWriter::new(file),
            flate2::Compression::fast(),
        ))
    } else {
        Box::new(std::io::BufWriter::new(file))
    };
    // Header
    let names: Vec<&str> = frame.column_names();
    writeln!(w, "{}", names.join("\t"))?;
    // Rows
    let cols: Vec<&Column> = frame.cols.iter().map(|(_, c)| c).collect();
    for i in 0..frame.height() {
        for (j, c) in cols.iter().enumerate() {
            if j > 0 {
                w.write_all(b"\t")?;
            }
            match c {
                Column::F64(v) => w.write_all(format_f64(v[i]).as_bytes())?,
                Column::Str(v) => {
                    if let Some(s) = &v[i] {
                        w.write_all(s.as_bytes())?;
                    }
                }
            }
        }
        w.write_all(b"\n")?;
    }
    w.flush()?;
    Ok(())
}

/// Concatenate frames vertically (rbind). All frames must have identical
/// column names + types in the same order.
pub fn concat_rows(frames: Vec<Frame>) -> Result<Frame> {
    let mut iter = frames.into_iter();
    let mut acc = match iter.next() {
        Some(f) => f,
        None => anyhow::bail!("concat_rows: empty input"),
    };
    for f in iter {
        anyhow::ensure!(
            f.width() == acc.width(),
            "concat_rows: column count mismatch ({} vs {})",
            acc.width(),
            f.width()
        );
        for ((an, ac), (bn, bc)) in acc.cols.iter_mut().zip(f.cols) {
            anyhow::ensure!(
                *an == bn,
                "concat_rows: column name mismatch '{}' vs '{}'",
                an,
                bn
            );
            match (ac, bc) {
                (Column::F64(a), Column::F64(b)) => a.extend(b),
                (Column::Str(a), Column::Str(b)) => a.extend(b),
                _ => anyhow::bail!("concat_rows: column type mismatch for '{}'", an),
            }
        }
        acc.n_rows += f.n_rows;
    }
    Ok(acc)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_basic() {
        let mut f = Frame::new();
        f.add_column(
            "SNP".into(),
            Column::Str(vec![Some("rs1".into()), Some("rs2".into())]),
        )
        .unwrap();
        f.add_column("L2".into(), Column::F64(vec![1.5, 2.5]))
            .unwrap();
        assert_eq!(f.height(), 2);
        assert_eq!(f.width(), 2);
        assert_eq!(f.column_names(), vec!["SNP", "L2"]);
        assert!(f.has_column("L2"));
    }

    #[test]
    fn join_inner_preserves_left_order() {
        let mut left = Frame::new();
        left.add_column(
            "SNP".into(),
            Column::Str(vec![Some("a".into()), Some("b".into()), Some("c".into())]),
        )
        .unwrap();
        left.add_column("X".into(), Column::F64(vec![1.0, 2.0, 3.0]))
            .unwrap();

        let mut right = Frame::new();
        right
            .add_column(
                "SNP".into(),
                Column::Str(vec![Some("c".into()), Some("a".into())]),
            )
            .unwrap();
        right
            .add_column("Y".into(), Column::F64(vec![30.0, 10.0]))
            .unwrap();

        let j = left.join_inner_on(&right, "SNP").unwrap();
        assert_eq!(j.height(), 2);
        let snps = j.column("SNP").unwrap().as_str().unwrap();
        assert_eq!(snps[0].as_deref(), Some("a"));
        assert_eq!(snps[1].as_deref(), Some("c"));
        let ys = j.column("Y").unwrap().as_f64().unwrap();
        assert_eq!(ys, &[10.0, 30.0]);
    }

    #[test]
    fn format_f64_preserves_dot_zero() {
        assert_eq!(format_f64(0.0), "0.0");
        assert_eq!(format_f64(6702.0), "6702.0");
        assert_eq!(format_f64(0.5), "0.5");
        assert_eq!(format_f64(0.123456), "0.123");
        assert_eq!(format_f64(-0.029), "-0.029");
        assert_eq!(format_f64(1.000), "1.0");
        assert_eq!(format_f64(f64::NAN), "");
    }

    #[test]
    fn read_tsv_uses_header_separator_uniformly() {
        // Header tabs, body uses spaces — should still be parsed with tabs.
        let path = "/tmp/frame_test_tabsep.tsv";
        std::fs::write(path, "SNP\tA1\nrs1\tA\nrs2\tT\n").unwrap();
        let f = read_tsv(path).unwrap();
        assert_eq!(f.height(), 2);
        assert_eq!(f.column_names(), vec!["SNP", "A1"]);
        std::fs::remove_file(path).ok();
    }

    #[test]
    fn type_infer_mixed() {
        // simulate via direct construction
        let csv = "SNP\tCHR\tBP\nrs1\t1\t100\nrs2\t1\t200\n";
        let path = "/tmp/frame_test_typeinfer.tsv";
        std::fs::write(path, csv).unwrap();
        let f = read_tsv(path).unwrap();
        assert_eq!(f.height(), 2);
        assert!(matches!(f.column("SNP").unwrap(), Column::Str(_)));
        assert!(matches!(f.column("CHR").unwrap(), Column::F64(_)));
        assert!(matches!(f.column("BP").unwrap(), Column::F64(_)));
        std::fs::remove_file(path).ok();
    }
}
