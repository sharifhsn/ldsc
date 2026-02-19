/// CLI argument definitions using clap derive macros.
use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[command(name = "ldsc", about = "LD Score Regression (Rust port)", version)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Subcommand)]
pub enum Command {
    /// Pre-process GWAS summary statistics (replaces munge_sumstats.py)
    MungeSumstats(MungeArgs),
    /// Compute LD scores from PLINK binary files (replaces --l2 in ldsc.py)
    Ldscore(LdscoreArgs),
    /// Estimate SNP heritability (replaces --h2 in ldsc.py)
    H2(H2Args),
    /// Estimate genetic correlation (replaces --rg in ldsc.py)
    Rg(RgArgs),
    /// Generate annotation files from BED regions or gene sets (replaces make_annot.py)
    MakeAnnot(MakeAnnotArgs),
}

// ---------------------------------------------------------------------------
// munge-sumstats
// ---------------------------------------------------------------------------

#[derive(Args)]
pub struct MungeArgs {
    /// Input summary statistics file (TSV/CSV, optionally gzipped)
    #[arg(long)]
    pub sumstats: String,

    /// Output file prefix
    #[arg(long)]
    pub out: String,

    /// Optional: merge alleles reference file
    #[arg(long)]
    pub merge_alleles: Option<String>,

    /// Minimum sample size to retain a SNP
    #[arg(long, default_value_t = 0)]
    pub n_min: u64,

    /// Minimum MAF to retain a SNP (0 to disable)
    #[arg(long, default_value_t = 0.01)]
    pub maf: f64,

    /// Minimum INFO score to retain a SNP (0 to disable)
    #[arg(long, default_value_t = 0.9)]
    pub info_min: f64,

    /// Remove strand-ambiguous SNPs (A/T, C/G)
    #[arg(long, default_value_t = true)]
    pub keep_mhc: bool,

    // --- Sample size overrides -----------------------------------------------

    /// Fix sample size for all SNPs (overrides any N column in the file)
    #[arg(long)]
    pub n: Option<f64>,

    /// Number of cases (N = N-cas + N-con; ignored if --n is set)
    #[arg(long)]
    pub n_cas: Option<f64>,

    /// Number of controls (N = N-cas + N-con; ignored if --n is set)
    #[arg(long)]
    pub n_con: Option<f64>,

    // --- Column name overrides (case-insensitive) ----------------------------
    // Use when the file uses a non-standard name not in the built-in synonym map.

    /// Name of the SNP ID column in the input file
    #[arg(long)]
    pub snp_col: Option<String>,

    /// Name of the sample size (N) column in the input file
    #[arg(long)]
    pub n_col: Option<String>,

    /// Name of the case sample size column (N = N-cas-col + N-con-col per row)
    #[arg(long)]
    pub n_cas_col: Option<String>,

    /// Name of the control sample size column (N = N-cas-col + N-con-col per row)
    #[arg(long)]
    pub n_con_col: Option<String>,

    /// Name of the effect allele (A1) column in the input file
    #[arg(long)]
    pub a1_col: Option<String>,

    /// Name of the other allele (A2) column in the input file
    #[arg(long)]
    pub a2_col: Option<String>,

    /// Name of the p-value column in the input file
    #[arg(long)]
    pub p_col: Option<String>,

    /// Name of the allele frequency column in the input file
    #[arg(long)]
    pub frq_col: Option<String>,

    /// Name of the imputation INFO column in the input file
    #[arg(long)]
    pub info_col: Option<String>,

    /// Signed summary statistic column and its null value: COLNAME,null
    /// (e.g., "Z,0" or "OR,1"). Used for P→Z conversion sign when the
    /// column name is not in the built-in map.
    #[arg(long)]
    pub signed_sumstats: Option<String>,

    /// Comma-separated list of column names to ignore (drop before any other processing).
    /// Useful when a column would otherwise be misidentified as a signed stat.
    #[arg(long)]
    pub ignore: Option<String>,

    /// Keep the allele frequency (MAF) column in the output file.
    /// By default only SNP, A1, A2, Z, N are written.
    #[arg(long, default_value_t = false)]
    pub keep_maf: bool,

    /// A1 is always the increasing allele — treat all Z-scores as positive when
    /// deriving Z from P-values. Eliminates the need for a signed summary stat.
    #[arg(long, default_value_t = false)]
    pub a1_inc: bool,

    /// Allow summary statistics files without allele columns (A1/A2).
    /// Skips strand-ambiguity filtering; output will not include A1/A2.
    #[arg(long, default_value_t = false)]
    pub no_alleles: bool,

    /// Comma-separated list of INFO column names (e.g. "INFO_EUR,INFO_EAS").
    /// When provided, the per-SNP INFO filter uses the MEAN of all listed columns.
    /// Takes precedence over the single --info-col / synonym-mapped INFO column.
    #[arg(long)]
    pub info_list: Option<String>,

    /// Name of the study-count column (number of studies a SNP was genotyped in).
    /// Used together with --nstudy-min to filter low-coverage SNPs.
    #[arg(long)]
    pub nstudy: Option<String>,

    /// Minimum number of studies a SNP must be genotyped in (requires --nstudy).
    #[arg(long)]
    pub nstudy_min: Option<u64>,
}

// ---------------------------------------------------------------------------
// ldscore (--l2)
// ---------------------------------------------------------------------------

#[derive(Args)]
pub struct LdscoreArgs {
    /// PLINK binary file prefix
    #[arg(long)]
    pub bfile: String,

    /// Output prefix
    #[arg(long)]
    pub out: String,

    /// LD window in cM (default 1.0; mutually exclusive with --ld-wind-kb / --ld-wind-snp)
    #[arg(long, default_value_t = 1.0)]
    pub ld_wind_cm: f64,

    /// LD window in kilobases (overrides --ld-wind-cm when set)
    #[arg(long)]
    pub ld_wind_kb: Option<f64>,

    /// LD window as number of flanking SNPs (overrides --ld-wind-cm when set)
    #[arg(long)]
    pub ld_wind_snp: Option<usize>,

    /// Annotation file prefix for partitioned LD scores.
    /// LDSC appends .annot or .annot.gz to this prefix.
    /// The annot file must have the same SNPs in the same order as the .bim file.
    #[arg(long)]
    pub annot: Option<String>,

    /// The annot file contains only annotation columns (no CHR/SNP/BP/CM columns).
    /// Use when working with "thin" annotation files.
    #[arg(long, default_value_t = false)]
    pub thin_annot: bool,

    /// File containing SNP IDs (one per line) to include in LD score computation.
    /// Only listed SNPs are used as targets and in windows; all others are excluded.
    #[arg(long)]
    pub extract: Option<String>,

    /// Exclude SNPs with minor allele frequency below FLOAT from LD computation.
    /// Applied after --extract; MAF is computed from the genotype data.
    #[arg(long)]
    pub maf: Option<f64>,

    /// File containing SNP IDs (one per line) to print LD scores for.
    /// Unlike --extract, all SNPs are still used in LD windows; only output is filtered.
    #[arg(long)]
    pub print_snps: Option<String>,

    /// File of individual IDs to restrict the reference panel to (one FID IID per line).
    /// Only individuals listed in this file are used for LD computation.
    #[arg(long)]
    pub keep: Option<String>,

    /// Compute per-allele LD scores: weight each r² by 2·p·(1−p) of the source SNP.
    /// L2[i] = Σ_j r²_unbiased(i,j) × 2·p_j·(1−p_j).
    #[arg(long, default_value_t = false)]
    pub per_allele: bool,

    /// Number of SNPs to process per BLAS chunk (default 50).
    /// Larger values change the window approximation slightly.
    #[arg(long, default_value_t = 50)]
    pub chunk_size: usize,

    /// Number of OpenBLAS threads for DGEMM calls (default 4).
    /// 4 is optimal for 1000G-scale data (n~2500); increase for biobank-scale data.
    #[arg(long, default_value_t = 4)]
    pub blas_threads: usize,
}

// ---------------------------------------------------------------------------
// h2 (heritability)
// ---------------------------------------------------------------------------

#[derive(Args)]
pub struct H2Args {
    /// Munged summary statistics file (.sumstats.gz)
    #[arg(long)]
    pub h2: String,

    /// LD score file prefix, per-chromosome (e.g. eas_ldscores/chr).
    /// LDSC appends .l2.ldscore.gz and the chromosome number.
    /// Provide exactly one of --ref-ld-chr or --ref-ld.
    #[arg(long)]
    pub ref_ld_chr: Option<String>,

    /// Single LD score file (non-chr-split alternative to --ref-ld-chr).
    /// When using this option, supply --m-snps explicitly (no .l2.M files to read).
    #[arg(long)]
    pub ref_ld: Option<String>,

    /// Regression weight LD score prefix (per-chromosome).
    /// Provide exactly one of --w-ld-chr or --w-ld.
    #[arg(long)]
    pub w_ld_chr: Option<String>,

    /// Single regression weight LD score file (non-chr-split alternative to --w-ld-chr).
    #[arg(long)]
    pub w_ld: Option<String>,

    /// Output prefix
    #[arg(long)]
    pub out: String,

    /// Number of SNPs in the LD score reference panel (e.g. 1190321 for HapMap3 EUR).
    /// Used as M in h2 = slope × M / N.  If omitted, reads from .l2.M_5_50 files
    /// alongside --ref-ld-chr; falls back to the number of regression SNPs if not found.
    #[arg(long)]
    pub m_snps: Option<f64>,

    /// Use .l2.M instead of .l2.M_5_50 as the M denominator.
    /// Python LDSC defaults to .l2.M_5_50; set this flag to replicate --not-M-5-50.
    #[arg(long, default_value_t = false)]
    pub not_m_5_50: bool,

    /// Sample prevalence (for liability-scale conversion; requires --pop-prev)
    #[arg(long)]
    pub samp_prev: Option<f64>,

    /// Population prevalence (for liability-scale conversion; requires --samp-prev)
    #[arg(long)]
    pub pop_prev: Option<f64>,

    /// Constrain the LD score regression intercept to 1 (for h2) rather than estimating it.
    /// Equivalent to Python's --no-intercept.
    #[arg(long, default_value_t = false)]
    pub no_intercept: bool,

    /// Fix the LD score regression intercept to a specific value (e.g. 1.02).
    /// Takes precedence over --no-intercept when both are given.
    #[arg(long)]
    pub intercept_h2: Option<f64>,

    /// Two-step estimator: use SNPs with chi² ≤ CHI2_MAX to estimate the intercept
    /// in step 1, then fix the intercept and estimate h2 on all SNPs in step 2.
    /// Recommended value: 30.
    #[arg(long)]
    pub two_step: Option<f64>,

    /// Remove SNPs with chi² > CHISQ_MAX before regression (default: no filter).
    #[arg(long)]
    pub chisq_max: Option<f64>,

    /// Number of jackknife blocks for SE estimation (default: 200)
    #[arg(long, default_value_t = 200)]
    pub n_blocks: usize,

    /// Print per-annotation τ coefficients and enrichment alongside total h2.
    /// Only meaningful when --ref-ld-chr points to partitioned (multi-column) LD score files.
    #[arg(long, default_value_t = false)]
    pub print_coefficients: bool,

    /// Print jackknife covariance matrix of regression estimates.
    #[arg(long, default_value_t = false)]
    pub print_cov: bool,

    /// Print per-block jackknife delete values for regression estimates.
    #[arg(long, default_value_t = false)]
    pub print_delete_vals: bool,

    /// Allow h2 < 0 or |rg| > 1 in output (Rust never clips; accepted for CLI parity).
    #[arg(long, default_value_t = false)]
    pub return_silly_things: bool,

    /// Force matrix inversion even when ill-conditioned (Rust uses least-squares solver;
    /// accepted for CLI parity).
    #[arg(long, default_value_t = false)]
    pub invert_anyway: bool,
}

// ---------------------------------------------------------------------------
// rg (genetic correlation)
// ---------------------------------------------------------------------------

#[derive(Args)]
pub struct RgArgs {
    /// Comma-separated list of .sumstats.gz files
    #[arg(long, value_delimiter = ',')]
    pub rg: Vec<String>,

    /// LD score file prefix (per-chromosome).
    /// Provide exactly one of --ref-ld-chr or --ref-ld.
    #[arg(long)]
    pub ref_ld_chr: Option<String>,

    /// Single LD score file (non-chr-split alternative to --ref-ld-chr).
    #[arg(long)]
    pub ref_ld: Option<String>,

    /// Regression weight LD score prefix (per-chromosome).
    /// Provide exactly one of --w-ld-chr or --w-ld.
    #[arg(long)]
    pub w_ld_chr: Option<String>,

    /// Single regression weight LD score file (non-chr-split alternative to --w-ld-chr).
    #[arg(long)]
    pub w_ld: Option<String>,

    /// Output prefix
    #[arg(long)]
    pub out: String,

    /// Number of SNPs in the LD score reference panel.
    /// If omitted, reads from .l2.M_5_50 files alongside --ref-ld-chr.
    #[arg(long)]
    pub m_snps: Option<f64>,

    /// Use .l2.M instead of .l2.M_5_50 as the M denominator.
    #[arg(long, default_value_t = false)]
    pub not_m_5_50: bool,

    /// Constrain the genetic covariance intercept to 0 rather than estimating it.
    #[arg(long, default_value_t = false)]
    pub no_intercept: bool,

    /// Two-step estimator for gencov: estimate intercept on SNPs with |Z1·Z2| ≤ CHI2_MAX,
    /// then fix it and re-estimate slope on all SNPs.
    #[arg(long)]
    pub two_step: Option<f64>,

    /// Remove SNPs with chi²₁·chi²₂ products above CHISQ_MAX before regression.
    #[arg(long)]
    pub chisq_max: Option<f64>,

    /// Number of jackknife blocks for SE estimation (default: 200)
    #[arg(long, default_value_t = 200)]
    pub n_blocks: usize,

    /// Sample prevalence for each trait (comma-separated), for liability-scale conversion.
    /// Requires --pop-prev.  Use one value per file in --rg.
    #[arg(long, value_delimiter = ',')]
    pub samp_prev: Vec<f64>,

    /// Population prevalence for each trait (comma-separated).
    /// Requires --samp-prev.  Use one value per file in --rg.
    #[arg(long, value_delimiter = ',')]
    pub pop_prev: Vec<f64>,

    /// Fix the gencov intercept to a specific value per trait-pair (comma-separated).
    /// One value per pair in order: for traits A,B,C the pairs are A-B, B-C.
    /// If fewer values than pairs are given, remaining pairs are freely estimated.
    #[arg(long, value_delimiter = ',')]
    pub intercept_gencov: Vec<f64>,

    /// Skip allele consistency checking between summary statistic files.
    /// Currently a no-op (allele checking not yet implemented); included for API parity.
    #[arg(long, default_value_t = false)]
    pub no_check_alleles: bool,

    /// Print jackknife covariance matrix of regression estimates.
    #[arg(long, default_value_t = false)]
    pub print_cov: bool,

    /// Print per-block jackknife delete values for regression estimates.
    #[arg(long, default_value_t = false)]
    pub print_delete_vals: bool,

    /// Allow |rg| > 1 or negative h2 in output (Rust never clips; accepted for CLI parity).
    #[arg(long, default_value_t = false)]
    pub return_silly_things: bool,

    /// Force matrix inversion even when ill-conditioned (accepted for CLI parity).
    #[arg(long, default_value_t = false)]
    pub invert_anyway: bool,
}

// ---------------------------------------------------------------------------
// make-annot
// ---------------------------------------------------------------------------

#[derive(Args)]
pub struct MakeAnnotArgs {
    /// PLINK .bim file listing the SNPs to annotate
    #[arg(long)]
    pub bimfile: String,

    /// Output annotation file path (e.g., prefix.annot or prefix.annot.gz).
    /// If the path ends in .gz, the output is gzip-compressed.
    #[arg(long)]
    pub annot_file: String,

    /// UCSC BED file defining annotation regions (mutually exclusive with --gene-set-file).
    /// Columns: chrom  start(0-based)  end(exclusive)  [name …]
    #[arg(long)]
    pub bed_file: Option<String>,

    /// Gene set file: one gene symbol per line.
    /// Requires --gene-coord-file.  SNPs within windowsize bp of each gene are annotated.
    #[arg(long)]
    pub gene_set_file: Option<String>,

    /// Gene coordinate file: tab-separated with columns GENE CHR START END.
    /// Requires --gene-set-file.
    #[arg(long)]
    pub gene_coord_file: Option<String>,

    /// Window size in base pairs to add around each region (symmetric extension).
    #[arg(long, default_value_t = 0)]
    pub windowsize: u32,

    /// Do not merge overlapping BED intervals before annotation.
    /// By default adjacent/overlapping intervals are merged.
    #[arg(long, default_value_t = false)]
    pub nomerge: bool,
}
