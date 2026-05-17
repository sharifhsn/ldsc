//! Plotters chart rendered into a `<canvas>` via the
//! `plotters-canvas` `CanvasBackend`. The L2-vs-MAF scatter is the
//! one plot LDSC researchers actually use for QC (positive r ≈ 0.27
//! expected in EUR panels per the LDSC wiki). The histogram was
//! dropped in workstream K — never published, percentiles in the
//! stat grid carry the same information.

use anyhow::{Context, Result};
use plotters::prelude::*;
use plotters_canvas::CanvasBackend;

/// Render an L2-vs-MAF scatter (one point per SNP).
pub fn render_l2_vs_maf(canvas_id: &str, maf: &[f64], l2: &[f64]) -> Result<()> {
    let backend =
        CanvasBackend::new(canvas_id).context("CanvasBackend::new — canvas element not found")?;
    let root = backend.into_drawing_area();
    root.fill(&WHITE).map_err(to_anyhow)?;

    if maf.is_empty() || l2.is_empty() {
        return Ok(());
    }

    let p99_l2 = percentile(l2, 0.99).unwrap_or(1.0).max(1.0);
    let blue = RGBColor(42, 113, 165).mix(0.35);

    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .caption("L2 vs MAF (one point per SNP)", ("Arial", 16, &BLACK))
        .build_cartesian_2d(0f64..0.5f64, 0f64..p99_l2)
        .map_err(to_anyhow)?;

    chart
        .configure_mesh()
        .x_desc("MAF")
        .y_desc("LD score (clipped at P99)")
        .axis_desc_style(("Arial", 13, &BLACK))
        .label_style(("Arial", 11, &BLACK))
        .draw()
        .map_err(to_anyhow)?;

    chart
        .draw_series(
            maf.iter()
                .zip(l2.iter())
                .filter(|(m, l)| m.is_finite() && l.is_finite())
                .map(|(&m, &l)| Circle::new((m, l.min(p99_l2)), 1.5, blue.filled())),
        )
        .map_err(to_anyhow)?;

    root.present().map_err(to_anyhow)?;
    Ok(())
}

/// q-th percentile via linear interpolation on a sorted copy. None
/// when the input is empty or all-NaN.
pub fn percentile(values: &[f64], q: f64) -> Option<f64> {
    let mut v: Vec<f64> = values.iter().copied().filter(|x| x.is_finite()).collect();
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = v.len();
    let pos = q * (n - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = (lo + 1).min(n - 1);
    let frac = pos - lo as f64;
    Some(v[lo] * (1.0 - frac) + v[hi] * frac)
}

/// plotters errors carry trait-object lifetimes that don't satisfy
/// anyhow's `Send + Sync + 'static` bound. Stringify them.
fn to_anyhow<E: std::fmt::Display>(e: E) -> anyhow::Error {
    anyhow::anyhow!("plotters: {}", e)
}
