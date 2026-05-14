//! Plotters charts rendered into a `<canvas>` via the
//! `plotters-canvas` `CanvasBackend`.
//!
//! Two charts: the LD-score histogram and the L2-vs-MAF scatter.
//! Both take pre-computed Vec<f64> buffers and render synchronously.

use anyhow::{Context, Result};
use plotters::prelude::*;
use plotters_canvas::CanvasBackend;

/// Render an LD-score histogram into the canvas with id `canvas_id`.
/// 50 evenly-spaced bins from 0 to the 99th percentile (clipping the
/// long tail so the bulk of the distribution is readable).
pub fn render_l2_histogram(canvas_id: &str, l2: &[f64]) -> Result<()> {
    let backend =
        CanvasBackend::new(canvas_id).context("CanvasBackend::new — canvas element not found")?;
    let root = backend.into_drawing_area();
    root.fill(&WHITE).map_err(to_anyhow)?;

    if l2.is_empty() {
        return Ok(());
    }

    let p99 = percentile(l2, 0.99).unwrap_or(0.0).max(1.0);
    let bins = 50usize;
    let mut counts = vec![0u32; bins];
    for &v in l2 {
        if v < 0.0 || !v.is_finite() {
            continue;
        }
        let bin = ((v / p99) * bins as f64).floor() as usize;
        let bin = bin.min(bins - 1);
        counts[bin] += 1;
    }
    let max_count = counts.iter().copied().max().unwrap_or(1).max(1);

    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .caption(
            format!("LD-score histogram (50 bins, 0 → P99 = {:.1})", p99),
            ("Arial", 16, &BLACK),
        )
        .build_cartesian_2d(0f64..p99, 0u32..(max_count + max_count / 8 + 1))
        .map_err(to_anyhow)?;

    chart
        .configure_mesh()
        .x_desc("LD score (ℓ_j)")
        .y_desc("count")
        .axis_desc_style(("Arial", 13, &BLACK))
        .label_style(("Arial", 11, &BLACK))
        .draw()
        .map_err(to_anyhow)?;

    let bin_w = p99 / bins as f64;
    let blue = RGBColor(42, 113, 165); // matches LDLink primary
    chart
        .draw_series(counts.iter().enumerate().map(|(i, &c)| {
            let x0 = i as f64 * bin_w;
            let x1 = x0 + bin_w;
            Rectangle::new([(x0, 0u32), (x1, c)], blue.filled())
        }))
        .map_err(to_anyhow)?;

    root.present().map_err(to_anyhow)?;
    Ok(())
}

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
