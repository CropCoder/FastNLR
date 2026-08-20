//! nlr-plot — statistics charts (plotters, outputs PNG).
//!
//! Generates from `RunStats`: a motif hit count bar chart and a per-chromosome NLR count bar chart.

use nlr_report::RunStats;
use plotters::prelude::*;
use std::path::Path;

/// Plot the motif hit count bar chart.
pub fn plot_motif_counts(path: &Path, stats: &RunStats) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(path, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    let max_count = stats.motif_counts.values().copied().max().unwrap_or(1).max(1) as i32;
    let mut chart = ChartBuilder::on(&root)
        .caption("Motif hit counts", ("sans-serif", 24))
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(1i32..21i32, 0i32..(max_count + 1))?;

    chart
        .configure_mesh()
        .x_desc("Motif")
        .y_desc("Count")
        .axis_desc_style(("sans-serif", 14))
        .draw()?;

    let data: Vec<(i32, i32)> = stats
        .motif_counts
        .iter()
        .map(|(id, c)| (*id as i32, *c as i32))
        .collect();

    chart.draw_series(
        Histogram::vertical(&chart)
            .style(BLUE.filled())
            .data(data.iter().map(|(id, c)| (*id, *c))),
    )?;

    root.present()?;
    Ok(())
}

/// Plot the per-chromosome NLR count bar chart.
pub fn plot_chromosome_nlrs(
    path: &Path,
    stats: &RunStats,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(path, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    let max_nlrs = stats
        .per_chromosome
        .values()
        .map(|(_, n, _)| *n as i32)
        .max()
        .unwrap_or(1)
        .max(1);

    let n = stats.per_chromosome.len().max(1) as i32;

    let mut chart = ChartBuilder::on(&root)
        .caption("NLR loci per chromosome", ("sans-serif", 24))
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(0i32..n, 0i32..(max_nlrs + 1))?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .x_labels(n as usize)
        .y_desc("NLR count")
        .axis_desc_style(("sans-serif", 14))
        .draw()?;

    let data: Vec<(i32, i32)> = stats
        .per_chromosome
        .values()
        .enumerate()
        .map(|(i, (_, nlr, _))| (i as i32, *nlr as i32))
        .collect();

    chart.draw_series(
        Histogram::vertical(&chart)
            .style(GREEN.filled())
            .data(data.iter().map(|(i, c)| (*i, *c))),
    )?;

    root.present()?;
    Ok(())
}
