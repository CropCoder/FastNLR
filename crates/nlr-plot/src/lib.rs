//! nlr-plot — statistics charts (plotters, outputs PNG).
//!
//! 从 `RunStats` 生成两张图：
//! - motif 命中计数柱状图（x 轴为 motif id，动态范围，支持扩展库 21-28）；
//! - 按染色体 NLR 计数柱状图（x 轴显示染色体名）。
//!
//! 字体说明：plotters 的 ab_glyph 后端默认不带任何字体，文字绘制会报
//! `FontUnavailable`。这里通过 `register_font` 注册内嵌的 Liberation Sans
//! （SIL Open Font License 1.1，与 GPL-3.0 兼容）。

use nlr_report::RunStats;
use plotters::coord::ranged1d::SegmentValue;
use plotters::prelude::*;
use std::path::Path;
use std::sync::Once;

/// 只注册一次内嵌字体。
static FONT_INIT: Once = Once::new();

fn ensure_font_registered() {
    FONT_INIT.call_once(|| {
        let _ = plotters::style::register_font(
            "sans-serif",
            plotters::style::FontStyle::Normal,
            include_bytes!("../assets/LiberationSans-Regular.ttf"),
        );
    });
}

/// motif 命中计数柱状图。
pub fn plot_motif_counts(path: &Path, stats: &RunStats) -> Result<(), Box<dyn std::error::Error>> {
    ensure_font_registered();
    let root = BitMapBackend::new(path, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    // 动态 motif id 范围（默认 1..=20，扩展库可达 28）。
    let ids: Vec<u8> = stats.motif_counts.keys().copied().collect();
    let min_id = ids.iter().copied().min().unwrap_or(1);
    let max_id = ids.iter().copied().max().unwrap_or(1);
    let x_hi = max_id as i32 + 1;
    let n_segments = (max_id - min_id + 1) as usize;
    let max_count = stats.motif_counts.values().copied().max().unwrap_or(1).max(1) as i32;

    let mut chart = ChartBuilder::on(&root)
        .caption("Motif hit counts", ("sans-serif", 24))
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d((min_id as i32..x_hi).into_segmented(), 0i32..(max_count + 1))?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .x_labels(n_segments)
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
            .margin(10)
            .data(data.iter().map(|(id, c)| (*id, *c))),
    )?;

    root.present()?;
    Ok(())
}

/// 按染色体 NLR 计数柱状图（x 轴显示染色体名）。
pub fn plot_chromosome_nlrs(
    path: &Path,
    stats: &RunStats,
) -> Result<(), Box<dyn std::error::Error>> {
    ensure_font_registered();
    let root = BitMapBackend::new(path, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    // BTreeMap 已按染色体名排序，keys() 顺序即输出顺序。
    let names: Vec<&String> = stats.per_chromosome.keys().collect();
    let n = names.len().max(1);
    let max_nlrs = stats
        .per_chromosome
        .values()
        .map(|(_, nlr, _)| *nlr as i32)
        .max()
        .unwrap_or(1)
        .max(1);

    let mut chart = ChartBuilder::on(&root)
        .caption("NLR loci per chromosome", ("sans-serif", 24))
        .margin(20)
        .x_label_area_size(70)
        .y_label_area_size(60)
        .build_cartesian_2d((0..n as i32).into_segmented(), 0i32..(max_nlrs + 1))?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .x_labels(n.min(60))
        .x_label_formatter(&|v: &SegmentValue<i32>| {
            let idx = match v {
                SegmentValue::Exact(i) | SegmentValue::CenterOf(i) => *i as usize,
                SegmentValue::Last => usize::MAX,
            };
            if idx < names.len() {
                let name = names[idx].as_str();
                if name.len() > 14 {
                    format!("{}...", &name[..14])
                } else {
                    name.to_string()
                }
            } else {
                String::new()
            }
        })
        .y_desc("NLR count")
        .axis_desc_style(("sans-serif", 14))
        .x_label_style(("sans-serif", 9))
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
            .margin(2)
            .data(data.iter().map(|(i, c)| (*i, *c))),
    )?;

    root.present()?;
    Ok(())
}

/// NLR 结构域类型计数柱状图（按数量降序）。
pub fn plot_nlr_type_counts(
    path: &Path,
    stats: &RunStats,
) -> Result<(), Box<dyn std::error::Error>> {
    ensure_font_registered();
    let root = BitMapBackend::new(path, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    // 按计数降序、类型名升序排序。
    let mut types: Vec<(&String, &usize)> = stats.nlr_type_counts.iter().collect();
    types.sort_by(|a, b| b.1.cmp(a.1).then_with(|| a.0.cmp(b.0)));
    let names: Vec<String> = types.iter().map(|(name, _)| name.to_string()).collect();
    let counts: Vec<i32> = types.iter().map(|(_, c)| **c as i32).collect();
    let n = names.len().max(1);
    let max_count = counts.iter().copied().max().unwrap_or(1).max(1);

    let mut chart = ChartBuilder::on(&root)
        .caption("NLR types", ("sans-serif", 24))
        .margin(20)
        .x_label_area_size(90)
        .y_label_area_size(60)
        .build_cartesian_2d((0..n as i32).into_segmented(), 0i32..(max_count + 1))?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .x_labels(n)
        .x_label_formatter(&|v: &SegmentValue<i32>| {
            let idx = match v {
                SegmentValue::Exact(i) | SegmentValue::CenterOf(i) => *i as usize,
                SegmentValue::Last => usize::MAX,
            };
            if idx < names.len() {
                names[idx].clone()
            } else {
                String::new()
            }
        })
        .y_desc("NLR count")
        .axis_desc_style(("sans-serif", 14))
        .x_label_style(("sans-serif", 9))
        .draw()?;

    chart.draw_series(
        Histogram::vertical(&chart)
            .style(BLUE.filled())
            .margin(10)
            .data((0..n as i32).zip(counts.iter().copied())),
    )?;

    root.present()?;
    Ok(())
}
