#!/usr/bin/env python3
"""Plot SKYPE benchmark error and high-copy-number nClose counts together.

The bars show the arithmetic mean across cell lines.  The dots show the
individual cell-line values used to calculate each mean.  The two metrics are
saved as panels in one combined figure.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, MaxNLocator


CONDITION_ORDER = [
    "skype",
    "ont_severus",
    "ont_savana",
    "ont_nanomonsv",
    "ont_sniffles2",
    "hifi_severus",
    "hifi_savana",
    "hifi_nanomonsv",
    "hifi_sniffles2",
    "illumina_svaba",
    "illumina_gripss",
    "illumina_manta",
]

CONDITION_INFO = {
    "skype": ("SKYPE", "SKYPE"),
    "ont_severus": ("ONT", "Severus"),
    "ont_savana": ("ONT", "SAVANA"),
    "ont_nanomonsv": ("ONT", "Nanomonsv"),
    "ont_sniffles2": ("ONT", "Sniffles2"),
    "hifi_severus": ("PacBio HiFi", "Severus"),
    "hifi_savana": ("PacBio HiFi", "SAVANA"),
    "hifi_nanomonsv": ("PacBio HiFi", "Nanomonsv"),
    "hifi_sniffles2": ("PacBio HiFi", "Sniffles2"),
    "illumina_svaba": ("Illumina", "SvABA"),
    "illumina_gripss": ("Illumina", "GRIPSS"),
    "illumina_manta": ("Illumina", "Manta"),
}

TECHNOLOGY_ORDER = ["SKYPE", "ONT", "PacBio HiFi", "Illumina"]
TECHNOLOGY_COLORS = {
    "SKYPE": "#178F83",
    "ONT": "#F6A51C",
    "PacBio HiFi": "#C4457A",
    "Illumina": "#30446F",
}
POINT_COLOR = "#67A2CC"

REQUIRED_COLUMNS = {
    "cell_line",
    "karyotype_type",
    "denoised_relative_error",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path("skype_bench_results/skype_bench.csv"),
        help="Benchmark CSV (default: %(default)s)",
    )
    parser.add_argument(
        "--artifacts-dir",
        type=Path,
        default=Path("skype_bench_results/artifacts"),
        help="Directory containing CELL_LINE/CONDITION/nclose_report.tsv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("skype_bench_results"),
        help="Output directory (default: %(default)s)",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.1,
        help="Count nClose rows with nclose_cn strictly greater than this value",
    )
    return parser.parse_args()


def load_plot_data(
    csv_path: Path, artifacts_dir: Path, threshold: float
) -> pd.DataFrame:
    data = pd.read_csv(csv_path)
    missing_columns = REQUIRED_COLUMNS - set(data.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"{csv_path} is missing required columns: {missing}")

    duplicates = data.duplicated(["cell_line", "karyotype_type"], keep=False)
    if duplicates.any():
        duplicate_rows = data.loc[
            duplicates, ["cell_line", "karyotype_type"]
        ].to_dict("records")
        raise ValueError(f"Duplicate benchmark rows found: {duplicate_rows}")

    unknown_conditions = sorted(set(data["karyotype_type"]) - set(CONDITION_INFO))
    if unknown_conditions:
        unknown = ", ".join(unknown_conditions)
        raise ValueError(f"Unknown karyotype_type value(s): {unknown}")

    missing_conditions = sorted(set(CONDITION_INFO) - set(data["karyotype_type"]))
    if missing_conditions:
        missing = ", ".join(missing_conditions)
        raise ValueError(f"Benchmark CSV is missing condition(s): {missing}")

    counts: list[int] = []
    report_rows: list[int] = []
    for row in data.itertuples(index=False):
        report_path = (
            artifacts_dir
            / str(row.cell_line)
            / str(row.karyotype_type)
            / "nclose_report.tsv"
        )
        if not report_path.is_file():
            raise FileNotFoundError(f"Missing nClose report: {report_path}")

        report = pd.read_csv(report_path, sep="\t")
        if "nclose_cn" not in report.columns:
            raise ValueError(f"{report_path} has no nclose_cn column")

        nclose_cn = pd.to_numeric(report["nclose_cn"], errors="coerce")
        invalid = report["nclose_cn"].notna() & nclose_cn.isna()
        if invalid.any():
            examples = report.loc[invalid, "nclose_cn"].head(5).tolist()
            raise ValueError(
                f"{report_path} contains non-numeric nclose_cn values: {examples}"
            )

        counts.append(int(nclose_cn.gt(threshold).sum()))
        report_rows.append(len(report))

    data = data.copy()
    data["nclose_report_rows"] = report_rows
    data["nclose_cn_gt_threshold_count"] = counts
    data["technology"] = data["karyotype_type"].map(
        lambda value: CONDITION_INFO[value][0]
    )
    data["caller"] = data["karyotype_type"].map(
        lambda value: CONDITION_INFO[value][1]
    )
    return data


def summarize(data: pd.DataFrame) -> pd.DataFrame:
    summary = (
        data.groupby("karyotype_type", as_index=False)
        .agg(
            denoised_relative_error_mean=("denoised_relative_error", "mean"),
            nclose_cn_gt_threshold_count_mean=(
                "nclose_cn_gt_threshold_count",
                "mean",
            ),
            cell_line_count=("cell_line", "nunique"),
        )
        .set_index("karyotype_type")
        .loc[CONDITION_ORDER]
        .reset_index()
    )
    summary["technology"] = summary["karyotype_type"].map(
        lambda value: CONDITION_INFO[value][0]
    )
    summary["caller"] = summary["karyotype_type"].map(
        lambda value: CONDITION_INFO[value][1]
    )
    summary["denoised_relative_error_mean"] = summary[
        "denoised_relative_error_mean"
    ].round(5)
    summary["nclose_cn_gt_threshold_count_mean"] = summary[
        "nclose_cn_gt_threshold_count_mean"
    ].round(1)
    return summary[
        [
            "karyotype_type",
            "technology",
            "caller",
            "cell_line_count",
            "denoised_relative_error_mean",
            "nclose_cn_gt_threshold_count_mean",
        ]
    ]


def text_color(background: str) -> str:
    red, green, blue = to_rgb(background)
    luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
    return "#1A1A1A" if luminance > 0.57 else "white"


def add_group_separators(axis: plt.Axes, summary: pd.DataFrame) -> None:
    technologies = summary["technology"].tolist()
    for index in range(len(technologies) - 1):
        if technologies[index] != technologies[index + 1]:
            axis.axhline(
                index + 0.5,
                color="white",
                linewidth=3.0,
                zorder=2.1,
            )


def draw_panel(
    axis: plt.Axes,
    data: pd.DataFrame,
    summary: pd.DataFrame,
    value_column: str,
    mean_column: str,
    title: str,
    xlabel: str,
    value_format,
) -> None:
    y_positions = np.arange(len(CONDITION_ORDER))
    means = summary[mean_column].to_numpy(dtype=float)
    bar_colors = [
        TECHNOLOGY_COLORS[technology] for technology in summary["technology"]
    ]

    axis.barh(
        y_positions,
        means,
        height=0.88,
        color=bar_colors,
        edgecolor="white",
        linewidth=0.8,
        zorder=2,
    )

    cell_lines = sorted(data["cell_line"].unique())
    offsets = np.linspace(-0.18, 0.18, len(cell_lines))
    offset_by_cell_line = dict(zip(cell_lines, offsets))
    for row_index, condition in enumerate(CONDITION_ORDER):
        condition_data = data[data["karyotype_type"] == condition]
        point_y = [
            row_index + offset_by_cell_line[cell_line]
            for cell_line in condition_data["cell_line"]
        ]
        axis.scatter(
            condition_data[value_column],
            point_y,
            s=31,
            color=POINT_COLOR,
            edgecolor="white",
            linewidth=0.45,
            alpha=0.95,
            zorder=3,
        )

    axis.set_yticks(y_positions, summary["caller"])
    axis.invert_yaxis()
    axis.set_title(title, loc="left", fontsize=15, pad=10)
    axis.set_xlabel(xlabel, fontsize=11)
    axis.tick_params(axis="both", labelsize=10)
    axis.xaxis.grid(True, color="#D9D9D9", linewidth=0.8, alpha=0.75)
    axis.set_axisbelow(True)
    axis.margins(y=0.015)
    add_group_separators(axis, summary)

    observed_max = float(data[value_column].max())
    axis.set_xlim(0, observed_max * 1.08)

    label_x = observed_max * 0.012
    for y_position, mean, color in zip(y_positions, means, bar_colors):
        axis.text(
            label_x,
            y_position,
            value_format(mean),
            va="center",
            ha="left",
            fontsize=10,
            color=text_color(color),
            zorder=4,
        )

    for spine in axis.spines.values():
        spine.set_color("#444444")
        spine.set_linewidth(0.8)


def technology_legend() -> list[Patch]:
    return [
        Patch(
            facecolor=TECHNOLOGY_COLORS[technology],
            edgecolor="none",
            label=technology,
        )
        for technology in TECHNOLOGY_ORDER
    ]


def add_figure_annotations(
    figure: plt.Figure, threshold: float, legend_columns: int = 4
) -> None:
    figure.legend(
        handles=technology_legend(),
        title="Input / technology",
        loc="upper center",
        bbox_to_anchor=(0.5, 0.995),
        ncol=legend_columns,
        frameon=False,
        fontsize=10,
        title_fontsize=11,
        handlelength=1.1,
        columnspacing=1.8,
    )
    figure.text(
        0.5,
        0.012,
        (
            "Bars show arithmetic means; dots show individual cell lines "
            f"(nClose threshold: nclose_cn > {threshold:g})."
        ),
        ha="center",
        va="bottom",
        fontsize=9,
        color="#555555",
    )


def save_all_formats(figure: plt.Figure, output_stem: Path) -> None:
    figure.savefig(
        output_stem.with_suffix(".png"),
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )
    figure.savefig(
        output_stem.with_suffix(".svg"),
        bbox_inches="tight",
        facecolor="white",
    )


def make_combined_plot(
    data: pd.DataFrame,
    summary: pd.DataFrame,
    output_dir: Path,
    threshold: float,
) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(14.5, 8.2))
    draw_panel(
        axes[0],
        data,
        summary,
        value_column="denoised_relative_error",
        mean_column="denoised_relative_error_mean",
        title="Mean denoised relative error",
        xlabel="Denoised relative error (lower is better)",
        value_format=lambda value: f"{value:.4f}",
    )
    axes[0].xaxis.set_major_locator(MaxNLocator(nbins=5))
    axes[0].xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:.2f}"))

    draw_panel(
        axes[1],
        data,
        summary,
        value_column="nclose_cn_gt_threshold_count",
        mean_column="nclose_cn_gt_threshold_count_mean",
        title=f"Mean nClose count (nclose_cn > {threshold:g})",
        xlabel=f"Number of nClose records with nclose_cn > {threshold:g}",
        value_format=lambda value: f"{value:,.1f}",
    )
    axes[1].xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
    axes[1].xaxis.set_major_formatter(
        FuncFormatter(lambda value, _: f"{value:,.0f}")
    )

    add_figure_annotations(figure, threshold)
    figure.subplots_adjust(left=0.105, right=0.985, top=0.89, bottom=0.105, wspace=0.34)
    save_all_formats(figure, output_dir / "skype_bench_metrics")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    data = load_plot_data(args.csv, args.artifacts_dir, args.threshold)
    summary = summarize(data)

    data.to_csv(args.output_dir / "skype_bench_plot_data.csv", index=False)
    summary.to_csv(args.output_dir / "skype_bench_means.csv", index=False)

    make_combined_plot(data, summary, args.output_dir, args.threshold)


if __name__ == "__main__":
    main()
