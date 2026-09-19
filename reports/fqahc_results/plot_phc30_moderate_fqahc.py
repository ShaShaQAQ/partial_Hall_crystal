#!/usr/bin/env python3
"""Plot and summarize the 30-site moderate-coupling FQAHC candidate."""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np


REPORT_DIR = Path(__file__).resolve().parent
DEFAULT_RESULT_DIR = (
    REPORT_DIR.parents[1] / "results" /
    "phc30_np12_v1_10_v2_2_v3_2_fqahc"
)
DEFAULT_OUTPUT_DIR = REPORT_DIR / "figures" / "optical_response"


def summarize_manifold(rows, manifold_size=15):
    """Summarize the lowest ``manifold_size`` levels in a spectrum table."""
    rows = np.asarray(rows, dtype=float)
    if rows.ndim != 2 or rows.shape[1] != 2:
        raise ValueError("spectrum must contain sector and energy columns")
    if len(rows) <= manifold_size:
        raise ValueError("spectrum has no level above the candidate manifold")

    order = np.argsort(rows[:, 1], kind="stable")
    selected = rows[order[:manifold_size]]
    ground_energy = selected[0, 1]
    top_energy = selected[-1, 1]
    next_energy = rows[order[manifold_size], 1]
    width = top_energy - ground_energy
    separation = next_energy - top_energy
    return {
        "sectors": tuple(sorted(selected[:, 0].astype(int).tolist())),
        "width": float(width),
        "next_energy": float(next_energy - ground_energy),
        "separation": float(separation),
        "isolated": bool(separation > width),
    }


def load_structure_factor(path):
    """Load per-sector and manifold-averaged structure-factor curves."""
    curves = {}
    average = None
    for line in Path(path).read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        if len(fields) != 5:
            raise ValueError("structure-factor rows must have five columns")
        label = fields[0]
        row = np.asarray([float(value) for value in fields[1:]], dtype=float)
        if label == "average":
            if average is None:
                average = []
            average.append(row)
        else:
            curves.setdefault(int(label), []).append(row)

    expected = set(range(15))
    if set(curves) != expected:
        raise ValueError("expected structure-factor curves for sectors 0 through 14")
    arrays = {sector: np.asarray(rows) for sector, rows in curves.items()}
    if any(array.shape != (15, 4) for array in arrays.values()):
        raise ValueError("each sector must contain 15 momentum samples")
    average_array = np.asarray(average if average is not None else [])
    if average_array.shape != (15, 4):
        raise ValueError("manifold average must contain 15 momentum samples")
    return arrays, average_array


def load_optical_curve(path):
    """Load one seven-column response-Lanczos conductivity table."""
    data = np.loadtxt(str(path), comments="#", dtype=float)
    data = np.atleast_2d(data)
    if data.ndim != 2 or data.shape[1] != 7:
        raise ValueError("optical curve must contain exactly seven columns")
    if data.shape[0] < 2:
        raise ValueError("optical curve must contain at least two frequencies")
    if not np.all(np.isfinite(data)):
        raise ValueError("optical curve contains non-finite values")
    if not np.all(np.diff(data[:, 0]) > 0):
        raise ValueError("frequencies must be strictly increasing")
    return data


def write_optical_manifold_average(optical_paths, output_path):
    """Write mean, standard deviation, minimum, and maximum over 15 sectors."""
    expected = set(range(15))
    if set(optical_paths) != expected:
        raise ValueError("optical average requires sectors 0 through 14")
    curves = [load_optical_curve(optical_paths[sector]) for sector in range(15)]
    frequencies = curves[0][:, 0]
    for sector, curve in enumerate(curves[1:], start=1):
        if not np.array_equal(curve[:, 0], frequencies):
            raise ValueError(
                "optical sector {} uses a different frequency mesh".format(sector)
            )
    components = np.stack([curve[:, 1:7] for curve in curves])
    statistics = np.column_stack(
        (
            frequencies,
            components.mean(axis=0),
            components.std(axis=0),
            components.min(axis=0),
            components.max(axis=0),
        )
    )
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    component_names = (
        "Re_total", "Im_total", "Re_regular", "Im_regular",
        "Re_drude", "Im_drude",
    )
    columns = ["omega"]
    for statistic in ("mean", "std", "min", "max"):
        columns.extend("{}_{}".format(statistic, name) for name in component_names)
    np.savetxt(
        str(output_path), statistics, fmt="%.16g",
        header=(
            "15 态等权光电导统计量；"
            "eta=0.065 M=600\n" + " ".join(columns)
        ),
    )
    return statistics


def load_spectrum(path):
    """Load the complete eight-level, fifteen-sector spectrum."""
    rows = np.loadtxt(str(path), comments="#", dtype=float)
    rows = np.atleast_2d(rows)
    if rows.shape != (120, 2):
        raise ValueError("expected 120 two-column spectrum rows")
    sectors = rows[:, 0].astype(int)
    if not np.allclose(rows[:, 0], sectors):
        raise ValueError("spectrum sector labels must be integers")
    counts = np.bincount(sectors, minlength=15)
    if len(counts) != 15 or not np.all(counts == 8):
        raise ValueError("expected eight levels in every sector 0 through 14")
    return rows


def _gnuplot_quote(path):
    text = str(path).replace("\\", "\\\\").replace('"', '\\"')
    return '"{}"'.format(text)


def generate_figures(spectrum_path, structure_path, optical_paths, output_dir):
    """Generate the spectrum, structure-factor, and four conductivity plots."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    spectrum = load_spectrum(spectrum_path)
    load_structure_factor(structure_path)
    average_path = (
        Path(optical_paths[0]).parent / "optical_manifold_average.dat"
    )
    statistics = write_optical_manifold_average(optical_paths, average_path)
    reference_frequencies = statistics[:, 0]

    outputs = {
        "spectrum": output_dir / "phc30_moderate_manybody_spectrum.png",
        "structure_factor": output_dir / "phc30_moderate_structure_factor.png",
        "regular_real": output_dir / "phc30_moderate_sigma_xx_regular_real.png",
        "regular_imag": output_dir / "phc30_moderate_sigma_xx_regular_imag.png",
        "total_real": output_dir / "phc30_moderate_sigma_xx_total_real.png",
        "total_imag": output_dir / "phc30_moderate_sigma_xx_total_imag.png",
    }
    gnuplot = shutil.which("gnuplot")
    if gnuplot is None:
        raise RuntimeError("gnuplot is required to generate result figures")

    summary = summarize_manifold(spectrum)
    spectrum_top = max(float(np.max(spectrum[:, 1])) * 1.05, 0.1)
    zoom_top = max(summary["next_energy"] * 1.25, summary["width"] * 2.0, 0.05)
    selected_sectors = ",".join(str(value) for value in summary["sectors"])
    spectrum_file = _gnuplot_quote(Path(spectrum_path).resolve())
    structure_file = _gnuplot_quote(Path(structure_path).resolve())

    lines = [
        "set datafile commentschars '#';",
        "set term pngcairo size 1800,760 enhanced font 'Noto Sans SC,14';",
        "set border linewidth 1.2; set grid ytics lc rgb '#d8d8d8' lw 1;",
        "set key opaque box top left;",
        "set output {};".format(_gnuplot_quote(outputs["spectrum"])),
        "set multiplot layout 1,2 title "
        "'倾斜 30-site ED：Np=12，t1=1，t3=0.2，V1=10，V2=V3=2';",
        "set xlabel '动量扇区 k'; set ylabel 'E-E0';",
        "set xrange [-0.5:14.5]; set xtics 0,1,14; set yrange [-0.02:{}];".format(spectrum_top),
        "set title '每个动量扇区保存的 8 个能级';",
        "plot {0} using (($2 <= {1:.16g}) ? $1 : 1/0):2 with points pt 7 ps 1.15 lc rgb '#c43b3b' title '最低 15 态', "
        "{0} using (($2 > {1:.16g}) ? $1 : 1/0):2 with points pt 7 ps 0.70 lc rgb '#607d9b' title '保存的高能态';".format(
            spectrum_file, summary["width"] + 1e-12
        ),
        "set ylabel ''; set yrange [{}:{}];".format(-0.005 * zoom_top, zoom_top),
        "set title '3x5 候选低能流形放大'; set key center right;",
        "set arrow 1 from graph 0, first {0:.16g} to graph 1, first {0:.16g} nohead dt 2 lw 1.5 lc rgb '#c43b3b';".format(summary["width"]),
        "set arrow 2 from graph 0, first {0:.16g} to graph 1, first {0:.16g} nohead dt 3 lw 1.5 lc rgb '#2b6ca3';".format(summary["next_energy"]),
        "set label 1 '最低态所在扇区：{}' at graph 0.03,0.94 front;".format(selected_sectors),
        "set label 2 '流形宽度 = {:.8g}' at graph 0.03,0.87 front;".format(summary["width"]),
        "set label 3 '第 15 至 16 态间隔 = {:.8g}' at graph 0.03,0.80 front;".format(summary["separation"]),
        "plot {0} using (($2 <= {1:.16g}) ? $1 : 1/0):2 with points pt 7 ps 1.15 lc rgb '#c43b3b' title '最低 15 态', "
        "{0} using (($2 > {1:.16g}) ? $1 : 1/0):2 with points pt 7 ps 0.70 lc rgb '#607d9b' title '保存的高能态';".format(
            spectrum_file, summary["width"] + 1e-12
        ),
        "unset multiplot; unset arrow 1; unset arrow 2; unset label 1; unset label 2; unset label 3; unset output;",
        "set term pngcairo size 1320,800 enhanced font 'Noto Sans SC,15';",
        "set output {};".format(_gnuplot_quote(outputs["structure_factor"])),
        "set title '温和参数 30-site 静态结构因子';",
        "set xlabel '动量 q'; set ylabel 'N(q)（原胞归一化）';",
        "set xrange [-0.3:14.3]; set xtics 0,1,14; set yrange [0:*]; set key opaque box top right;",
    ]
    structure_plots = [
        "{0} using (strcol(1) eq '{1}' ? $2 : 1/0):3 with linespoints pt 7 ps 0.35 lw 0.8 lc {2} notitle".format(
            structure_file, sector, sector + 1
        )
        for sector in range(15)
    ]
    structure_plots.append(
        "{} using (strcol(1) eq 'average' ? $2 : 1/0):3 "
        "with linespoints pt 6 ps 1.0 lw 2.7 lc rgb '#171717' title '15 态平均'".format(
            structure_file
        )
    )
    lines.append("plot " + ", ".join(structure_plots) + "; unset output;")

    components = (
        ("regular_real", 4, "Re sigma_xx^{reg}", "正则电导实部"),
        ("regular_imag", 5, "Im sigma_xx^{reg}", "正则电导虚部"),
        ("total_real", 2, "Re sigma_xx^{total}", "总电导实部"),
        ("total_imag", 3, "Im sigma_xx^{total}", "总电导虚部"),
    )
    for key, column, ylabel, title in components:
        curves = [_gnuplot_quote(Path(optical_paths[sector]).resolve()) for sector in (0, 5, 10)]
        average_file = _gnuplot_quote(average_path.resolve())
        minimum_column = column + 12
        maximum_column = column + 18
        lines.extend([
            "set output {};".format(_gnuplot_quote(outputs[key])),
            "set title '{}'; set xlabel '频率 omega'; set ylabel '{}';".format(title, ylabel),
            "set xrange [{}:{}]; set autoscale y; set key opaque box top right;".format(
                reference_frequencies[0], reference_frequencies[-1]
            ),
            "plot {4} using 1:{5}:{6} with filledcurves lc rgb '#d9d9d9' title '扇区极值范围', "
            "{0} using 1:{3} with lines lw 1.2 dt 3 lc rgb '#1f77b4' title '扇区 k=0', "
            "{1} using 1:{3} with lines lw 1.1 dt 2 lc rgb '#d62728' title '扇区 k=5', "
            "{2} using 1:{3} with lines lw 1.1 dt 4 lc rgb '#2f855a' title '扇区 k=10', "
            "{4} using 1:{3} with lines lw 2.8 dt 1 lc rgb '#171717' title '15 态平均'; unset output;".format(
                curves[0], curves[1], curves[2], column,
                average_file, minimum_column, maximum_column,
            ),
        ])

    with tempfile.NamedTemporaryFile(mode="w", suffix=".gnuplot", delete=False) as handle:
        script_path = Path(handle.name)
        handle.write("\n".join(lines) + "\n")
    try:
        environment = os.environ.copy()
        environment.update({"LC_ALL": "C.utf8", "LANG": "C.utf8"})
        subprocess.run(
            [gnuplot, str(script_path)], check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            env=environment,
        )
    finally:
        script_path.unlink()
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-dir", type=Path, default=DEFAULT_RESULT_DIR)
    parser.add_argument("--optical-dir", type=Path)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    arguments = parser.parse_args()
    data_dir = arguments.result_dir / "data"
    optical_dir = arguments.optical_dir or data_dir / "optical"
    optical_paths = {
        sector: optical_dir / "sector_{}_optical_response.dat".format(sector)
        for sector in range(15)
    }
    outputs = generate_figures(
        data_dir / "spectrum.dat",
        data_dir / "structure_factor_ground_manifold.dat",
        optical_paths,
        arguments.output_dir,
    )
    for key, path in outputs.items():
        print("{}: {}".format(key, path))


if __name__ == "__main__":
    main()
