#!/usr/bin/env python3

import argparse
import csv
import math
from pathlib import Path
from xml.sax.saxutils import escape


class PumpPoint:
    def __init__(
        self,
        series,
        step,
        phi,
        energy,
        q_left,
        q_mid,
        q_right,
        sweeps,
        attempts,
        converged,
        reason,
    ):
        self.series = series
        self.step = step
        self.phi = phi
        self.energy = energy
        self.q_left = q_left
        self.q_mid = q_mid
        self.q_right = q_right
        self.sweeps = sweeps
        self.attempts = attempts
        self.converged = converged
        self.reason = reason

    @property
    def x(self):
        return self.phi / (2.0 * math.pi)


SERIES = [
    ("forward_pi12", "forward, dphi=pi/12", "#1f77b4", "circle"),
    ("branch_from_fwd17_pi24", "fine branch, dphi=pi/24", "#d62728", "square"),
    ("branch_from_fwd17_pi36", "fine branch, dphi=pi/36", "#2ca02c", "triangle"),
]


def read_pumping(path, series):
    points = []
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            row = stripped.split()
            points.append(
                PumpPoint(
                    series=series,
                    step=int(row[0]),
                    phi=float(row[1]),
                    energy=float(row[2]),
                    q_left=float(row[3]),
                    q_mid=float(row[4]),
                    q_right=float(row[5]),
                    sweeps=int(row[6]),
                    attempts=int(row[7]),
                    converged=row[8].lower() == "true",
                    reason=row[9],
                )
            )
    return sorted(points, key=lambda p: (p.x, p.step))


def detect_switch_window(series_points):
    windows = []
    for points in series_points:
        ordered = sorted(points, key=lambda p: p.x)
        for prev, curr in zip(ordered, ordered[1:]):
            if prev.x > 0.70 and prev.q_mid > 0.15 and curr.q_mid < 0.10:
                windows.append((prev.x, curr.x))
    if not windows:
        return None
    return min(windows, key=lambda window: window[1] - window[0])


def write_plot_data(path, series_points):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "series",
                "step",
                "phi",
                "phi_over_2pi",
                "energy",
                "q_mid_raw",
                "sweeps",
                "attempts",
                "converged",
                "reason",
            ]
        )
        for point in sorted(series_points, key=lambda p: (p.series, p.x, p.step)):
            writer.writerow(
                [
                    point.series,
                    point.step,
                    "%.16g" % point.phi,
                    "%.16g" % point.x,
                    "%.16g" % point.energy,
                    "%.16g" % point.q_mid,
                    point.sweeps,
                    point.attempts,
                    str(point.converged).lower(),
                    point.reason,
                ]
            )


def fmt(value):
    if abs(value) < 1e-12:
        value = 0.0
    return "%.3g" % value


def polyline(points, xmap, ymap):
    return " ".join("%.2f,%.2f" % (xmap(x), ymap(y)) for x, y in points)


def marker_svg(kind, x, y, color):
    if kind == "circle":
        return '<circle cx="%.2f" cy="%.2f" r="4.1" fill="%s"/>' % (x, y, color)
    if kind == "square":
        return '<rect x="%.2f" y="%.2f" width="8" height="8" fill="%s"/>' % (x - 4, y - 4, color)
    return '<path d="M %.2f %.2f L %.2f %.2f L %.2f %.2f Z" fill="%s"/>' % (
        x,
        y - 4.8,
        x - 4.7,
        y + 4.1,
        x + 4.7,
        y + 4.1,
        color,
    )


def draw_panel(
    parts,
    panel,
    title,
    ylabel,
    ydomain,
    loaded,
    switch_window,
):
    left = 96
    top = panel["top"]
    width = 850
    height = 270
    bottom = top + height
    xdomain = (-0.02, 1.02)

    def xmap(x):
        return left + (x - xdomain[0]) / (xdomain[1] - xdomain[0]) * width

    def ymap(y):
        return bottom - (y - ydomain[0]) / (ydomain[1] - ydomain[0]) * height

    parts.append('<text x="%d" y="%d" class="title">%s</text>' % (left, top - 24, escape(title)))
    parts.append('<rect x="%d" y="%d" width="%d" height="%d" fill="white" stroke="#333" stroke-width="1"/>' % (left, top, width, height))

    if switch_window is not None:
        x0 = xmap(switch_window[0])
        x1 = xmap(switch_window[1])
        parts.append('<rect x="%.2f" y="%d" width="%.2f" height="%d" fill="#ddd" opacity="0.75"/>' % (x0, top, x1 - x0, height))

    for tick in [0.0, 0.25, 0.5, 0.75, 1.0]:
        x = xmap(tick)
        parts.append('<line x1="%.2f" y1="%d" x2="%.2f" y2="%d" class="grid"/>' % (x, top, x, bottom))
        parts.append('<line x1="%.2f" y1="%d" x2="%.2f" y2="%d" class="axis-tick"/>' % (x, bottom, x, bottom + 6))
        if panel.get("xlabel", False):
            parts.append('<text x="%.2f" y="%d" class="tick" text-anchor="middle">%s</text>' % (x, bottom + 22, fmt(tick)))

    ytick_start = math.ceil(ydomain[0] / 0.1) * 0.1
    ticks = []
    v = ytick_start
    while v <= ydomain[1] + 1e-9:
        ticks.append(round(v, 10))
        v += 0.1
    for tick in ticks:
        y = ymap(tick)
        parts.append('<line x1="%d" y1="%.2f" x2="%d" y2="%.2f" class="grid"/>' % (left, y, left + width, y))
        parts.append('<line x1="%d" y1="%.2f" x2="%d" y2="%.2f" class="axis-tick"/>' % (left - 6, y, left, y))
        parts.append('<text x="%d" y="%.2f" class="tick" text-anchor="end" dominant-baseline="middle">%s</text>' % (left - 10, y, fmt(tick)))

    ideal = [(i / 200.0, i / 200.0 / 3.0) for i in range(0, 201)]
    parts.append(
        '<polyline points="%s" fill="none" stroke="#333" stroke-width="1.6" stroke-dasharray="7,5"/>'
        % polyline(ideal, xmap, ymap)
    )
    parts.append('<line x1="%d" y1="%.2f" x2="%d" y2="%.2f" stroke="#aaa" stroke-width="1"/>' % (left, ymap(0.0), left + width, ymap(0.0)))
    parts.append('<line x1="%d" y1="%.2f" x2="%d" y2="%.2f" stroke="#aaa" stroke-width="1"/>' % (left, ymap(1.0 / 3.0), left + width, ymap(1.0 / 3.0)))

    for _rel, label, color, marker, points in loaded:
        values = [(p.x, p.q_mid) for p in points]
        parts.append(
            '<polyline points="%s" fill="none" stroke="%s" stroke-width="2.2"/>'
            % (polyline(values, xmap, ymap), color)
        )
        for x, y in values:
            parts.append(marker_svg(marker, xmap(x), ymap(y), color))

    parts.append('<text x="%d" y="%d" class="axis-label" text-anchor="middle">Phi_y / 2pi</text>' % (left + width / 2, bottom + 46))
    parts.append(
        '<text x="28" y="%.2f" class="axis-label" text-anchor="middle" transform="rotate(-90 28 %.2f)">%s</text>'
        % (top + height / 2, top + height / 2, escape(ylabel))
    )
    if switch_window is not None:
        sx = xmap(0.5 * (switch_window[0] + switch_window[1]))
        sy = ymap(0.02)
        parts.append('<text x="%.2f" y="%.2f" class="note" text-anchor="middle">sector jump</text>' % (sx, sy + 34))
        parts.append('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="#444" stroke-width="1.2" marker-end="url(#arrow)"/>' % (sx - 30, sy + 24, sx - 4, sy + 4))


def write_svg(path, loaded, switch_window):
    width = 1050
    height = 430
    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">' % (width, height, width, height),
        "<defs>",
        '<marker id="arrow" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto" markerUnits="strokeWidth">',
        '<path d="M0,0 L0,6 L6,3 z" fill="#444"/>',
        "</marker>",
        "</defs>",
        "<style>",
        "text{font-family:Arial,Helvetica,sans-serif;fill:#111}.title{font-size:18px;font-weight:700}.tick{font-size:12px}.axis-label{font-size:15px}.note{font-size:13px;fill:#444}.legend{font-size:13px}.grid{stroke:#e7e7e7;stroke-width:1}.axis-tick{stroke:#333;stroke-width:1}",
        "</style>",
        '<rect width="100%" height="100%" fill="white"/>',
    ]

    draw_panel(
        parts,
        {"top": 86, "xlabel": True},
        "Lx=15, Ly=6, nu=1/3: raw warm-start flux-insertion charge",
        "raw Q_mid",
        (-0.10, 0.36),
        loaded,
        switch_window,
    )

    legend_x = 680
    legend_y = 38
    parts.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#333" stroke-width="1.6" stroke-dasharray="7,5"/>' % (legend_x, legend_y, legend_x + 34, legend_y))
    parts.append('<text x="%d" y="%d" class="legend" dominant-baseline="middle">ideal nu=1/3</text>' % (legend_x + 44, legend_y))
    y = legend_y + 22
    for _rel, label, color, marker, _points in loaded:
        parts.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" stroke-width="2.2"/>' % (legend_x, y, legend_x + 34, y, color))
        parts.append(marker_svg(marker, legend_x + 17, y, color))
        parts.append('<text x="%d" y="%d" class="legend" dominant-baseline="middle">%s</text>' % (legend_x + 44, y, escape(label)))
        y += 22

    parts.append("</svg>")
    path.write_text("\n".join(parts))


def main():
    parser = argparse.ArgumentParser(description="Plot Lx=15 charge pump data as SVG.")
    parser.add_argument(
        "--campaign",
        type=Path,
        default=Path(__file__).resolve().parent / "flux_Lx15_pump_prod_20260710_170505",
    )
    parser.add_argument("--outdir", type=Path, default=None)
    args = parser.parse_args()

    campaign = args.campaign.resolve()
    outdir = (args.outdir or campaign / "plots").resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    loaded = []
    all_points = []
    for rel, label, color, marker in SERIES:
        path = campaign / rel / "pumping.dat"
        if not path.is_file():
            continue
        points = read_pumping(path, rel)
        loaded.append((rel, label, color, marker, points))
        all_points.extend(points)

    if not loaded:
        raise SystemExit("no pumping.dat files found under %s" % campaign)

    write_plot_data(outdir / "charge_pump_plot_data.csv", all_points)
    switch_window = detect_switch_window([points for _rel, _label, _color, _marker, points in loaded])
    svg = outdir / "charge_pump_raw.svg"
    write_svg(svg, loaded, switch_window)

    print("campaign=%s" % campaign)
    print("series=%d points=%d" % (len(loaded), len(all_points)))
    if switch_window is not None:
        print("switch_window_phi_over_2pi=%.6f,%.6f" % (switch_window[0], switch_window[1]))
    print("wrote=%s" % svg)
    print("wrote=%s" % (outdir / "charge_pump_plot_data.csv"))


if __name__ == "__main__":
    main()
