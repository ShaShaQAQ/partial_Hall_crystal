#!/usr/bin/env python3
"""Package traceable small files from the full 30-site optical run."""

import argparse
from datetime import datetime
import hashlib
from pathlib import Path
import shutil

import numpy as np

from plot_phc30_moderate_fqahc import (
    load_optical_curve,
    write_optical_manifold_average,
)


REPORT_DIR = Path(__file__).resolve().parent
DEFAULT_RESULT_DIR = (
    REPORT_DIR.parents[1] / "results" /
    "phc30_np12_v1_10_v2_2_v3_2_fqahc"
)
DEFAULT_SOURCE_DIR = Path(
    "/home/public/shajy/codex_runs/"
    "phc30_np12_v1_10_v2_2_v3_2_fqahc/optical"
)


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_response_header(path):
    first_line = Path(path).read_text().splitlines()[0]
    if not first_line.startswith("# "):
        raise ValueError("response file has no metadata header: {}".format(path))
    fields = {}
    for token in first_line[2:].split():
        if "=" in token:
            key, value = token.split("=", 1)
            fields[key] = value
    return fields


def validate_response_header(fields, sector, path):
    expected = {
        "sector": str(sector),
        "lattice": "corrected",
        "Np": "12",
        "t1": "1.0",
        "t3": "0.2",
        "V1": "10.0",
        "V2": "2.0",
        "V3": "2.0",
        "requested_mmax": "600",
        "lanczos_steps": "600",
        "lanczos_breakdown": "false",
    }
    for key, value in expected.items():
        if fields.get(key) != value:
            raise ValueError(
                "{}: expected {}={}, got {}".format(
                    path, key, value, fields.get(key)
                )
            )
    residual = float(fields.get("residual", "nan"))
    if not np.isfinite(residual) or residual >= 1e-7:
        raise ValueError("{}: invalid residual {}".format(path, residual))
    return residual


def _toml_string(value):
    return '"{}"'.format(str(value).replace("\\", "\\\\").replace('"', '\\"'))


def package_optical_results(source_dir, result_dir):
    source_dir = Path(source_dir).resolve()
    result_dir = Path(result_dir).resolve()
    output_dir = result_dir / "data" / "optical"
    output_dir.mkdir(parents=True, exist_ok=True)

    records = []
    copied_paths = {}
    curves = {}
    for sector in range(15):
        stem = "sector_{}_optical_response".format(sector)
        data_source = source_dir / (stem + ".dat")
        jld2_source = source_dir / (stem + ".jld2")
        if not data_source.is_file() or not jld2_source.is_file():
            raise FileNotFoundError("missing optical outputs for sector {}".format(sector))
        fields = parse_response_header(data_source)
        residual = validate_response_header(fields, sector, data_source)
        curve = load_optical_curve(data_source)
        data_destination = output_dir / data_source.name
        shutil.copy2(str(data_source), str(data_destination))
        copied_paths[sector] = data_destination
        curves[sector] = curve
        records.append(
            {
                "sector": sector,
                "data_path": data_source,
                "data_size": data_source.stat().st_size,
                "data_sha256": sha256_file(data_source),
                "jld2_path": jld2_source,
                "jld2_size": jld2_source.stat().st_size,
                "jld2_sha256": sha256_file(jld2_source),
                "ground_energy": float(fields["Eg"]),
                "residual": residual,
                "source_norm2": float(fields["source_norm2"]),
                "diamagnetic_expectation": float(fields["Kexp"]),
                "drude_weight": float(fields["drude_weight"]),
            }
        )

    average_path = output_dir / "optical_manifold_average.dat"
    statistics = write_optical_manifold_average(copied_paths, average_path)
    manifest_path = result_dir / "optical_manifest.toml"
    with manifest_path.open("w") as output:
        output.write("schema_version = 1\n")
        output.write('result_id = "phc30_np12_v1_10_v2_2_v3_2_fqahc"\n')
        output.write('storage_host = "W003"\n')
        output.write("source_directory = {}\n".format(_toml_string(source_dir)))
        output.write("generated_at = {}\n".format(_toml_string(datetime.now().isoformat())))
        for record in records:
            output.write("\n[[responses]]\n")
            output.write("sector = {}\n".format(record["sector"]))
            output.write("data_path = {}\n".format(_toml_string(record["data_path"])))
            output.write("data_size_bytes = {}\n".format(record["data_size"]))
            output.write("data_sha256 = {}\n".format(_toml_string(record["data_sha256"])))
            output.write("jld2_path = {}\n".format(_toml_string(record["jld2_path"])))
            output.write("jld2_size_bytes = {}\n".format(record["jld2_size"]))
            output.write("jld2_sha256 = {}\n".format(_toml_string(record["jld2_sha256"])))
            for key in (
                "ground_energy", "residual", "source_norm2",
                "diamagnetic_expectation", "drude_weight",
            ):
                output.write("{} = {:.16g}\n".format(key, record[key]))

    drude_weights = np.asarray([record["drude_weight"] for record in records])
    maximum_residual = max(record["residual"] for record in records)
    average_regular_real = statistics[:, 3]
    peak_index = int(np.argmax(average_regular_real))
    diagnostics_path = result_dir / "data" / "optical_diagnostics.txt"
    with diagnostics_path.open("w") as output:
        output.write("sectors={}\n".format(",".join(str(k) for k in range(15))))
        output.write("maximum_residual={:.16g}\n".format(maximum_residual))
        output.write("drude_weight_mean={:.16g}\n".format(drude_weights.mean()))
        output.write("drude_weight_std={:.16g}\n".format(drude_weights.std()))
        output.write("drude_weight_min={:.16g}\n".format(drude_weights.min()))
        output.write("drude_weight_max={:.16g}\n".format(drude_weights.max()))
        output.write("average_regular_peak_omega={:.16g}\n".format(statistics[peak_index, 0]))
        output.write("average_regular_peak_value={:.16g}\n".format(average_regular_real[peak_index]))

    return {
        "sectors": tuple(range(15)),
        "maximum_residual": maximum_residual,
        "average_path": average_path,
        "manifest_path": manifest_path,
        "diagnostics_path": diagnostics_path,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--result-dir", type=Path, default=DEFAULT_RESULT_DIR)
    arguments = parser.parse_args()
    summary = package_optical_results(arguments.source_dir, arguments.result_dir)
    print("packaged sectors: {}".format(summary["sectors"]))
    print("maximum residual: {:.6e}".format(summary["maximum_residual"]))
    print("manifest: {}".format(summary["manifest_path"]))


if __name__ == "__main__":
    main()
