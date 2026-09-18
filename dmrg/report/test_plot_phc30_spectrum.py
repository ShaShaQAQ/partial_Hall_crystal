from pathlib import Path
import sys

import numpy as np
import pytest


REPORT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(REPORT_DIR))

from plot_phc30_spectrum import (
    format_summary_text,
    load_spectrum,
    summarize_low_energy,
)


DATA_PATH = REPORT_DIR / "data" / "optical_response" / "phc30_spectrum_V1_100.dat"


def test_loads_eight_levels_in_each_of_fifteen_sectors():
    rows = load_spectrum(DATA_PATH)

    assert rows.shape == (120, 2)
    assert set(rows[:, 0].astype(int)) == set(range(15))
    assert all(np.count_nonzero(rows[:, 0] == k) == 8 for k in range(15))


def test_summarizes_fifteen_state_low_energy_manifold():
    summary = summarize_low_energy(load_spectrum(DATA_PATH))

    assert summary["manifold_size"] == 15
    assert summary["manifold_width"] == pytest.approx(0.0287121891)
    assert summary["sixteenth_energy"] == pytest.approx(0.2157939408)
    assert summary["separation"] == pytest.approx(0.1870817517)


def test_summary_annotation_uses_real_line_breaks():
    summary = summarize_low_energy(load_spectrum(DATA_PATH))

    text = format_summary_text(summary)

    assert "\\n" not in text
    assert text.count("\n") == 2
