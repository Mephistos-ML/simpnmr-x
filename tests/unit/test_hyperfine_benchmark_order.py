# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Suturina Group

"""Unit tests for functional ordering in hyperfine benchmarks."""

from simpnmr_x.core.benchmarks.hyperfine.summary import (
    sort_functionals_by_mean_absolute_signal_value,
)


def test_functionals_order_by_mean_absolute_signal_mean():
    """Opposite signed couplings contribute by magnitude to the ranking."""
    summary = {
        "large": {"H": {"one": {"mean": -8.0}, "two": {"mean": 8.0}}},
        "small": {"H": {"one": {"mean": -1.0}, "two": {"mean": 3.0}}},
        "middle": {"H": {"one": {"mean": -4.0}, "two": {"mean": 6.0}}},
    }

    assert sort_functionals_by_mean_absolute_signal_value(summary, "H") == [
        "small",
        "middle",
        "large",
    ]
