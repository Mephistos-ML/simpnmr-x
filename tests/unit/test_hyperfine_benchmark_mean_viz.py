# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Suturina Group

"""Unit tests for signal-mean hyperfine benchmark curves."""

import matplotlib.pyplot as plt
import numpy as np
import pytest

from simpnmr_x.viz.plots import benchmarks
from simpnmr_x.viz.style.theme import build_spec


@pytest.mark.parametrize(
    "y_label, title_metric",
    [
        (
            r"$\overline{A}_\mathregular{FC}$ (ppm Å$^\mathregular{-3}$)",
            r"$\overline{A}_\mathregular{FC}$",
        ),
        (
            r"$\overline{A}_{\mathregular{SD}}^{\mathregular{ax}}$ "
            r"(ppm Å$^\mathregular{-3}$)",
            r"$\overline{A}_{\mathregular{SD}}^{\mathregular{ax}}$",
        ),
    ],
)
@pytest.mark.parametrize(
    "nucleus_label, isotope, formatted_isotope",
    [
        ("H", "1H", r"$^\mathregular{1} \mathregular{H}$"),
        ("C", "13C", r"$^\mathregular{13} \mathregular{C}$"),
    ],
)
def test_functional_curves_plot_every_signal_mean_with_unique_math_labels(
    monkeypatch, y_label, title_metric, nucleus_label, isotope, formatted_isotope
):
    """One coloured curve per signal follows the source functional order."""
    monkeypatch.setattr(benchmarks, "render_figure", lambda fig, **kwargs: None)
    summary = {
        "B3LYP": {
            nucleus_label: {
                "H_BDI_γ": {
                    "mean": -3.0,
                    "signal_math_label": r"$\mathrm{H}_{\gamma}$",
                },
                "H_BDI_pMe": {
                    "mean": 2.0,
                    "signal_math_label": r"$\mathrm{H}_{\mathrm{pMe}}$",
                },
            }
        },
        "PBE": {
            nucleus_label: {
                "H_BDI_γ": {
                    "mean": -4.0,
                    "signal_math_label": r"$\mathrm{H}_{\gamma}$",
                },
                "H_BDI_pMe": {
                    "mean": 5.0,
                    "signal_math_label": r"$\mathrm{H}_{\mathrm{pMe}}$",
                },
            }
        },
    }
    spec = build_spec("paper")
    with spec.context():
        fig, ax = benchmarks.plot_hyperfine_metric_functional_mean_curves(
            nucleus_label,
            summary,
            isotope=isotope,
            functional_order=["PBE", "B3LYP"],
            spec=spec,
            y_label=y_label,
            title_metric=title_metric,
            save_name="mean_curves.pdf",
            window_title="Mean benchmark",
            save=False,
            show=False,
        )

    lines = [line for line in ax.lines if line.get_label().startswith("$")]
    assert ax.get_ylabel() == y_label
    assert ax.get_title() == (
        f"{formatted_isotope} {title_metric}\nFunctional dependence"
    )
    assert len(lines) == 2
    np.testing.assert_allclose(lines[0].get_ydata(), [-4.0, -3.0])
    np.testing.assert_allclose(lines[1].get_ydata(), [5.0, 2.0])
    assert lines[0].get_label() == r"$\mathrm{H}_{\gamma}$"
    assert lines[1].get_label() == r"$\mathrm{H}_{\mathrm{pMe}}$"
    assert len({tuple(line.get_color()) for line in lines}) == 2
    assert [tick.get_text() for tick in ax.get_xticklabels()] == ["PBE", "B3LYP"]
    assert all(tick.get_rotation() == 45 for tick in ax.get_xticklabels())
    assert all(
        tick.get_horizontalalignment() == "right" for tick in ax.get_xticklabels()
    )
    legend = ax.get_legend()
    assert legend is not None
    assert {text.get_text() for text in legend.get_texts()} == {
        r"$\mathrm{H}_{\gamma}$",
        r"$\mathrm{H}_{\mathrm{pMe}}$",
    }
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    axes_bounds = ax.get_window_extent(renderer)
    legend_bounds = legend.get_window_extent(renderer)
    assert legend_bounds.x0 >= axes_bounds.x1
    assert (legend_bounds.y0 + legend_bounds.y1) / 2 == pytest.approx(
        (axes_bounds.y0 + axes_bounds.y1) / 2
    )
    text_positions = [
        text.get_window_extent(renderer).x0 for text in legend.get_texts()
    ]
    assert text_positions[0] == pytest.approx(text_positions[1])
    plt.close(fig)
