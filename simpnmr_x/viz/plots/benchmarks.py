# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Suturina Group

"""Plot benchmark summaries."""

import matplotlib.pyplot as plt
import numpy as np

from simpnmr_x.viz.layout.axes import rotate_x_tick_labels
from simpnmr_x.viz.layout.canvas import create_canvas
from simpnmr_x.viz.layout.export import render_figure
from simpnmr_x.viz.layout.legend import place_legend_right
from simpnmr_x.viz.style.theme import PlotSpec
from simpnmr_x.viz.utils.fmt import isotope_format


def plot_hyperfine_metric_functional_mean_curves(
    nucleus_label: str,
    summary: dict[str, dict[str, dict[str, dict[str, object]]]],
    *,
    isotope: str,
    functional_order: list[str],
    spec: PlotSpec,
    y_label: str,
    title_metric: str,
    save: bool,
    show: bool,
    save_name: str,
    window_title: str,
) -> tuple[plt.Figure, plt.Axes]:
    """Plot one signal-mean series per label across configured functionals.

    Args:
        nucleus_label: Nucleus key used to select entries from the summary.
        summary: Functional, nucleus, and signal summaries with ``mean`` values.
        isotope: Domain isotope label displayed in the title, e.g. ``"13C"``.
        functional_order: Functional labels in the required x-axis order.
        spec: Resolved plotting style.
        y_label: Y-axis label, including the metric and units.
        title_metric: Metric label used in the figure title.
        save: Whether to save the figure.
        show: Whether to display the figure.
        save_name: Output image file name.
        window_title: Matplotlib window title.

    Returns:
        A tuple ``(fig, ax)``.
    """
    fig, ax = create_canvas(
        spec.profile,
        variant="vertical_extended",
        window_title=window_title,
        layout="constrained",
    )

    glyphs = spec.glyphs
    spec.skin_axes(ax)
    palette = spec.palette

    functionals = functional_order
    signals = dict.fromkeys(
        signal_label
        for functional_summary in summary.values()
        for signal_label in functional_summary.get(nucleus_label, {})
    )
    x_positions = np.arange(len(functionals))
    colours = plt.colormaps["turbo"](np.linspace(0.06, 0.94, len(signals)))
    for colour, signal_label in zip(colours, signals):
        values = [
            float(summary[functional][nucleus_label][signal_label]["mean"])
            if signal_label in summary[functional].get(nucleus_label, {})
            else np.nan
            for functional in functionals
        ]
        signal_entry = next(
            functional_summary[nucleus_label][signal_label]
            for functional_summary in summary.values()
            if signal_label in functional_summary.get(nucleus_label, {})
        )
        math_label = str(signal_entry["signal_math_label"])
        ax.plot(
            x_positions,
            values,
            color=colour,
            lw=glyphs.line_lw,
            marker=glyphs.marker,
            markersize=glyphs.ms,
            label=math_label,
        )

    ax.axhline(
        0.0,
        color=palette.primary,
        lw=glyphs.line_lw,
    )
    ax.grid(axis="y", color=palette.grid)

    ax.set_xticks(x_positions)
    ax.set_xticklabels(functionals)
    rotate_x_tick_labels(ax, angle=45)

    ax.set_ylabel(y_label)
    ax.set_xlabel("Functional")
    ax.set_title(f"{isotope_format(isotope)} {title_metric}\nFunctional dependence")
    place_legend_right(ax)

    render_figure(
        fig,
        save=save,
        show=show,
        save_name=save_name,
    )

    return fig, ax
