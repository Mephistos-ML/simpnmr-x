# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Suturina Group

"""Generic plotting orchestration for hyperfine benchmark pipelines."""

import logging
import os
import re

from simpnmr_x.core.benchmarks.hyperfine.summary import (
    sort_functionals_by_mean_absolute_signal_value,
)
from simpnmr_x.viz.plots.benchmarks import (
    plot_hyperfine_metric_functional_mean_curves,
)
from simpnmr_x.viz.style.theme import PlotSpec

logger = logging.getLogger(__name__)


def plot_hyperfine_functional_mean_summary(
    summary: dict[str, dict[str, dict[str, dict[str, object]]]],
    *,
    isotopes_by_nucleus: dict[str, str],
    output_dir: str,
    spec: PlotSpec,
    show: bool,
    y_label: str,
    title_metric: str,
    filename_metric: str,
    log_metric: str,
    window_metric: str,
) -> None:
    """Plot signal-mean hyperfine curves across functionals for each nucleus.

    Args:
        summary: Functional, nucleus, and signal summaries with ``mean`` values.
        isotopes_by_nucleus: Nucleus labels mapped to domain isotope labels.
        output_dir: Directory for generated plots.
        spec: Resolved plotting style.
        show: Whether to display the plots.
        y_label: Y-axis label, including the metric and units.
        title_metric: Metric label used in figure titles.
        filename_metric: Metric token used in output file names.
        log_metric: Metric label used in log messages.
        window_metric: Metric label used in Matplotlib window titles.
    """
    plotted_nuclei: list[str] = []
    nuclei = dict.fromkeys(
        nucleus_label
        for functional_summary in summary.values()
        for nucleus_label in functional_summary
    )
    for nucleus_label in nuclei:
        isotope = isotopes_by_nucleus[nucleus_label]
        safe_nucleus = safe_filename_token(nucleus_label)
        save_name = os.path.join(
            output_dir,
            f"{safe_nucleus}_{filename_metric}_benchmark_mean_curves",
        )
        plot_hyperfine_metric_functional_mean_curves(
            nucleus_label=nucleus_label,
            isotope=isotope,
            summary=summary,
            functional_order=sort_functionals_by_mean_absolute_signal_value(
                summary, nucleus_label
            ),
            spec=spec,
            y_label=y_label,
            title_metric=title_metric,
            save=True,
            show=show,
            save_name=save_name,
            window_title=f"{window_metric} mean benchmark: {isotope}",
        )
        plotted_nuclei.append(nucleus_label)

    logger.info(
        "%s signal-mean %s benchmark plots saved to %s",
        ", ".join(plotted_nuclei),
        log_metric,
        output_dir,
    )


def safe_filename_token(value: str) -> str:
    """Return a filesystem-safe token for generated benchmark plot names."""
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return token.strip("_") or "functional"
