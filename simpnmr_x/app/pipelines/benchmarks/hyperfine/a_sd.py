# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Suturina Group

"""Run A_sd benchmark workflows."""

import os

from simpnmr_x.app.params.options import BenchmarkAsdRunOptions
from simpnmr_x.app.pipelines.benchmarks.hyperfine.runner import (
    plot_hyperfine_functional_mean_summary,
)
from simpnmr_x.app.pipelines.benchmarks.hyperfine.sources import (
    group_loaded_sources_by_functional,
    load_hyperfine_benchmark_sources,
)
from simpnmr_x.core.benchmarks.hyperfine.a_sd import (
    summarize_a_sd_ranges_by_functional_and_nucleus,
)
from simpnmr_x.viz.style.theme import apply_profile


def run_benchmark_a_sd(config, options: BenchmarkAsdRunOptions | None = None) -> int:
    """Run the A_sd benchmark workflow from a YAML configuration."""
    if options is None:
        raise ValueError("BenchmarkAsdRunOptions is required")

    os.makedirs(config.project_name, exist_ok=True)
    spec = apply_profile(options.runtime.plot_profile)

    if options.dry_run:
        return 0

    signals = load_hyperfine_benchmark_sources(config)

    a_sd_summary = summarize_a_sd_ranges_by_functional_and_nucleus(
        group_loaded_sources_by_functional(signals)
    )
    plot_hyperfine_functional_mean_summary(
        a_sd_summary,
        isotopes_by_nucleus={
            nucleus.label_nn: nucleus.isotope
            for source in signals
            for nucleus in source["molecule"].nuclei
        },
        output_dir=config.project_name,
        spec=spec,
        show=options.runtime.show_plots,
        y_label=(
            r"$\overline{A}_{\mathregular{SD}}^{\mathregular{ax}}$ "
            r"(ppm Å$^\mathregular{-3}$)"
        ),
        title_metric=r"$\overline{A}_{\mathregular{SD}}^{\mathregular{ax}}$",
        filename_metric="A_SD",
        log_metric="A_sd",
        window_metric="A_sd",
    )

    return 0
