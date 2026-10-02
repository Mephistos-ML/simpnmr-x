# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Suturina Group

"""Run the temperature-dependent susceptibility plotting workflow."""

from __future__ import annotations

import numpy as np

from simpnmr_x.app.params.options import PlotChiTRunOptions
from simpnmr_x.cfg.plot_chit import ChiTSourceConfig, PlotChiTConfig
from simpnmr_x.core.conv.a3_to_cm3mol import a3_to_cm3mol
from simpnmr_x.core.domain.tensor import canonical_principal_axes
from simpnmr_x.core.fitting.variable_temperatures.components import (
    calculate_E_D_components,
    compute_analytic_component,
    compute_curie_prefactor,
    compute_g_components,
    compute_g_sq_components,
    rotate_tensors_to_frame,
    validate_common_principal_axes,
)
from simpnmr_x.io.qc import gateway as rdrs
from simpnmr_x.viz.plots.susc import plot_chit_comparison
from simpnmr_x.viz.style.theme import apply_profile

_G_FRAME_ALIGNMENT_TOLERANCE = 1.0e-2


def run_plot_chit(config: PlotChiTConfig, options: PlotChiTRunOptions) -> int:
    """Load configured ORCA sources and render a χT comparison plot.

    Args:
        config: Validated ``plot_chit`` configuration.
        options: Runtime plotting options from the CLI.

    Returns:
        Zero after the plot has been rendered.
    """

    series: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    if config.xrd is not None:
        series["XRD Geometry"] = _read_chi_t_series(config.xrd)
    if config.opt is not None:
        opt_series = _read_chi_t_series(config.opt)
        series["Opt. Geometry"] = opt_series
        maximum_temperature = float(np.max(opt_series[0]))
        if config.tip is not None:
            reference_temperature = config.tip.reference_temperature
            analytic_reference_temperature = (
                maximum_temperature
                if reference_temperature == "max"
                else float(reference_temperature)
            )
            analytic_reference = _analytic_chi_t_series(
                config.opt,
                opt_series[0],
                analytic_reference_temperature,
            )
            series["Opt. Geom (TIP excl.)"] = _remove_analytic_tip(
                opt_series,
                analytic_reference,
                reference_temperature,
            )

    temperature_limits = None
    if config.temperature is not None:
        temperature_limits = (
            config.temperature.minimum,
            config.temperature.maximum,
        )

    spec = apply_profile(options.runtime.plot_profile)
    with spec.context():
        plot_chit_comparison(
            series,
            spec,
            show=options.runtime.show_plots,
            save=True,
            save_name=config.output_file,
            temperature_limits=temperature_limits,
        )

    return 0


def _read_chi_t_series(
    source: ChiTSourceConfig,
) -> tuple[np.ndarray, np.ndarray]:
    values = rdrs.read_orca_chi_t(source.file, source.section)
    temperatures = np.asarray(list(values.keys()), dtype=float)
    chi_t = np.asarray(list(values.values()), dtype=float)
    return temperatures, chi_t


def _remove_analytic_tip(
    opt_series: tuple[np.ndarray, np.ndarray],
    analytic_series: tuple[np.ndarray, np.ndarray],
    reference_temperature: str | float,
) -> tuple[np.ndarray, np.ndarray]:
    """Remove the fitted TIP contribution from the optical-geometry series.

    The TIP slope is obtained from the difference between the experimental and
    analytic chiT values at the selected reference temperature.

    Args:
        opt_series: Temperatures and chiT values from the optical geometry.
        analytic_series: Temperatures and analytic chiT values on the same grid.
        reference_temperature: Temperature used to determine the TIP slope, or
            ``"max"`` to use the highest temperature in the series.

    Returns:
        The original temperatures and chiT values with the TIP slope removed.
    """
    temperatures, chi_t = opt_series
    _, analytic_chi_t = analytic_series
    reference_index = _reference_index(temperatures, reference_temperature)
    reference = float(temperatures[reference_index])
    tip_chi_t_reference = float(
        chi_t[reference_index] - analytic_chi_t[reference_index]
    )
    tip_chi = tip_chi_t_reference / reference

    return temperatures, chi_t - tip_chi * temperatures


def _analytic_chi_t_series(
    source: ChiTSourceConfig,
    temperatures: np.ndarray,
    reference_temperature: str | float,
) -> tuple[np.ndarray, np.ndarray]:
    """Calculate analytic chiT values over a temperature grid.

    The source tensors are transformed using the susceptibility frame at the
    selected reference temperature.

    Args:
        source: ORCA source containing susceptibility, g, and effective-H data.
        temperatures: Temperatures at which to evaluate analytic chiT.
        reference_temperature: Temperature selecting the susceptibility frame,
            or ``"max"`` to use the highest temperature in the grid.

    Returns:
        The input temperatures and corresponding analytic chiT values.
    """
    temperatures = np.asarray(temperatures, dtype=float)
    reference_index = _reference_index(temperatures, reference_temperature)
    reference = float(temperatures[reference_index])

    tensors = rdrs.read_orca_susceptibility(source.file, source.section)
    tensor = _get_temperature_value(tensors, reference)
    _, chi_frame = canonical_principal_axes(tensor)

    g_tensor = rdrs.read_g_tensor_ab_initio(source.file, source.section)
    eff_h = rdrs.read_eff_hamiltonian_tensor(source.file, source.section)
    if g_tensor is None or eff_h is None:
        raise ValueError(
            "Analytic chiT requires both the ORCA g-tensor and effective Hamiltonian"
        )

    eff_h_frame, g_frame = rotate_tensors_to_frame(eff_h, g_tensor, chi_frame)
    validate_common_principal_axes(
        eff_h_frame,
        g_frame,
        _G_FRAME_ALIGNMENT_TOLERANCE,
    )

    spin = rdrs.read_orca_spin(source.file)
    D_J, E_J = calculate_E_D_components(eff_h_frame)
    analytic_chi = compute_analytic_component(
        "iso",
        temperatures,
        compute_g_sq_components(g_frame),
        compute_g_components(g_frame),
        D_J,
        E_J,
        spin,
    )
    analytic_chi_t = a3_to_cm3mol(
        analytic_chi * temperatures * compute_curie_prefactor(spin)
    )
    return temperatures, analytic_chi_t


def _reference_index(
    temperatures: np.ndarray,
    reference_temperature: str | float,
) -> int:
    if reference_temperature == "max":
        return int(np.argmax(temperatures))

    matches = np.flatnonzero(
        np.isclose(temperatures, float(reference_temperature), rtol=0.0, atol=1.0e-8)
    )
    if matches.size == 0:
        raise ValueError(
            "TIP reference_temperature is not present in the OPT temperature grid"
        )
    return int(matches[0])


def _get_temperature_value(values: dict[float, np.ndarray], temperature: float):
    for value_temperature, value in values.items():
        if np.isclose(value_temperature, temperature, rtol=0.0, atol=1.0e-8):
            return value
    raise ValueError("TIP reference_temperature is not present in the OPT tensor grid")
