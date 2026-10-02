"""Verify susceptibility output scaling and lossless reduced CSV import."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from simpnmr_x.core.conv.a3.a3_to_cm3mol import A3_TO_CM3MOL
from simpnmr_x.core.conv.a3.a3_to_reduced import a3_to_reduced
from simpnmr_x.core.conv.reduced.reduced_to_a3 import reduced_to_a3
from simpnmr_x.core.fitting.variable_temperatures.components import (
    compute_curie_prefactor,
)
from simpnmr_x.io.csv.susc import read_susceptibilities_csv, save_susc


@pytest.mark.unit
@pytest.mark.parametrize("spin,total_J", [(1.0, None), (1.0, 2.5)])
def test_reduced_csv_round_trip_and_uncertainties(tmp_path, spin, total_J):
    tensor = np.array([[1.0, 0.2, 0.3], [0.2, 2.0, 0.4], [0.3, 0.4, 3.0]])
    molecules = []
    models = []
    for temperature in [200.0, 300.0]:
        susceptibility = SimpleNamespace(
            temperature=temperature,
            iso=2.0,
            axiality=1.5,
            rhombicity=0.5,
            tensor=tensor,
            dtensor=tensor - 2.0 * np.eye(3),
            eigvals=np.linalg.eigvalsh(tensor),
            alpha=10.0,
            beta=20.0,
            gamma=30.0,
        )
        molecules.append(
            SimpleNamespace(
                susc=susceptibility,
                electronic=SimpleNamespace(spin_S=spin, total_J=total_J),
            )
        )
        models.append(
            SimpleNamespace(
                r2=0.9,
                adj_r2=0.8,
                mae=0.1,
                rmse=0.2,
                fit_stdev={"iso": 0.1, "ax": 0.2, "alpha": 1.0},
            )
        )
    path = tmp_path / "reduced.csv"
    save_susc(molecules, str(path), susc_models=models, susc_units="reduced")
    table = pd.read_csv(path, comment="#")
    restored = read_susceptibilities_csv(str(path))
    for index, temperature in enumerate([200.0, 300.0]):
        scale = float(
            a3_to_reduced(1.0, temperature, compute_curie_prefactor(spin, total_J))
        )
        assert table.loc[index, "chi_iso (reduced)"] == pytest.approx(2.0 * scale)
        assert table.loc[index, "chi_iso-s-dev (reduced)"] == pytest.approx(0.1 * scale)
        assert table.loc[index, "alpha (degrees)"] == 10.0
        assert table.loc[index, "alpha-s-dev (degrees)"] == 1.0
        assert table.loc[index, "RMSE (ppm)"] == 0.2
        np.testing.assert_allclose(restored[index][0], tensor)
        assert restored[index][1] == temperature
        assert restored[index][2] == pytest.approx(2.0)
    assert "chi_alpha-s-dev (reduced)" not in table


@pytest.mark.unit
@pytest.mark.parametrize("convert", [a3_to_reduced, reduced_to_a3])
@pytest.mark.parametrize(
    "temperature,prefactor",
    [
        (0.0, 1.0),
        (float("nan"), 1.0),
        (300.0, 0.0),
        (300.0, -1.0),
        (300.0, float("inf")),
    ],
)
def test_reduced_conversion_rejects_invalid_inputs(convert, temperature, prefactor):
    with pytest.raises(ValueError):
        convert(1.0, temperature, prefactor)


@pytest.mark.unit
def test_reduced_conversions_support_arrays_and_preserve_inputs():
    values = np.array([1.0, 2.0])
    temperatures = np.array([200.0, 300.0])
    prefactor = compute_curie_prefactor(1.0, total_J=2.5)
    reduced = a3_to_reduced(values, temperatures, prefactor)
    np.testing.assert_allclose(reduced, values * temperatures / prefactor)
    np.testing.assert_allclose(reduced_to_a3(reduced, temperatures, prefactor), values)
    np.testing.assert_array_equal(values, [1.0, 2.0])


@pytest.mark.unit
@pytest.mark.parametrize("units", ["A3", "cm3 mol-1", "reduced"])
def test_plot_units_scale_chi_but_preserve_ratio_and_angles(monkeypatch, units):
    import matplotlib.pyplot as plt

    from simpnmr_x.viz.plots import fitted_shifts
    from simpnmr_x.viz.style.theme import build_spec

    nucleus = SimpleNamespace(
        signal_label="H1",
        signal_math_label="",
        label_nn="H",
        shift=SimpleNamespace(avg=1.0),
    )
    molecule = SimpleNamespace(
        nuclei=[nucleus],
        susc=SimpleNamespace(temperature=300.0, alpha=30.0, beta=20.0, gamma=10.0),
        electronic=SimpleNamespace(spin_S=1.0, total_J=2.5),
    )
    model = SimpleNamespace(
        VARNAMES=["iso", "rho_over_ax", "alpha"],
        final_var_values={"iso": 2.0, "rho_over_ax": 0.25, "alpha": 30.0},
        fit_vars={},
        fit_stdev={},
        adj_r2=0.9,
        mae=0.1,
        rmse=0.2,
    )
    captured = []

    def capture_table(ax, blocks, spec, **kwargs):
        captured.extend(blocks)

    monkeypatch.setattr(fitted_shifts, "render_compact_table", capture_table)
    monkeypatch.setattr(fitted_shifts, "render_figure", lambda fig, **kwargs: None)
    fig, ax = fitted_shifts.plot_fitted_shifts(
        molecule,
        {"H1": SimpleNamespace(shift=1.1)},
        model,
        spec=build_spec("paper"),
        susc_units=units,
        save=False,
        show=False,
        show_point_labels=False,
    )
    try:
        scale = {
            "A3": 1.0,
            "cm3 mol-1": A3_TO_CM3MOL,
            "reduced": float(
                a3_to_reduced(1.0, 300.0, compute_curie_prefactor(1.0, total_J=2.5))
            ),
        }[units]
        lines = captured[1][1]
        assert lines[0].endswith(f": {2.0 * scale:.3f}")
        assert lines[1].endswith(": 0.250")
        assert lines[2].endswith(": 30.000")
        assert "ppm" in ax.get_xlabel()
    finally:
        plt.close(fig)
