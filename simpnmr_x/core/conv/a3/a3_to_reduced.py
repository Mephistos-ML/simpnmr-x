# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Suturina Group

"""Convert susceptibility values from Å³ to reduced χT/C."""

import numpy as np
from numpy.typing import NDArray


def a3_to_reduced(
    values: NDArray | float,
    temperature: NDArray | float,
    curie_prefactor: NDArray | float,
) -> NDArray:
    """Convert susceptibility values from Å³ to reduced χT/C.

    Args:
        values: Scalar or array in Å³.
        temperature: Positive finite temperature in kelvin.
        curie_prefactor: Positive finite Curie prefactor in Å³ K, computed
            using J when available, otherwise S.

    Returns:
        Values in reduced χT/C, with NumPy broadcasting of the supplied inputs.

    Raises:
        ValueError: If temperature or Curie prefactor is non-positive or non-finite.
    """
    temperatures = np.asarray(temperature, dtype=float)
    prefactors = np.asarray(curie_prefactor, dtype=float)
    if not np.all(np.isfinite(temperatures) & (temperatures > 0)):
        raise ValueError("Reduced conversion requires positive finite temperature")
    if not np.all(np.isfinite(prefactors) & (prefactors > 0)):
        raise ValueError("Reduced conversion requires positive finite Curie prefactor")
    return np.asarray(values, dtype=float) * (temperatures / prefactors)
