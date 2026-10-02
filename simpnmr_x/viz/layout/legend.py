# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Suturina Group

"""Legend layout helpers for visualization."""

from matplotlib.axes import Axes
from matplotlib.legend import Legend


def place_legend_right(ax: Axes) -> Legend:
    """Place a single-column legend outside the axes, centred on their right.

    The legend participates in constrained layout so space is reserved beside
    the axes. Appearance is inherited from the active plotting style.

    Args:
        ax: Axes containing the labelled artists.

    Returns:
        The legend positioned to the right of the axes.
    """
    return ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), ncols=1)
