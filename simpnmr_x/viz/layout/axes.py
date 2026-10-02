# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Suturina Group

"""Axis layout helpers for visualization."""

from matplotlib.axes import Axes


def rotate_x_tick_labels(ax: Axes, *, angle: float) -> None:
    """Rotate x-axis tick labels around their right-aligned anchors.

    Args:
        ax: Axes whose x-axis tick labels will be rotated.
        angle: Rotation angle in degrees.
    """
    for label in ax.get_xticklabels():
        label.set_rotation(angle)
        label.set_horizontalalignment("right")
        label.set_rotation_mode("anchor")
