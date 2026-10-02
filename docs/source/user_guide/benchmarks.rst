Hyperfine Benchmarks
=====================

``simpnmr_x`` provides dedicated benchmark workflows for comparing hyperfine
sources. The supported benchmark subcommands are:

- ``benchmark a_fc``: compare isotropic Fermi-contact values ``A_fc``.
- ``benchmark a_sd``: compare axial spin-dipolar values ``A_sd``.

Benchmark workflows use a standalone YAML configuration file with an explicit
list of hyperfine sources, and they save diagnostic plots into the configured
project directory.

Configuration
-------------

The minimal configuration for ``benchmark a_fc`` looks like this:

.. code-block:: yaml

   project:
     name: A_FC_Benchmark

   signal_labels:
     file: signal_labels.csv

   nuclei:
     include: H

   benchmark:
     max_label_tolerance: 0.05

   hyperfine:
     - functional: B3LYP
       method: dft
       file: hfc/B3LYP.out

     - functional: PBE0
       method: csv
       file: hfc/PBE0_hyperfine.csv


Required configuration blocks:

- ``project``: the output directory for benchmark results.
- ``signal_labels``: a signal-label mapping file for grouping nuclei.
- ``nuclei``: the nuclei to include in the benchmark.
- ``hyperfine``: a list of hyperfine source blocks.

Each ``hyperfine`` entry contains:

- ``functional``: the functional name associated with the source.
- ``file``: the hyperfine source file path.
- ``method``: ``dft`` or ``csv``. Defaults to ``dft``.

The optional ``benchmark`` block supports:

- ``max_label_tolerance``: relative tolerance used when selecting a majority
  signal label for the max benchmark curve.

Running benchmarks
------------------

Run the ``A_fc`` benchmark:

.. code-block:: console

   simpnmr-x benchmark a_fc benchmark_hyperfine.yml

Run the ``A_sd`` benchmark:

.. code-block:: console

   simpnmr-x benchmark a_sd benchmark_hyperfine.yml

Plots are always saved to the project directory. To suppress interactive display,
use the ``--hide`` flag:

.. code-block:: console

   simpnmr-x --hide benchmark a_fc benchmark_hyperfine.yml

Outputs
-------

Both workflows create the directory specified in ``project:name``.

Files produced by ``benchmark a_fc``:

- ``A_FC_benchmark_max.csv``: maximum ``A_fc`` values by functional and nucleus.
- ``<NUCLEUS>_A_FC_benchmark_mean_curves.pdf``: one colour-coded mean ``A_fc``
  curve per signal across functionals, with math-formatted signal labels.

Files produced by ``benchmark a_sd``:

- ``<NUCLEUS>_A_SD_benchmark_mean_curves.pdf``: one colour-coded mean ``A_sd``
  curve per signal across functionals, with math-formatted signal labels.

Plot titles and y-axis labels denote the mean with an overline:
:math:`\overline{A}_{\mathrm{FC}}` and
:math:`\overline{A}_{\mathrm{SD}}^{\mathrm{ax}}`.
Each plotted value is the arithmetic mean of all nucleus values sharing the
same signal label and nucleus type across sources for that functional.

In each mean-curves plot, functionals are ordered from lower to higher mean
absolute coupling across the signals of that nucleus and hyperfine component.
The ordering is therefore determined separately for the H/C and A_fc/A_sd plots.

Notes
-----

- The ``A_fc`` benchmark compares Fermi-contact hyperfine components.
- The ``A_sd`` benchmark compares axial spin-dipolar hyperfine components.
- Unknown top-level YAML keys cause a configuration validation error.
