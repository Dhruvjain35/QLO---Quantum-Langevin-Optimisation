"""Offline analysis utilities: barren-plateau statistics, shot-noise theory, fits.

These may use the analytic reference gradient; the estimators in
``qlo.gradients.finite_shot`` may not.
"""

from qlo.analysis.bp_variance import (
    bp_variance_per_parameter,
    bp_variance_summary,
    barren_plateau_variance,
    gradient_samples_across_inits,
)
from qlo.analysis.shot_theory import (
    theoretical_parameter_shift_variance,
    theoretical_parameter_shift_variance_vector,
)
from qlo.analysis.statistics import fit_loglog_slope, summarize_shot_noise

__all__ = [
    "gradient_samples_across_inits",
    "bp_variance_per_parameter",
    "bp_variance_summary",
    "barren_plateau_variance",
    "theoretical_parameter_shift_variance",
    "theoretical_parameter_shift_variance_vector",
    "fit_loglog_slope",
    "summarize_shot_noise",
]
