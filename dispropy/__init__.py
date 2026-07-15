"""Disproportionality metrics for pharmacovigilance."""

from .diagnostics import GPSFitWarning
from .metrics import (
    calculate_disproportionality,
    calculate_ebgm,
    calculate_ic,
    calculate_prr,
    calculate_ror,
)

ror = calculate_ror
prr = calculate_prr
ic = calculate_ic
ebgm = calculate_ebgm
disproportionality = calculate_disproportionality

__all__ = [
    "calculate_ror",
    "calculate_prr",
    "calculate_ic",
    "calculate_ebgm",
    "calculate_disproportionality",
    "ror",
    "prr",
    "ic",
    "ebgm",
    "disproportionality",
    "GPSFitWarning",
]

__version__ = "0.4.0"
