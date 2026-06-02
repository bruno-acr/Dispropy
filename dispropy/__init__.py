"""Public API for dispropy."""

from .metrics import (
    calculate_disproportionality,
    calculate_ic,
    calculate_prr,
    calculate_ror,
)

ror = calculate_ror
prr = calculate_prr
ic = calculate_ic
disproportionality = calculate_disproportionality

__all__ = [
    "calculate_ror",
    "calculate_prr",
    "calculate_ic",
    "calculate_disproportionality",
    "ror",
    "prr",
    "ic",
    "disproportionality",
]
