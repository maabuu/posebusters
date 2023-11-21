"""Provides functions for calculating errors."""

from __future__ import annotations


def ae(val_pred: float, val_true: float) -> float:
    """Calculate absolute error."""
    return abs(val_pred - val_true)


def pe(val_pred: float, val_true: float) -> float:
    """Calculate percentage error."""
    return (val_pred - val_true) / val_true


def ape(val_pred: float, val_true: float) -> float:
    """Calculate absolute percentage error."""
    return abs(pe(val_pred, val_true))


def bae(val: float, lb: float, ub: float) -> float:
    """Calculate out of bounds absolute error."""
    if val < lb:
        return ae(val, lb)
    if val > ub:
        return ae(val, ub)
    return 0.0


def bpe(val: float, lb: float, ub: float) -> float:
    """Calculate out of bounds percentage error."""
    if val < lb:
        return pe(val, lb)
    if val > ub:
        return pe(val, ub)
    return 0.0


def bape(val: float, lb: float, ub: float) -> float:
    """Calculate out of bounds absolute percentage error."""
    if val < lb:
        return ape(val, lb)
    if val > ub:
        return ape(val, ub)
    return 0.0
