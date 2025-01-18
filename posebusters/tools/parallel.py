"""Provides some functions for parallel processing."""

from __future__ import annotations

from collections.abc import Iterable
from multiprocessing import Pool
from typing import Callable

from tqdm import tqdm


def parallel_map(function: Callable, iterable: Iterable, max_workers: int | None = None, chunk_size: int = 1, **kwargs):
    """Parallel map function."""
    with Pool(processes=max_workers, **kwargs) as pool:
        results = tqdm(
            pool.imap_unordered(function, iterable, chunk_size),
        )
    return results
