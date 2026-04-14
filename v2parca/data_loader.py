"""Convenience loader for the pre-computed ParCa state shipped in ``data/``.

The repo ships a gzipped pickle of the full v2parca composite state
(``data/parca_state.pkl.gz``) so downstream repos can skip the 70-minute
ParCa run. This module is the one-liner to load it.
"""

from __future__ import annotations

import gzip
import pickle
from pathlib import Path
from typing import Any


_DEFAULT_PATH = Path(__file__).resolve().parent.parent / 'data' / 'parca_state.pkl.gz'


def load_parca_state(path: str | Path | None = None) -> dict[str, Any]:
    """Return the gzipped, pre-computed v2parca composite state.

    Args:
        path: optional path to a ``parca_state.pkl.gz`` file. Defaults to
            the one shipped at ``data/parca_state.pkl.gz``.

    Returns:
        Dict of top-level stores (``process``, ``cell_specs``, ``mass``,
        ``constants``, ...). See ``data/README.md`` for the full layout.
    """
    p = Path(path) if path is not None else _DEFAULT_PATH
    with gzip.open(p, 'rb') as f:
        return pickle.load(f)
