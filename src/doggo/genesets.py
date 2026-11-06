"""Loading and filtering gene set collections (GMT format)."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable


def read_gmt(path: str | Path) -> dict[str, list[str]]:
    """Read a GMT file into ``{term: [genes, ...]}``.

    Each line is ``term<TAB>description<TAB>gene1<TAB>gene2...``. Duplicate
    genes within a set are dropped (first occurrence kept); empty lines and
    lines starting with ``#`` are ignored. If a term appears more than once,
    the gene lists are merged.
    """
    sets: dict[str, list[str]] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{lineno}: expected at least 2 tab-separated fields")
            term = fields[0].strip()
            genes = [g.strip() for g in fields[2:] if g.strip()]
            existing = sets.setdefault(term, [])
            seen = set(existing)
            for g in genes:
                if g not in seen:
                    existing.append(g)
                    seen.add(g)
    return sets


def filter_gene_sets(
    sets: dict[str, list[str]],
    prefixes: Iterable[str] | None = None,
    terms: Iterable[str] | None = None,
) -> dict[str, list[str]]:
    """Keep sets whose name starts with any of ``prefixes`` and/or is in ``terms``.

    With neither given, all sets are returned.
    """
    prefixes = tuple(prefixes or ())
    terms = set(terms or ())
    if not prefixes and not terms:
        return dict(sets)
    return {
        t: g
        for t, g in sets.items()
        if (prefixes and t.startswith(prefixes)) or t in terms
    }
