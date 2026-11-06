"""Competitive gene set tests for a per-gene score.

For each gene set, the rows of the input table whose gene is in the set are
compared with all rows (the background). The tested quantity is the set's mean
of a per-row value: |score| (``abs``), score (``greater``/``two-sided``) or
-score (``less``).

Two methods compute the p-value:

``zscore`` (default)
    Closed-form test of the set mean against the rest, following limma's
    ``cameraPR`` (Wu & Smyth, 2012). The variance of the set mean is
    inflated by ``1 + (n - 1) * inter_gene_cor`` to account for genes in a
    set being correlated; with ``inter_gene_cor = 0`` it reduces to a
    two-sample t-test that assumes independent genes. p-values come from a
    t distribution with N - 2 degrees of freedom: instant, and no resolution
    floor.

``resample``
    The original doggo Monte Carlo test. The statistic is computed on random
    samples of the same size drawn from all rows, and the p-value is the
    fraction at least as extreme as the observed one. It makes no
    distributional assumption, but assumes independent genes (no
    correlation adjustment) and its smallest p-value is 1/(max_reps+1).

Sampling uses sequential early stopping (Besag & Clifford, 1991): drawing
stops once ``min_hits`` extreme samples have been seen, giving p = h / l.
Sets that never reach ``min_hits`` within ``max_reps`` draws get the
add-one estimate p = (k + 1) / (max_reps + 1), so p is never 0.
"""

from __future__ import annotations

import logging
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

log = logging.getLogger(__name__)

ALTERNATIVES = ("abs", "greater", "less", "two-sided")
METHODS = ("zscore", "resample")

# Cap on the number of random indices held in memory per batch.
_MAX_BATCH_ELEMENTS = 4_000_000


@dataclass(frozen=True)
class EnrichmentConfig:
    """Parameters of the enrichment test.

    method:        ``zscore`` (closed form, correlation-adjusted) or ``resample``.
    inter_gene_cor: assumed average correlation between genes in a set
                   (``zscore`` only; cameraPR default 0.01; 0 = independent).

    alternative:
        ``abs``       mean |score| of the set is larger than chance
                      (the original doggo "2-tailed" test; direction-agnostic
                      "is this set perturbed at all").
        ``greater``   mean score of the set is larger than chance.
        ``less``      mean score of the set is smaller than chance.
        ``two-sided`` mean score of the set is further from the overall mean
                      than chance, in either direction.
    max_reps:      maximum number of random samples per set.
    min_hits:      stop early after this many extreme samples (0 disables).
    replace:       sample rows with replacement (original behaviour) or not.
    min_size:      skip sets with fewer matched rows than this.
    max_size:      skip sets with more matched rows than this (None = no cap).
    correction:    any ``statsmodels.stats.multitest.multipletests`` method.
    alpha:         family-wise / FDR level used for the ``significant`` column.
    """

    alternative: str = "abs"
    method: str = "zscore"
    inter_gene_cor: float = 0.01
    max_reps: int = 1000
    min_hits: int = 21
    replace: bool = False
    min_size: int = 3
    max_size: int | None = None
    correction: str = "fdr_bh"
    alpha: float = 0.05
    batch_size: int = 1000

    def __post_init__(self) -> None:
        if self.alternative not in ALTERNATIVES:
            raise ValueError(f"alternative must be one of {ALTERNATIVES}, got {self.alternative!r}")
        if self.method not in METHODS:
            raise ValueError(f"method must be one of {METHODS}, got {self.method!r}")
        if not -1.0 < self.inter_gene_cor < 1.0:
            raise ValueError("inter_gene_cor must be between -1 and 1")
        if self.max_reps < 1:
            raise ValueError("max_reps must be >= 1")
        if self.min_hits < 0:
            raise ValueError("min_hits must be >= 0")
        if self.min_size < 1:
            raise ValueError("min_size must be >= 1")
        if self.max_size is not None and self.max_size < self.min_size:
            raise ValueError("max_size must be >= min_size")
        if self.batch_size < 1:
            raise ValueError("batch_size must be >= 1")


def _transform(scores: np.ndarray, alternative: str) -> np.ndarray:
    """Per-row values whose sample mean is the test statistic (before centering)."""
    if alternative == "abs":
        return np.abs(scores)
    if alternative == "less":
        return -scores
    return scores  # greater, two-sided


def _statistic(sample_means: np.ndarray, alternative: str, pop_mean: float) -> np.ndarray:
    if alternative == "two-sided":
        return np.abs(sample_means - pop_mean)
    return sample_means


def _draw(rng: np.random.Generator, n_total: int, n: int, size: int, replace: bool) -> np.ndarray:
    if replace or n == n_total:
        if not replace:  # sampling all rows without replacement: every draw is the full set
            return np.broadcast_to(np.arange(n_total), (size, n_total))
        return rng.integers(0, n_total, size=(size, n))
    # Without replacement: the n smallest of N uniform keys form a uniform random subset.
    keys = rng.random((size, n_total))
    return np.argpartition(keys, n - 1, axis=1)[:, :n]


def resample_test(
    scores: np.ndarray,
    member_idx: np.ndarray,
    config: EnrichmentConfig,
    rng: np.random.Generator,
) -> dict:
    """Run the resampling test for one set.

    ``scores`` are all rows' scores (the background); ``member_idx`` indexes
    the rows that belong to the set.
    """
    scores = np.asarray(scores, dtype=float)
    n_total = scores.shape[0]
    n = int(len(member_idx))
    if n == 0 or n > n_total:
        raise ValueError("member_idx must be non-empty and no larger than scores")

    alt = config.alternative
    values = _transform(scores, alt)
    pop_mean = float(values.mean())
    observed = float(_statistic(np.array([values[member_idx].mean()]), alt, pop_mean)[0])

    per_batch_cap = max(1, _MAX_BATCH_ELEMENTS // (n if config.replace else n_total))
    batch = min(config.batch_size, per_batch_cap)

    hits = 0
    draws = 0
    stopped_early = False
    while draws < config.max_reps:
        size = min(batch, config.max_reps - draws)
        idx = _draw(rng, n_total, n, size, config.replace)
        null = _statistic(values[idx].mean(axis=1), alt, pop_mean)
        extreme = null >= observed - 1e-12 * max(1.0, abs(observed))
        if config.min_hits > 0:
            cum = hits + np.cumsum(extreme)
            reached = np.flatnonzero(cum >= config.min_hits)
            if reached.size:
                draws += int(reached[0]) + 1
                hits = config.min_hits
                stopped_early = True
                break
        hits += int(extreme.sum())
        draws += size

    if stopped_early:
        p = hits / draws
    else:
        p = (hits + 1) / (draws + 1)

    return {"statistic": observed, "n_draws": draws, "n_extreme": hits, "p_value": float(p)}


def zscore_test(
    scores: np.ndarray,
    member_idx: np.ndarray,
    alternative: str = "abs",
    inter_gene_cor: float = 0.01,
) -> dict:
    """Closed-form, correlation-adjusted test of one set (cameraPR formula).

    Returns the moderated t statistic (positive = set mean above the rest,
    for the transformed values) and its p-value.
    """
    values = _transform(np.asarray(scores, dtype=float), alternative)
    N = values.shape[0]
    n = int(len(member_idx))
    if n == 0 or n >= N - 1:
        return {"statistic": np.nan, "p_value": np.nan}
    mean_all = values.mean()
    var_all = values.var(ddof=1)
    delta = N / (N - n) * (values[member_idx].mean() - mean_all)  # set mean minus mean of the rest
    var_pooled = ((N - 1) * var_all - delta**2 * n * (N - n) / N) / (N - 2)
    vif = 1.0 + (n - 1) * inter_gene_cor
    denom = var_pooled * (vif / n + 1.0 / (N - n))
    if denom <= 0:
        return {"statistic": np.nan, "p_value": np.nan}
    t = delta / np.sqrt(denom)
    df = N - 2
    if alternative == "two-sided":
        p = 2.0 * stats.t.sf(abs(t), df)
    else:  # abs / greater / less are upper-tail tests on the transformed values
        p = stats.t.sf(t, df)
    return {"statistic": float(t), "p_value": float(min(1.0, p))}


def _match(genes: pd.Series, gene_sets: dict[str, list[str]], ignore_case: bool):
    if pd.api.types.is_float_dtype(genes) and (genes.dropna() % 1 == 0).all():
        genes = genes.astype("Int64")  # e.g. Entrez IDs read as floats: 1234.0 -> "1234"
    keys = genes.astype(str).str.strip()
    if ignore_case:
        keys = keys.str.upper()
    positions: dict[str, list[int]] = {}
    for i, g in enumerate(keys.to_numpy()):
        positions.setdefault(g, []).append(i)
    out = {}
    for term, members in gene_sets.items():
        idx: list[int] = []
        for g in dict.fromkeys(m.upper() if ignore_case else m for m in members):
            idx.extend(positions.get(g, ()))
        out[term] = (len(members), np.asarray(sorted(idx), dtype=np.intp))
    return out


_SHARED: dict = {}


def _init_worker(scores: np.ndarray, config: EnrichmentConfig) -> None:
    _SHARED["scores"] = scores
    _SHARED["config"] = config


def _describe(scores: np.ndarray, idx: np.ndarray) -> dict:
    s = scores[idx]
    return {
        "n_matched": len(idx),
        "set_mean": float(s.mean()),
        "set_mean_abs": float(np.abs(s).mean()),
        "frac_up": float((s > 0).mean()),
        "background_mean": float(scores.mean()),
        "background_mean_abs": float(np.abs(scores).mean()),
        "effect": float(s.mean() - scores.mean()),
    }


def _worker(args):
    term, n_in_set, idx, seed = args
    scores, config = _SHARED["scores"], _SHARED["config"]
    res = resample_test(scores, idx, config, np.random.default_rng(seed))
    return {"term": term, "set_size": n_in_set, **_describe(scores, idx), **res}


DESCRIPTIVE_COLUMNS = [
    "term", "set_size", "n_matched", "set_mean", "set_mean_abs", "frac_up",
    "background_mean", "background_mean_abs", "effect",
]
RESULT_COLUMNS = {
    "zscore": DESCRIPTIVE_COLUMNS + ["statistic", "p_value", "p_adj", "significant"],
    "resample": DESCRIPTIVE_COLUMNS + ["statistic", "n_draws", "n_extreme", "p_value", "p_adj", "significant"],
}


def run_enrichment(
    data: pd.DataFrame,
    gene_sets: dict[str, list[str]],
    gene_col: str,
    score_col: str,
    config: EnrichmentConfig = EnrichmentConfig(),
    seed: int | None = None,
    jobs: int = 1,
    ignore_case: bool = False,
) -> pd.DataFrame:
    """Test every gene set against one table of per-gene scores.

    Rows with a missing gene or score are dropped. Each row is one draw
    unit, so a gene present on several rows counts several times (as in the
    original notebooks). Results are reproducible for a given ``seed``
    regardless of ``jobs``.
    """
    for col in (gene_col, score_col):
        if col not in data.columns:
            raise KeyError(f"column {col!r} not in input; available: {list(data.columns)}")

    df = data[[gene_col, score_col]].copy()
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    before = len(df)
    df = df.dropna().reset_index(drop=True)
    if len(df) < before:
        log.warning("dropped %d rows with missing gene or non-numeric score", before - len(df))
    if df.empty:
        raise ValueError("no usable rows after dropping missing values")

    scores = df[score_col].to_numpy(dtype=float)
    matched = _match(df[gene_col], gene_sets, ignore_case)

    tasks = []
    skipped = 0
    for term, (n_in_set, idx) in matched.items():
        n = len(idx)
        if n < config.min_size or (config.max_size is not None and n > config.max_size) or n > len(scores):
            skipped += 1
            continue
        tasks.append((term, n_in_set, idx))
    log.info("testing %d sets (%d skipped by size filters) against %d rows", len(tasks), skipped, len(scores))
    if not tasks:
        n_hit = sum(1 for _, idx in matched.values() if len(idx))
        log.warning(
            "no gene set had >= %d matching rows (%d sets matched any row). Check --gene-col, "
            "gene identifier type/species, or try --ignore-case.", config.min_size, n_hit,
        )

    if config.method == "zscore":
        rows = [
            {"term": t, "set_size": s, **_describe(scores, i),
             **zscore_test(scores, i, config.alternative, config.inter_gene_cor)}
            for t, s, i in tasks
        ]
        return _finish(rows, config)

    seeds = np.random.SeedSequence(seed).spawn(len(tasks))
    payload = [(t, s, i, ss) for (t, s, i), ss in zip(tasks, seeds)]

    if jobs == 0:
        jobs = os.cpu_count() or 1
    if jobs > 1 and len(payload) > 1:
        with ProcessPoolExecutor(max_workers=jobs, initializer=_init_worker, initargs=(scores, config)) as ex:
            rows = list(ex.map(_worker, payload, chunksize=max(1, len(payload) // (jobs * 8))))
    else:
        _init_worker(scores, config)
        rows = []
        for k, p in enumerate(payload, 1):
            rows.append(_worker(p))
            if k % 250 == 0:
                log.info("  %d / %d sets", k, len(payload))

    res = _finish(rows, config)
    floor = 1.0 / (config.max_reps + 1)
    n_floor = int((res["p_value"] <= floor * (1 + 1e-9)).sum())
    if n_floor:
        log.warning(
            "%d of %d sets reached the resolution limit p = 1/(max_reps+1) = %.2g; "
            "increase --max-reps or use --method zscore.",
            n_floor, len(res), floor,
        )
    return res


def _finish(rows: list[dict], config: EnrichmentConfig) -> pd.DataFrame:
    res = pd.DataFrame(rows, columns=RESULT_COLUMNS[config.method][:-2])
    if res.empty:
        res["p_adj"] = pd.Series(dtype=float)
        res["significant"] = pd.Series(dtype=bool)
        return res
    p = res["p_value"].to_numpy(dtype=float)
    ok = ~np.isnan(p)
    p_adj = np.full_like(p, np.nan)
    reject = np.zeros(len(p), dtype=bool)
    if ok.any():
        reject[ok], p_adj[ok], _, _ = multipletests(p[ok], alpha=config.alpha, method=config.correction)
    res["p_adj"] = p_adj
    res["significant"] = reject
    return res.sort_values(["p_value", "term"], kind="mergesort").reset_index(drop=True)
