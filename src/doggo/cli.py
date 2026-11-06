"""Command-line interface: ``doggo``."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

from doggo import __version__
from doggo.enrichment import ALTERNATIVES, METHODS, EnrichmentConfig, run_enrichment
from doggo.genesets import filter_gene_sets, read_gmt
from doggo import msigdb

log = logging.getLogger("doggo")


def _read_table(path: str, sep: str | None) -> pd.DataFrame:
    if sep is None:
        suffix = "".join(Path(path.removesuffix(".gz")).suffixes[-1:]).lower()
        sep = "\t" if suffix in (".tsv", ".txt", ".tab") else ","
    return pd.read_csv(path, sep=sep)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="doggo",
        description="Resampling-based gene set enrichment for any per-gene score "
        "(log2FC, model parameters, statistics, ...).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    io = p.add_argument_group("input / output")
    io.add_argument("-i", "--input", required=True, help="table with one row per gene (CSV/TSV, optionally .gz)")
    src = io.add_mutually_exclusive_group(required=True)
    src.add_argument("-g", "--gmt", help="gene set collection in GMT format (any source)")
    src.add_argument("--org", "-org", dest="org",
                     help="use the full MSigDB collection for this organism: human or mouse "
                          "(downloaded once and cached; aliases hs, mm)")
    io.add_argument("-o", "--output", required=True, help="output CSV path")
    io.add_argument("--sep", default=None, help="input delimiter (default: tab for .tsv/.txt/.tab, else comma)")
    io.add_argument("--gene-col", required=True, help="column holding gene identifiers matching the GMT")
    io.add_argument("--score-col", required=True, help="numeric column to test")
    io.add_argument("--group-col", default=None,
                    help="run independently for each value of this column (e.g. cluster); "
                         "multiple-testing correction is applied within each group")
    io.add_argument("--filter", action="append", default=[], metavar="COL=VALUE",
                    help="keep only rows where COL equals VALUE (repeatable), e.g. --filter Param=Max")
    io.add_argument("--ignore-case", action="store_true", help="match gene symbols case-insensitively")

    gs = p.add_argument_group("gene sets")
    gs.add_argument("--msigdb-version", default=msigdb.DEFAULT_VERSION,
                    help="MSigDB release used with --org, e.g. 2026.2 or 2022.1")
    gs.add_argument("--id-type", choices=msigdb.ID_TYPES, default="symbols",
                    help="gene identifiers used with --org: symbols or entrez (NCBI Gene IDs)")
    gs.add_argument("--prefix", action="append", default=[],
                    help="only test sets whose name starts with this (repeatable), e.g. REACTOME_ or GOBP_")
    gs.add_argument("--term", action="append", default=[], help="test this exact set name (repeatable)")
    gs.add_argument("--min-size", type=int, default=3, help="minimum number of matched rows to test a set")
    gs.add_argument("--max-size", type=int, default=None, help="maximum number of matched rows to test a set")

    st = p.add_argument_group("statistics")
    st.add_argument("--alternative", choices=ALTERNATIVES, default="abs",
                    help="abs: mean |score| larger than chance (original doggo 2-tailed test); "
                         "greater / less: mean score shifted up / down; "
                         "two-sided: mean score shifted either way")
    st.add_argument("--method", choices=METHODS, default="zscore",
                    help="zscore: closed-form test with correlation adjustment (fast, no p floor); "
                         "resample: Monte Carlo test of the original doggo (assumes independent genes)")
    st.add_argument("--inter-gene-cor", type=float, default=0.01,
                    help="zscore only: assumed average correlation between genes in a set "
                         "(cameraPR default; 0 assumes independence and gives smaller p-values)")
    st.add_argument("--max-reps", type=int, default=1000, help="resample only: maximum random samples per set")
    st.add_argument("--min-hits", type=int, default=21,
                    help="resample only: stop sampling a set after this many extreme draws (0 = always run max-reps)")
    st.add_argument("--replace", action="store_true",
                    help="resample only: sample background rows with replacement (as the v0.3 notebooks did)")
    st.add_argument("--correction", default="fdr_bh",
                    help="statsmodels multipletests method: fdr_bh, bonferroni, holm, hs, fdr_by, ...")
    st.add_argument("--alpha", type=float, default=0.05, help="level for the 'significant' column")

    run = p.add_argument_group("execution")
    run.add_argument("--seed", type=int, default=None, help="random seed for reproducible results")
    run.add_argument("-j", "--jobs", type=int, default=1, help="worker processes (0 = all CPUs)")
    run.add_argument("-q", "--quiet", action="store_true", help="only print warnings and errors")
    return p


def _apply_filters(df: pd.DataFrame, filters: list[str]) -> pd.DataFrame:
    for f in filters:
        if "=" not in f:
            raise SystemExit(f"--filter expects COL=VALUE, got {f!r}")
        col, val = f.split("=", 1)
        if col not in df.columns:
            raise SystemExit(f"--filter column {col!r} not in input")
        df = df[df[col].astype(str) == val]
    return df


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.WARNING if args.quiet else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stderr,
    )

    try:
        config = EnrichmentConfig(
            alternative=args.alternative,
            method=args.method,
            inter_gene_cor=args.inter_gene_cor,
            max_reps=args.max_reps,
            min_hits=args.min_hits,
            replace=args.replace,
            min_size=args.min_size,
            max_size=args.max_size,
            correction=args.correction,
            alpha=args.alpha,
        )
    except ValueError as e:
        raise SystemExit(f"error: {e}")

    if args.org:
        try:
            gmt = msigdb.fetch_gmt(args.org, args.msigdb_version, args.id_type)
        except (ValueError, RuntimeError) as e:
            raise SystemExit(f"error: {e}")
    else:
        gmt = args.gmt
    sets = filter_gene_sets(read_gmt(gmt), args.prefix, args.term)
    log.info("loaded %d gene sets from %s", len(sets), gmt)
    if not sets:
        raise SystemExit("error: no gene sets left after --prefix/--term filtering")

    data = _apply_filters(_read_table(args.input, args.sep), args.filter)
    log.info("loaded %d rows from %s", len(data), args.input)

    common = dict(gene_col=args.gene_col, score_col=args.score_col, config=config,
                  jobs=args.jobs, ignore_case=args.ignore_case)
    try:
        if args.group_col is None:
            result = run_enrichment(data, sets, seed=args.seed, **common)
        else:
            if args.group_col not in data.columns:
                raise KeyError(f"--group-col {args.group_col!r} not in input")
            parts = []
            for k, (group, sub) in enumerate(data.groupby(args.group_col, sort=True)):
                log.info("group %s=%s (%d rows)", args.group_col, group, len(sub))
                seed = None if args.seed is None else args.seed + k
                r = run_enrichment(sub, sets, seed=seed, **common)
                r.insert(0, args.group_col, group)
                parts.append(r)
            result = pd.concat(parts, ignore_index=True)
    except (KeyError, ValueError) as e:
        raise SystemExit(f"error: {e.args[0] if e.args else e}")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out, index=False)
    log.info("wrote %d rows to %s (%d significant)", len(result), out, int(result["significant"].sum()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
