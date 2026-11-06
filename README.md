# doggo

Resampling-based gene set enrichment for **any per-gene score**: log2 fold changes,
fitted model parameters, test statistics, etc. For every gene set, doggo compares
the scores of the genes in the set with random gene samples of the same size from
your table, and reports an empirical p-value plus a multiple-testing correction.

## Install

### Conda (recommended)

```bash
git clone https://github.com/qveu/doggo.git && cd doggo
conda env create -f environment.yml     # mamba / micromamba work too
conda activate doggo
pip install --no-deps -e .
doggo --help
```

`environment.yml` pins Python and every direct dependency. For an exact rebuild
of every transitive package, generate a lock file once with
[conda-lock](https://github.com/conda/conda-lock) and commit it:

```bash
conda-lock lock -f environment.yml -p linux-64 -p osx-64 -p osx-arm64 -p win-64
conda-lock install -n doggo conda-lock.yml
```

### pip / venv

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements-lock.txt
pip install --no-deps -e .
```

Run the tests with `pytest`.

## Usage

Inputs:

* **A table** (CSV or TSV, optionally gzipped) with one row per gene, containing a
  gene identifier column and a numeric score column. Extra columns are fine.
* **Gene sets** in GMT format, e.g. from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb).
  Identifiers must match your gene column (use mouse sets for mouse symbols).
  `msigdb.v2022.1.Mm.symbols.gmt.txt` (mouse) and `msigdb.v7.1.symbols.gmt` (human) are included.

```bash
# Bulk table, Reactome sets only
doggo -i results.csv -g msigdb.v2022.1.Mm.symbols.gmt.txt -o enrichment.csv \
      --gene-col ext_gene --score-col b --prefix REACTOME_ --seed 1

# 10x / Seurat FindAllMarkers output: one analysis per cluster, 4 CPUs
doggo -i 10x_doggo/Myeloid_cells_updated_0.4res.all_markers_genes.csv \
      -g msigdb.v2022.1.Mm.symbols.gmt.txt -o markers_enrichment.csv \
      --gene-col gene --score-col avg_log2FC --group-col cluster \
      --prefix REACTOME_ --seed 1 -j 4

# Long table with several parameters per gene: pick one, uppercase symbols vs. mouse sets
doggo -i allParameters_PV-CVB_22920.csv -g msigdb.v2022.1.Mm.symbols.gmt.txt -o pv_max.csv \
      --gene-col target --score-col SigFactor_PV --filter Param=Max \
      --ignore-case --prefix REACTOME_ --prefix GOBP_ --max-reps 10000
```

Key options (see `doggo --help` for all):

| Option | Default | Meaning |
|---|---|---|
| `--prefix` / `--term` | all sets | Restrict to set names starting with a prefix (`REACTOME_`, `GOBP_`, `HALLMARK_`, ...) or exact names. Repeatable. |
| `--group-col` | none | Run separately per value (e.g. cluster). Correction is applied within each group. |
| `--filter COL=VALUE` | none | Keep only matching rows. Repeatable. |
| `--alternative` | `abs` | `abs`: mean \|score\| of the set larger than chance (the original doggo test). `greater` / `less`: mean score shifted up / down. `two-sided`: shifted either way. |
| `--max-reps` | 1000 | Maximum random samples per set. The smallest attainable p is 1/(max-reps+1); raise it when testing many sets. |
| `--min-hits` | 21 | Stop sampling a set early once this many random samples are as extreme as the observed one. |
| `--min-size` / `--max-size` | 3 / none | Matched-row limits for testing a set. |
| `--correction` | `fdr_bh` | Any statsmodels `multipletests` method (`bonferroni`, `holm`, `hs`, ...). |
| `--replace` | off | Sample with replacement, as the v0.3 notebooks did. |
| `--seed`, `-j/--jobs` | none, 1 | Results are reproducible for a given seed, whatever the number of jobs. |

Output columns: `term`, `set_size` (genes in the GMT set), `n_matched` (rows in your
table that matched), `set_mean`, `set_mean_abs`, `background_mean`, `effect`
(set mean − background mean), `statistic` (the tested quantity), `n_draws`,
`n_extreme`, `p_value`, `p_adj`, `significant`.

From Python:

```python
import pandas as pd
from doggo import EnrichmentConfig, filter_gene_sets, read_gmt, run_enrichment

sets = filter_gene_sets(read_gmt("msigdb.v2022.1.Mm.symbols.gmt.txt"), ["REACTOME_"])
res = run_enrichment(pd.read_csv("results.csv"), sets, gene_col="ext_gene", score_col="b",
                     config=EnrichmentConfig(max_reps=5000), seed=1, jobs=4)
```

## Method

For a set whose matched rows have scores *x*, the statistic is the mean of
*|x|* (`abs`), *x* (`greater`), −*x* (`less`) or |mean(*x*) − background mean|
(`two-sided`). The same statistic is computed on random samples of
`n_matched` rows from the whole table. Sampling stops after `min-hits` samples
at least as extreme as the observed one (Besag & Clifford, 1991, p = h/l);
otherwise it runs `max-reps` samples and p = (k+1)/(max-reps+1).

### Changes from the v0.3 notebooks

The notebooks are kept for reference. The packaged version differs in a few
deliberate ways:

* Sets that never reach `min-hits` get p = (k+1)/(n+1) instead of k/n, so p is never 0.
* Sets below `min-size` are left out of the output and the correction instead of being given p = 1.
* Default correction is Benjamini–Hochberg (`fdr_bh`); the notebooks used the statsmodels default, Holm–Šidák (`--correction hs`).
* Default sampling is without replacement (`--replace` restores the old behaviour).
* The notebooks' 1-tailed mode chose the direction from the data; use `--alternative greater`, `less` or `two-sided` instead.
* In `10x_doggo/Doggo_v0.3_10x.ipynb` the per-cluster loop passed each cluster to `f()`, but `f()` used the full table, so every cluster was scored against all clusters' rows. `--group-col` actually subsets.
