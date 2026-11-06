# doggo

**doggo finds the biological pathways that change in your experiment.**

Give it a table with one number per gene (for example, a fold change from a
differential expression analysis). It checks hundreds or thousands of known gene
groups ("gene sets", such as Reactome pathways or GO terms) and tells you
which groups changed more than you'd expect by chance.

It works with any per-gene number: bulk RNA-seq fold changes, single-cell
cluster markers, model parameters, or test statistics. It covers human and
mouse.

---

## Contents

1. [Quick start](#1-quick-start)
2. [Installation](#2-installation)
3. [Preparing your table](#3-preparing-your-table)
4. [Running doggo](#4-running-doggo)
5. [Choosing gene sets](#5-choosing-gene-sets)
6. [Reading the results](#6-reading-the-results)
7. [Choosing settings](#7-choosing-settings)
8. [All options](#8-all-options)
9. [Troubleshooting](#9-troubleshooting)
10. [How it works](#10-how-it-works)
11. [Using doggo from Python](#11-using-doggo-from-python)
12. [Changes from the original notebooks](#12-changes-from-the-original-notebooks)

---

## 1. Quick start

```bash
# one-time setup (see section 2 for details)
git clone https://github.com/qveu/doggo.git
cd doggo
conda env create -f environment.yml
conda activate doggo
pip install --no-deps -e .

# run it
doggo -i my_results.csv -org mouse -o pathways.csv \
      --gene-col gene --score-col log2FC --prefix REACTOME_
```

This tests every Reactome pathway against the `log2FC` column of
`my_results.csv` and writes the results to `pathways.csv`, with the most
significant pathways first. The first run downloads the mouse gene sets
(about 16 MB). Later runs reuse that copy.

---

## 2. Installation

You need either **conda** (Anaconda, Miniconda, Mamba or Micromamba) or a
plain **Python 3.10+**. Both routes install exact, tested versions of every
package, so results don't drift when packages update.

### Option A: conda (recommended)

```bash
git clone https://github.com/qveu/doggo.git
cd doggo
conda env create -f environment.yml
conda activate doggo
pip install --no-deps -e .
```

Check that it worked:

```bash
doggo --version
```

Each time you open a new terminal, run `conda activate doggo` before using doggo.

### Option B: Python virtual environment

```bash
git clone https://github.com/qveu/doggo.git
cd doggo
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements-lock.txt
pip install --no-deps -e .
```

### Optional: check the installation

```bash
pytest
```

All tests should pass.

<details>
<summary>For maintainers: freezing the conda environment completely</summary>

`environment.yml` pins every package doggo uses directly. To also pin
every package those depend on, for exact rebuilds years later, create a
lock file with [conda-lock](https://github.com/conda/conda-lock) and commit it:

```bash
conda-lock lock -f environment.yml -p linux-64 -p osx-64 -p osx-arm64 -p win-64
conda-lock install -n doggo conda-lock.yml
```
</details>

---

## 3. Preparing your table

doggo reads a **CSV** (comma-separated) or **TSV** (tab-separated) file. It
can be gzipped (`.gz`). You need:

| | What it is | Example |
|---|---|---|
| **A gene column** | Gene names (symbols, like `Cd14` or `CD14`) or NCBI Gene IDs | `gene` |
| **A score column** | One number per gene | `log2FC` |

Other columns are fine; doggo ignores them unless you ask it to use them.

```
gene,log2FC,pvalue
Cd14,2.31,1e-12
Lyz2,1.85,3e-9
Actb,0.02,0.91
...
```

**Include all the genes you measured, not just the significant ones.** doggo
compares each pathway with the rest of your table, so the table must
represent the whole background.

**Which score to use:** a plain fold change works. If you have one, a
*shrunken* fold change (for example DESeq2's `lfcShrink`) or a test
statistic is better: raw fold changes of rarely detected genes are very
noisy and can make random pathways look changed (see
[section 7](#which-score-column)).

**Single-cell or multi-condition tables** can hold several results in one
file, e.g. one block of rows per cluster. Use `--group-col` to analyze each
block separately, or `--filter` to keep just one (see the next section).

---

## 4. Running doggo

The basic pattern:

```bash
doggo -i INPUT -org ORGANISM -o OUTPUT --gene-col GENE_COLUMN --score-col SCORE_COLUMN
```

### Example: bulk RNA-seq results

```bash
doggo -i deseq_results.csv -org human -o hallmark.csv \
      --gene-col symbol --score-col log2FoldChange --prefix HALLMARK_
```

### Example: single-cell markers, one analysis per cluster

A Seurat `FindAllMarkers` table has a `cluster` column. `--group-col`
analyzes each cluster on its own and stacks the results in one file:

```bash
doggo -i 10x_doggo/Myeloid_cells_updated_0.4res.all_markers_genes.csv -org mouse \
      -o cluster_pathways.csv --gene-col gene --score-col avg_log2FC \
      --group-col cluster --prefix REACTOME_
```

### Example: keep only some rows

`--filter COLUMN=VALUE` keeps matching rows only. This table has two rows
per gene (`Param` is `MP` or `Max`), and its gene names are in capitals while
the mouse gene sets use `Cd14`-style names, so `--ignore-case` is needed:

```bash
doggo -i allParameters_PV-CVB_22920.csv -org mouse -o pv_max.csv \
      --gene-col target --score-col SigFactor_PV --filter Param=Max \
      --ignore-case --prefix REACTOME_ --prefix GOBP_
```

Progress messages print as it runs. Add `-q` to hide them. Run `doggo --help`
for the full list of options.

---

## 5. Choosing gene sets

### Built-in: MSigDB for human or mouse

`-org human` or `-org mouse` uses the full
[MSigDB](https://www.gsea-msigdb.org/gsea/msigdb) collection (release
2026.2), the most widely used gene set database. MSigDB only covers these two
organisms.

The collection is large, so pick the families you care about with `--prefix`
(you can repeat it):

| Prefix | What it contains | Human | Mouse |
|---|---|---:|---:|
| `HALLMARK_` | 50 well-defined core processes; a good first look | 50 | 50 |
| `REACTOME_` | Curated pathways | 1,867 | 1,350 |
| `KEGG_` | KEGG pathways | 976 | – |
| `WP_` | WikiPathways | 927 | 200 |
| `BIOCARTA_` | BioCarta pathways | 292 | 252 |
| `GOBP_` | Gene Ontology: biological processes | 6,905 | 7,738 |
| `GOMF_` | Gene Ontology: molecular functions | 1,918 | 1,945 |
| `GOCC_` | Gene Ontology: cellular components | 1,094 | 1,088 |
| `HP_` | Human phenotypes | 5,932 | – |
| `MP_` | Mouse phenotypes | – | 6,757 |

To test specific sets by name, use `--term REACTOME_APOPTOSIS` (repeatable).

**Gene names vs IDs:** by default doggo matches gene symbols. If your gene
column holds NCBI Gene IDs (numbers), add `--id-type entrez`.

**Older releases:** to reproduce a past analysis, add e.g.
`--msigdb-version 2022.1`.

### Downloading ahead of time

Gene sets download automatically the first time you use `-org`. To download in
advance, for example before working on a computer without internet access:

```bash
doggo-fetch                                     # human and mouse
doggo-fetch --org mouse --id-type symbols entrez
```

Files are stored in `~/.cache/doggo`. Set the `DOGGO_CACHE` environment
variable to use a different folder, e.g. a shared lab folder. Downloads of the
default release are checked against known fingerprints, so a damaged
download is caught. The files are not stored in this repository.

### Your own gene sets

Use `-g my_sets.gmt` instead of `-org`. GMT is a simple text format, one set
per line, with tabs between fields:

```
MY_SET_NAME<TAB>description<TAB>Gene1<TAB>Gene2<TAB>Gene3 ...
```

Single collections can also be downloaded by hand from
<https://data.broadinstitute.org/gsea-msigdb/msigdb/release/>.

---

## 6. Reading the results

The output is a CSV, one row per gene set, sorted from most to least
significant.

**Start with these columns:**

| Column | Meaning |
|---|---|
| `term` | Gene set name |
| `p_adj` | **Main result.** Adjusted p-value (accounts for testing many sets at once). Below 0.05 is the usual cut-off. |
| `significant` | `True` if `p_adj` is below `--alpha` (0.05 by default) |
| `n_matched` | How many of the set's genes were found in your table |
| `frac_up` | Fraction of those genes with a positive score. Near 1 = mostly up, near 0 = mostly down, around 0.5 = mixed |
| `set_mean_abs` vs `background_mean_abs` | Average size of change in the set vs in all genes |

**Other columns:**

| Column | Meaning |
|---|---|
| `set_size` | Number of genes in the set's definition |
| `set_mean`, `background_mean` | Average score (with sign) in the set and in all genes |
| `effect` | `set_mean` − `background_mean` |
| `statistic` | The test score. Larger = stronger evidence |
| `p_value` | p-value before adjusting for multiple tests. Don't use it to pick hits. |
| `n_draws`, `n_extreme` | Only with `--method resample`: random samples drawn, and how many were at least as extreme |
| *group column* | Only with `--group-col`: which group the row belongs to |

**Why use `p_adj` and not `p_value`?** When you test 1,000 gene sets, about
50 will have `p_value` < 0.05 by pure chance, even when nothing is going on.
`p_adj` corrects for this.

**Interpreting hits:** a significant set with `frac_up` near 0.5 has genes
changing in *both* directions. Look at the individual genes to see what is
going on.

**Small sets:** sets matching only a few genes (`n_matched` of 3–5) can be
driven by one or two genes. Consider `--min-size 10` for a more robust list.

---

## 7. Choosing settings

The defaults suit most analyses. These are the choices that matter most.

### What counts as "changed"? (`--alternative`)

| Setting | Asks | Use when |
|---|---|---|
| `abs` *(default)* | Are the set's genes changing more than usual, in **either** direction? | You want any perturbed pathway, including those with some genes up and others down |
| `greater` | Are the set's genes shifted **up**? | You only care about activation |
| `less` | Are the set's genes shifted **down**? | You only care about repression |
| `two-sided` | Are the set's genes shifted consistently up **or** consistently down? | You want coordinated changes and their direction |

`abs` catches the most kinds of change, but it is less sensitive to a
modest, consistent shift in one direction. Use `greater`, `less` or
`two-sided` if that's what you expect.

### Which score column?

With `abs`, a gene counts as "changed" if its score is far from zero for
any reason, including noise. Raw fold changes of rarely detected genes are
noisy, so pathways full of such genes can look changed when they aren't.
Prefer, in order:

1. a shrunken fold change (DESeq2 `lfcShrink`, edgeR `predFC`, or similar),
2. a test statistic (e.g. DESeq2's `stat`, limma's `t`),
3. a plain fold change.

### How doggo computes p-values (`--method`)

| Setting | Notes |
|---|---|
| `zscore` *(default)* | Uses a formula. Takes seconds even for thousands of sets, and reports any p-value however small. Allows for genes in a pathway tending to change together (see next setting). |
| `resample` | The original doggo approach: compares each set with thousands of random gene sets of the same size. Makes no assumptions about the shape of your data, but is much slower and can't produce p-values smaller than 1/(`--max-reps`+1). Use it to double-check results or for very unusual score distributions. |

On the example single-cell data, the two methods give closely matching
p-values. `zscore` ran all 8 clusters in about 3 seconds; `resample` with
20,000 draws per set took about 6 minutes for a single cluster.

### Genes that change together (`--inter-gene-cor`)

Genes in a pathway are often controlled together, so they rise and fall
together. A group of such genes swings more than a random group, and
ignoring this makes pathways look more significant than they are. doggo
allows for a small typical amount of this by default (0.01, the value used
by the widely used CAMERA method). Use `0` to turn it off (more hits, but
less trustworthy) or a larger value like `0.05` to be more cautious.

### Multiple-testing correction (`--correction`)

The default, `fdr_bh` (Benjamini–Hochberg), controls the share of false hits
among your significant sets, and is standard for this kind of screen. Use
`bonferroni` or `holm` for a stricter rule.

---

## 8. All options

**Input and output**

| Option | Default | Description |
|---|---|---|
| `-i`, `--input` | *required* | Your table (CSV/TSV, optionally `.gz`) |
| `-o`, `--output` | *required* | Where to write results (CSV) |
| `--gene-col` | *required* | Column with gene names or IDs |
| `--score-col` | *required* | Column with the numbers to test |
| `--group-col` | – | Analyze each value of this column separately (e.g. `cluster`) |
| `--filter COL=VALUE` | – | Keep only rows where `COL` equals `VALUE`. Repeatable. |
| `--ignore-case` | off | Match gene names regardless of upper/lower case |
| `--sep` | auto | Column separator. Tab for `.tsv`/`.txt`, comma otherwise |

**Gene sets** (use either `-org` or `-g`)

| Option | Default | Description |
|---|---|---|
| `-org`, `--org` | – | `human` or `mouse`: use MSigDB |
| `-g`, `--gmt` | – | Use your own GMT file |
| `--prefix` | all | Only sets whose name starts with this, e.g. `REACTOME_`. Repeatable. |
| `--term` | – | Only this exact set. Repeatable. |
| `--min-size` / `--max-size` | 3 / none | Skip sets with fewer/more matched genes |
| `--msigdb-version` | 2026.2 | MSigDB release |
| `--id-type` | `symbols` | `symbols` or `entrez` (NCBI Gene IDs) |

**Statistics**

| Option | Default | Description |
|---|---|---|
| `--alternative` | `abs` | `abs`, `greater`, `less` or `two-sided` (section 7) |
| `--method` | `zscore` | `zscore` or `resample` (section 7) |
| `--inter-gene-cor` | 0.01 | How much genes in a set change together (`zscore` only) |
| `--correction` | `fdr_bh` | Multiple-testing correction: `fdr_bh`, `fdr_by`, `bonferroni`, `holm`, `hs` |
| `--alpha` | 0.05 | Cut-off for the `significant` column |

**Resampling only** (`--method resample`)

| Option | Default | Description |
|---|---|---|
| `--max-reps` | 1000 | Maximum random gene sets drawn per set |
| `--min-hits` | 21 | Stop early once this many random sets beat the real one |
| `--replace` | off | Allow the same gene twice in a random draw (old behavior) |
| `--seed` | – | Fix the random numbers so reruns give identical results |
| `-j`, `--jobs` | 1 | CPU cores to use (`0` = all) |

**Other:** `-q` (quiet), `--version`, `--help`.

---

## 9. Troubleshooting

**"no gene set had >= 3 matching rows"**: doggo didn't recognize your gene
names. Check that:
- `--gene-col` points at the gene names, not e.g. a row-number column;
- the organism matches (`-org mouse` for mouse genes);
- capitalization matches. Mouse sets use `Cd14`, human sets use `CD14`. Add `--ignore-case`;
- you are using symbols vs IDs correctly (`--id-type entrez` for numeric IDs).

**"dropped N rows with missing gene or non-numeric score"**: some rows had
empty cells or text such as `NA` in the score column. Those rows are skipped.
This is usually harmless.

**"could not download …"**: the computer can't reach the Broad Institute
server. Run `doggo-fetch` on a computer with internet, copy `~/.cache/doggo`
across, or download the GMT file by hand and use `-g`.

**"N sets reached the resolution limit"** (resample method): many p-values
hit the smallest value the chosen number of draws allows. Use the default
`--method zscore`, or raise `--max-reps` (around 20 × the number of sets).

**Nothing is significant**: try a smaller collection (`--prefix HALLMARK_`),
check that your table includes all measured genes, and check that the score
column is the one you meant.

**`doggo: command not found`**: activate the environment first
(`conda activate doggo` or `source .venv/bin/activate`).

---

## 10. How it works

For each gene set, doggo finds the set's genes in your table, takes the
average of their scores (or of the size of their scores, for `abs`), and asks
how unusual that average is compared with the rest of the genes:

- **`zscore` method:** calculates how much the average of a randomly chosen
  group of that size would normally vary, allows for genes in a set moving
  together, and converts the difference into a p-value. This is the method
  of `cameraPR` in the limma R package (Wu & Smyth, 2012, *Nucleic Acids
  Research* 40:e133).
- **`resample` method:** draws random groups of genes of the same size,
  computes the same average for each, and reports the fraction that are at
  least as extreme as the real set. Sampling stops early once the answer is
  clearly unremarkable (Besag & Clifford, 1991, *Biometrika* 78:301).

The p-values for all sets are then adjusted for multiple testing.

This is a *competitive* test: it asks whether a pathway changes **more than
other genes** in the same experiment, not whether it changes at all.

---

## 11. Using doggo from Python

```python
import pandas as pd
from doggo import EnrichmentConfig, filter_gene_sets, read_gmt, run_enrichment
from doggo.msigdb import fetch_gmt

sets = filter_gene_sets(read_gmt(fetch_gmt("mouse")), prefixes=["REACTOME_"])
table = pd.read_csv("my_results.csv")

results = run_enrichment(
    table, sets, gene_col="gene", score_col="log2FC",
    config=EnrichmentConfig(alternative="abs", method="zscore"),
)
print(results.head())
```

---

## 12. Changes from the original notebooks

The original notebooks (`Doggo_v0.*.ipynb`, `10x_doggo/`) are kept for
reference. The command-line tool differs from them in these ways:

- **New default method.** `zscore` replaces resampling as the default (faster,
  no p-value floor, allows for genes changing together). `--method resample
  --inter-gene-cor 0` is the closest match to the old behavior.
- **Fixed per-cluster bug.** In `10x_doggo/Doggo_v0.3_10x.ipynb`, each cluster was
  meant to be analyzed separately, but the code used the full table every time.
  `--group-col` analyzes each group correctly. The per-cluster CSVs in
  `10x_doggo/` were affected.
- **p-values are never exactly 0.** With resampling, sets that never met a
  random set as extreme now get 1/(draws+1) instead of 0.
- **Too-small sets are left out** rather than given p = 1, so they no longer
  affect the multiple-testing correction.
- **Default correction is Benjamini–Hochberg.** The notebooks used Holm–Šidák
  (`--correction hs`).
- **Resampling draws without repeats** by default (`--replace` for the old way).
- **The old "1-tailed" mode** picked its direction from the data. Use
  `--alternative greater`, `less` or `two-sided` instead.
- **Gene sets are downloaded on demand** (`-org`) instead of being read from
  pickle files with hard-coded paths.

---

## License

MIT. See [LICENSE](LICENSE).
