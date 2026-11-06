import numpy as np
import pandas as pd
import pytest

from doggo import EnrichmentConfig, filter_gene_sets, read_gmt, resample_test, run_enrichment
from doggo.cli import main


@pytest.fixture
def toy(tmp_path):
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(500)]
    scores = rng.normal(0, 1, 500)
    scores[:20] += 3.0  # UP set strongly shifted up
    scores[20:40] -= 3.0  # DOWN set strongly shifted down
    df = pd.DataFrame({"gene": genes, "score": scores})
    gmt = tmp_path / "sets.gmt"
    gmt.write_text(
        "TEST_UP\tna\t" + "\t".join(genes[:20]) + "\n"
        "TEST_DOWN\tna\t" + "\t".join(genes[20:40]) + "\n"
        "TEST_NULL\tna\t" + "\t".join(genes[200:260]) + "\n"
        "TEST_TINY\tna\tG300\tG301\n"
        "OTHER_SET\tna\t" + "\t".join(genes[:20]) + "\n"
    )
    return df, gmt


def test_read_and_filter_gmt(toy):
    _, gmt = toy
    sets = read_gmt(gmt)
    assert set(sets) == {"TEST_UP", "TEST_DOWN", "TEST_NULL", "TEST_TINY", "OTHER_SET"}
    assert len(sets["TEST_UP"]) == 20
    assert set(filter_gene_sets(sets, ["TEST_"])) == {"TEST_UP", "TEST_DOWN", "TEST_NULL", "TEST_TINY"}
    assert set(filter_gene_sets(sets, terms=["OTHER_SET"])) == {"OTHER_SET"}


@pytest.mark.parametrize("replace", [False, True])
def test_alternatives(toy, replace):
    df, gmt = toy
    sets = filter_gene_sets(read_gmt(gmt), ["TEST_"])

    def p(alt):
        cfg = EnrichmentConfig(alternative=alt, replace=replace, max_reps=2000)
        r = run_enrichment(df, sets, "gene", "score", cfg, seed=1).set_index("term")
        assert "TEST_TINY" not in r.index  # below min_size
        return r["p_value"]

    for alt in ("abs", "two-sided"):
        pv = p(alt)
        assert pv["TEST_UP"] < 0.01 and pv["TEST_DOWN"] < 0.01 and pv["TEST_NULL"] > 0.05
    pv = p("greater")
    assert pv["TEST_UP"] < 0.01 and pv["TEST_DOWN"] > 0.9
    pv = p("less")
    assert pv["TEST_DOWN"] < 0.01 and pv["TEST_UP"] > 0.9


def test_pvalue_never_zero_and_early_stop():
    scores = np.r_[np.full(10, 10.0), np.zeros(990)]
    cfg = EnrichmentConfig(max_reps=200, min_hits=5)
    res = resample_test(scores, np.arange(10), cfg, np.random.default_rng(0))
    assert res["p_value"] == pytest.approx(1 / 201)
    res = resample_test(scores, np.arange(500, 510), cfg, np.random.default_rng(0))
    assert res["n_extreme"] == 5 and res["n_draws"] == 5 and res["p_value"] == 1.0


def test_reproducible_across_jobs(toy):
    df, gmt = toy
    sets = read_gmt(gmt)
    a = run_enrichment(df, sets, "gene", "score", seed=42, jobs=1)
    b = run_enrichment(df, sets, "gene", "score", seed=42, jobs=2)
    pd.testing.assert_frame_equal(a, b)


def test_cli_groups_and_filters(toy, tmp_path):
    df, gmt = toy
    long = pd.concat([df.assign(cluster=0, kind="a"), df.assign(cluster=1, kind="a"), df.assign(cluster=1, kind="b")])
    inp = tmp_path / "in.tsv"
    long.to_csv(inp, sep="\t", index=False)
    out = tmp_path / "out" / "res.csv"
    rc = main(["-i", str(inp), "-g", str(gmt), "-o", str(out), "--gene-col", "gene", "--score-col", "score",
               "--group-col", "cluster", "--filter", "kind=a", "--prefix", "TEST_", "--seed", "3", "-q"])
    assert rc == 0
    res = pd.read_csv(out)
    assert sorted(res["cluster"].unique()) == [0, 1]
    assert (res.groupby("cluster")["n_matched"].max() == 60).all()  # filter removed duplicate kind=b rows
    assert set(res.columns) >= {"term", "p_value", "p_adj", "significant"}


def test_missing_column_errors(toy, tmp_path):
    df, gmt = toy
    inp = tmp_path / "in.csv"
    df.to_csv(inp, index=False)
    with pytest.raises(SystemExit):
        main(["-i", str(inp), "-g", str(gmt), "-o", str(tmp_path / "o.csv"), "--gene-col", "nope",
              "--score-col", "score", "-q"])


def test_org_aliases_and_urls():
    from doggo import msigdb

    assert msigdb.normalize_org("Mouse") == msigdb.normalize_org("mm") == "Mm"
    assert msigdb.normalize_org("human") == msigdb.normalize_org("HS") == "Hs"
    with pytest.raises(ValueError):
        msigdb.normalize_org("zebrafish")
    assert msigdb.gmt_url("mouse", "2022.1").endswith("/2022.1.Mm/msigdb.v2022.1.Mm.symbols.gmt")


def test_cli_org_uses_cache(toy, tmp_path, monkeypatch):
    df, gmt = toy
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "msigdb.v9999.1.Mm.symbols.gmt").write_text(gmt.read_text())
    monkeypatch.setenv("DOGGO_CACHE", str(cache))
    inp = tmp_path / "in.csv"
    df.to_csv(inp, index=False)
    out = tmp_path / "res.csv"
    rc = main(["-i", str(inp), "-org", "mouse", "--msigdb-version", "9999.1", "-o", str(out),
               "--gene-col", "gene", "--score-col", "score", "--prefix", "TEST_", "--seed", "1", "-q"])
    assert rc == 0
    assert set(pd.read_csv(out)["term"]) == {"TEST_UP", "TEST_DOWN", "TEST_NULL"}


def test_cli_requires_one_gene_set_source(toy, tmp_path):
    df, gmt = toy
    inp = tmp_path / "in.csv"
    df.to_csv(inp, index=False)
    base = ["-i", str(inp), "-o", str(tmp_path / "o.csv"), "--gene-col", "gene", "--score-col", "score"]
    with pytest.raises(SystemExit):
        main(base)
    with pytest.raises(SystemExit):
        main(base + ["-g", str(gmt), "--org", "mouse"])
