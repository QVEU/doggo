"""Download and cache MSigDB gene set collections.

MSigDB publishes collections for two organisms: human (Hs) and mouse (Mm).
Files are fetched once from the Broad Institute and cached locally; the
pinned default release is verified against known SHA-256 checksums.

Cache location: ``$DOGGO_CACHE`` if set, else ``$XDG_CACHE_HOME/doggo``,
else ``~/.cache/doggo``.
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import os
import shutil
import sys
import urllib.request
from pathlib import Path

log = logging.getLogger(__name__)

BASE_URL = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release"
DEFAULT_VERSION = "2026.2"
ID_TYPES = ("symbols", "entrez")

_ORG_ALIASES = {
    "hs": "Hs", "human": "Hs", "homo_sapiens": "Hs", "hsa": "Hs",
    "mm": "Mm", "mouse": "Mm", "mus_musculus": "Mm", "mmu": "Mm",
}
ORGANISMS = ("human", "mouse")

# Checksums of the pinned default release (downloaded 2026-11-06).
_SHA256 = {
    "msigdb.v2026.2.Hs.symbols.gmt": "dfaa03d6061cb060ba6430fbf3f19fe5ff74bb58bdfaf7f3ea767e6b70f17c09",
    "msigdb.v2026.2.Hs.entrez.gmt": "afa81be6e493413fcd99237363c8b2e276e7b7d6a6bac468d208c446c680414a",
    "msigdb.v2026.2.Mm.symbols.gmt": "10ef5fa078473fba9d17196cec8b08e4221ddedf30b4d14da80bd3ce5f01c691",
    "msigdb.v2026.2.Mm.entrez.gmt": "20fb6a0d71cf208d40d71d77cc65f93898bfae3a3b5490601c530bdfa4db7ad7",
}


def normalize_org(org: str) -> str:
    """Map ``human``/``mouse``/``hs``/``mm``/... to MSigDB's ``Hs``/``Mm``."""
    key = org.strip().lower().replace(" ", "_").replace("-", "_")
    if key not in _ORG_ALIASES:
        raise ValueError(f"unknown organism {org!r}; MSigDB provides: human (Hs), mouse (Mm)")
    return _ORG_ALIASES[key]


def cache_dir() -> Path:
    if os.environ.get("DOGGO_CACHE"):
        return Path(os.environ["DOGGO_CACHE"]).expanduser()
    xdg = os.environ.get("XDG_CACHE_HOME")
    return (Path(xdg).expanduser() if xdg else Path.home() / ".cache") / "doggo"


def gmt_filename(org: str, version: str = DEFAULT_VERSION, id_type: str = "symbols") -> str:
    if id_type not in ID_TYPES:
        raise ValueError(f"id_type must be one of {ID_TYPES}")
    return f"msigdb.v{version}.{normalize_org(org)}.{id_type}.gmt"


def gmt_url(org: str, version: str = DEFAULT_VERSION, id_type: str = "symbols") -> str:
    o = normalize_org(org)
    return f"{BASE_URL}/{version}.{o}/{gmt_filename(o, version, id_type)}"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_gmt(
    org: str,
    version: str = DEFAULT_VERSION,
    id_type: str = "symbols",
    directory: str | Path | None = None,
    force: bool = False,
) -> Path:
    """Return a local path to the full MSigDB GMT for ``org``, downloading if needed."""
    name = gmt_filename(org, version, id_type)
    dest_dir = Path(directory) if directory else cache_dir()
    dest = dest_dir / name
    expected = _SHA256.get(name)

    if dest.exists() and not force:
        if expected is None or _sha256(dest) == expected:
            return dest
        log.warning("cached %s failed checksum; downloading again", dest)

    dest_dir.mkdir(parents=True, exist_ok=True)
    url = gmt_url(org, version, id_type)
    log.info("downloading %s", url)
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "doggo"})
        with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as out:
            shutil.copyfileobj(r, out)
    except Exception as e:  # noqa: BLE001 - surface a clear message for any network failure
        tmp.unlink(missing_ok=True)
        raise RuntimeError(
            f"could not download {url}: {e}. Download it manually and pass it with -g/--gmt."
        ) from e

    if expected is not None:
        got = _sha256(tmp)
        if got != expected:
            tmp.unlink(missing_ok=True)
            raise RuntimeError(f"checksum mismatch for {name}: expected {expected}, got {got}")
    else:
        log.warning("no pinned checksum for %s (not the default release); not verified", name)
    tmp.replace(dest)
    return dest


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="doggo-fetch",
        description="Download MSigDB gene set collections into the doggo cache and print their paths.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--org", nargs="+", default=list(ORGANISMS),
                   help="organisms to fetch: human, mouse (aliases: hs, mm)")
    p.add_argument("--msigdb-version", default=DEFAULT_VERSION, help="MSigDB release, e.g. 2026.2 or 2022.1")
    p.add_argument("--id-type", nargs="+", choices=ID_TYPES, default=["symbols"], help="gene identifier type(s)")
    p.add_argument("--dir", default=None, help=f"download directory (default: {cache_dir()})")
    p.add_argument("--force", action="store_true", help="download again even if cached")
    args = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s", stream=sys.stderr)
    try:
        for org in args.org:
            for t in args.id_type:
                print(fetch_gmt(org, args.msigdb_version, t, args.dir, args.force))
    except (ValueError, RuntimeError) as e:
        raise SystemExit(f"error: {e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
