"""Whole-pipeline CLI determinism: same seed => byte-identical output trees.

Ported from ``thekswenson/Zombi`` ``tests/test_randomization.py`` (which ran two same-seed
projects through T -> G -> S and asserted identical output directories). zombi2 already has
API-level ``*_reproducible`` tests; this is the end-to-end **CLI / output-file** analogue across
``species -> genomes -> sequence``, so it also guards the *writers* (file naming, row ordering,
serialization) — not just the RNG stream.

Gzipped outputs embed an mtime in the gzip header, so files are compared by their *decompressed*
payload rather than raw bytes.
"""

import gzip
import hashlib

from zombi2.cli import main


def _content_hash(path):
    data = path.read_bytes()
    if path.suffix == ".gz":
        try:
            data = gzip.decompress(data)  # ignore the mtime carried in the gzip header
        except OSError:
            pass
    return hashlib.sha256(data).hexdigest()


def _fingerprint(root):
    """relpath -> content hash for every file under ``root`` (order-independent).

    ``.log`` files are skipped: they record wall-clock timing (``... in 0.01 s``), which is
    provenance, not simulated content, and so legitimately differs run to run.
    """
    return {
        str(p.relative_to(root)): _content_hash(p)
        for p in root.rglob("*")
        if p.is_file() and p.suffix != ".log"
    }


def _assert_identical(d1, d2, label):
    f1, f2 = _fingerprint(d1), _fingerprint(d2)
    assert f1, f"{label}: no output files were produced — nothing to compare"
    assert set(f1) == set(f2), f"{label}: file sets differ: {set(f1) ^ set(f2)}"
    differing = [k for k in f1 if f1[k] != f2[k]]
    assert not differing, f"{label}: same seed produced differing files: {differing}"


def _run_pipeline(base):
    sp, gen, seq = base / "S", base / "G", base / "Q"
    assert main(["species", "--birth", "1", "--death", "0.3", "--tips", "12", "--age", "4",
                 "--seed", "7", "-o", str(sp)]) == 0
    assert main(["genomes", "--tree", str(sp / "species_tree.nwk"),
                 "--write", "trace", "profiles",  # trace needed so 'sequence' can replay
                 "--seed", "7", "-o", str(gen)]) == 0
    assert main(["sequences", "--genomes", str(gen), "--family-speed", "0.5",
                 "--seed", "7", "-o", str(seq)]) == 0
    return {"species": sp, "genomes": gen, "sequence": seq}


def test_pipeline_is_deterministic_under_fixed_seed(tmp_path):
    a = _run_pipeline(tmp_path / "proj1")
    b = _run_pipeline(tmp_path / "proj2")
    for stage in ("species", "genomes", "sequence"):
        _assert_identical(a[stage], b[stage], stage)
