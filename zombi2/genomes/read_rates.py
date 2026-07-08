"""Read user-supplied rate tables from disk.

Two small tab/whitespace-separated formats, both with an optional header row and ``#`` comments:

* **per-family** (:func:`read_family_rates`) — one row per gene family, giving that family's own
  duplication / transfer / loss rates. Feeds :class:`~zombi2.genomes.rates.FamilySampledRates`
  (``rates=``): listed families use their tabulated rates, unlisted ones fall back to the model's
  distributions. Columns ``family duplication transfer loss`` (aliases ``dup``/``trans``/``loss``)::

      family  duplication  transfer  loss
      1       3            2         1
      2       4            0         1

* **per-branch** (:func:`read_branch_rates`) — one row per species-tree branch, giving that branch's
  transfer **emission** factor (how often it donates; scales the transfer rate via
  :class:`~zombi2.genomes.rates.BranchRates`) and **receptivity** weight (how likely it is to
  receive; biases recipient choice via :class:`~zombi2.genomes.transfers.TransferModel`). Columns
  ``branch emission receptivity`` — either optional::

      branch  emission  receptivity
      n3      5.0       1.0
      n7      1.0       10.0
"""

from __future__ import annotations

from pathlib import Path

_FAMILY_ALIASES = {"duplication": "duplication", "dup": "duplication",
                   "transfer": "transfer", "trans": "transfer",
                   "loss": "loss"}


def _rows(path):
    """Yield the split fields of each non-blank, non-comment line."""
    for raw in Path(path).read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if line:
            yield line.split()


def _looks_numeric(fields) -> bool:
    try:
        [float(x) for x in fields[1:]]
        return True
    except ValueError:
        return False


def read_family_rates(path) -> dict[str, tuple[float, float, float]]:
    """Read a per-family rate table into ``{family_id: (dup, transfer, loss)}``.

    A header row (``family``/``dup``/``trans``/``loss``, case-insensitive, any order) is honoured if
    present; otherwise columns are taken positionally as ``family dup transfer loss``.
    """
    rows = list(_rows(path))
    if not rows:
        return {}
    cols = {"family": 0, "duplication": 1, "transfer": 2, "loss": 3}
    start = 0
    if not _looks_numeric(rows[0]):  # header row: remap columns by name
        header = [h.lower() for h in rows[0]]
        cols = {"family": header.index("family") if "family" in header else 0}
        for i, name in enumerate(header):
            if name in _FAMILY_ALIASES:
                cols[_FAMILY_ALIASES[name]] = i
        for need in ("duplication", "transfer", "loss"):
            if need not in cols:
                raise ValueError(f"{path}: header is missing a '{need}' column")
        start = 1
    out: dict[str, tuple[float, float, float]] = {}
    for r in rows[start:]:
        fam = str(r[cols["family"]])
        out[fam] = (float(r[cols["duplication"]]), float(r[cols["transfer"]]),
                    float(r[cols["loss"]]))
    return out


def read_branch_rates(path) -> tuple[dict[str, float], dict[str, float]]:
    """Read a per-branch table into ``(emission_factors, receptivity_weights)``.

    ``emission_factors`` scales each listed branch's transfer *donation* rate; ``receptivity_weights``
    biases how likely each listed branch is to *receive*. Either column may be omitted (the
    corresponding map is then empty). A header row (``branch``/``emission``/``receptivity``) is
    honoured if present; otherwise columns are ``branch emission receptivity`` positionally.
    """
    rows = list(_rows(path))
    if not rows:
        return {}, {}
    emit_col, recept_col, branch_col = 1, 2, 0
    start = 0
    if not _looks_numeric(rows[0]):
        header = [h.lower() for h in rows[0]]
        branch_col = header.index("branch") if "branch" in header else 0
        emit_col = header.index("emission") if "emission" in header else None
        recept_col = header.index("receptivity") if "receptivity" in header else None
        if emit_col is None and recept_col is None:
            raise ValueError(f"{path}: header needs an 'emission' and/or 'receptivity' column")
        start = 1
    emission: dict[str, float] = {}
    receptivity: dict[str, float] = {}
    for r in rows[start:]:
        branch = str(r[branch_col])
        if emit_col is not None and emit_col < len(r):
            emission[branch] = float(r[emit_col])
        if recept_col is not None and recept_col < len(r):
            receptivity[branch] = float(r[recept_col])
    return emission, receptivity
