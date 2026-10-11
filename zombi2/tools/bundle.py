"""Convert a written run between its two layouts: one file per family, or one file per run.

A run written with ``bundle=True`` (``--bundle``) holds each per-family output as one file
(`zombi2._runtime.outputs.BUNDLED`). `bundle_run` turns a run written without it into that layout,
and `unbundle_run` turns it back, so a run already on disk needs no new simulation either way.

Both conversions are exact: the files one writes are the files the run would have written in the
other layout, byte for byte. Each checks that before it touches anything else, by reading back what
it wrote and comparing it with what it read, and stops with an error if the two differ. The files it
converted from stay unless ``remove=True``.
"""
from __future__ import annotations

import os
import pathlib
import re
from dataclasses import dataclass

from .._runtime.outputs import BUNDLED, tree_table_header, tree_table_row


@dataclass(frozen=True)
class _Output:
    """One per-family output: its directory, the pattern its files are named by, and whether each
    file is a tree (a ``complete`` and ``extant`` pair) or a FASTA."""

    directory: str
    pattern: re.Pattern
    trees: bool

    @property
    def bundled(self) -> str:
        return BUNDLED[self.directory]


#: Every per-family output, in the order a conversion reports them. A file's ``stem`` is ``fam``, or
#: ``block`` on a nucleotide run's sequences, and ``which`` is ``complete`` or ``extant``.
_OUTPUTS = (
    _Output("gene_trees", re.compile(r"gene_tree_(?P<stem>fam)(?P<id>\d+)_(?P<which>complete|extant)\.nwk"),
            trees=True),
    _Output("phylograms",
            re.compile(r"phylogram_(?P<stem>fam|block)(?P<id>\d+)_(?P<which>complete|extant)\.nwk"),
            trees=True),
    _Output("alignments", re.compile(r"(?P<stem>fam|block)(?P<id>\d+)\.fasta"), trees=False),
    _Output("ancestral", re.compile(r"sequences_ancestral_(?P<stem>fam|block)(?P<id>\d+)\.fasta"),
            trees=False),
)

#: What a tree table's first column is called, by the stem of the files it came from.
_KEY = {"fam": "family", "block": "block"}


@dataclass(frozen=True)
class Converted:
    """What one conversion did, for one output in one directory."""

    #: the bundled file, or the per-family directory, that was written
    wrote: str
    #: how many families (or blocks) it holds
    units: int
    #: ``family``, or ``block`` for a nucleotide run's sequences
    unit: str
    #: how many files it was converted from
    files: int
    #: whether those files were removed
    removed: bool


# --- the per-family files, read and written -------------------------------------------------------

def _per_family_files(place: pathlib.Path, out: _Output, flat: bool) -> dict[str, pathlib.Path]:
    """``{file name: path}`` for every per-family file of ``out`` at ``place``: in its directory,
    or straight in ``place`` under ``flat``."""
    folder = place if flat else place / out.directory
    if not folder.is_dir():
        return {}
    return {p.name: p for p in sorted(folder.iterdir()) if out.pattern.fullmatch(p.name)}


def _bundled_text(out: _Output, files: dict[str, pathlib.Path]) -> str:
    """The bundled file built from the per-family ``files``, in family order."""
    by_unit: dict[tuple[str, int], dict[str, str]] = {}
    for name, path in files.items():
        m = out.pattern.fullmatch(name)
        assert m is not None
        by_unit.setdefault((m["stem"], int(m["id"])), {})[m["which"] if out.trees else "fasta"] = (
            path.read_text(encoding="utf-8"))
    units = sorted(by_unit, key=lambda unit: unit[1])
    if len({stem for stem, _ in units}) > 1:
        raise ValueError(f"{out.directory}/ mixes families and blocks; it is not one run's output")
    if out.trees:
        key = _KEY[units[0][0]]
        rows = [tree_table_header(key)]
        for unit in units:
            pair = by_unit[unit]
            if "complete" not in pair:
                raise ValueError(f"{out.directory}/ has an extant tree for {unit[0]}{unit[1]} but no "
                                 f"complete one, which no run writes")
            rows.append(tree_table_row(unit[1], _one_line(pair["complete"], out, unit),
                                       _one_line(pair["extant"], out, unit)
                                       if "extant" in pair else None))
        return "".join(row + "\n" for row in rows)
    parts: list[str] = []
    for stem, uid in units:
        lines = by_unit[(stem, uid)]["fasta"].splitlines(keepends=True)
        parts.extend(line.rstrip("\n") + f" {stem}{uid}\n" if line.startswith(">") else line
                     for line in lines)
    return "".join(parts)


def _one_line(text: str, out: _Output, unit: tuple[str, int]) -> str:
    """A Newick file's one tree, without its newline. A file holding anything else has no row it
    could become, and saying so beats a table that quietly loses it."""
    tree = text[:-1] if text.endswith("\n") else text
    if "\n" in tree or "\t" in tree or not tree:
        raise ValueError(f"{out.directory}/ file for {unit[0]}{unit[1]} is not one Newick tree on "
                         f"one line, so it has no exact row in {out.bundled}")
    return tree


# --- the bundled file, read and written ------------------------------------------------------------

def _per_family_texts(out: _Output, text: str) -> dict[str, str]:
    """``{file name: contents}``: the per-family files a bundled file stands for."""
    files: dict[str, str] = {}
    if out.trees:
        lines = text.splitlines()
        if not lines or lines[0].split("\t")[1:] != ["complete", "extant"]:
            raise ValueError(f"{out.bundled} does not start with the header "
                             f"'family<TAB>complete<TAB>extant'")
        stem = {v: k for k, v in _KEY.items()}.get(lines[0].split("\t")[0])
        if stem is None:
            raise ValueError(f"{out.bundled}'s first column is neither 'family' nor 'block'")
        for n, line in enumerate(lines[1:], start=2):
            cells = line.split("\t")
            if len(cells) != 3 or not cells[0].isdigit() or not cells[1]:
                raise ValueError(f"{out.bundled} line {n} is not '<id><TAB>complete<TAB>extant'")
            uid, complete, extant = cells
            prefix = "gene_tree_" if out.directory == "gene_trees" else "phylogram_"
            files[f"{prefix}{stem}{uid}_complete.nwk"] = complete + "\n"
            if extant:
                files[f"{prefix}{stem}{uid}_extant.nwk"] = extant + "\n"
        return files
    name = "{}.fasta" if out.directory == "alignments" else "sequences_ancestral_{}.fasta"
    current = None
    for n, line in enumerate(text.splitlines(keepends=True), start=1):
        if line.startswith(">"):
            header, _, tag = line.rstrip("\n").rpartition(" ")
            if not header or not re.fullmatch(r"(fam|block)\d+", tag):
                raise ValueError(f"{out.bundled} line {n} names no family: a header here is "
                                 f"'>name fam<N>'")
            current = name.format(tag)
            files[current] = files.get(current, "") + header + "\n"
        elif current is None:
            raise ValueError(f"{out.bundled} line {n} comes before the first header")
        else:
            files[current] += line
    return files


# --- the two conversions --------------------------------------------------------------------------

def _places(directory) -> list[pathlib.Path]:
    """Where a run keeps per-family outputs: its ``genomes/`` and ``sequences/`` directories, and the
    directory itself, which is where a ``--flat`` run or a single level's ``write`` puts them."""
    d = pathlib.Path(directory)
    if not d.is_dir():
        raise FileNotFoundError(f"{d} is not a directory")
    return [p for p in (d, d / "genomes", d / "sequences") if p.is_dir()]


def bundle_run(directory, *, remove: bool = False) -> list[Converted]:
    """Bundle every per-family output under ``directory``: each output's per-family files become its
    one file (`BUNDLED`), beside them. ``directory`` is a run, or one level of it.

    The bundled file is read back and compared with the files it came from before anything else
    happens; ``remove=True`` then removes those files, and their directory once it is empty. An
    output already in both layouts is left as it is when the two agree, and is an error when they
    do not: which of the two describes the run cannot be told from the files."""
    done = []
    for place in _places(directory):
        for out in _OUTPUTS:
            in_folder = _per_family_files(place, out, flat=False)
            files = in_folder or _per_family_files(place, out, flat=True)
            if not files:
                continue
            originals = {name: path.read_text(encoding="utf-8") for name, path in files.items()}
            target = place / out.bundled
            if target.exists():
                _check_agree(out, place, target, originals)
            else:
                target.write_text(_bundled_text(out, files), encoding="utf-8")
                if _per_family_texts(out, target.read_text(encoding="utf-8")) != originals:
                    target.unlink()
                    raise ValueError(f"{out.directory} at {place} did not survive bundling exactly, "
                                     f"so nothing was changed")
            if remove:
                for path in files.values():
                    path.unlink()
                if in_folder and not any((place / out.directory).iterdir()):
                    (place / out.directory).rmdir()
            done.append(_converted(out, str(target), originals, len(files), remove))
    return done


def unbundle_run(directory, *, flat: bool = False, remove: bool = False) -> list[Converted]:
    """Unbundle every bundled output under ``directory``: each bundled file becomes its per-family
    files, in its directory, or beside it under ``flat``.

    The files written are read back and compared with the bundled file before anything else
    happens; ``remove=True`` then removes the bundled file. An output already in both layouts is
    left as it is when the two agree, and is an error when they do not, as in `bundle_run`."""
    done = []
    for place in _places(directory):
        for out in _OUTPUTS:
            source = place / out.bundled
            if not source.is_file():
                continue
            texts = _per_family_texts(out, source.read_text(encoding="utf-8"))
            folder = place if flat else place / out.directory
            existing = (_per_family_files(place, out, flat=False)
                        or _per_family_files(place, out, flat=True))
            if existing:
                _check_agree(out, place, source,
                             {name: path.read_text(encoding="utf-8")
                              for name, path in existing.items()})
                folder = next(iter(existing.values())).parent
            else:
                folder.mkdir(exist_ok=True)
                for name, text in texts.items():
                    (folder / name).write_text(text, encoding="utf-8")
                written = {name: folder / name for name in texts}
                # compared as files rather than as one text: a streamed run bundles its families
                # in the order they finished, and `_bundled_text` puts them in family order
                if _per_family_texts(out, _bundled_text(out, written)) != texts:
                    for path in written.values():
                        path.unlink()
                    raise ValueError(f"{source} did not survive unbundling exactly, so nothing was "
                                     f"changed")
            if remove:
                source.unlink()
            done.append(_converted(out, str(folder) + os.sep, texts, 1, remove))
    return done


def _check_agree(out: _Output, place: pathlib.Path, bundled: pathlib.Path,
                 per_family: dict[str, str]) -> None:
    """Refuse an output found in both layouts unless the two hold the same files."""
    if _per_family_texts(out, bundled.read_text(encoding="utf-8")) != per_family:
        raise ValueError(f"{place} holds {out.directory} in both layouts, {bundled.name} and "
                         f"per-family files, and they differ; remove the one that is not this run's")


def _converted(out: _Output, wrote: str, files: dict[str, str], n_files: int,
               removed: bool) -> Converted:
    """What a conversion reports, read off the per-family file names it handled."""
    matches = [out.pattern.fullmatch(name) for name in files]
    return Converted(wrote, len({m["id"] for m in matches if m}),
                     _KEY[matches[0]["stem"]] if matches and matches[0] else "family",
                     n_files, removed)


__all__ = ["Converted", "bundle_run", "unbundle_run"]
