"""``run.zombi2`` — one human-readable report for a whole run directory.

A run accumulates level by level in one directory, each command dropping its own machine records there:
a ``<level>.log`` (versions, timestamp, the resolved parameters, the input files by content hash, the
result line) and a ``<level>_summary.json`` (the headline numbers). Those are complete but scattered —
to see what a run *is* you open four files across two sub-directories, and none of them ties the levels
together or says what ``profiles.tsv`` holds.

This builds the missing single view. It is a **derived projection** of those per-level records, never a
source of truth: ``write_run_report`` reads whatever records the directory currently holds and
rewrites ``run.zombi2`` from them, so running a new level, or re-running one, just refreshes it — there
is no document the levels take turns appending to. ``zombi2 report DIR`` rebuilds it on demand.

It reports, per level: the parameters, the inputs it was computed on, the result, and **every file it
wrote** (a per-family directory is summarised by count, not enumerated). Because it records the input
files by hash, it can also see when a downstream level was computed on an upstream one that has since
changed, and say so.
"""
from __future__ import annotations

import hashlib
import os

from zombi2._runtime.summary import read_summary

RUN_REPORT_NAME = "run.zombi2"

#: The levels, in pipeline order, each grown into ``<run>/<level>/``. A ``joint`` run is not here: it
#: writes one report spanning ``species/`` and its driver (see ``_joint_section``).
_LEVELS = ("species", "genomes", "sequences", "traits")

#: The record files a level drops beside its data — not listed one-per-line among the outputs, but named
#: together on one line so nothing is undocumented.
_RECORD_SUFFIXES = (".log", "_summary.json")

#: Parameters that are plumbing, not model: how the run was driven, not what was simulated. Dropped from
#: the parameters line and the reproduce command (the ``.log`` keeps the complete set either way).
_SKIP_PARAMS = frozenset({"command", "run", "seed", "params", "flat", "force", "quiet", "write",
                          "parallel", "stream", "max_lineages", "source", "fasta", "tip_fates"})

#: A readable order for the parameters that do show — model choice first, then rates, then the rest
#: (anything unlisted follows, alphabetically). Only affects display order.
_PARAM_ORDER = ("resolution", "kind", "model", "birth", "death", "n_extant", "total_time", "sampling",
                "fossils", "mass_extinction", "initial_families", "origination", "duplication",
                "transfer", "loss", "transfer_to", "self_transfer", "replacement", "max_family_size",
                "inversion", "chromosomes", "length", "divergence", "kappa", "frequencies",
                "substitution", "family_speed", "rate", "reverts_to", "pull", "states", "switch")

#: Prettier labels for a few summary keys / parameters; everything else has its underscores turned to
#: spaces. Kept small on purpose — a key that reads fine as-is does not need an entry.
_LABELS = {
    "genes_per_genome": "genes/genome", "copies_per_family_per_genome": "copies/family",
    "total_branch_length": "total branch length", "stem_length": "stem", "realised_rates": "realised",
    "value_at_root_node": "value at root", "families_with_sequences": "families with sequences",
    "ancestral_sequences": "ancestral sequences", "assembled_genomes": "assembled genomes",
    "event_rows": "event rows", "extant_genomes": "extant genomes", "n_extant": "n-extant",
    "total_time": "total-time", "max_family_size": "max-family-size", "transfer_to": "transfer-to",
}

#: What each output file is, by basename within a level directory (per-family directories keyed by their
#: directory name). Data files must all appear here — ``test_report`` runs a full pipeline and fails on
#: any output that is neither listed nor a record, so this cannot silently fall behind the outputs.
_GLOSS = {
    # species/
    "species_extant.nwk": "Newick — the sampled extant tips only, the tree most analyses want",
    "species_complete.nwk": "Newick — the whole tree, all nodes including extinct lineages",
    "species_events.tsv": "one row per speciation/extinction (time, kind, lineage, children)",
    "species_fates.tsv": "each tip's fate: extant / extinct / unsampled",
    # genomes/
    "genomes.tsv": "gene content of every genome — one row per gene copy (lineage, family, copy)",
    "profiles.tsv": "family × extant-genome copy-number matrix",
    "initial_genome.tsv": "the genome present at t=0, that the initial families started from",
    "genome_events.tsv": "every duplication / transfer / loss / origination, per gene copy",
    "gene_trees": "gene trees, one per family — complete (all copies) and extant (survivors), Newick",
    # sequences/
    "clock_species_tree_complete.nwk": "the species tree rescaled to substitutions by the clock, all nodes",
    "clock_species_tree_extant.nwk": "the clock-rescaled tree, sampled tips only",
    "alignments": "one FASTA alignment per family (the extant sequences)",
    "phylograms": "gene trees with branch lengths in substitutions/site, one per family",
    # traits/
    "trait_tree.nwk": "the tree with each node's trait annotated ([&trait=…])",
    "trait_values.tsv": "the trait value (or state) at every node",
    "trait_events.tsv": "each discrete state change along a branch (time, lineage, from→to)",
}

_RULE = "─" * 80
_HEAVY = "═" * 80


# ── reading the per-level records ───────────────────────────────────────────────────────────────

def _read_log(path: str) -> dict | None:
    """Parse a ``<level>.log`` into its parts, or ``None`` if it is not there. The file is the header
    lines (version/timestamp/command), then ``input<TAB>sha<TAB>path`` rows, then the resolved
    parameters as ``key<TAB>value``, then ``result``."""
    if not os.path.isfile(path):
        return None
    header_keys = {"zombi2_version", "python_version", "numpy_version", "platform", "timestamp",
                   "command_line", "result"}
    out: dict = {"inputs": [], "params": {}}
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            key = parts[0]
            if key == "input" and len(parts) == 3:
                out["inputs"].append((parts[2], parts[1]))         # (path, sha256)
            elif key in header_keys:
                out[key] = parts[1] if len(parts) > 1 else ""
            elif len(parts) >= 2:
                out["params"][key] = parts[1]
    return out


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


# ── formatting helpers ──────────────────────────────────────────────────────────────────────────

def _label(key: str) -> str:
    return _LABELS.get(key, key.replace("_", " "))


def _num(x) -> str:
    """A number for a human: integers plain, floats to four significant figures without trailing zeros."""
    if isinstance(x, bool):
        return str(x)
    if isinstance(x, int):
        return str(x)
    if isinstance(x, float):
        return f"{x:.4g}"
    return str(x)


def _stat_lines(summary: dict) -> list[str]:
    """The result block for a level, from its ``_summary.json``: one readable line per top-level entry.

    Counts read ``value noun`` (``30 extant``); measured quantities read ``noun value`` (``height 20.1``);
    a ``{min,mean,max}`` triple collapses to ``min – max (mean …)``. A family-size cap that was reached is
    flagged, because a truncated genome is the one summary number a reader most needs to see."""
    lines = []
    for key, value in summary.items():
        if key in ("level", "seed"):
            continue
        if key == "mean_pairwise_identity":
            lines.append(f"mean pairwise identity {float(value) * 100:.1f}%")
        elif key == "family_size_cap":
            hit = value.get("families_at_cap", 0)
            if hit:
                lines.append(f"⚠ {_label('max_family_size')} {value['cap']}: {hit} "
                             f"famil{'y' if hit == 1 else 'ies'} reached it, "
                             f"{value['cells_at_cap']} genome cell(s) truncated")
        elif isinstance(value, dict):
            lines.append(f"{_label(key)} {_dict_phrase(value)}")
        elif isinstance(value, list):
            lines.append(f"{_label(key)} {', '.join(map(str, value))}")
        else:
            lines.append(f"{_label(key)} {_num(value)}")
    return lines


def _dict_phrase(d: dict) -> str:
    """One ``·``-joined phrase for a sub-dict of a summary. A ``{min,mean,max}`` triple becomes a range;
    a dict of integer counts reads ``value noun``; anything else (rates, lengths) reads ``noun value``."""
    if set(d) <= {"min", "mean", "max"}:
        if d.get("min") == d.get("max"):
            return _num(d.get("min"))
        span = f"{_num(d['min'])} – {_num(d['max'])}"
        return f"{span} (mean {_num(d['mean'])})" if "mean" in d else span
    counts = all(isinstance(v, int) and not isinstance(v, bool) for v in d.values())
    if counts:
        return " · ".join(f"{_num(v)} {_label(k)}" for k, v in d.items())
    return " · ".join(f"{_label(k)} {_num(v)}" for k, v in d.items())


def _visible_params(params: dict) -> list[tuple[str, str]]:
    """The parameters worth showing, in reading order — the model parameters, dropping plumbing
    (`_SKIP_PARAMS`), anything left unset (``None`` / ``False`` / blank), and ``substitution`` when
    ``divergence`` is set (the two are mutually exclusive on the command line, and divergence is the one
    the user chose — the log keeps the substitution rate it solved to)."""
    drop = set(_SKIP_PARAMS)
    if params.get("divergence") not in (None, "None", ""):
        drop.add("substitution")
    items = [(k, v) for k, v in params.items() if k not in drop and v not in ("None", "False", "")]
    items.sort(key=lambda kv: (_PARAM_ORDER.index(kv[0]) if kv[0] in _PARAM_ORDER else len(_PARAM_ORDER),
                               kv[0]))
    return items


def _params_phrase(params: dict) -> str:
    return " · ".join(f"{_label(k)} {v}" for k, v in _visible_params(params))


def _reproduce(command: str, run: str, params: dict, seed: str | None, source: str | None) -> str:
    """A clean, copy-pasteable command that reproduces this level — rebuilt from the resolved parameters
    (not the raw command line, whose shell quoting is already lost), so a rate expression comes back
    correctly quoted."""
    parts = [f"zombi2 {command} {run}"]
    if source:
        parts.append(f"--from {source}")
    for k, v in _visible_params(params):
        flag, val = f"--{k.replace('_', '-')}", _flag_value(v)
        parts.append(f"{flag} {val}" if val else flag)       # a store-true flag takes no value
    if seed not in (None, "None"):
        parts.append(f"--seed {seed}")
    return " ".join(parts)


def _flag_value(value: str) -> str:
    """One parameter as it goes back on the command line. A list the log rendered as ``[a, b, c]`` (a
    multi-value flag like ``--frequencies``) becomes space-separated bare values; ``True`` (a store-true
    flag) drops its value; anything with shell-significant characters is quoted."""
    if value == "True":
        return ""
    if value.startswith("[") and value.endswith("]"):
        import ast
        try:
            return " ".join(_flatten(ast.literal_eval(value)))
        except (ValueError, SyntaxError):
            pass
    return f'"{value}"' if any(c in value for c in " *(){},'") else value


def _flatten(seq) -> list[str]:
    out = []
    for x in seq if isinstance(seq, (list, tuple)) else [seq]:
        out.extend(_flatten(x)) if isinstance(x, (list, tuple)) else out.append(str(x))
    return out


# ── output-file listing ─────────────────────────────────────────────────────────────────────────

def _list_outputs(level_dir: str) -> tuple[list[tuple[str, str]], list[str]]:
    """``(data files, record files)`` for a level directory. Each data entry is ``(name, gloss)``; a
    sub-directory of per-family files is one entry, its count folded into the gloss. Record files (the
    ``.log`` and ``_summary.json``) are returned by name only."""
    data, records = [], []
    for name in sorted(os.listdir(level_dir)):
        full = os.path.join(level_dir, name)
        if name == RUN_REPORT_NAME:
            continue
        if name.endswith(_RECORD_SUFFIXES):
            records.append(name)
        elif os.path.isdir(full):
            n = len([f for f in os.listdir(full) if not f.startswith(".")])
            gloss = _GLOSS.get(name, "per-family output files")
            data.append((f"{name}/", f"{gloss} ({n} files)"))
        else:
            data.append((name, _GLOSS.get(name, "")))
    return data, records


# ── staleness ───────────────────────────────────────────────────────────────────────────────────

def _stale_notes(run: str, sections: list[dict]) -> list[str]:
    """Warnings for any level computed on a file inside this run that has since changed — the input's
    recorded hash no longer matches the file on disk, so an upstream level was re-run without this one."""
    run_abs = os.path.abspath(run)
    notes = []
    for sec in sections:
        for path, recorded in sec["log"].get("inputs", []):
            abspath = os.path.abspath(path)
            if os.path.commonpath([run_abs, abspath]) != run_abs or not os.path.isfile(abspath):
                continue                                    # an input from outside the run, or now gone
            if _sha256(abspath) != recorded:
                rel = os.path.relpath(abspath, run_abs)
                notes.append(f"⚠ {sec['title']} was computed on {rel}, which has changed since — "
                             f"re-run {sec['command']} so the run agrees with itself")
    return notes


# ── assembling the report ───────────────────────────────────────────────────────────────────────

def _section(run: str, title: str, command: str, level_dir: str, log: dict, summary: dict) -> dict:
    return {"title": title, "command": command, "dir": level_dir, "log": log, "summary": summary}


def _collect_sections(run: str) -> list[dict]:
    """Every command that wrote into this run, as sections in pipeline order. A ``joint`` run contributes
    one section spanning ``species/`` and its driver; each ordinary level with its own ``.log`` another."""
    sections = []
    joint_summary = os.path.join(run, "joint_summary.json")
    if os.path.isfile(joint_summary):
        sections.append(_joint_section(run, read_summary(joint_summary)))
    for level in _LEVELS:
        level_dir = os.path.join(run, level)
        log = _read_log(os.path.join(level_dir, f"{level}.log"))
        summary_path = os.path.join(level_dir, _SUMMARY_FILE.get(level, f"{level}_summary.json"))
        if log is not None and os.path.isfile(summary_path):
            sections.append(_section(run, level.upper(), level, level_dir, log,
                                     read_summary(summary_path)))
    return sections


#: the summary file each level writes (all ``<level>_summary.json`` except genomes' ``genome_summary``).
_SUMMARY_FILE = {"genomes": "genome_summary.json", "traits": "trait_summary.json"}


def _joint_section(run: str, summary: dict) -> dict:
    """The section for a ``joint`` run: its one ``joint.log`` lives in ``species/``, and it wrote the
    species tree plus its driver (a trait, or a named genomes family)."""
    driver = summary.get("driver", "trait")
    driver_dir = "traits" if driver == "trait" else "genomes"
    log = _read_log(os.path.join(run, "species", "joint.log")) or {"params": {}, "inputs": []}
    section = _section(run, "JOINT (species + driver)", "joint", os.path.join(run, "species"),
                       log, summary)
    section["extra_dirs"] = [os.path.join(run, driver_dir)]
    return section


def _one_liner(sections: list[dict]) -> str:
    """A single sentence for the top of the report, synthesised from the sections present."""
    by = {s["command"]: s["summary"] for s in sections}
    if "joint" in by:
        sp = by["joint"].get("species", {}).get("tips", {})
        return (f"A species tree of {sp.get('extant', '?')} tips, its diversification driven by "
                f"a {by['joint'].get('driver', 'trait')}.")
    bits = []
    if "species" in by:
        bits.append(f"A species tree of {by['species']['tips']['extant']} sampled tips")
    if "genomes" in by:
        bits.append(f"{by['genomes']['families']['surviving']} gene families evolved along it"
                    if bits else "gene families along a species tree")
    if "sequences" in by:
        bits.append("sequences down each gene tree")
    if "traits" in by:
        bits.append(f"a {by['traits'].get('kind', '')} trait".strip())
    return (", on which ".join(bits[:2]) + ("; " + "; ".join(bits[2:]) if bits[2:] else "") + "."
            if bits else "An empty run.")


def build_run_report(run: str) -> str | None:
    """The text of ``run.zombi2`` for a run directory, or ``None`` if it holds no level records yet
    (e.g. a ``--flat`` run, whose levels share one directory and cannot be told apart)."""
    sections = _collect_sections(run)
    if not sections:
        return None

    from zombi2 import __version__
    env = next((s["log"] for s in sections if s["log"].get("zombi2_version")), {})
    lines = [f"ZOMBI2 {env.get('zombi2_version', __version__)}  ·  run report", _HEAVY,
             f"  {_one_liner(sections)}"]
    built = " · ".join(p for p in (
        f"zombi2 {env['zombi2_version']}" if env.get("zombi2_version") else "",
        f"numpy {env['numpy_version']}" if env.get("numpy_version") else "",
        f"python {env['python_version']}" if env.get("python_version") else "",
        env.get("platform", "")) if p)
    if built:
        lines.append(f"  built with {built}")
    for note in _stale_notes(run, sections):
        lines.append(f"  {note}")
    lines.append("")

    for i, sec in enumerate(sections, 1):
        lines += _render_section(run, sec, i, len(sections))

    lines += ["TO REPRODUCE", _RULE]
    for sec in sections:
        seed = sec["log"].get("params", {}).get("seed") or sec["summary"].get("seed")
        lines.append("  " + _reproduce(sec["command"], run, sec["log"].get("params", {}),
                                       seed, _source_of(sec)))
    lines.append("")
    return "\n".join(lines)


def _source_of(section: dict) -> str | None:
    src = section["log"].get("params", {}).get("source")
    return src if src not in (None, "None", "") else None


def _field(label: str, contents: list[str]) -> list[str]:
    """A labelled, aligned block: the label on the first content line, continuations blank, everything
    lined up at the same column."""
    return [f"  {(label if k == 0 else ''):<11}  {c}" for k, c in enumerate(contents)]


def _render_section(run: str, sec: dict, i: int, n: int) -> list[str]:
    log, summary = sec["log"], sec["summary"]
    seed = summary.get("seed", log.get("params", {}).get("seed", "?"))
    tag = f"level {i}/{n}"
    lines = [f"{sec['title']}{' ' * max(1, 80 - len(sec['title']) - len(tag))}{tag}", _RULE]

    stamp = " · ".join(p for p in (f"seed {seed}", log.get("timestamp", "").replace("T", " ")) if p)
    if stamp:
        lines.append(f"  {stamp}")
    params = _params_phrase(log.get("params", {}))
    if params:
        lines += _field("parameters", _wrap(params))
    lines += _field("computed on", [os.path.relpath(p, run) if _within(run, p) else p
                                    for p, _ in log.get("inputs", [])])

    stats = _stat_lines(summary if sec["command"] != "joint" else _flatten_joint(summary))
    lines += _field("result", stats)

    data, records = _gather_outputs(sec)
    file_lines = [f"{name:<28}  {gloss}".rstrip() for name, gloss in data]
    if records:
        file_lines.append(f"records: {', '.join(records)}")
    lines += _field("files", file_lines)
    lines.append("")
    return lines


def _gather_outputs(sec: dict) -> tuple[list[tuple[str, str]], list[str]]:
    """A section's outputs, across its own directory and any extra directories (a joint driver)."""
    data, records = _list_outputs(sec["dir"])
    for extra in sec.get("extra_dirs", []):
        if os.path.isdir(extra):
            name = os.path.basename(extra)
            d, r = _list_outputs(extra)
            data += [(f"{name}/{n}", g) for n, g in d]
            records += [f"{name}/{x}" for x in r]
    return data, records


def _flatten_joint(summary: dict) -> dict:
    """A joint summary nests the species and driver summaries; show the species numbers as the section's
    result (the driver has its own files listed), prefixed so they read in context."""
    sp = dict(summary.get("species", {}))
    sp.pop("level", None)
    sp.pop("seed", None)
    return {"driver": summary.get("driver", "trait"), **sp}


def _within(run: str, path: str) -> bool:
    run_abs, p_abs = os.path.abspath(run), os.path.abspath(path)
    return os.path.commonpath([run_abs, p_abs]) == run_abs


def _wrap(phrase: str, width: int = 64) -> list[str]:
    """Break a ``·``-joined phrase into lines no wider than ``width`` (never mid-item)."""
    lines, cur = [], ""
    for item in phrase.split(" · "):
        piece = (cur + " · " + item) if cur else item
        if len(piece) > width and cur:
            lines.append(cur)
            cur = item
        else:
            cur = piece
    if cur:
        lines.append(cur)
    return lines


def write_run_report(run: str) -> str | None:
    """(Re)write ``run/run.zombi2`` from the level records the directory now holds; return its path, or
    ``None`` if there is nothing to report yet. Safe to call after every level run — it is a projection,
    so it always reflects the directory's current state."""
    text = build_run_report(run)
    if text is None:
        return None
    path = os.path.join(run, RUN_REPORT_NAME)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)
    return path


__all__ = ["RUN_REPORT_NAME", "build_run_report", "write_run_report"]
