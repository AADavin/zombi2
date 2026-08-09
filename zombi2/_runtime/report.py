"""``run.zombi2`` — one human-readable report for a whole run directory.

A run accumulates level by level in one directory, each command dropping its own machine records there:
a ``<level>.log`` (versions, timestamp, the resolved parameters, the input files by content hash, the
result line) and a ``<level>_summary.json`` (the headline numbers). Those are complete but scattered —
to see what a run *is* you open four files across two sub-directories, and none of them ties the levels
together or says what ``profiles.tsv`` holds.

This builds the missing single view. It is a **derived projection** of those per-level records, never a
source of truth: every ``zombi2 <level>`` run calls ``write_run_report``, which reads whatever records
the directory currently holds and rewrites ``run.zombi2`` from them — so running a new level, or
re-running one, just refreshes it, and there is no document the levels take turns appending to.

It reports, per level: the parameters, the inputs it was computed on, the result, and **every file it
wrote** (a per-family directory is summarised by count, not enumerated). Because it records the input
files by hash, it can also see when a downstream level was computed on an upstream one that has since
changed, and say so.
"""
from __future__ import annotations

import hashlib
import os
import shlex

from zombi2._runtime.summary import read_summary

RUN_REPORT_NAME = "run.zombi2"

#: The levels, in pipeline order, each grown into ``<run>/<level>/``. A ``joint`` run is not here: it
#: writes one report spanning ``species/`` and its driver (see ``_joint_section``).
_LEVELS = ("species", "genomes", "sequences", "traits")

#: The record files a level drops beside its data — not listed one-per-line among the outputs, but named
#: together on one line so nothing is undocumented.
_RECORD_SUFFIXES = (".log", "_summary.json")

#: Plumbing, not model — dropped from both the parameters line and the reproduce command: how the run
#: was driven or laid out (never what was simulated), or handled elsewhere (``seed`` is appended,
#: ``source`` becomes ``--from``). The ``.log`` keeps the complete set regardless.
_PLUMBING = frozenset({"command", "run", "seed", "params", "flat", "force", "quiet", "write",
                       "parallel", "stream", "source"})

#: Additionally hidden from the *display* only — inputs shown on the "computed on" line (``tip_fates``,
#: ``fasta``) and the ``max_lineages`` guard rail. The reproduce command keeps them, so a run that used
#: ``--tip-fates`` or a raised ``--max-lineages`` still reproduces from the printed command alone.
_DISPLAY_EXTRA = frozenset({"max_lineages", "fasta", "tip_fates"})

#: A readable order for the parameters that do show — model choice first, then rates, then the rest
#: (anything unlisted follows, alphabetically). Only affects display order.
_PARAM_ORDER = ("resolution", "kind", "model", "birth", "death", "n_extant", "total_time", "sampling",
                "fossils", "mass_extinction", "initial_families", "origination", "duplication",
                "transfer", "loss", "transfer_to", "self_transfer", "replacement", "max_family_size",
                "inversion", "chromosomes", "length", "divergence", "kappa", "frequencies",
                "substitution", "rate", "reverts_to", "pull", "states", "switch")

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
    "species_events.tsv": "one row per speciation/extinction (time, kind, parents, children)",
    "species_fates.tsv": "each tip's fate: extant / extinct / unsampled",
    # genomes/ (family resolution)
    "genomes.tsv": "gene content of every genome — one row per gene copy (lineage, family, copy)",
    "profiles.tsv": "family × extant-genome copy-number matrix",
    "initial_genome.tsv": "the genome present at t=0, that the initial families started from",
    "genome_events.tsv": "every duplication / transfer / loss / origination, per gene copy",
    "gene_trees": "gene trees, one per family — complete (all copies) and extant (survivors), Newick",
    # genomes/ (ordered resolution adds gene order and chromosome-level events)
    "gene_order.tsv": "the gene arrangement of every genome (lineage, chromosome, topology, position, strand, family, copy)",
    "chromosome_events.tsv": "the chromosome network's edges — initial / origination / speciation / fission / fusion / loss (time, kind, parents, children)",
    "rearrangement_events.tsv": "every inversion / transposition / translocation — the segment moved and where it went",
    # genomes/ (nucleotide resolution: the genome as ancestry blocks along a sequence)
    "blocks.tsv": "the ancestry blocks of every genome (chromosome, position, source, start, end, strand, copy, gene)",
    "block_events.tsv": "every block-level event along the genomes (duplication / transfer / loss / origination, per ancestral interval)",
    "genes.tsv": "the genes annotated on the initial genome (family, name, source, start, end, strand)",
    "bed": "per-genome gene annotations in BED format",
    "gff": "per-genome gene annotations in GFF format",
    # sequences/
    "clock_species_tree_complete.nwk": "the species tree rescaled to substitutions by the clock, all nodes",
    "clock_species_tree_extant.nwk": "the clock-rescaled tree, sampled tips only",
    "alignments": "one FASTA alignment per family (the extant sequences)",
    "phylograms": "gene trees with branch lengths in substitutions/site, one per family",
    # traits/
    "trait_tree.nwk": "the tree with every node's trait annotated ([&trait=…])",
    "trait_values.tsv": "the trait value (or state) at every node (tips, extinct lineages, internal)",
    "trait_events.tsv": "each discrete state change along a branch (time, lineage, from→to)",
    # written by genomes / traits when the input tree came from elsewhere and carries its own labels
    "names.tsv": "your tree's tip labels, mapped to ZOMBI's n<id> node ids",
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
        if key in ("level", "seed", "unit"):        # level/seed are in the header; unit is plumbing
            continue
        if key == "assembled_genomes" and value == 0:
            continue                                 # assembly is nucleotide-only; 0 here reads as failure
        if key == "empty_genomes" and value == 0:
            continue                                 # the healthy case; the line is only news when it fires
        if key == "mean_pairwise_identity":
            # `None` is the declared value when there are no alignments to compare (a run whose
            # genomes emptied, or one started with --initial-families 0): mean pairwise identity is
            # undefined, not zero, so there is no line to write. This used to render it regardless
            # and raise, which turned a legitimately empty run into a crash *after* every file had
            # been written.
            if value is not None:
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


#: How the report words each scope a rate can be written from. `Global` is not a "per" at all — its
#: total is one shared budget that scales with nothing — so it gets a phrase instead. The order here
#: is the order the rows print in, coarsest unit first.
_SCOPE_UNITS = {
    "Global": "for the whole tree",
    "PerLineage": "per lineage",
    "PerCopy": "per gene copy",
    "PerChromosome": "per chromosome",
    "PerSite": "per site",
}

#: What each rate slot is counted per when this run left the scope unwritten — the level's default,
#: which a bare number takes (SPEC §18). A written scope wins over it, because since scope overrides
#: landed the unit is a property of the **value** wherever one was written: ``--loss PerLineage(0.4)``
#: is a per-lineage budget however this table reads. Time is time⁻¹ throughout — tree time, in
#: whatever unit the tree is in, which is the run's own and is not calibrated to anything.
_DEFAULT_SCOPES = {
    "birth": "PerLineage", "death": "PerLineage", "fossils": "PerLineage",
    "origination": "PerLineage", "chromosome_origination": "PerLineage", "switch": "PerLineage",
    "duplication": "PerCopy", "transfer": "PerCopy", "loss": "PerCopy", "inversion": "PerCopy",
    "transposition": "PerCopy", "translocation": "PerCopy",
    "fission": "PerChromosome", "fusion": "PerChromosome", "chromosome_loss": "PerChromosome",
    "substitution": "PerSite",
}

#: The nucleotide resolution counts a gene event per lineage rather than per copy: the rate says how
#: often a lineage does the event and the extent says how much DNA it moves, so the number means the
#: same whatever the genome holds (SPEC §18, and the scopes ``genomes/nucleotide.py`` fills in).
_NUCLEOTIDE_SCOPES = dict.fromkeys(("duplication", "transfer", "loss", "inversion", "transposition",
                                    "translocation"), "PerLineage")


def _scope_of(slot: str, value: str, resolution: str | None) -> str | None:
    """The scope one rate was counted per in this run, as a scope class name — or ``None`` for a
    parameter that is not a rate, which has no unit to state.

    Every written rate starts from its scope (``PerLineage(0.4).scaled_by(…)``), and the log records
    rates in that written form, so the head of the value answers the question whenever the run
    answered it. Anything else — a bare number, a ``{'a->b': rate}`` dict, a k×k matrix — took the
    level's default, which is the slot's entry here."""
    head = value.split("(", 1)[0].strip()
    if head in _SCOPE_UNITS:
        return head
    if resolution == "nucleotide" and slot in _NUCLEOTIDE_SCOPES:
        return _NUCLEOTIDE_SCOPES[slot]
    return _DEFAULT_SCOPES.get(slot)


def _units_line(params: dict) -> list[str]:
    """One line naming what the rates in *this* run are counted per — only the slots it used, each
    under the scope this run actually gave it.

    A run directory is read long after the CLI help that answers this is to hand: "duplication 0.15"
    says nothing about per what, and a reviewer handed a folder could not tell whether a tree height
    of 5.391 was time or something else. The variance of a continuous trait is not a rate and is not
    listed here; nor is anything dimensionless (a probability, a fraction, a count)."""
    resolution = params.get("resolution")
    slots: dict[str, list[str]] = {}
    for slot, value in _visible_params(params, drop_zeros=True):
        scope = _scope_of(slot, value, resolution)
        if scope is not None:
            slots.setdefault(_SCOPE_UNITS[scope], []).append(_label(slot))
    rows = [(per, ", ".join(slots[per])) for per in _SCOPE_UNITS.values() if per in slots]
    if not rows:
        return []
    # An aligned block rather than one wrapped phrase: run the groups together and the closing
    # "per unit of tree time" reads as if it qualified only the last of them.
    width = max(len(per) for per, _ in rows)
    closing = "and all per unit of tree time (time⁻¹)" if len(rows) > 1 else \
              "per unit of tree time (time⁻¹)"        # "all" of one group reads oddly
    return _field("rate units", [f"{per:<{width}}   {names}" for per, names in rows] + [closing])


def _visible_params(params: dict, *, for_reproduce: bool = False,
                    drop_zeros: bool = False) -> list[tuple[str, str]]:
    """The parameters worth showing, in reading order — the model parameters, dropping plumbing
    (`_PLUMBING`) and anything left unset (``None`` / ``False`` / blank).

    ``substitution`` under ``--divergence`` is a special case: the log records the rate the run
    *solved to*, and that number will not go back into the flag, so it is shown as the clock's shape
    instead (`_clock_shape`), which is what the two flags compose as. A strict clock has no shape and
    is dropped — the solved number is what ``--divergence`` computes, so the flag alone rebuilds it.

    ``for_reproduce`` keeps the inputs and guards that display hides (`_DISPLAY_EXTRA`), so the printed
    command pins the run completely. ``drop_zeros`` hides rates left at zero — a process that is not
    happening — to declutter the display (the nucleotide model has a dozen such knobs); reproduce keeps
    them."""
    drop = set(_PLUMBING) if for_reproduce else set(_PLUMBING | _DISPLAY_EXTRA)
    shape = None
    if params.get("divergence") not in (None, "None", ""):
        shape = _clock_shape(params.get("substitution", ""))
        if shape is None:
            drop.add("substitution")
    items = [(k, v) for k, v in params.items() if k not in drop and v not in ("None", "False", "")
             and not (drop_zeros and v in ("0", "0.0"))]
    if shape is not None:
        items = [(k, shape if k == "substitution" else v) for k, v in items]
    items.sort(key=lambda kv: (_PARAM_ORDER.index(kv[0]) if kv[0] in _PARAM_ORDER else len(_PARAM_ORDER),
                               kv[0]))
    return items


def _clock_shape(written: str) -> str | None:
    """A substitution rate with its base taken off — ``PerSite().varying_among('lineages',
    LogNormal(0.0, 0.3))`` — or ``None`` when it carries no shape to keep.

    ``--divergence`` sets the *scale* of a clock and refuses a substitution rate that names a base of
    its own, so the solved rate the log records is the one form of this run's clock that cannot go
    back into the flag it came from. What was written was the shape, and the shape plus the
    divergence is what rebuilds the run: without it the printed command silently rebuilds a strict
    clock. A rate with no modifiers *is* only its base, so there is nothing to carry across.

    Rendering goes through `Rate.__repr__`, the one renderer for the written form, rather than by
    cutting the number out of the string. A value that does not parse as a rate is left alone: this
    runs while the report is being written, after every output file is already on disk, so an
    unreadable parameter must not be what turns a finished run into a traceback."""
    import dataclasses

    from zombi2.rates.parse import parse_rate
    from zombi2.rates.rate import Rate
    try:
        rate = parse_rate(written)
    except ValueError:            # every parser refusal, including an old log's retired spelling
        return None
    if not isinstance(rate, Rate) or rate.scope is None or not rate.modifiers:
        return None
    return repr(dataclasses.replace(rate, base=None))


def _params_phrase(params: dict) -> str:
    return " · ".join(f"{_label(k)} {v}" for k, v in _visible_params(params, drop_zeros=True))


def _given_flags(log: dict) -> set | None:
    """The long-option flags the user actually typed, from the recorded command line — or ``None`` when
    it is unavailable (an old log) or a ``--params`` file was used (then the reproduce lists every
    resolved parameter, so it stands alone without the file)."""
    cl = log.get("command_line")
    if not cl:
        return None
    tokens = cl.split()
    if any(t == "--params" or t.startswith("--params=") for t in tokens):
        return None
    return {t.split("=", 1)[0] for t in tokens if t.startswith("--")}


def _reproduce(command: str, run: str, params: dict, seed: str | None, source: str | None,
               given: set | None = None) -> str:
    """A clean, copy-pasteable command that reproduces this level — rebuilt from the resolved parameters
    (not the raw command line, whose shell quoting is already lost), so a rate expression comes back
    correctly quoted. ``given`` is the flags the user actually typed: when set, only those are emitted,
    so a run left on its defaults reproduces as the short bare command it was — re-running it fills the
    same defaults. When ``None``, every resolved parameter is shown (an old log, or a --params run)."""
    parts = [f"zombi2 {command} {shlex.quote(run)}"]
    if source:
        parts.append(f"--from {shlex.quote(source)}")
    for k, v in _visible_params(params, for_reproduce=True):
        if given is not None and _flag_name(k) not in given:
            continue                                         # a default the user did not type
        parts.append(_flag(k, v))
    if seed not in (None, "None"):
        parts.append(f"--seed {seed}")
    return " ".join(parts)


#: The options that take several values at once, so a recorded list goes back as several bare
#: arguments: ``--frequencies 0.25 0.25 0.25 0.25``. These are the options argparse declares with
#: ``nargs`` and the reproduce block prints — ``--write`` is plumbing (`_PLUMBING`) and never
#: reaches here. Every other parameter is one argument to one flag **even when its value looks like
#: a list**, and is quoted whole. Naming them is the point: a k×k ``--switch`` matrix is a single
#: expression, and flattening it produced four bare numbers argparse rejects as unrecognised.
_MULTI_VALUE = frozenset({"frequencies", "exchangeabilities"})

#: The options that are *repeated* rather than given several values — ``--mass-extinction TIME
#: FRACTION`` takes two values and can be given again for the next pulse, so the log's list of pairs
#: goes back as one flag per pair. Run together on one flag, argparse takes the first pair and calls
#: the rest unrecognised.
_REPEATED = frozenset({"mass_extinction"})


def _flag_name(key: str) -> str:
    """The long option a logged parameter came from: ``max_family_size`` → ``--max-family-size``."""
    return f"--{key.replace('_', '-')}"


def _flag(key: str, value: str) -> str:
    """One parameter as it goes back on the command line, flag and value together.

    ``True`` (a store-true flag) is the flag alone. A `_MULTI_VALUE` or `_REPEATED` option spreads its
    recorded list back over the command line as argparse wants to read it. Everything else is one
    shell-quoted argument, so a rate expression, a matrix, a path with spaces, or any shell
    metacharacter survives a copy-paste intact."""
    flag = _flag_name(key)
    if value == "True":
        return flag                                          # a store-true flag takes no value
    if key in _MULTI_VALUE or key in _REPEATED:
        groups = _value_groups(value)
        if groups is not None:
            joined = [" ".join(shlex.quote(v) for v in group) for group in groups]
            return " ".join(f"{flag} {g}" for g in joined) if key in _REPEATED else \
                f"{flag} {' '.join(joined)}"
    return f"{flag} {shlex.quote(value)}"


def _value_groups(value: str) -> list[list[str]] | None:
    """A recorded list as the groups of values it goes back as — one group per repetition of the flag
    — or ``None`` when there are none to spread: it is not a list, or an empty one. ``[0.25, 0.25]``
    is one group of two; ``[[2.0, 0.5], [4.0, 0.4]]`` is two groups of two."""
    if not (value.startswith("[") and value.endswith("]")):
        return None
    import ast
    try:
        parsed = ast.literal_eval(value)
    except (ValueError, SyntaxError):
        return None
    if not isinstance(parsed, list) or not parsed:
        return None
    if all(isinstance(x, (list, tuple)) for x in parsed):
        return [_flatten(x) for x in parsed]
    return [_flatten(parsed)]


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


def output_signposts(level_dir: str) -> list[tuple[str, str]]:
    """The data files a level wrote, as ``(name, one-line description)`` — a per-family directory folded
    to a count, record files (``.log`` / ``_summary.json``) left out. The end-of-command signpost and
    the report's OUTPUT FILES section both read this, so the terminal and ``run.zombi2`` never disagree."""
    return _list_outputs(level_dir)[0]


# ── staleness ───────────────────────────────────────────────────────────────────────────────────

def _stale_notes(run: str, sections: list[dict]) -> list[str]:
    """Warnings for any level computed on a file inside this run that has since changed — the input's
    recorded hash no longer matches the file on disk, so an upstream level was re-run without this one.

    Each input is located by its path *relative to the run directory as the recording command named it*
    and then re-resolved against ``run`` as it is now — so the check holds regardless of the working
    directory the report is regenerated from, or a move/rename of the whole run directory. (Resolving
    the recorded path against the current CWD instead would silently miss both.)"""
    notes = []
    for sec in sections:
        recorded_run = sec["log"].get("params", {}).get("run")
        for path, recorded_hash in sec["log"].get("inputs", []):
            rel = _intra_run_rel(path, recorded_run)
            if rel is None:
                continue                                    # an input from outside this run
            current = os.path.join(run, rel)
            if os.path.isfile(current) and _sha256(current) != recorded_hash:
                notes.append(f"⚠ {sec['title']} was computed on {rel}, which has changed since — "
                             f"re-run {sec['command']} so the run agrees with itself")
    return notes


def _intra_run_rel(input_path: str, recorded_run: str | None) -> str | None:
    """An input's path relative to the run directory as the recording command named it, or ``None`` if it
    lies outside the run (an external tree). ``os.path.relpath`` cancels the working directory shared by
    the two recorded paths, so the result is the same string whether they were recorded relative or
    absolute — which is what makes the staleness check independent of where it is later run."""
    if not recorded_run:
        return None
    try:
        rel = os.path.relpath(input_path, recorded_run)
    except ValueError:                                      # different drives on Windows
        return None
    return None if rel.startswith("..") or os.path.isabs(rel) else rel


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
        for level_dir, label in _level_dirs(run, level):
            log = _read_log(os.path.join(level_dir, f"{level}.log"))
            if log is None:
                continue                                 # this level was not run into this directory
        # the summary is the numbers; a level that writes none (nucleotide genomes today) still gets a
        # section from its log, falling back to the log's own result line.
            summary_path = os.path.join(level_dir, _SUMMARY_FILE.get(level, f"{level}_summary.json"))
            summary = read_summary(summary_path) if os.path.isfile(summary_path) else {}
            sections.append(_section(run, label, level, level_dir, log, summary))
    return sections


def _level_dirs(run: str, level: str):
    """``(directory, section title)`` for each run of ``level`` in this run directory.

    One for every level but traits, which a run may hold several of — ``--name`` puts each in
    ``traits/<name>/`` so that one can drive another. Each named trait gets its own section, titled
    with its name, so the report says which is which; an unnamed one keeps ``traits/``."""
    top = os.path.join(run, level)
    if level != "traits":
        yield top, level.upper()
        return
    if os.path.isfile(os.path.join(top, "traits.log")):
        yield top, "TRAITS"
    for name in sorted(os.listdir(top)) if os.path.isdir(top) else []:
        nested = os.path.join(top, name)
        if os.path.isfile(os.path.join(nested, "traits.log")):
            yield nested, f"TRAITS ({name})"


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
    """A single sentence for the top of the report, synthesised from the sections present. Every summary
    lookup is guarded — a level that wrote no summary (nucleotide genomes) still contributes a phrase."""
    by = {s["command"]: s["summary"] for s in sections}
    if "joint" in by:
        n = by["joint"].get("species", {}).get("tips", {}).get("extant", "?")
        return (f"A species tree of {n} tips, its diversification driven by "
                f"a {by['joint'].get('driver', 'trait')}.")
    bits = []
    if "species" in by:
        n = by["species"].get("tips", {}).get("extant")
        bits.append(f"A species tree of {n} sampled tips" if n is not None else "A species tree")
    if "genomes" in by:
        surv = by["genomes"].get("families", {}).get("surviving")
        bits.append(f"{surv} gene families evolved along it" if (surv is not None and bits)
                    else ("gene families evolved along it" if bits else "gene families"))
    if "sequences" in by:
        bits.append("sequences down each gene tree")
    if "traits" in by:
        bits.append(f"a {by['traits'].get('kind', '')} trait".strip())
    return (", on which ".join(bits[:2]) + ("; " + "; ".join(bits[2:]) if bits[2:] else "") + "."
            if bits else "An empty run.")


#: Where a reader of a written run finds the software. In the report because the report is what
#: travels with the data — a run directory is read years later, on a machine that never had ZOMBI2.
_HOME_URL = "https://github.com/AADavin/zombi2"


def build_run_report(run: str) -> str | None:
    """The text of ``run.zombi2`` for a run directory, or ``None`` if it holds no level records yet
    (e.g. a ``--flat`` run, whose levels share one directory and cannot be told apart)."""
    sections = _collect_sections(run)
    if not sections:
        return None

    from zombi2 import __version__
    env: dict = next((s["log"] for s in sections if s["log"].get("zombi2_version")), {})
    lines = [f"ZOMBI2 {env.get('zombi2_version', __version__)}  ·  run report", _HEAVY,
             f"  {_one_liner(sections)}"]
    built = " · ".join(p for p in (
        f"zombi2 {env['zombi2_version']}" if env.get("zombi2_version") else "",
        f"numpy {env['numpy_version']}" if env.get("numpy_version") else "",
        f"python {env['python_version']}" if env.get("python_version") else "",
        env.get("platform", "")) if p)
    if built:
        lines.append(f"  built with {built}")
    # Where the tool that made this came from. A folder handed on — deposited, emailed, inherited —
    # outlives the environment that produced it, and the version string alone is only a head start:
    # a research assistant who reconstructed a whole run from this file could not say where ZOMBI2
    # lived or what to install. Two lines, and the run report answers it without anyone asking.
    version = env.get("zombi2_version", __version__)
    lines.append(f"  {_HOME_URL}")
    lines.append(f"  reinstall it with:  pip install zombi2=={version}")
    for note in _stale_notes(run, sections):
        lines.append(f"  {note}")
    lines.append("")

    for i, sec in enumerate(sections, 1):
        lines += _render_section(run, sec, i, len(sections))

    lines += ["TO REPRODUCE", _RULE]
    if built:
        # the RNG stream is numpy's, so byte-identical output is only guaranteed on the same versions —
        # say so, rather than let a run silently fail to reproduce on a different machine.
        lines.append(f"  # recorded under {built} — byte-identical output requires the same environment")
    for sec in _runnable_order(run, sections):
        seed = sec["log"].get("params", {}).get("seed") or sec["summary"].get("seed")
        lines.append("  " + _reproduce(sec["command"], run, sec["log"].get("params", {}),
                                       seed, _source_of(sec), _given_flags(sec["log"])))
    lines.append("")
    return "\n".join(lines)


def _runnable_order(run: str, sections: list[dict]) -> list[dict]:
    """The sections in an order that actually **runs**, which is not the order they are read in.

    The report lists levels in pipeline order — species, genomes, sequences, traits — because that is
    how a reader wants them. TO REPRODUCE cannot use that order, for two reasons at once. A level
    whose rate is conditioned on a trait must run *after* the trait that writes the file it reads,
    and traits come last; but promoting that trait also drags the level below it, and a sequences run
    replays the genomes run it must still come after. Copy-pasted, the block failed on the line that
    read a file nothing had written yet — in the block the command line tells you to open first.

    Both dependencies are already recorded, in one place and at file granularity: every level's log
    lists the files it read (the same ``input`` rows the staleness guard hashes). So a section
    depends on whichever section wrote each of its inputs, and this is a topological sort over that.
    Reading the inputs rather than the ``conditioned_on`` marker is what makes it complete — the
    marker names levels, so it cannot tell one named trait from another, and it says nothing at all
    about the ordinary pipeline edges.

    It is a stable sort, not a reshuffle: whenever pipeline order already runs it is what comes out,
    and the ordering moves only what would otherwise not run.
    """
    # Every directory each section wrote into, relative to the run — a joint section owns its driver's
    # directory as well as species/, having written both.
    owners = [(k, os.path.relpath(d, run)) for k, sec in enumerate(sections)
              for d in (sec["dir"], *sec.get("extra_dirs", [])) if d]
    needs: dict[int, set[int]] = {}
    for k, sec in enumerate(sections):
        recorded_run = sec["log"].get("params", {}).get("run")
        for path, _digest in sec["log"].get("inputs", []):
            rel = _intra_run_rel(path, recorded_run)
            if rel is None:
                continue                    # an input from outside the run: nothing here writes it
            # the longest matching directory, so a named trait's own beats the traits/ it sits under
            wrote = [o for o in owners if _inside(rel, o[1])]
            if wrote:
                j = max(wrote, key=lambda o: len(o[1]))[0]
                if j != k:                  # a level never blocks on its own directory
                    needs.setdefault(k, set()).add(j)

    ordered: list[dict] = []
    placed: set[int] = set()
    while len(placed) < len(sections):
        # the first section, in reading order, whose inputs are all written by sections already out
        for k, sec in enumerate(sections):
            if k not in placed and needs.get(k, set()) <= placed:
                ordered.append(sec)
                placed.add(k)
                break
        else:                               # a cycle: emit the rest as they came rather than hang
            ordered.extend(s for k, s in enumerate(sections) if k not in placed)
            break
    return ordered


def _inside(rel: str, directory: str) -> bool:
    """Whether the run-relative path ``rel`` is that directory or something under it."""
    return rel == directory or rel.startswith(directory + os.sep)


def _source_of(section: dict) -> str | None:
    src = section["log"].get("params", {}).get("source")
    return src if src not in (None, "None", "") else None


def _field(label: str, contents: list[str]) -> list[str]:
    """A labelled, aligned block: the label on the first content line, continuations blank, everything
    lined up at the same column."""
    return [f"  {(label if k == 0 else ''):<11}  {c}" for k, c in enumerate(contents)]


def _written(log: dict) -> set:
    """The output tokens a level's log records having written (its effective ``--write`` list)."""
    raw = log.get("params", {}).get("write", "")
    if isinstance(raw, str) and raw.startswith("["):
        import ast
        try:
            return set(ast.literal_eval(raw))
        except (ValueError, SyntaxError):
            pass
    return set()


def _render_section(run: str, sec: dict, i: int, n: int) -> list[str]:
    log, summary = sec["log"], sec["summary"]
    # a count for a file the run did not write reads as a dangling reference: the ancestral-sequence
    # reconstruction is computed always but written only on --write ancestral, so drop its count unless
    # the ancestral output is actually on disk.
    if summary.get("ancestral_sequences") and "ancestral" not in _written(log):
        summary = {k: v for k, v in summary.items() if k != "ancestral_sequences"}
    seed = summary.get("seed", log.get("params", {}).get("seed", "?"))
    tag = f"level {i}/{n}"
    lines = [f"{sec['title']}{' ' * max(1, 80 - len(sec['title']) - len(tag))}{tag}", _RULE]

    stamp = " · ".join(p for p in (f"seed {seed}", log.get("timestamp", "").replace("T", " ")) if p)
    if stamp:
        lines.append(f"  {stamp}")
    params = _params_phrase(log.get("params", {}))
    if params:
        lines += _field("parameters", _wrap(params))
        lines += _units_line(log.get("params", {}))
    recorded_run = log.get("params", {}).get("run")
    lines += _field("computed on", [f"{_intra_run_rel(p, recorded_run) or p}  sha256 {sha[:8]}…"
                                    for p, sha in log.get("inputs", [])])

    if sec["command"] == "joint":
        stats = _stat_lines(_flatten_joint(summary))
    elif summary:
        stats = _stat_lines(summary)
    else:                                                # a level that wrote no summary: its result line
        stats = [s for s in (log.get("result", ""),) if s]
    lines += _field("result", stats)

    data, records = _gather_outputs(sec)
    file_lines = [f"{name:<28}  {gloss}".rstrip() for name, gloss in data]
    if records:
        # the .log is always here; only some resolutions also drop a _summary.json — don't advertise
        # machine-readable stats for one (ordered/nucleotide genomes) that wrote only the log.
        has_stats = any(r.endswith("_summary.json") for r in records)
        gloss = "run parameters + machine-readable stats" if has_stats else "run parameters"
        file_lines.append(f"records: {', '.join(records)}  ({gloss})")
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


__all__ = ["RUN_REPORT_NAME", "build_run_report", "write_run_report", "output_signposts"]
