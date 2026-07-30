"""``zombi2 genomes`` — evolve gene families along a species tree.

``--resolution`` picks the model: ``family`` (the D/T/L/O gene-family core,
`zombi2.genomes.simulate_genomes_family()`), ``ordered`` (genes with a position and
orientation on chromosomes — segmental rearrangements and the chromosome tier,
`simulate_genomes_ordered()`), or ``nucleotide`` (the genome as a nucleotide
sequence of ancestry blocks, with declared indivisible genes and intergenic spacer,
`simulate_genomes_nucleotide()`). Long options are the API keyword names, and
every rate takes the written form (SPEC §5): a bare number on its natural scope, or the same
``scope(base) × modifiers`` expression the Python API takes — ``--loss "0.25 * OnTime({0: 1.0, 3:
2.0})"``. The nucleotide engine wires ``OnTime`` and ``DrivenBy``, so any other modifier is rejected
there rather than silently ignored."""
from __future__ import annotations

import argparse
import os
import time

from zombi2.genomes import (WIRED_MODIFIERS, simulate_genomes_nucleotide, simulate_genomes_ordered,
                            simulate_genomes_family)
from zombi2.genomes.nucleotide import WIRED_MODIFIERS as _NUC_WIRED
from zombi2.rates.parse import parse_rate
from zombi2.rates.scope import Global, PerLineage
from zombi2.tree import node_label, read_newick
from zombi2.cli.framework import (resolve_seed, _add_flat_arg, _add_force_arg, _add_quiet_arg, _add_parallel_arg,
                                  _add_from_arg, _add_params_arg, _add_run_arg, _rate, _rates_help,
                                  _read_tip_fates, _write_params_log, check_stale_downstream,
                                  clear_stale_downstream, conditioned_levels, default_outputs,
                                  defaults_used, guidance, input_digests, level_dir,
                                  parallel_from_args, record_conditioning, resolve_tree,
                                  sibling_fates, warn, warn_if_fates_were_inferred)

#: the RATES block for ``zombi2 genomes -h``, built from the level's own declaration
RATES_HELP = _rates_help(
    WIRED_MODIFIERS, "--loss",
    note="Each rate keeps its natural scope here (D/T/L per copy, origination per lineage), so "
         "there is no scope wrapper to write. DrivenBy is wired for all four gene-family rates "
         "(on --transfer it drives how often a lineage DONATES); --transfer-to takes the same "
         "DrivenBy, on its own, as a recipient weight. --resolution ordered wires OnTime and "
         "ByFamily (the weight lands on the segment an event covers, not on the gene it started "
         "from); --resolution nucleotide wires OnTime and DrivenBy — every rate there, rearrangements "
         "included.")

# the write vocabularies, mirroring each Result.write (there is no exported constant to import)
_FAMILY_OUTPUTS = ("events", "profiles", "genomes", "initial_genome", "gene_trees",
                   "summary")
_ORDERED_OUTPUTS = ("events", "profiles", "gene_order", "initial_genome", "gene_trees",
                    "chromosome_events", "summary")
_NUCLEOTIDE_OUTPUTS = ("events", "genes", "blocks", "initial_genome", "gene_trees",
                       "chromosome_events", "gff", "bed")
_OUTPUTS = {"family": _FAMILY_OUTPUTS, "ordered": _ORDERED_OUTPUTS,
            "nucleotide": _NUCLEOTIDE_OUTPUTS}

# knobs that need a *structured* genome — (attribute, default) pairs — rejected under the family resolution
_STRUCTURED_ONLY = (
    ("inversion", 0.0), ("transposition", 0.0), ("translocation", 0.0),
    ("chromosomes", 1), ("topology", "circular"),
    ("fission", 0.0), ("fusion", 0.0),
    ("chromosome_origination", 0.0), ("chromosome_loss", 0.0),
    ("inversion_probability", 0.0),
)

# knobs only the nucleotide engine has — rejected under the family and ordered resolutions
_NUCLEOTIDE_ONLY = (
    ("root_length", 1000), ("genes", 0), ("gene_length", 100), ("gff", None),
    ("trim_overlaps", False),
    ("inversion_extent", 50.0), ("transposition_extent", 50.0), ("translocation_extent", 50.0),
    ("loss_extent", 50.0), ("duplication_extent", 50.0), ("transfer_extent", 50.0),
    ("origination_extent", 50.0),
)

# A genome the command starts with. The library function defaults to 0 — an explicit caller says what
# it wants — but a bare `zombi2 genomes -t tree.nwk` should hand back a genome rather than 100 empty
# ones, and origination stays 0 so nothing arrives that was not asked for. The run log records the
# resolved value, so a run is never ambiguous about which it used.
_DEFAULT_INITIAL_FAMILIES = 100

# knobs the nucleotide engine does not have — it starts from a sequence, not from a family count,
# and its transfers are always additive. Paired with the default, so leaving the flag alone is not
# mistaken for setting it.
#: The per-genome family cap, the same default the two engines carry: how many copies of ONE family
#: one genome may hold. Duplication compounds, so a run is bounded unless you ask otherwise;
#: `--max-family-size none` is that ask. Ten is a genome-shaped number — most families in a real
#: genome are single-copy and the big ones run to tens — so the default says what an ordinary genome
#: looks like rather than merely catching an explosion. Raise it, or lift it with `none`, when the
#: model is about the large families. (It reads as `PerLineage(10)` did, which is the point: that
#: spelling resolved to 10 × the node count of the species tree — 1470 on a 147-node tree — so the
#: number here was never the number that applied.)
_DEFAULT_MAX_FAMILY_SIZE = 10

_NOT_IN_NUCLEOTIDE = (("initial_families", None), ("replacement", False),
                      ("max_family_size", _DEFAULT_MAX_FAMILY_SIZE), ("family_speed", None))


def _add_genomes_args(p: argparse.ArgumentParser) -> None:
    _add_run_arg(p, "genomes evolve along the species tree it already holds")
    g = p.add_argument_group("general")
    _add_params_arg(g)
    _add_from_arg(g, "the species tree — a Newick file, or another run's directory")
    g.add_argument("--resolution", choices=("family", "ordered", "nucleotide"),
                   default="family", metavar="RESOLUTION",
                   help="family (gene-family counts, default), ordered (genes positioned on "
                        "chromosomes, with rearrangements), or nucleotide (the genome as a "
                        "sequence of ancestry blocks, with indivisible genes)")
    g.add_argument("--seed", type=int, default=None, metavar="N",
                   help="RNG seed for reproducibility")
    g.add_argument("--tip-fates", metavar="FILE", dest="tip_fates",
                   help="[external non-ultrametric trees] a TSV 'tip_name<TAB>extant|extinct|unsampled' "
                        "declaring each tip's fate; required when the input tree is not ultrametric "
                        "(ZOMBI won't guess extinct lineages from early-sampled tips). A species run's "
                        "own species_fates.tsv is in this format and can be passed directly.")

    g = p.add_argument_group("gene-family events (D/T/L/O)", "rates on their natural scope — see RATES below")
    g.add_argument("--duplication", type=_rate, default=None, metavar="RATE",
                   help="gene duplication rate (per copy)")
    g.add_argument("--transfer", type=_rate, default=None, metavar="RATE",
                   help="horizontal transfer rate (per copy)")
    g.add_argument("--loss", type=_rate, default=None, metavar="RATE",
                   help="gene loss rate (per copy)")
    g.add_argument("--origination", type=_rate, default=None, metavar="RATE",
                   help="new-family origination rate (per lineage)")

    g = p.add_argument_group("transfer & content")
    g.add_argument("--transfer-to", type=_transfer_to, default="uniform",
                   metavar="RULE", dest="transfer_to",
                   help="recipient rule for a transfer: uniform (any contemporaneous lineage, "
                        "default), distance (closer relatives likelier), or a DrivenBy weight — "
                        "\"DrivenBy('trait_events.tsv', {'competent': 2.0, 'normal': 1.0})\" — "
                        "which redistributes transfers without changing how many there are "
                        "(family resolution only)")
    g.add_argument("--replacement", action="store_true",
                   help="a transfer overwrites a homologous copy in the recipient (replacing HGT)")
    g.add_argument("--self-transfer", action="store_true", dest="self_transfer",
                   help="allow a lineage to transfer to itself")
    g.add_argument("--initial-families", type=int, default=None, metavar="N",
                   dest="initial_families",
                   help=f"number of gene families the root genome starts with (default "
                        f"{_DEFAULT_INITIAL_FAMILIES}); 0 starts empty, so every family must then "
                        f"arrive by --origination")
    g.add_argument("--max-family-size", type=_family_cap, default=_DEFAULT_MAX_FAMILY_SIZE,
                   metavar="CAP", dest="max_family_size",
                   help=f"cap on how many copies of one family a single genome may hold — a whole "
                        f"number (default {_DEFAULT_MAX_FAMILY_SIZE}), or 'none' to remove it. "
                        f"Duplication compounds, so a run is bounded unless you ask otherwise. When "
                        f"the cap does bite, the duplication or arriving transfer is dropped rather "
                        f"than retried, so realised rates in a family at the cap sit below the ones "
                        f"you declared — lower it deliberately, and use 'none' to measure rates")
    g.add_argument("--family-speed", type=_family_speed, default=None, metavar="DRAW",
                   dest="family_speed",
                   help="one per-family tempo scaling every rate that family has, as a ByFamily draw "
                        "— \"ByFamily(spread=0.5)\" — so a fast family is fast at everything (a "
                        "ByFamily on a single rate varies that rate alone)")

    g = p.add_argument_group("structured genome", "only with --resolution ordered or nucleotide")
    g.add_argument("--inversion", type=_rate, default=0.0, metavar="RATE",
                   help="segmental inversion rate (per copy)")
    g.add_argument("--transposition", type=_rate, default=0.0, metavar="RATE",
                   help="segmental transposition rate — move a run within a chromosome "
                        "(per copy)")
    g.add_argument("--translocation", type=_rate, default=0.0, metavar="RATE",
                   help="segmental translocation rate — move a run to another chromosome (per copy)")
    g.add_argument("--chromosomes", type=int, default=1, metavar="N",
                   help="number of chromosomes at the origin (default 1)")
    g.add_argument("--topology", choices=("circular", "linear"), default="circular", metavar="TOPO",
                   help="chromosome topology (default circular) — a segmental run wraps past the "
                        "origin on a circular chromosome, stops at the end on a linear one")
    g.add_argument("--fission", type=_rate, default=0.0, metavar="RATE",
                   help="chromosome fission rate — split one in two (per chromosome)")
    g.add_argument("--fusion", type=_rate, default=0.0, metavar="RATE",
                   help="chromosome fusion rate — merge two into one (per chromosome)")
    g.add_argument("--chromosome-origination", type=_rate, default=0.0, metavar="RATE",
                   dest="chromosome_origination",
                   help="new-chromosome origination rate — a de-novo plasmid (per lineage)")
    g.add_argument("--chromosome-loss", type=_rate, default=0.0, metavar="RATE",
                   dest="chromosome_loss",
                   help="whole-chromosome loss rate, never the last one (per chromosome)")
    g.add_argument("--inversion-probability", type=float, default=0.0, metavar="P",
                   dest="inversion_probability",
                   help="probability a transposed/translocated block lands inverted (default 0)")

    g = p.add_argument_group("nucleotide genome", "only with --resolution nucleotide")
    g.add_argument("--root-length", type=int, default=1000, metavar="BP", dest="root_length",
                   help="length in bp of each initial replicon (default 1000)")
    g.add_argument("--genes", type=int, default=0, metavar="N",
                   help="number of evenly-spaced genes to declare on each initial replicon (default 0 "
                        "— an all-intergenic genome). Use --gff instead to declare real ones")
    g.add_argument("--gene-length", type=int, default=100, metavar="BP", dest="gene_length",
                   help="length in bp of each evenly-spaced gene (default 100)")
    g.add_argument("--gff", metavar="FILE",
                   help="a GFF3 declaring the initial genome's replicons and genes at exact "
                        "coordinates — the 'start from a real genome' path (excludes --genes)")
    g.add_argument("--trim-overlaps", action="store_true", dest="trim_overlaps",
                   help="[--gff] shorten overlapping gene annotations instead of refusing the file")
    g.add_argument("--fasta", metavar="FILE",
                   help="[--gff] the initial genome's DNA — one >seqid record per GFF ##sequence-region, "
                        "each exactly its declared length. The sequence level then founds every block "
                        "from this DNA, so an assembled genome descends from exactly what you supply; "
                        "without it the run is pure ancestry and letters come from the model")
    for knob, what in (("inversion", "inverted"), ("transposition", "moved within a chromosome"),
                       ("translocation", "moved to another chromosome"), ("loss", "deleted"),
                       ("duplication", "copied in tandem"), ("transfer", "copied to a recipient"),
                       ("origination", "laid down as new material")):
        g.add_argument(f"--{knob}-extent", type=float, default=50.0, metavar="BP",
                       dest=f"{knob}_extent",
                       help=f"mean bp {what} per event — the extent (geometric, default 50)")

    g = p.add_argument_group("outputs")
    g.add_argument("--write", nargs="+", choices=sorted({o for v in _OUTPUTS.values() for o in v}),
                   default=None, metavar="PART",
                   help="which outputs to write (default: each resolution's own, which is all of "
                        "them). family: events, profiles, genomes, initial_genome, gene_trees. "
                        "ordered: those with gene_order for genomes, plus chromosome_events. "
                        "nucleotide: events, genes, blocks, initial_genome, gene_trees, "
                        "chromosome_events, gff, bed. "
                        "'events' is the whole history in one time-ordered "
                        "table — the genealogy, where each event happened, and the rearrangements, "
                        "which used to be three files; 'genomes' is every node's gene content, "
                        "ancestors included, where "
                        "'profiles' counts only the extant tips; 'initial_genome' is the genome the "
                        "run started with, at the start of the root branch, which belongs to no node "
                        "and so gets its own file; 'gene_trees' writes one Newick per family, "
                        "complete and extant; 'gff' and 'bed' annotate every genome in its own "
                        "coordinates — the genes and the blocks respectively — named to join the "
                        "FASTA the sequence level writes.")
    _add_flat_arg(g)
    _add_parallel_arg(g)
    g.add_argument("--stream", action="store_true",
                   help="[family] write each gene family straight to disk instead of building the "
                        "whole run in memory — for a very large number of families, where the "
                        "in-memory result would not fit. Composes with --parallel and --write; the "
                        "same files are written and the disk is the handoff to the sequence level "
                        "(gene trees are grouped under gene_trees/ regardless of --flat). Like "
                        "--parallel it is a separate engine: it draws families in a different order, "
                        "so for the same seed the run it produces DIFFERS from a serial one (both "
                        "valid samples). Fix the mode alongside the seed to reproduce a run")
    _add_quiet_arg(g)
    _add_force_arg(g)


def _family_cap(text: str):
    """The argparse ``type`` for ``--max-family-size``: how many copies of one family a genome may
    hold, as a plain whole number, or ``none`` for no cap.

    No scope wrapper: the cap is counted inside one genome, so "per what?" has only one answer and
    saying it out loud only created room to say it wrong. ``PerLineage(N)`` used to be accepted and
    multiplied N by the size of the *species tree*, which made the real cap depend on something the
    user had not chosen; it is refused now, with the arithmetic, so a stale command line fails loudly
    instead of running at a different cap than it reads.

    There is a cap by default because duplication **compounds**: a family whose duplication rate sits
    above its loss rate multiplies without bound. ``none`` is how you ask for that on purpose.
    """
    s = text.strip()
    if s.lower() in ("none", "off"):
        return None
    try:
        return _positive_int(s)
    except ValueError:
        pass
    try:                                        # a scope wrapper: the old spelling, worth naming
        cap = parse_rate(s)
    except ValueError:
        cap = None
    if isinstance(cap, PerLineage):
        raise argparse.ArgumentTypeError(
            f"--max-family-size is a plain number of copies in one genome now, so write "
            f"--max-family-size {cap.base:g} (or 'none'). {text!r} used to mean {cap.base:g} × the "
            f"number of nodes in the species tree, which is why the cap you got depended on how big "
            f"the tree was")
    if isinstance(cap, Global):
        raise argparse.ArgumentTypeError(
            f"--max-family-size is a plain number of copies in one genome now, so write "
            f"--max-family-size {cap.base:g} (or 'none') — {text!r} meant exactly that")
    raise argparse.ArgumentTypeError(
        f"--max-family-size takes a whole number of copies (e.g. 20) or 'none' for no cap; "
        f"got {text!r}")


def _positive_int(s: str) -> int:
    """``s`` as a whole number ≥ 1, or ``ValueError``. A float spelling is refused rather than
    rounded: ``10`` against ``10.0`` is exactly the ambiguity this option has finished with."""
    n = int(s)                                  # raises on "10.0", which is the point
    if n < 1:
        raise ValueError(n)
    return n


def _family_speed(text: str):
    """The argparse ``type`` for ``--family-speed``: one per-family tempo scaling every rate a family
    has, written as a ``ByFamily`` draw — ``--family-speed "ByFamily(spread=0.5)"``.

    Parsed by the same ast-whitelist parser every rate flag uses, so the expression is the one you
    would write in Python and nothing is evaluated. It differs from a ``ByFamily`` on a single rate:
    there each rate varies on its own, here one draw moves them together.
    """
    from zombi2.rates.modifiers import ByFamily
    from zombi2.rates.parse import parse_rate

    try:
        value = parse_rate(text)
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"--family-speed: {e}") from None
    if not isinstance(value, ByFamily):
        raise argparse.ArgumentTypeError(
            f"--family-speed takes a ByFamily draw, e.g. \"ByFamily(spread=0.5)\"; got {text!r}")
    return value


def _transfer_to(text: str):
    """The argparse ``type`` for ``--transfer-to``: the recipient rule of a transfer.

    ``uniform`` and ``distance`` are the two named rules; anything else is read as the written form
    of a ``DrivenBy`` — ``--transfer-to "DrivenBy('trait_events.tsv', {'competent': 2.0})"`` — the
    **choice slot** of SPEC §5, where the mapping's numbers are per-candidate weights rather than
    rate multipliers. Parsed by the same ast-whitelist parser every rate flag uses, so the expression
    is the one you would write in Python and nothing is evaluated.
    """
    from zombi2.rates.modifiers import DrivenBy
    from zombi2.rates.parse import parse_rate

    if text in ("uniform", "distance"):
        return text
    value, detail = None, ""
    try:
        value = parse_rate(text)
    except ValueError as e:
        # only quote the parser when the text was meant as an expression; for a plain misspelt rule
        # ("uniforn") its "unknown name" reading is noise, and the flag's own list is the answer
        detail = f"\n{e}" if "(" in text else ""
    if not isinstance(value, DrivenBy):
        raise argparse.ArgumentTypeError(
            f"--transfer-to takes 'uniform', 'distance', or a DrivenBy recipient weight written on "
            f"its own — e.g. \"DrivenBy('trait_events.tsv', {{'competent': 2.0}})\" — got {text!r}. "
            f"The numbers there are weights over the candidate recipients, not a rate, so there is "
            f"no base number in front of it.{detail}")
    return value


#: Every rate flag with the value it holds when unset — what "the caller described no model at all"
#: means, and so when a bare run may fill in demonstration rates. The structural rates sit at 0.0
#: rather than None because absent genuinely means off for them: a genome with no inversion rate is
#: a genome that does not invert, whereas an unset D/T/L is a question this level has to answer.
_RATE_FLAGS = (("duplication", None), ("transfer", None), ("loss", None), ("origination", None),
               ("inversion", 0.0), ("transposition", 0.0), ("translocation", 0.0),
               ("fission", 0.0), ("fusion", 0.0),
               ("chromosome_origination", 0.0), ("chromosome_loss", 0.0))


def _stray(args, knobs) -> list[str]:
    """The flags in ``knobs`` the user actually set (their value differs from the default)."""
    return [f"--{attr.replace('_', '-')}" for attr, default in knobs
            if getattr(args, attr) != default]


def run(args, parser):
    # Fill the rate defaults first, so everything below validates the run that will actually happen
    # rather than a half-specified one. A bare `zombi2 genomes out/` runs, and shows what this level
    # is for: with every rate left at zero it would only inherit, which demonstrates nothing. But
    # defaulting a rate the caller left out *beside* ones they set would be a surprise —
    # `--duplication 0.3` alone plainly means no transfer — so this applies only when none was given.
    # Filled after the resolution checks above, so the nucleotide engine can still tell "not given"
    # from "given as the default" — passing --initial-families 100 to a nucleotide run used to slip
    # past the stray check and be silently ignored while the log recorded it as if it applied.
    if args.initial_families is None and args.resolution != "nucleotide":
        args.initial_families = _DEFAULT_INITIAL_FAMILIES

    _CORE = ("duplication", "transfer", "loss", "origination")
    # "gave no rate" means no rate of *any* kind, structural ones included: a run given --inversion
    # has had its model described, and silently adding gene turnover to it would be the surprise.
    if not _stray(args, _RATE_FLAGS):
        warn(defaults_used(args, duplication=0.2, transfer=0.1, loss=0.25, origination=0.5))
    else:
        for r in _CORE:                        # unset beside a set one means off, not defaulted
            if getattr(args, r) is None:
                setattr(args, r, 0.0)

    # a flag a resolution does not have is an error, never silently ignored — otherwise
    # `--inversion` under the family resolution, or `--initial-families` under nucleotide, would quietly
    # produce a run that is not the one asked for
    if args.resolution == "family":
        if stray := _stray(args, _STRUCTURED_ONLY):
            parser.error(f"these options need --resolution ordered or nucleotide: "
                         f"{', '.join(stray)} (the gene-family core has no chromosomes or positions)")
    else:
        for flag, given in (("--parallel", args.parallel is not None), ("--stream", args.stream)):
            if given:
                parser.error(f"{flag} applies to --resolution family only, where gene families are "
                             f"independent and evolve one per worker; the {args.resolution} resolution "
                             f"couples families by position (inversions, translocations), so it has no "
                             f"per-family engine")
    if args.resolution != "nucleotide":
        if stray := _stray(args, _NUCLEOTIDE_ONLY):
            parser.error(f"these options need --resolution nucleotide: {', '.join(stray)} "
                         f"(the {args.resolution} resolution counts genes, not base pairs)")
    else:
        if stray := _stray(args, _NOT_IN_NUCLEOTIDE):
            # each of these is absent for its own reason, so say which rather than give one blanket
            # explanation that fits some of them and not the rest
            why = {"--initial-families": "the genome is founded from a sequence, not from a family "
                                         "count — see --root-length / --genes / --gff",
                   "--replacement": "a nucleotide transfer is always additive",
                   "--max-family-size": "a quota counts copies of a family, and here an event takes "
                                        "an arc of DNA that may cover several families or none",
                   "--family-speed": "a per-family tempo has to reach the arc an event covers, which "
                                     "is per-family weighting this resolution does not wire"}
            parser.error("; ".join(f"the nucleotide resolution has no {f} ({why[f]})"
                                   for f in stray))
        if args.gff and args.genes:
            parser.error("pass either --gff or --genes, not both — a GFF already declares the genes")
        if args.fasta and not args.gff:
            parser.error("--fasta needs --gff: the FASTA's records are matched to the GFF's "
                         "replicons by id, so there is nothing to lay down without one")
        # the nucleotide engine holds each rate constant, so a modifier expression would be
        # accepted and then dropped; refuse it instead
        modulated = [f"--{n}" for n in ("duplication", "transfer", "loss", "origination", "inversion",
                                        "transposition", "translocation", "fission", "fusion")
                     if not isinstance(getattr(args, n), float)
                     and any(not isinstance(m, _NUC_WIRED) for m in getattr(args, n).modifiers)]
        if modulated:
            parser.error(f"--resolution nucleotide wires only "
                         f"{', '.join(w.__name__ for w in _NUC_WIRED)}, but "
                         f"{', '.join(modulated)} carries another modifier")

    vocab = _OUTPUTS[args.resolution]
    if args.write:
        bad = [o for o in args.write if o not in vocab]
        if bad:
            parser.error(f"--write {' '.join(bad)} not available for --resolution "
                         f"{args.resolution}; choose from: {', '.join(vocab)}")

    # refuse up front if re-running would orphan a later level already in the run (unless --force)
    resolve_seed(args)                      # a run must be reproducible from its own log
    check_stale_downstream(args, "genomes")

    tree_path = resolve_tree(args.source or args.run, is_run_dir=args.source is None)
    # an explicit --tip-fates wins; otherwise pick up the run's own species_fates.tsv so extinct and
    # unsampled tips are read from the record rather than guessed from tip depth
    tip_fates = _read_tip_fates(args.tip_fates) if args.tip_fates else sibling_fates(tree_path)
    try:
        with open(tree_path, encoding="utf-8") as f:
            tree, names = read_newick(f.read(), tip_fates=tip_fates)
        warn_if_fates_were_inferred(tree, args)
    except FileNotFoundError:
        raise FileNotFoundError(f"tree file not found: {tree_path}") from None

    common = dict(duplication=args.duplication, transfer=args.transfer, loss=args.loss,
                  origination=args.origination, transfer_to=args.transfer_to,
                  self_transfer=args.self_transfer, seed=args.seed)
    # the two knobs only the family-tier engines have — kept out of `common`, which the nucleotide
    # engine shares and which has neither
    family_knobs = dict(max_family_size=args.max_family_size, family_speed=args.family_speed)
    structured = dict(inversion=args.inversion, transposition=args.transposition,
                      translocation=args.translocation, chromosomes=args.chromosomes,
                      topology=args.topology, fission=args.fission, fusion=args.fusion,
                      chromosome_origination=args.chromosome_origination,
                      chromosome_loss=args.chromosome_loss,
                      inversion_probability=args.inversion_probability)

    clear_stale_downstream(args, "genomes")   # --force: drop the now-stale downstream (run succeeded)
    os.makedirs(args.run, exist_ok=True)
    out = level_dir(args.run, "genomes", args.flat)
    streaming = args.stream and args.resolution == "family"

    t0 = time.perf_counter()
    if args.resolution == "ordered":
        result = simulate_genomes_ordered(
            tree, replacement=args.replacement, initial_families=args.initial_families,
            progress=not args.quiet, **structured, **family_knobs, **common)
    elif args.resolution == "nucleotide":
        result = simulate_genomes_nucleotide(
            tree, root_length=args.root_length, genes=args.genes, gene_length=args.gene_length,
            gff=args.gff, fasta=args.fasta, trim_overlaps=args.trim_overlaps,
            inversion_extent=args.inversion_extent,
            transposition_extent=args.transposition_extent,
            translocation_extent=args.translocation_extent, loss_extent=args.loss_extent,
            duplication_extent=args.duplication_extent, transfer_extent=args.transfer_extent,
            origination_extent=args.origination_extent, progress=not args.quiet,
            **structured, **common)
    elif streaming:
        # each family written straight to disk (no whole run in memory) — the engine writes `out`
        # itself, so there is no result.write below; a StreamedRun handle comes back.
        result = simulate_genomes_family(
            tree, replacement=args.replacement, initial_families=args.initial_families,
            parallel=parallel_from_args(args, parser), stream_to=out,
            outputs=tuple(args.write) if args.write else None, progress=not args.quiet,
            **family_knobs, **common)
    else:
        result = simulate_genomes_family(
            tree, replacement=args.replacement, initial_families=args.initial_families,
            parallel=parallel_from_args(args, parser), progress=not args.quiet,
            **family_knobs, **common)
    dt = time.perf_counter() - t0

    # a genome run is on a fixed tree, so its complete tree is the input; a StreamedRun does not carry
    # one, so read it from `tree` there. The rest of the CLI's bookkeeping is identical either way.
    complete_tree = tree if streaming else result.complete_tree
    if streaming:
        # A streamed run wrote its own files as it went, so there is nothing to write here — but the
        # log still has to say which, in the same words as an in-memory run, and the engine already
        # resolved that (a bare --stream leaves `outputs=None` for it to default).
        wanted = result.outputs
    else:
        # The many-files-per-run outputs (gene trees, gff, bed) get a subdirectory apiece, and
        # `write` is where that is decided — so a run written from Python has this layout too, and
        # --flat is simply passed through rather than being a second layout the CLI knows about.
        # ...minus the species tree. A Result writes it so a directory written from Python is a
        # dataset on its own; a *run* already keeps one canonical copy under species/, shared by
        # every level, and a second under genomes/ would be two files for one fact.
        wanted = tuple(args.write) if args.write else tuple(
            o for o in default_outputs(result) if o != "species_tree")
        result.write(out, outputs=wanted, flat=args.flat)
    # The events index against the tree canonicalised to n<id> labels, so the run needs that exact
    # tree to be replayable. A run grown here already has it — `zombi2 species` wrote the identical
    # file — so only a run reading its tree from elsewhere (--from) needs a copy, and it goes where
    # a species run keeps one rather than under a second name for the same thing.
    species_dir = level_dir(args.run, "species", args.flat)
    canonical = os.path.join(species_dir, "species_complete.nwk")
    if not os.path.exists(canonical):
        with open(canonical, "w", encoding="utf-8") as f:
            f.write(complete_tree.to_newick() + "\n")
        # write the fate table beside the canonical tree too, so a later level on this run reads each
        # tip's fate from the record instead of guessing it from depth (matching a species run's output)
        fate_rows = ["lineage\tfate"] + [f"{node_label(n.id)}\t{n.fate}"
                                         for n in sorted(complete_tree.leaves(), key=lambda x: x.id)]
        with open(os.path.join(species_dir, "species_fates.tsv"), "w", encoding="utf-8") as f:
            f.write("\n".join(fate_rows) + "\n")
    if names:  # an external tree: map ZOMBI's n<id> back to the user's labels (join on profiles cols)
        rows = ["node\tname"] + [f"{node_label(i)}\t{lbl}" for i, lbl in sorted(names.items())]
        with open(os.path.join(out, "names.tsv"), "w", encoding="utf-8") as f:
            f.write("\n".join(rows) + "\n")

    if streaming:                               # a StreamedRun carries counts, not the run in memory
        summary = (f"{result.n_families} gene families, {result.n_events} events, streamed to disk "
                   f"(family)")
    elif args.resolution == "nucleotide":       # no phyletic profiles here: the unit is a base pair
        extant = [n.id for n in result.complete_tree.extant()]
        bp = sum(result.genomes[s].length for s in extant)
        summary = (f"{len(result.gene_spans)} genes and {bp} bp across {len(extant)} extant "
                   f"genomes (nucleotide)")
    else:
        n_families, n_species = result.profiles.shape
        summary = f"{n_families} gene families across {n_species} extant genomes ({args.resolution})"
    print(f"wrote {args.run}/ ({summary}) in {dt:.3g} s")
    guidance(args, f"genomes and gene trees under {out}/")
    if names:
        guidance(args, f"your tree's tip labels, mapped to ZOMBI's n<id>: {os.path.join(out, 'names.tsv')}")
    if not args.flat:                             # record which same-run levels drove a rate (if any),
        record_conditioning(out, conditioned_levels(   # so re-running one of them knows it orphans this
            args.run, (args.duplication, args.transfer, args.loss, args.origination, args.transfer_to)))
    # The log is this run's parameters, not the parser's: a family run has no --root-length and no
    # --inversion, and recording them at their defaults reads as though it had them and chose those
    # values. Each resolution's own gates already say which options belong to which.
    other = {"family": (*_STRUCTURED_ONLY, *_NUCLEOTIDE_ONLY),
             "ordered": _NUCLEOTIDE_ONLY,
             "nucleotide": _NOT_IN_NUCLEOTIDE}[args.resolution]
    _write_params_log(os.path.join(out, "genomes.log"), args, summary,
                      omit={attr for attr, _default in other},
                      effective={"write": list(wanted)},
                      inputs=input_digests(tree_path, args.tip_fates,
                                           os.path.join(os.path.dirname(tree_path),
                                                        "species_fates.tsv"),
                                           args.duplication, args.transfer, args.loss,
                                           args.origination, args.transfer_to))
    return 0
