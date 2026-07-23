"""``zombi2 genomes`` — evolve gene families along a species tree.

``--resolution`` picks the model: ``unordered`` (the D/T/L/O gene-family core,
:func:`zombi2.genomes.simulate_genomes_unordered`), ``ordered`` (genes with a position and
orientation on chromosomes — segmental rearrangements and the chromosome tier,
:func:`~zombi2.genomes.simulate_genomes_ordered`), or ``nucleotide`` (the genome as a nucleotide
sequence of ancestry blocks, with declared indivisible genes and intergenic spacer,
:func:`~zombi2.genomes.simulate_genomes_nucleotide`). Long options are the API keyword names, and
every rate takes the written form (SPEC §5): a bare number on its natural scope, or the same
``scope(base) × modifiers`` expression the Python API takes — ``--loss "0.25 * OnTime({0: 1.0, 3:
2.0})"``. The nucleotide engine takes **constant rates only**, so a modifier expression is rejected
there rather than silently ignored."""
from __future__ import annotations

import argparse
import os
import time

from zombi2.genomes import (WIRED_MODIFIERS, simulate_genomes_nucleotide, simulate_genomes_ordered,
                            simulate_genomes_unordered)
from zombi2.species import read_newick
from zombi2.cli.framework import (_add_flat_arg, _add_quiet_arg, _add_from_arg, _add_params_arg, _add_run_arg,
                                  _rate, _rates_help, _write_params_log, default_outputs, level_dir,
                                  resolve_tree)

#: the RATES block for ``zombi2 genomes -h``, built from the level's own declaration
RATES_HELP = _rates_help(
    WIRED_MODIFIERS, "--loss",
    note="Each rate keeps its natural scope here (D/T/L per copy, origination per lineage), so "
         "there is no scope wrapper to write. DrivenBy is wired for all four gene-family rates "
         "(on --transfer it drives how often a lineage DONATES); --transfer-to takes the same "
         "DrivenBy, on its own, as a recipient weight. --resolution ordered wires OnTime only; "
         "--resolution nucleotide takes constant rates only.")

# the write vocabularies, mirroring each Result.write (there is no exported constant to import)
_UNORDERED_OUTPUTS = ("events", "profiles", "genomes", "initial_genome", "gene_trees")
_ORDERED_OUTPUTS = ("events", "profiles", "gene_order", "initial_genome", "gene_trees",
                    "chromosome_events")
_NUCLEOTIDE_OUTPUTS = ("events", "genes", "blocks", "initial_genome", "gene_trees",
                       "chromosome_events", "gff", "bed")
_OUTPUTS = {"unordered": _UNORDERED_OUTPUTS, "ordered": _ORDERED_OUTPUTS,
            "nucleotide": _NUCLEOTIDE_OUTPUTS}

# knobs that need a *structured* genome — (attribute, default) pairs — rejected under unordered
_STRUCTURED_ONLY = (
    ("inversion", 0.0), ("transposition", 0.0), ("translocation", 0.0),
    ("chromosomes", 1), ("topology", "circular"),
    ("fission", 0.0), ("fusion", 0.0),
    ("chromosome_origination", 0.0), ("chromosome_loss", 0.0),
    ("inversion_probability", 0.0),
)

# knobs only the nucleotide engine has — rejected under unordered and ordered
_NUCLEOTIDE_ONLY = (
    ("root_length", 1000), ("genes", 0), ("gene_length", 100), ("gff", None),
    ("trim_overlaps", False),
    ("inversion_length", 50.0), ("transposition_length", 50.0), ("translocation_length", 50.0),
    ("loss_length", 50.0), ("duplication_length", 50.0), ("transfer_length", 50.0),
    ("origination_length", 50.0),
)

# A genome the command starts with. The library function defaults to 0 — an explicit caller says what
# it wants — but a bare `zombi2 genomes -t tree.nwk` should hand back a genome rather than 100 empty
# ones, and origination stays 0 so nothing arrives that was not asked for. The run log records the
# resolved value, so a run is never ambiguous about which it used.
_DEFAULT_INITIAL_FAMILIES = 100

# knobs the nucleotide engine does not have — it seeds from a sequence, not from a family count,
# and its transfers are always additive. Paired with the default, so leaving the flag alone is not
# mistaken for setting it.
_NOT_IN_NUCLEOTIDE = (("initial_families", _DEFAULT_INITIAL_FAMILIES), ("replacement", False))


def _add_genomes_args(p: argparse.ArgumentParser) -> None:
    _add_run_arg(p, "genomes evolve along the species tree it already holds")
    g = p.add_argument_group("general")
    _add_params_arg(g)
    _add_from_arg(g, "the species tree — a Newick file, or another run's directory")
    g.add_argument("--resolution", choices=("unordered", "ordered", "nucleotide"),
                   default="unordered", metavar="RESOLUTION",
                   help="unordered (gene-family counts, default), ordered (genes positioned on "
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
    g.add_argument("--duplication", type=_rate, default=0.0, metavar="RATE",
                   help="gene duplication rate (per copy)")
    g.add_argument("--transfer", type=_rate, default=0.0, metavar="RATE",
                   help="horizontal transfer rate (per copy)")
    g.add_argument("--loss", type=_rate, default=0.0, metavar="RATE",
                   help="gene loss rate (per copy)")
    g.add_argument("--origination", type=_rate, default=0.0, metavar="RATE",
                   help="new-family origination rate (per lineage)")

    g = p.add_argument_group("transfer & content")
    g.add_argument("--transfer-to", type=_transfer_to, default="uniform",
                   metavar="RULE", dest="transfer_to",
                   help="recipient rule for a transfer: uniform (any contemporaneous lineage, "
                        "default), distance (closer relatives likelier), or a DrivenBy weight — "
                        "\"DrivenBy('trait_events.tsv', {'competent': 2.0, 'normal': 1.0})\" — "
                        "which redistributes transfers without changing how many there are "
                        "(unordered only)")
    g.add_argument("--replacement", action="store_true",
                   help="a transfer overwrites a homologous copy in the recipient (replacing HGT)")
    g.add_argument("--self-transfer", action="store_true", dest="self_transfer",
                   help="allow a lineage to transfer to itself")
    g.add_argument("--initial-families", type=int, default=_DEFAULT_INITIAL_FAMILIES, metavar="N",
                   dest="initial_families",
                   help=f"number of gene families the root genome starts with (default "
                        f"{_DEFAULT_INITIAL_FAMILIES}); 0 starts empty, so every family must then "
                        f"arrive by --origination")

    g = p.add_argument_group("structured genome", "only with --resolution ordered or nucleotide")
    g.add_argument("--inversion", type=_rate, default=0.0, metavar="RATE",
                   help="segmental inversion rate (per copy)")
    g.add_argument("--transposition", type=_rate, default=0.0, metavar="RATE",
                   help="segmental transposition rate — move a run within a chromosome "
                        "(per copy)")
    g.add_argument("--translocation", type=_rate, default=0.0, metavar="RATE",
                   help="segmental translocation rate — move a run to another chromosome (per copy)")
    g.add_argument("--chromosomes", type=int, default=1, metavar="N",
                   help="number of chromosomes at the crown (default 1)")
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
                   help="length in bp of each seed replicon (default 1000)")
    g.add_argument("--genes", type=int, default=0, metavar="N",
                   help="number of evenly-spaced genes to declare on each seed replicon (default 0 "
                        "— an all-intergenic genome). Use --gff instead to declare real ones")
    g.add_argument("--gene-length", type=int, default=100, metavar="BP", dest="gene_length",
                   help="length in bp of each evenly-spaced gene (default 100)")
    g.add_argument("--gff", metavar="FILE",
                   help="a GFF3 declaring the seed genome's replicons and genes at exact "
                        "coordinates — the 'start from a real genome' path (excludes --genes)")
    g.add_argument("--trim-overlaps", action="store_true", dest="trim_overlaps",
                   help="[--gff] shorten overlapping gene annotations instead of refusing the file")
    for knob, what in (("inversion", "inverted"), ("transposition", "moved within a chromosome"),
                       ("translocation", "moved to another chromosome"), ("loss", "deleted"),
                       ("duplication", "copied in tandem"), ("transfer", "copied to a recipient"),
                       ("origination", "laid down as new material")):
        g.add_argument(f"--{knob}-length", type=float, default=50.0, metavar="BP",
                       dest=f"{knob}_length",
                       help=f"mean bp {what} per event (geometric, default 50)")

    g = p.add_argument_group("outputs")
    g.add_argument("--write", nargs="+", choices=sorted({o for v in _OUTPUTS.values() for o in v}),
                   default=None, metavar="PART",
                   help="which outputs to write (default: each resolution's own, which is all of "
                        "them). unordered: events, profiles, genomes, initial_genome, gene_trees. "
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
    _add_quiet_arg(g)


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


def _read_tip_fates(path: str) -> dict:
    """Parse a ``--tip-fates`` file into ``{tip_name: fate}``: one
    ``tip_name<TAB>extant|extinct|unsampled`` row per tip (whitespace also accepted; blank lines and
    ``#`` comments skipped). This is the same shape ``species_fates.tsv`` is written in, so that output
    feeds straight back in — its ``lineage<TAB>fate`` header row is recognised and skipped. The values
    are checked against the tree by :func:`~zombi2.species.read_newick`."""
    fates = {}
    try:
        with open(path) as f:
            for lineno, raw in enumerate(f, 1):
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t") if "\t" in line else line.split()
                if parts == ["lineage", "fate"]:
                    continue  # the species_fates.tsv header, so that file is a valid --tip-fates input
                if len(parts) != 2:
                    raise ValueError(f"{path}:{lineno}: expected 'tip_name<TAB>extant|extinct|unsampled', "
                                     f"got {raw.rstrip()!r}")
                fates[parts[0]] = parts[1]
    except FileNotFoundError:
        raise FileNotFoundError(f"tip-fates file not found: {path}") from None
    return fates


def _stray(args, knobs) -> list[str]:
    """The flags in ``knobs`` the user actually set (their value differs from the default)."""
    return [f"--{attr.replace('_', '-')}" for attr, default in knobs
            if getattr(args, attr) != default]


def run(args, parser):
    # a flag a resolution does not have is an error, never silently ignored — otherwise
    # `--inversion` under unordered, or `--initial-families` under nucleotide, would quietly
    # produce a run that is not the one asked for
    if args.resolution == "unordered":
        if stray := _stray(args, _STRUCTURED_ONLY):
            parser.error(f"these options need --resolution ordered or nucleotide: "
                         f"{', '.join(stray)} (the unordered core has no chromosomes or positions)")
    if args.resolution != "nucleotide":
        if stray := _stray(args, _NUCLEOTIDE_ONLY):
            parser.error(f"these options need --resolution nucleotide: {', '.join(stray)} "
                         f"(the {args.resolution} resolution counts genes, not base pairs)")
    else:
        if stray := _stray(args, _NOT_IN_NUCLEOTIDE):
            parser.error(f"the nucleotide resolution has no {', '.join(stray)} (it is seeded from "
                         "a sequence — see --root-length / --genes / --gff — and its transfers are "
                         "additive)")
        if args.gff and args.genes:
            parser.error("pass either --gff or --genes, not both — a GFF already declares the genes")
        # the nucleotide engine holds each rate constant, so a modifier expression would be
        # accepted and then dropped; refuse it instead
        modulated = [f"--{n}" for n in ("duplication", "transfer", "loss", "origination", "inversion",
                                        "transposition", "translocation", "fission", "fusion")
                     if not isinstance(getattr(args, n), float)]
        if modulated:
            parser.error(f"--resolution nucleotide takes constant rates only, but "
                         f"{', '.join(modulated)} carries a modifier")

    vocab = _OUTPUTS[args.resolution]
    if args.write:
        bad = [o for o in args.write if o not in vocab]
        if bad:
            parser.error(f"--write {' '.join(bad)} not available for --resolution "
                         f"{args.resolution}; choose from: {', '.join(vocab)}")

    tip_fates = _read_tip_fates(args.tip_fates) if args.tip_fates else None
    tree_path = resolve_tree(args.source or args.run)
    try:
        with open(tree_path) as f:
            tree, names = read_newick(f.read(), tip_fates=tip_fates)
    except FileNotFoundError:
        raise FileNotFoundError(f"tree file not found: {tree_path}") from None

    common = dict(duplication=args.duplication, transfer=args.transfer, loss=args.loss,
                  origination=args.origination, transfer_to=args.transfer_to,
                  self_transfer=args.self_transfer, seed=args.seed)
    structured = dict(inversion=args.inversion, transposition=args.transposition,
                      translocation=args.translocation, chromosomes=args.chromosomes,
                      topology=args.topology, fission=args.fission, fusion=args.fusion,
                      chromosome_origination=args.chromosome_origination,
                      chromosome_loss=args.chromosome_loss,
                      inversion_probability=args.inversion_probability)

    t0 = time.perf_counter()
    if args.resolution == "ordered":
        result = simulate_genomes_ordered(
            tree, replacement=args.replacement, initial_families=args.initial_families,
            progress=not args.quiet, **structured, **common)
    elif args.resolution == "nucleotide":
        result = simulate_genomes_nucleotide(
            tree, root_length=args.root_length, genes=args.genes, gene_length=args.gene_length,
            gff=args.gff, trim_overlaps=args.trim_overlaps,
            inversion_length=args.inversion_length,
            transposition_length=args.transposition_length,
            translocation_length=args.translocation_length, loss_length=args.loss_length,
            duplication_length=args.duplication_length, transfer_length=args.transfer_length,
            origination_length=args.origination_length, progress=not args.quiet,
            **structured, **common)
    else:
        result = simulate_genomes_unordered(
            tree, replacement=args.replacement, initial_families=args.initial_families,
            progress=not args.quiet, **common)
    dt = time.perf_counter() - t0

    os.makedirs(args.run, exist_ok=True)
    out = level_dir(args.run, "genomes", args.flat)
    # gene trees are one file per family per view, so a hundred families is hundreds of files —
    # they get their own directory unless --flat says otherwise, whether asked for or defaulted
    wanted = tuple(args.write) if args.write else default_outputs(result)
    if rest := [o for o in wanted if o != "gene_trees"]:
        result.write(out, outputs=rest)
    if "gene_trees" in wanted:
        result.write(level_dir(out, "gene_trees", args.flat), outputs=("gene_trees",))
    # The events index against the tree canonicalised to n<id> labels, so the run needs that exact
    # tree to be replayable. A run grown here already has it — `zombi2 species` wrote the identical
    # file — so only a run reading its tree from elsewhere (--from) needs a copy, and it goes where
    # a species run keeps one rather than under a second name for the same thing.
    species_dir = level_dir(args.run, "species", args.flat)
    canonical = os.path.join(species_dir, "species_complete.nwk")
    if not os.path.exists(canonical):
        with open(canonical, "w") as f:
            f.write(result.complete_tree.to_newick() + "\n")
    if names:  # an external tree: map ZOMBI's n<id> back to the user's labels (join on profiles cols)
        rows = ["node\tname"] + [f"n{i}\t{lbl}" for i, lbl in sorted(names.items())]
        with open(os.path.join(out, "names.tsv"), "w") as f:
            f.write("\n".join(rows) + "\n")

    if args.resolution == "nucleotide":         # no phyletic profiles here: the unit is a base pair
        extant = [n.id for n in result.complete_tree.extant()]
        bp = sum(result.genomes[s].length for s in extant)
        summary = (f"{len(result.gene_spans)} genes and {bp} bp across {len(extant)} extant "
                   f"genomes (nucleotide)")
    else:
        n_families, n_species = result.profiles.shape
        summary = f"{n_families} gene families across {n_species} extant genomes ({args.resolution})"
    print(f"wrote {args.run}/ ({summary}) in {dt:.3g} s")
    _write_params_log(os.path.join(out, "genomes.log"),
                      args, summary)
    return 0
