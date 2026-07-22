"""``zombi2 genomes`` — evolve gene families along a species tree.

``--resolution`` picks the model: ``unordered`` (the D/T/L/O gene-family core,
:func:`zombi2.genomes.simulate_genomes_unordered`) or ``ordered`` (genes with a position and
orientation on chromosomes — segmental rearrangements and the chromosome tier,
:func:`~zombi2.genomes.simulate_genomes_ordered`). Long options are the API keyword names, and every
rate takes the written form (SPEC §5): a bare number on its natural scope, or the same ``scope(base)
× modifiers`` expression the Python API takes — ``--loss "0.25 * OnTime({0: 1.0, 3: 2.0})"``."""
from __future__ import annotations

import argparse
import os
import time

from zombi2.genomes import WIRED_MODIFIERS, simulate_genomes_ordered, simulate_genomes_unordered
from zombi2.species import read_newick
from zombi2.cli.framework import _add_params_arg, _rate, _rates_help, _write_params_log

#: the RATES block for ``zombi2 genomes -h``, built from the level's own declaration
RATES_HELP = _rates_help(
    WIRED_MODIFIERS, "--loss",
    note="Each rate keeps its natural scope here (D/T/L per copy, origination per lineage), so "
         "there is no scope wrapper to write. DrivenBy is wired for --loss, --duplication and "
         "--origination; --resolution ordered wires OnTime only.")

# the write vocabularies, mirroring each Result.write (there is no exported constant to import)
_UNORDERED_OUTPUTS = ("events", "profiles")
_ORDERED_OUTPUTS = ("events", "profiles", "gene_order", "rearrangements", "chromosome_events")

# ordered-only knobs — (attribute, default) pairs — rejected under --resolution unordered
_ORDERED_ONLY = (
    ("inversion", 0.0), ("transposition", 0.0), ("translocation", 0.0),
    ("chromosomes", 1), ("topology", "circular"),
    ("fission", 0.0), ("fusion", 0.0),
    ("chromosome_origination", 0.0), ("chromosome_loss", 0.0),
    ("inversion_probability", 0.0),
)


def _add_genomes_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("general")
    _add_params_arg(g)
    g.add_argument("-t", "--tree", required=True, metavar="FILE",
                   help="the species tree, Newick (a ZOMBI 'species_complete.nwk' or any external "
                        "tree; genomes evolve on the complete tree, extinct lineages included)")
    g.add_argument("-o", "--output", required=True, metavar="DIR", dest="output",
                   help="output directory (created if needed)")
    g.add_argument("--resolution", choices=("unordered", "ordered"), default="unordered",
                   metavar="RESOLUTION",
                   help="unordered (gene-family counts, default) or ordered (genes positioned on "
                        "chromosomes, with rearrangements)")
    g.add_argument("--seed", type=int, default=None, metavar="N",
                   help="RNG seed for reproducibility")
    g.add_argument("--tip-fates", metavar="FILE", dest="tip_fates",
                   help="[external non-ultrametric trees] a TSV 'tip_name<TAB>extant|extinct' "
                        "declaring each tip's fate; required when the input tree is not ultrametric "
                        "(ZOMBI won't guess extinct lineages from early-sampled tips)")

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
    g.add_argument("--transfer-to", choices=("uniform", "distance"), default="uniform",
                   metavar="RULE", dest="transfer_to",
                   help="recipient rule for a transfer: uniform (any contemporaneous lineage, "
                        "default) or distance (closer relatives likelier)")
    g.add_argument("--replacement", action="store_true",
                   help="a transfer overwrites a homologous copy in the recipient (replacing HGT)")
    g.add_argument("--self-transfer", action="store_true", dest="self_transfer",
                   help="allow a lineage to transfer to itself")
    g.add_argument("--initial-families", type=int, default=0, metavar="N", dest="initial_families",
                   help="number of gene families present at the crown (default 0)")

    g = p.add_argument_group("structured genome", "only with --resolution ordered")
    g.add_argument("--inversion", type=_rate, default=0.0, metavar="RATE",
                   help="segmental inversion rate (per chromosome)")
    g.add_argument("--transposition", type=_rate, default=0.0, metavar="RATE",
                   help="segmental transposition rate — move a run within a chromosome "
                        "(per chromosome)")
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

    g = p.add_argument_group("outputs")
    g.add_argument("--write", nargs="+", choices=_ORDERED_OUTPUTS, default=None, metavar="PART",
                   help="which outputs to write (default: events, profiles [+ gene_order when "
                        "ordered]). unordered: events, profiles. ordered adds: gene_order, "
                        "rearrangements, chromosome_events.")


def _read_tip_fates(path: str) -> dict:
    """Parse a ``--tip-fates`` file into ``{tip_name: fate}``: one ``tip_name<TAB>extant|extinct``
    row per tip (whitespace also accepted; blank lines and ``#`` comments skipped). The values are
    checked against the tree by :func:`~zombi2.species.read_newick`."""
    fates = {}
    try:
        with open(path) as f:
            for lineno, raw in enumerate(f, 1):
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t") if "\t" in line else line.split()
                if len(parts) != 2:
                    raise ValueError(f"{path}:{lineno}: expected 'tip_name<TAB>extant|extinct', "
                                     f"got {raw.rstrip()!r}")
                fates[parts[0]] = parts[1]
    except FileNotFoundError:
        raise FileNotFoundError(f"tip-fates file not found: {path}") from None
    return fates


def run(args, parser):
    # reject ordered-only knobs under the unordered resolution, so a silently-ignored flag can't
    # give a misleading run (e.g. --inversion with --resolution unordered)
    if args.resolution == "unordered":
        stray = [f"--{attr.replace('_', '-')}" for attr, default in _ORDERED_ONLY
                 if getattr(args, attr) != default]
        if stray:
            parser.error(f"these options need --resolution ordered: {', '.join(stray)} "
                         "(the unordered core has no chromosomes or positions)")

    vocab = _ORDERED_OUTPUTS if args.resolution == "ordered" else _UNORDERED_OUTPUTS
    if args.write:
        bad = [o for o in args.write if o not in vocab]
        if bad:
            parser.error(f"--write {' '.join(bad)} not available for --resolution "
                         f"{args.resolution}; choose from: {', '.join(vocab)}")

    tip_fates = _read_tip_fates(args.tip_fates) if args.tip_fates else None
    try:
        with open(args.tree) as f:
            tree, names = read_newick(f.read(), tip_fates=tip_fates)
    except FileNotFoundError:
        raise FileNotFoundError(f"tree file not found: {args.tree}") from None

    common = dict(duplication=args.duplication, transfer=args.transfer, loss=args.loss,
                  origination=args.origination, transfer_to=args.transfer_to,
                  replacement=args.replacement, self_transfer=args.self_transfer,
                  initial_families=args.initial_families, seed=args.seed)

    t0 = time.perf_counter()
    if args.resolution == "ordered":
        result = simulate_genomes_ordered(
            tree, inversion=args.inversion, transposition=args.transposition,
            translocation=args.translocation, chromosomes=args.chromosomes, topology=args.topology,
            fission=args.fission, fusion=args.fusion,
            chromosome_origination=args.chromosome_origination,
            chromosome_loss=args.chromosome_loss,
            inversion_probability=args.inversion_probability, **common)
    else:
        result = simulate_genomes_unordered(tree, **common)
    dt = time.perf_counter() - t0

    os.makedirs(args.output, exist_ok=True)
    if args.write:
        result.write(args.output, outputs=args.write)
    else:
        result.write(args.output)               # each Result.write's own default
    # the tree the events are indexed against, canonicalised to n<id> labels so its ids match the
    # event log's `lineage` column — this makes the run self-describing and lets `zombi2 sequences
    # --genomes DIR` rebuild the gene trees (from genome_events.tsv + this tree) with no other input.
    with open(os.path.join(args.output, "genome_species_tree.nwk"), "w") as f:
        f.write(result.complete_tree.to_newick() + "\n")
    if names:  # an external tree: map ZOMBI's n<id> back to the user's labels (join on profiles cols)
        rows = ["node\tname"] + [f"n{i}\t{lbl}" for i, lbl in sorted(names.items())]
        with open(os.path.join(args.output, "names.tsv"), "w") as f:
            f.write("\n".join(rows) + "\n")

    n_families, n_species = result.profiles.shape
    summary = f"{n_families} gene families across {n_species} extant genomes ({args.resolution})"
    print(f"wrote {args.output}/ ({summary}) in {dt:.3g} s")
    _write_params_log(os.path.join(args.output, "genomes.log"), args, summary)
    return 0
