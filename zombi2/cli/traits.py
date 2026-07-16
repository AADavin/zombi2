"""zombi2 traits command."""
from __future__ import annotations

import argparse
import os

import numpy as np

from zombi2.traits.biogeography import DEC, simulate_biogeography
from zombi2.traits.models import (
    BrownianMotion, OrnsteinUhlenbeck, EarlyBurst, Mk, ThresholdModel, simulate_traits,
)
from zombi2.tree import read_newick
from zombi2.cli.framework import _write_params_log


def _add_trait_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("general")
    g.add_argument("-t", "--tree", required=True, metavar="FILE",
                   help="input species tree in Newick format (e.g. species_tree.nwk)")
    g.add_argument("--model", choices=("bm", "ou", "eb", "mk", "threshold", "dec"), default="bm",
                   metavar="MODEL",
                   help="trait model: bm=Brownian motion, ou=Ornstein-Uhlenbeck, "
                        "eb=early burst/ACDC, mk=discrete k-state, threshold, "
                        "dec=geographic-range evolution (default: bm)")
    g.add_argument("--replicates", type=int, default=1, metavar="N",
                   help="simulate the trait this many times with the same parameters; writes "
                        "traits.tsv with one column per replicate (default: 1)")
    g.add_argument("--seed", type=int, default=None, metavar="N",
                   help="RNG seed for reproducibility")
    g.add_argument("-o", "--out", required=True, metavar="DIR", help="output directory")

    g = p.add_argument_group("continuous traits", "bm / ou / eb / threshold")
    g.add_argument("--sigma2", type=float, default=1.0, metavar="S2",
                   help="diffusion rate (default: 1.0)")
    g.add_argument("--x0", type=float, default=None, metavar="X0",
                   help="root value [bm/eb/threshold]; OU defaults it to --theta")
    g.add_argument("--trend", type=float, default=0.0, metavar="MU", help="directional drift [bm/eb]")

    g = p.add_argument_group("ornstein-uhlenbeck", "--model ou")
    g.add_argument("--alpha", type=float, default=1.0, metavar="A",
                   help="mean-reversion strength (default: 1.0)")
    g.add_argument("--theta", type=float, default=0.0, metavar="T", help="optimum (default: 0.0)")

    g = p.add_argument_group("early burst & Mk rate", "--model eb / --model mk")
    g.add_argument("--rate", type=float, default=1.0, metavar="R",
                   help="EB rate-of-change (negative = early burst) [eb], or the per-transition "
                        "rate [mk] (default: 1.0)")

    g = p.add_argument_group("discrete Mk", "--model mk")
    g.add_argument("--states", type=int, default=2, metavar="K",
                   help="number of states for the mk model (default: 2)")
    g.add_argument("--ordered", action="store_true",
                   help="only allow transitions between adjacent states (i <-> i±1)")
    g.add_argument("--q-matrix", default=None, metavar="FILE",
                   help="path to a whitespace/comma-separated k x k rate matrix (an arbitrary "
                        "Markov chain); overrides --states/--rate/--ordered")

    g = p.add_argument_group("threshold", "--model threshold")
    g.add_argument("--thresholds", default="0.0", metavar="CUTS",
                   help="comma-separated liability cut points (default: 0.0)")

    g = p.add_argument_group("DEC biogeography", "--model dec")
    g.add_argument("--areas", default="3", metavar="SPEC",
                   help="number of areas (e.g. 3) or comma-separated area labels (e.g. A,B,C) "
                        "(default: 3)")
    g.add_argument("--dispersal", type=float, default=0.1, metavar="RATE",
                   help="rate of gaining an area (dispersal) (default: 0.1)")
    g.add_argument("--extinction", type=float, default=0.1, metavar="RATE",
                   help="rate of losing an area (local extinction) (default: 0.1)")
    g.add_argument("--max-range-size", type=int, default=None, metavar="N",
                   help="maximum number of areas a range may span (default: all)")
    g.add_argument("--root-range", default=None, metavar="AREAS",
                   help="comma-separated area labels for the root range (e.g. A); default: random")


def _read_q_matrix(path: str):
    """Read a ``k x k`` rate matrix from a whitespace/comma-separated file.

    Blank lines and ``#`` comments are skipped; the diagonal is ignored (recomputed by
    :class:`~zombi2.Mk`). Each row is the *from*-state, each column the *to*-state.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            rows.append([float(x) for x in line.replace(",", " ").split()])
    if rows and any(len(r) != len(rows[0]) for r in rows):
        raise ValueError("q-matrix rows must all have the same length (a square k x k matrix)")
    return rows


def _build_trait_model(args):
    x0 = args.x0
    if args.model == "bm":
        return BrownianMotion(sigma2=args.sigma2, x0=(0.0 if x0 is None else x0), trend=args.trend)
    if args.model == "ou":
        return OrnsteinUhlenbeck(sigma2=args.sigma2, alpha=args.alpha, theta=args.theta, x0=x0)
    if args.model == "eb":
        return EarlyBurst(sigma2=args.sigma2, rate=args.rate,
                          x0=(0.0 if x0 is None else x0), trend=args.trend)
    if args.model == "mk":
        if args.q_matrix:                                    # arbitrary user-supplied Markov chain
            return Mk(_read_q_matrix(args.q_matrix))
        if args.ordered:                                     # adjacent-only (meristic) character
            return Mk.ordered(args.states, args.rate)
        return Mk.equal_rates(args.states, args.rate)        # equal rates (all-to-all)
    # threshold
    thresholds = [float(t) for t in str(args.thresholds).split(",")]
    return ThresholdModel(thresholds=thresholds, sigma2=args.sigma2,
                          x0=(0.0 if x0 is None else x0))


def _parse_areas(text: str):
    """``--areas``: an integer count, or a comma-separated list of area labels."""
    text = str(text)
    if "," in text:
        return [a.strip() for a in text.split(",")]
    try:
        return int(text)
    except ValueError:
        return [text]


def _build_dec_model(args) -> DEC:
    return DEC(areas=_parse_areas(args.areas), dispersal=args.dispersal,
               extinction=args.extinction, max_range_size=args.max_range_size)


def _dec_root(args):
    """The root range from ``--root-range`` (a set of area labels), or ``None``."""
    if args.root_range is None:
        return None
    return {a.strip() for a in str(args.root_range).split(",")}


def _fmt_cell(value) -> str:
    """Format one trait value for a TSV cell (range tuple -> {A,B}, float -> 6 sig figs)."""
    if isinstance(value, tuple):
        return "{" + ",".join(str(v) for v in value) + "}"
    return f"{value:.6g}" if isinstance(value, float) else str(value)


def _replicate_table(results) -> str:
    """A wide ``node``-by-replicate table: one column (rep_1..rep_N) per simulation, one row per
    node (tips and ancestral)."""
    tree = results[0].tree
    header = "node\t" + "\t".join(f"rep_{i + 1}" for i in range(len(results)))
    lines = [header]
    for node in tree.nodes():
        cells = [_fmt_cell(res.label(res.node_values[node])) for res in results]
        lines.append(node.name + "\t" + "\t".join(cells))
    return "\n".join(lines) + "\n"


def _run_trait(args) -> str:
    """Evolve a trait along the supplied species tree and write the output folder.

    DEC (geographic ranges) runs the biogeography driver (it splits ranges at speciations);
    every other model overlays a trait with the standard driver. Both return a ``TraitResult``,
    so the output writing is shared.
    """
    if args.replicates < 1:
        raise ValueError("--replicates must be >= 1")
    with open(args.tree) as f:
        tree = read_newick(f.read())
    os.makedirs(args.out, exist_ok=True)

    if args.model == "dec":
        dec = _build_dec_model(args)
        root = _dec_root(args)
        simulate = lambda rng: simulate_biogeography(tree, dec, root_state=root, rng=rng)  # noqa: E731
    else:
        model = _build_trait_model(args)
        simulate = lambda rng: simulate_traits(tree, model, rng=rng)  # noqa: E731

    if args.replicates > 1:
        rng = np.random.default_rng(args.seed)
        results = [simulate(rng) for _ in range(args.replicates)]
        with open(os.path.join(args.out, "traits.tsv"), "w") as f:
            f.write(_replicate_table(results))       # one column per replicate, all nodes
        return (f"wrote {args.replicates} trait replicates to {args.out}/traits.tsv "
                f"(model={args.model}; {len(tree.nodes())} nodes x {args.replicates} columns)")

    res = simulate(np.random.default_rng(args.seed))
    with open(os.path.join(args.out, "traits.tsv"), "w") as f:
        f.write(res.to_tsv(nodes="all"))             # every node: tips AND ancestral states
    with open(os.path.join(args.out, "trait_tree.nwk"), "w") as f:
        f.write(res.to_newick() + "\n")              # values annotated on every node too
    return (f"wrote traits to {args.out}/ (model={args.model}; "
            f"{len(tree.extant_leaves())} tip + {len(tree.internal_nodes())} ancestral values)")


# ═══════════════════════════════════════════════════════════════════════════════
# coevolve --couple traits:genes: trait-conditioned gene-family dynamics
# ═══════════════════════════════════════════════════════════════════════════════
def _add_trait_model_args(g) -> None:
    """Scalar trait-model flags (DEC ranges do not apply), added to the argument group ``g``.
    Used by the ``coevolve --couple traits:genes`` edge for the trait it simulates.

    ``--trait-model`` stores into ``args.model`` so :func:`_build_trait_model` is reused as-is.
    The Mk rate matrix reuses the command's shared ``--q-matrix`` (not re-added here).
    """
    g.add_argument("--trait-model", dest="model", default="bm", metavar="MODEL",
                   choices=("bm", "ou", "eb", "mk", "threshold"),
                   help="trait to evolve then couple to gene families: bm=Brownian motion, "
                        "ou=Ornstein-Uhlenbeck, eb=early burst, mk=discrete k-state, threshold "
                        "(default: bm). Use --trait-file to supply a precomputed trait instead")
    g.add_argument("--sigma2", type=float, default=1.0, metavar="S2",
                   help="diffusion rate [bm/ou/eb/threshold] (default: 1.0)")
    g.add_argument("--x0", type=float, default=None, metavar="X0",
                   help="root value [bm/eb/threshold]; OU defaults it to --theta")
    g.add_argument("--trend", type=float, default=0.0, metavar="MU", help="directional drift [bm/eb]")
    g.add_argument("--alpha", type=float, default=1.0, metavar="A",
                   help="OU mean-reversion strength [ou]")
    g.add_argument("--theta", type=float, default=0.0, metavar="T", help="OU optimum [ou]")
    g.add_argument("--rate", type=float, default=1.0, metavar="R",
                   help="EB rate-of-change (negative = early burst) [eb], or per-transition rate [mk]")
    g.add_argument("--states", type=int, default=2, metavar="K",
                   help="number of states for the mk model (default: 2)")
    g.add_argument("--ordered", action="store_true",
                   help="[mk] only allow transitions between adjacent states (i <-> i±1)")
    g.add_argument("--thresholds", default="0.0", metavar="CUTS",
                   help="comma-separated liability cut points [threshold] (default: 0.0)")


def run(args, parser):
    summary = _run_trait(args)
    print(summary)
    _write_params_log(os.path.join(args.out, "traits.log"), args, summary)
    return 0
