"""zombi2 coevolve command."""
from __future__ import annotations

import argparse
import os
import time

import numpy as np


from zombi2.genomes.simulation import Genomes, simulate_genomes
from zombi2.coevolve.sse import BiSSE, MuSSE, QuaSSE, HiSSE, simulate_sse
from zombi2.coevolve.gene_diversification import (
    GeneDiversification, simulate_gene_diversification, simulate_co_diversification,
)
from zombi2.coevolve.cladogenetic_genome import (
    CladogeneticGenome, simulate_cladogenetic_genome, _branch_count_and_length,
)
from zombi2.coevolve.gene_conditioned_trait import GeneConditionedTrait, simulate_gene_conditioned_trait
from zombi2.coevolve.trait_coupling import TraitGeneCoupling, simulate_trait_linked_genomes
from zombi2.coevolve.trait_gene_feedback import TraitGeneFeedback, simulate_trait_gene_feedback
from zombi2.traits.models import (
    BrownianMotion, Mk, TraitResult,
    Cladogenesis, simulate_traits,
)
from zombi2.tree import Tree, prune, read_newick

from zombi2.cli.framework import _write_params_log
from zombi2.cli.genomes import _write_profiles_only
from zombi2.cli.trait import _add_trait_model_args, _build_trait_model, _read_q_matrix

def _add_traits_genes_args(p: argparse.ArgumentParser) -> None:
    """The ``coevolve --couple traits:genes`` flags — a trait conditions a gene-family panel's
    loss/gain (formerly the standalone ``coevolve-genetrait`` command). Runs on a GIVEN -t tree.
    The trait's Mk matrix reuses the shared ``--q-matrix``; the panel writes with ``--write``."""
    g = p.add_argument_group("traits:genes trait model",
                             "--couple traits:genes: the trait to simulate (on a GIVEN tree; needs "
                             "-t). mk uses the shared --q-matrix")
    _add_trait_model_args(g)
    g.add_argument("--trait-file", default=None, metavar="TSV",
                   help="use a precomputed trait instead of simulating one: a node<TAB>value table "
                        "over ALL nodes (tips and ancestors), as 'zombi2 trait' writes with "
                        "nodes=all; values must be numeric (encode discrete states as numbers). "
                        "Overrides --trait-model")
    g.add_argument("--trait-center", action="store_true",
                   help="[discrete trait] center the state values around their mean so the trait "
                        "pushes retention both up and down — recommended for a binary "
                        "aerobic/anaerobic trait; by default states are 0,1,..,k-1")
    g.add_argument("--trait-steps", type=int, default=16, metavar="K",
                   help="[continuous trait] within-branch resolution: sub-segment each branch into "
                        "K pieces (default 16; ignored for discrete traits, which use their exact "
                        "stochastic map)")

    g = p.add_argument_group("traits:genes gene panel & coupling",
                             "the panel, its trait-neutral base rates, and which families respond")
    g.add_argument("--panel", type=int, default=50, metavar="N",
                   help="number of gene families in the panel (default 50)")
    g.add_argument("--loss", type=float, default=0.5, metavar="RATE",
                   help="baseline per-copy loss rate — the loss where the trait is neutral (default 0.5)")
    g.add_argument("--trans", type=float, default=1.0, metavar="RATE",
                   help="per-copy transfer (HGT) rate — the field-blind gain channel (default 1.0)")
    g.add_argument("--dup", type=float, default=0.0, metavar="RATE",
                   help="per-copy duplication rate, trait-independent (default 0)")
    g.add_argument("--orig", type=float, default=0.0, metavar="RATE",
                   help="background origination rate of brand-new, uncoupled families (default 0)")
    g.add_argument("--responsive", default="0.3", metavar="SPEC",
                   help="which families respond to the trait: an integer count, a fraction "
                        "(e.g. 0.3), a comma-separated id/index list (e.g. F3,F7,12), or @FILE of "
                        "ids/indices (default: 0.3 = 30%% of the panel, chosen at random)")
    g.add_argument("--weight", type=float, default=1.0, metavar="W",
                   help="coupling weight of each responsive family (default 1.0)")
    g.add_argument("--signed", action="store_true",
                   help="randomise the sign of each responsive weight (some families favoured by a "
                        "high trait value, some by a low one); by default all favour a high value")
    g.add_argument("--effect-loss", type=float, default=2.0, metavar="B",
                   help="retention coupling strength: a responsive family's loss scales by "
                        "exp(-effect_loss * weight * trait) (default 2.0; 0 = no coupling)")
    g.add_argument("--effect-gain", type=float, default=0.0, metavar="B",
                   help="optional HGT-activity coupling: a lineage's transfer rate scales by "
                        "exp(effect_gain * trait) (default 0 = trait-blind gain)")
    g.add_argument("--panel-root-fraction", dest="panel_root_fraction", type=float, default=0.5,
                   metavar="F",
                   help="[traits:genes + genes:traits JOINT model only] fraction of the panel present "
                        "at the root (default 0.5). In the joint model --theta-absent/--theta-present "
                        "are the trait's optima at an empty/full panel, and --trait-alpha/--trait-sigma2 "
                        "its OU dynamics")

    g = p.add_argument_group("traits:genes output")
    g.add_argument("--write", dest="write", nargs="+", metavar="PART",
                   choices=(*Genomes.WRITE_PARTS, "all"), default=["profiles", "trees"],
                   help="which gene-family outputs to write — any of {profiles, trace, trees, "
                        "events, transfers, summary} or 'all' (default: profiles trees). "
                        "traits.tsv, trait_tree.nwk and coupling.tsv (the responsive-family "
                        "manifest) are always written alongside")
    g.add_argument("--sparse", action="store_true",
                   help="write the profile as a sparse long table (needs 'profiles' in --write)")
    g.add_argument("--annotate-species", action="store_true",
                   help="label internal gene-tree nodes <gid>|<species-branch> (e.g. g570|i5)")

def _parse_responsive(text: str):
    """``--responsive``: a count (int), a fraction (float), or an id/index list (``@FILE`` or CSV)."""
    text = str(text).strip()
    if text.startswith("@"):
        with open(text[1:]) as f:
            content = f.read()
        return [tok for tok in content.replace(",", " ").split() if tok]
    if "," in text:
        return [tok.strip() for tok in text.split(",") if tok.strip()]
    if "." in text:
        return float(text)
    return int(text)

def _load_trait_result(tree: Tree, path: str) -> TraitResult:
    """Load a precomputed trait: a ``node<TAB>value`` table over every node (numeric values).

    Returns a continuous-kind :class:`~zombi2.traits.TraitResult` (values used at node
    resolution; a supplied trait carries no within-branch stochastic map). Every node — tips and
    ancestors — must be present, since the gene simulation reads the trait on every branch.
    """
    name2node = {n.name: n for n in tree.nodes()}
    values: dict = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2 or parts[0] == "node":  # skip a header row
                continue
            node = name2node.get(parts[0])
            if node is None:
                continue
            try:
                values[node] = float(parts[1])
            except ValueError:
                raise ValueError(
                    f"--trait-file needs numeric trait values; got {parts[1]!r} for node "
                    f"{parts[0]!r}. Encode discrete states as numbers (e.g. 0 / 1).") from None
    missing = [n.name for n in tree.nodes() if n not in values]
    if missing:
        raise ValueError(
            f"--trait-file is missing values for {len(missing)} node(s) (e.g. {missing[:3]}); it "
            "must cover every node — tips AND ancestors — like the traits.tsv that 'zombi2 trait' "
            "writes (its nodes=all output)")
    return TraitResult(tree=tree, model=None, node_values=values, history=None, kind="continuous")

def _write_coupling_manifest(out: str, coupling: TraitGeneCoupling) -> None:
    """Write ``coupling.tsv`` — the per-family coupling weights plus the effect sizes, so the
    trait↔gene linkage that generated the profiles is recorded for downstream inference."""
    lines = [f"# effect_loss\t{coupling.effect_loss}",
             f"# effect_gain\t{coupling.effect_gain}",
             f"# base_loss\t{coupling.base_loss}",
             f"# transfer\t{coupling.transfer}",
             f"# duplication\t{coupling.duplication}",
             f"# origination\t{coupling.origination}",
             "family\tweight"]
    for i, fam in enumerate(coupling.panel_ids):
        lines.append(f"{fam}\t{coupling.weights[i]:.6g}")
    with open(os.path.join(out, "coupling.tsv"), "w") as f:
        f.write("\n".join(lines) + "\n")

def _run_traits_genes(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """``coevolve --couple traits:genes``: simulate a trait, then evolve a gene-family panel whose
    loss/gain is conditioned on it. An overlay edge — runs along a GIVEN -t tree."""
    if not args.tree:
        parser.error("traits:genes runs on a GIVEN tree — pass -t/--tree (a trait conditions gene "
                     "content along it; there is nothing to grow)")
    if args.age is not None or args.tips is not None:
        parser.error("traits:genes uses the given -t tree; --age/--tips only apply to the "
                     "into-species edges that grow a tree")
    with open(args.tree) as f:
        tree = read_newick(f.read())
    parts = set(Genomes.WRITE_PARTS) if "all" in args.write else set(args.write)
    if args.sparse and "profiles" not in parts:
        raise ValueError("--sparse affects the profile output; add 'profiles' to --write")
    if args.panel < 1:
        raise ValueError("--panel must be >= 1")
    rng = np.random.default_rng(args.seed)

    # 1) the trait: simulate one (--trait-model) or load a precomputed one (--trait-file)
    if args.trait_file:
        result = _load_trait_result(tree, args.trait_file)
        trait_desc = f"file:{os.path.basename(args.trait_file)}"
    else:
        result = simulate_traits(tree, _build_trait_model(args), rng=rng)
        trait_desc = args.model

    # optional: center discrete states so the coupling is two-sided (recommended for binary)
    state_values = None
    if args.trait_center and result.kind == "discrete":
        k = len(result.model.states)
        state_values = [i - (k - 1) / 2.0 for i in range(k)]

    # 2) the coupling: choose the responsive families and the effect sizes
    coupling = TraitGeneCoupling.build(
        args.panel, _parse_responsive(args.responsive), weight=args.weight, signed=args.signed,
        effect_loss=args.effect_loss, effect_gain=args.effect_gain, base_loss=args.loss,
        transfer=args.trans, duplication=args.dup, origination=args.orig,
        state_values=state_values, rng=rng)

    if args.null == "cid":
        return _run_traits_genes_cid_null(args, parser)
    if args.null == "timing":
        parser.error("traits:genes has no 'timing' null; use --null neutral")
    if args.null == "neutral":
        coupling = coupling.null("neutral")

    # 3) run
    t0 = time.perf_counter()
    res = simulate_trait_linked_genomes(tree, result, coupling, trait_steps=args.trait_steps, rng=rng)
    dt = time.perf_counter() - t0

    os.makedirs(args.out, exist_ok=True)
    res.genomes().write(args.out, include=parts, sparse=args.sparse,
                        annotate_species=args.annotate_species)
    with open(os.path.join(args.out, "traits.tsv"), "w") as f:
        f.write(res.trait.to_tsv(nodes="all"))
    with open(os.path.join(args.out, "trait_tree.nwk"), "w") as f:
        f.write(res.trait.to_newick() + "\n")
    _write_coupling_manifest(args.out, coupling)
    if args.null == "neutral":
        _write_null_manifest(args.out, "traits:genes", "neutral", cut="trait -> panel loss/gain",
                             preserved="removed (effect_loss = effect_gain = 0; uncoupled panel)")
    tag = " [neutral null]" if args.null == "neutral" else ""
    return (f"wrote{tag} [{' '.join(sorted(parts))}] + trait to {args.out}/ (trait={trait_desc}, "
            f"panel {coupling.n_families} families, {coupling.n_responsive} responsive, "
            f"{len(tree.extant_leaves())} tips) in {dt:.3g} s")

_COEVOLVE_NODES = ("species", "traits", "genes")

_COEVOLVE_EDGES = {
    "traits:species", "genes:species", "species:traits",
    "species:genes", "traits:genes", "genes:traits",
}

def _add_coevolve_mode_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("general")
    g.add_argument("--couple", action="append", nargs="+", metavar="DRIVER:TARGET", default=None,
                   help="a directed coupling edge 'driver:target' over {species, traits, genes} — "
                        "the driver's state modulates the target's rates. Implemented: "
                        "'traits:species' (SSE), 'species:traits' (cladogenetic), their "
                        "combination = ClaSSE, 'genes:species' (key innovations), 'species:genes' "
                        "(cladogenetic genome), 'genes:traits' (a modifier gene switches a trait "
                        "optimum) and 'traits:genes' (a trait conditions a gene-family panel). "
                        "Repeatable; default traits:species. See docs/guide/coevolution.md for "
                        "the full edge set")
    g.add_argument("-t", "--tree", default=None, metavar="FILE",
                   help="input species tree (Newick) — required for the on-a-given-tree edges "
                        "(species:traits, species:genes, genes:traits, traits:genes). Omit for the "
                        "into-species edges (traits:species / ClaSSE / genes:species), which GROW "
                        "the tree via --age/--tips")
    g.add_argument("--age", type=float, default=None, metavar="T",
                   help="[into-species] crown age to grow for (the extant tip count is random)")
    g.add_argument("--tips", type=int, default=None, metavar="N",
                   help="[into-species] stop when this many extant tips first coexist (age random)")
    g.add_argument("--seed", type=int, default=None, metavar="N",
                   help="RNG seed for reproducibility")
    g.add_argument("--null", choices=("none", "neutral", "cid", "timing"), default="none",
                   help="generate the matched DECOUPLED null instead of the coupled model — cut the "
                        "driver→target arrow while keeping the target's variance (for calibrating a "
                        "detector's false-positive rate). 'neutral' (all edges): the driver stops "
                        "setting the rates. 'cid': the variance comes from a HIDDEN, uncorrelated "
                        "driver (traits:species natively; the gene/trait edges via a neutral "
                        "observed channel). 'timing': an at-speciation burst is spread along "
                        "branches (species:traits, species:genes). See "
                        "docs/guide/coevolution_nulls.md")
    g.add_argument("--hidden", type=int, default=2, metavar="H",
                   help="[--null cid, traits:species] number of hidden rate classes (2 = CID-2, "
                        "4 = CID-4)")
    g.add_argument("-o", "--out", required=True, metavar="DIR", help="output directory")

    g = p.add_argument_group("SSE model", "--couple traits:species (trait drives diversification)")
    g.add_argument("--sse-model", dest="sse_model", choices=("bisse", "musse", "quasse", "hisse"),
                   default="bisse", metavar="MODEL",
                   help="which state-dependent model drives diversification: bisse (binary trait), "
                        "musse (k-state), quasse (continuous trait), hisse (binary trait + hidden "
                        "diversification classes) (default: bisse)")
    g.add_argument("--root-state", type=int, default=None, metavar="I",
                   help="[bisse/musse] root state index (default: the character's stationary "
                        "distribution)")

    g = p.add_argument_group("BiSSE", "--sse-model bisse (binary trait)")
    g.add_argument("--lambda0", type=float, default=1.0, metavar="RATE", help="speciation in state 0")
    g.add_argument("--lambda1", type=float, default=2.0, metavar="RATE", help="speciation in state 1")
    g.add_argument("--mu0", type=float, default=0.3, metavar="RATE", help="extinction in state 0")
    g.add_argument("--mu1", type=float, default=0.3, metavar="RATE", help="extinction in state 1")
    g.add_argument("--q01", type=float, default=0.1, metavar="RATE", help="transition rate 0 -> 1")
    g.add_argument("--q10", type=float, default=0.1, metavar="RATE", help="transition rate 1 -> 0")

    g = p.add_argument_group("MuSSE", "--sse-model musse (k-state trait)")
    g.add_argument("--birth", type=float, nargs="+", default=None, metavar="RATE",
                   help="per-state speciation rates (k values)")
    g.add_argument("--death", type=float, nargs="+", default=None, metavar="RATE",
                   help="per-state extinction rates (k values)")
    g.add_argument("--q-matrix", default=None, metavar="FILE",
                   help="path to a k x k anagenetic transition-rate matrix (same format as "
                        "'zombi2 trait --q-matrix')")

    g = p.add_argument_group("QuaSSE", "--sse-model quasse (continuous trait)")
    g.add_argument("--spec-low", type=float, default=0.5, metavar="RATE",
                   help="speciation rate at low trait values")
    g.add_argument("--spec-high", type=float, default=2.0, metavar="RATE",
                   help="speciation rate at high trait values")
    g.add_argument("--spec-center", type=float, default=0.0, metavar="X",
                   help="trait value at the middle of the speciation sigmoid")
    g.add_argument("--spec-slope", type=float, default=1.0, metavar="S",
                   help="steepness of the speciation sigmoid")
    g.add_argument("--qmu", type=float, default=0.1, metavar="RATE", help="constant extinction rate")
    g.add_argument("--diffusion", type=float, default=1.0, metavar="S2",
                   help="trait diffusion rate sigma^2 (Brownian motion)")
    g.add_argument("--root-value", type=float, default=0.0, metavar="X0", help="root trait value x0")

    g = p.add_argument_group("HiSSE", "--sse-model hisse (binary trait + hidden diversification "
                                      "classes)")
    g.add_argument("--hidden-classes", dest="hidden_classes", type=int, default=2, metavar="H",
                   help="number of hidden diversification classes (>= 2; default 2)")
    g.add_argument("--hidden-scale", dest="hidden_scale", type=float, default=3.0, metavar="S",
                   help="speciation spread across the hidden classes: the classes span the base "
                        "--lambda0/--lambda1 up to S times faster (geometric; default 3.0). The "
                        "observed transitions --q01/--q10 and extinction --mu0/--mu1 are shared")
    g.add_argument("--hidden-switch", dest="hidden_switch", type=float, default=0.1, metavar="RATE",
                   help="rate of switching between hidden classes (symmetric; default 0.1)")

    g = p.add_argument_group("cladogenetic kernel",
                             "--couple species:traits (speciation jumps the trait)")
    g.add_argument("--clado-shift", dest="clado_shift", type=float, default=0.3, metavar="P",
                   help="[discrete trait] probability a daughter hops to another state AT each "
                        "speciation (cladogenetic change; default 0.3)")
    g.add_argument("--clado-jump", dest="clado_jump", type=float, default=1.0, metavar="S2",
                   help="[continuous trait] variance of the Gaussian jump added to each daughter's "
                        "value AT each speciation (default 1.0)")

    g = p.add_argument_group(
        "gene-driven diversification",
        "--couple genes:species (key-innovation gene families; base rates reuse --lambda0/--mu0)")
    g.add_argument("--drivers", type=int, default=2, metavar="N",
                   help="number of binary 'driver' (key-innovation) gene families")
    g.add_argument("--driver-speciation", dest="driver_speciation", type=float, default=1.0,
                   metavar="B",
                   help="per-driver effect on log speciation: a present driver scales lambda by "
                        "exp(this) (>0 = a key innovation; default 1.0)")
    g.add_argument("--driver-extinction", dest="driver_extinction", type=float, default=0.0,
                   metavar="B",
                   help="per-driver effect on log extinction: a present driver scales mu by "
                        "exp(this) (default 0)")
    g.add_argument("--driver-loss", dest="driver_loss", type=float, default=0.1, metavar="RATE",
                   help="rate a present driver is lost/deleted (default 0.1)")
    g.add_argument("--driver-origination", dest="driver_origination", type=float, default=0.05,
                   metavar="RATE", help="rate an absent driver appears de novo (default 0.05)")
    g.add_argument("--driver-transfer", dest="driver_transfer", type=float, default=0.5,
                   metavar="RATE",
                   help="per-donor HGT rate of a driver — frequency-dependent gain: a driver in "
                        "more live genomes spreads faster (default 0.5)")
    g.add_argument("--root-drivers", dest="root_drivers", type=int, default=0, metavar="M",
                   help="number of drivers present at the root (the first m; default 0 = drivers "
                        "enter by origination)")
    g.add_argument("--driver-clado-loss", dest="driver_clado_loss", type=float, default=0.0,
                   metavar="P",
                   help="probability a daughter drops each driver it carries AT each speciation. "
                        ">0 adds the species:genes burst, making this the species<->genes JOINT "
                        "model (--couple genes:species --couple species:genes; default 0)")
    g.add_argument("--driver-clado-gain", dest="driver_clado_gain", type=float, default=0.0,
                   metavar="P",
                   help="probability a daughter gains each absent driver AT each speciation "
                        "(the species:genes burst; default 0)")

    g = p.add_argument_group("cladogenetic genome",
                             "--couple species:genes (on a GIVEN tree; needs -t)")
    g.add_argument("--genome-size", dest="genome_size", type=int, default=30, metavar="N",
                   help="number of families in the root genome (default 30)")
    g.add_argument("--gene-loss", dest="gene_loss", type=float, default=0.0, metavar="RATE",
                   help="anagenetic per-family loss rate along a branch (default 0)")
    g.add_argument("--gene-origination", dest="gene_origination", type=float, default=0.0,
                   metavar="RATE",
                   help="anagenetic origination rate of new families, per lineage (default 0). "
                        "With both anagenetic rates 0 the change is purely cladogenetic")
    g.add_argument("--clado-gene-loss", dest="clado_gene_loss", type=float, default=0.1, metavar="P",
                   help="probability a daughter drops each family AT each speciation (the "
                        "founder-effect burst; default 0.1)")
    g.add_argument("--clado-gene-gain", dest="clado_gene_gain", type=float, default=2.0,
                   metavar="MEAN",
                   help="mean number of new families a daughter gains AT each speciation "
                        "(Poisson; default 2.0)")

    g = p.add_argument_group("gene-conditioned trait",
                             "--couple genes:traits (on a GIVEN tree; needs -t)")
    g.add_argument("--modifier-gain", dest="modifier_gain", type=float, default=0.5, metavar="RATE",
                   help="rate the modifier gene is gained (absent -> present; default 0.5)")
    g.add_argument("--modifier-loss", dest="modifier_loss", type=float, default=0.5, metavar="RATE",
                   help="rate the modifier gene is lost (present -> absent; default 0.5)")
    g.add_argument("--root-modifier", dest="root_modifier", action="store_true",
                   help="start with the modifier gene present at the root")
    g.add_argument("--theta-absent", dest="theta_absent", type=float, default=0.0, metavar="T",
                   help="the trait's OU optimum while the modifier is absent (default 0)")
    g.add_argument("--theta-present", dest="theta_present", type=float, default=5.0, metavar="T",
                   help="the trait's OU optimum while the modifier is present (default 5) — "
                        "acquiring the gene pulls the trait toward this peak")
    g.add_argument("--trait-alpha", dest="trait_alpha", type=float, default=1.0, metavar="A",
                   help="OU mean-reversion strength of the trait (0 = Brownian; default 1.0)")
    g.add_argument("--trait-sigma2", dest="trait_sigma2", type=float, default=1.0, metavar="S2",
                   help="trait diffusion rate sigma^2 (default 1.0)")
    g.add_argument("--trait-x0", dest="trait_x0", type=float, default=None, metavar="X0",
                   help="root trait value (default: the optimum of the root modifier state)")

    _add_traits_genes_args(p)   # --couple traits:genes (a trait conditions a gene-family panel)

def _build_anagenetic_trait(args: argparse.Namespace, parser: argparse.ArgumentParser):
    """The along-branch trait model for the ``species:traits`` edge on a **given** tree, taken from
    ``--sse-model``. No diversification happens on a fixed tree, so only the transition/diffusion
    structure is used (bisse/musse -> the Q as an :class:`Mk`; quasse -> Brownian ``--diffusion``);
    the speciation/extinction rates are inactive here. Returns ``(model, kind_label)``."""
    if args.sse_model == "quasse":
        return BrownianMotion(sigma2=args.diffusion, x0=args.root_value), "continuous"
    if args.sse_model == "musse":
        if args.q_matrix is None:
            parser.error("species:traits with --sse-model musse needs --q-matrix (the k-state "
                         "anagenetic transition matrix)")
        return Mk(_read_q_matrix(args.q_matrix)), "discrete"
    # bisse -> a binary Mk from the q01/q10 rates
    return Mk([[0.0, args.q01], [args.q10, 0.0]]), "discrete"

def _build_sse_model(args: argparse.Namespace, parser: argparse.ArgumentParser):
    """Construct the traits:species (SSE) model selected by ``--sse-model`` from the CLI args."""
    if args.sse_model == "bisse":
        return BiSSE(args.lambda0, args.lambda1, args.mu0, args.mu1, args.q01, args.q10)
    if args.sse_model == "musse":
        if args.birth is None or args.death is None or args.q_matrix is None:
            parser.error("--sse-model musse needs --birth and --death (k rates each) plus "
                         "--q-matrix (a k x k transition-rate matrix file)")
        return MuSSE(birth=args.birth, death=args.death, Q=_read_q_matrix(args.q_matrix))
    if args.sse_model == "hisse":
        H = args.hidden_classes
        if H < 2:
            parser.error("--sse-model hisse needs --hidden-classes >= 2")
        if args.hidden_scale <= 0:
            parser.error("--hidden-scale must be > 0")
        # H hidden classes spanning the base rates up to hidden_scale x faster (geometric); the
        # observed trait (q01/q10) and extinction (mu0/mu1) are shared across classes.
        factors = np.geomspace(1.0, args.hidden_scale, H)
        classes = [BiSSE(args.lambda0 * f, args.lambda1 * f, args.mu0, args.mu1, args.q01, args.q10)
                   for f in factors]
        return HiSSE(classes, hidden_transition=args.hidden_switch)
    # quasse: sigmoidal speciation in the trait + constant extinction (bounded for exact thinning)
    spec = QuaSSE.sigmoid(args.spec_low, args.spec_high, args.spec_center, args.spec_slope)
    bound = max(args.spec_low, args.spec_high) + args.qmu
    return QuaSSE(spec, lambda x: args.qmu, sigma2=args.diffusion,
                  rate_bound=bound, x0=args.root_value)

def _sse_tip_signal(res: TraitResult) -> str:
    """A short summary of the tip-state distribution — the diversification signal, for the log."""
    vals = list(res.labeled_values().values())
    if not vals:
        return ""
    if res.kind == "continuous":
        return f", tip trait mean {sum(vals) / len(vals):.3g}"
    from collections import Counter
    counts = Counter(vals)
    total = len(vals)
    frac = " ".join(f"{k}:{100 * n / total:.0f}%"
                    for k, n in sorted(counts.items(), key=lambda kv: str(kv[0])))
    return f", tip states {frac}"

def _write_coevolve_outputs(out: str, tree: Tree, res: TraitResult) -> None:
    """Write the shared coevolve outputs: the tree, the trait at every node, and the annotated
    trait tree."""
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick() + "\n")
    with open(os.path.join(out, "traits.tsv"), "w") as f:
        f.write(res.to_tsv(nodes="all"))              # every node: tips AND ancestral states
    with open(os.path.join(out, "trait_tree.nwk"), "w") as f:
        f.write(res.to_newick() + "\n")               # trait annotated on every node

def _write_null_manifest(out: str, edge: str, kind: str, *, cut: str, preserved: str,
                         extra: dict | None = None) -> None:
    """Record the null's provenance — which arrow was cut and how the target's variance was kept —
    so a downstream calibration is self-documenting. See docs/guide/coevolution_nulls.md."""
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "null_manifest.tsv"), "w") as f:
        f.write("field\tvalue\n")
        f.write(f"null\t{kind}\n")
        f.write(f"edge\t{edge}\n")
        f.write(f"cut_arrow\t{cut}\n")
        f.write(f"variance_preserved_by\t{preserved}\n")
        for k, v in (extra or {}).items():
            f.write(f"{k}\t{v}\n")

def _null_sse_model(model, args: argparse.Namespace, parser: argparse.ArgumentParser):
    """Apply ``--null`` to a traits:species SSE model; returns ``(null_model, manifest_kwargs)``."""
    if args.null == "neutral":
        return model.null("neutral"), dict(cut="trait -> (lambda, mu)",
                                            preserved="removed (rates no longer depend on the trait)")
    if args.null == "cid":
        if isinstance(model, QuaSSE):
            parser.error("--null cid needs a discrete character (bisse/musse); QuaSSE is "
                         "continuous — use --null neutral")
        if isinstance(model, HiSSE):
            parser.error("--sse-model hisse is already a hidden-state model; its null is "
                         "--null neutral (a constant-rate tree)")
        # bisse -> the binary CID; musse -> the k-state CID (both hide the hidden class in output)
        return (model.null("cid", n_hidden=args.hidden),
                dict(cut="trait -> (lambda, mu)",
                     preserved=f"a hidden CID class (CID-{args.hidden})",
                     extra={"hidden_classes": args.hidden}))
    parser.error("traits:species has no 'timing' null; use --null neutral or --null cid")

def _null_species_traits_timing(args: argparse.Namespace, tree: Tree,
                                parser: argparse.ArgumentParser):
    """The species:traits ``timing`` null: drop cladogenesis, add a matched anagenetic rate so the
    same expected change is spread **along branches** (analytic, from the tree's branch stats).
    Returns ``(anagenetic_model, manifest_extra)``."""
    n_branches, total_len = _branch_count_and_length(tree)
    if total_len <= 0.0:
        parser.error("tree has zero total branch length; cannot spread the cladogenetic change")
    per_len = n_branches / total_len
    if args.sse_model == "quasse":                         # continuous: boost Brownian diffusion
        extra = args.clado_jump * per_len
        return (BrownianMotion(sigma2=args.diffusion + extra, x0=args.root_value),
                {"matched_sigma2_extra": f"{extra:.6g}"})
    if args.sse_model == "bisse":                          # binary: boost the Mk transition rates
        s = args.clado_shift * per_len
        return Mk([[0.0, args.q01 + s], [args.q10 + s, 0.0]]), {"matched_rate_extra": f"{s:.6g}"}
    parser.error("--null timing for species:traits supports --sse-model bisse or quasse (a matched "
                 "anagenetic rate); musse timing is not implemented — use --null neutral")

_NULL_OVERLAY = dict(duplication=0.0, transfer=1.0, loss=0.5, origination=0.0)   # neutral genome

def _run_genes_species_cid_null(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """genes:species CID null: the drivers shape a genuinely heterogeneous tree; a NEUTRAL overlay
    genome is the decoupled observed channel; the drivers are withheld as ground-truth."""
    if args.tree:
        parser.error("genes:species grows the tree (it is an OUTPUT); give --age/--tips, not -t")
    if (args.age is None) == (args.tips is None):
        parser.error("genes:species grows the tree — give exactly one of --age or --tips")
    model = GeneDiversification(
        args.drivers, lambda0=args.lambda0, mu0=args.mu0,
        driver_speciation=args.driver_speciation, driver_extinction=args.driver_extinction,
        loss=args.driver_loss, origination=args.driver_origination,
        transfer=args.driver_transfer, root_drivers=args.root_drivers)
    t0 = time.perf_counter()
    res = simulate_gene_diversification(model, age=args.age, n_tips=args.tips, seed=args.seed)
    tree = prune(res.tree)                                  # extant-only, for a clean null dataset
    seed2 = None if args.seed is None else args.seed + 1
    profiles = simulate_genomes(tree, output="profiles", initial_families=args.genome_size,
                                seed=seed2, **_NULL_OVERLAY)
    dt = time.perf_counter() - t0

    _write_profiles_only(args.out, tree, profiles)         # tree + OBSERVED neutral genome
    with open(os.path.join(args.out, "drivers_ground_truth.tsv"), "w") as f:
        f.write(res.to_tsv(nodes="all"))                   # the HIDDEN drivers (withheld from analysis)
    _write_drivers_manifest(args.out, model)
    _write_null_manifest(
        args.out, "genes:species", "cid", cut="gene content -> (lambda, mu)",
        preserved="a hidden driver panel (drivers_ground_truth.tsv) shaped the tree; Profiles.tsv is "
                  "a neutral overlay genome, decoupled from diversification",
        extra={"observed": f"neutral overlay: {args.genome_size} families, transfer=1.0 loss=0.5",
               "hidden_drivers": model.n_drivers})
    return (f"wrote genes:species [cid null] to {args.out}/ ({len(tree.leaves())} tips; "
            f"{len(profiles.families)} neutral observed families, {model.n_drivers} hidden drivers) "
            f"in {dt:.3g} s")

def _run_genes_traits_cid_null(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """genes:traits CID null: a modifier gene shapes a trait with real optimum shifts; a NEUTRAL
    overlay genome is the observed gene content; the modifier is withheld as ground-truth."""
    if not args.tree:
        parser.error("genes:traits runs on a GIVEN tree — pass -t/--tree")
    if args.age is not None or args.tips is not None:
        parser.error("genes:traits uses the given -t tree; --age/--tips only apply to into-species edges")
    with open(args.tree) as f:
        tree = read_newick(f.read())
    model = GeneConditionedTrait(
        gene_gain=args.modifier_gain, gene_loss=args.modifier_loss, root_gene=args.root_modifier,
        theta_absent=args.theta_absent, theta_present=args.theta_present,
        alpha=args.trait_alpha, sigma2=args.trait_sigma2, x0=args.trait_x0)
    t0 = time.perf_counter()
    res = simulate_gene_conditioned_trait(tree, model, seed=args.seed)     # modifier shapes the trait
    seed2 = None if args.seed is None else args.seed + 1
    profiles = simulate_genomes(tree, output="profiles", initial_families=args.genome_size,
                                seed=seed2, **_NULL_OVERLAY)
    dt = time.perf_counter() - t0

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick() + "\n")
    with open(os.path.join(args.out, "Profiles.tsv"), "w") as f:
        f.write(profiles.to_tsv())                         # OBSERVED neutral genome
    with open(os.path.join(args.out, "Presence.tsv"), "w") as f:
        f.write(profiles.to_tsv(presence=True))
    with open(os.path.join(args.out, "traits.tsv"), "w") as f:  # the trait (real optimum shifts), observed
        f.write("node\ttrait\n")
        for n in tree.nodes():
            f.write(f"{n.name}\t{res.node_trait[n]:.6g}\n")
    with open(os.path.join(args.out, "modifier_ground_truth.tsv"), "w") as f:   # the HIDDEN modifier
        f.write("node\tmodifier\n")
        for n in tree.nodes():
            f.write(f"{n.name}\t{int(res.gene.node_values[n])}\n")
    _write_null_manifest(
        args.out, "genes:traits", "cid", cut="gene presence -> trait optimum",
        preserved="a hidden modifier (modifier_ground_truth.tsv) shifted the trait's optimum; "
                  "Profiles.tsv is a neutral overlay genome, decoupled from the trait",
        extra={"observed": f"neutral overlay: {args.genome_size} families, transfer=1.0 loss=0.5"})
    return (f"wrote genes:traits [cid null] to {args.out}/ ({len(tree.extant_leaves())} tips; "
            f"{len(profiles.families)} neutral observed families, hidden modifier) in {dt:.3g} s")

def _run_traits_genes_cid_null(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """traits:genes CID null: a HIDDEN trait drives the panel's retention; a SECOND, independent
    neutral trait is the observed trait; the driving trait is withheld as ground-truth."""
    if not args.tree:
        parser.error("traits:genes runs on a GIVEN tree — pass -t/--tree")
    if args.age is not None or args.tips is not None:
        parser.error("traits:genes uses the given -t tree; --age/--tips only apply to into-species edges")
    if args.trait_file:
        parser.error("--null cid for traits:genes simulates its own hidden + observed traits; "
                     "--trait-file is not supported here")
    with open(args.tree) as f:
        tree = read_newick(f.read())
    parts = set(Genomes.WRITE_PARTS) if "all" in args.write else set(args.write)
    if args.panel < 1:
        raise ValueError("--panel must be >= 1")
    rng = np.random.default_rng(args.seed)

    hidden_trait = simulate_traits(tree, _build_trait_model(args), rng=rng)   # drives the panel
    state_values = None
    if args.trait_center and hidden_trait.kind == "discrete":
        k = len(hidden_trait.model.states)
        state_values = [i - (k - 1) / 2.0 for i in range(k)]
    coupling = TraitGeneCoupling.build(
        args.panel, _parse_responsive(args.responsive), weight=args.weight, signed=args.signed,
        effect_loss=args.effect_loss, effect_gain=args.effect_gain, base_loss=args.loss,
        transfer=args.trans, duplication=args.dup, origination=args.orig,
        state_values=state_values, rng=rng)
    t0 = time.perf_counter()
    res = simulate_trait_linked_genomes(tree, hidden_trait, coupling, trait_steps=args.trait_steps,
                                        rng=rng)
    observed_trait = simulate_traits(tree, _build_trait_model(args), rng=rng)   # independent, decoupled
    dt = time.perf_counter() - t0

    os.makedirs(args.out, exist_ok=True)
    res.genomes().write(args.out, include=parts, sparse=args.sparse,
                        annotate_species=args.annotate_species)   # the panel (shaped by hidden trait)
    with open(os.path.join(args.out, "traits.tsv"), "w") as f:
        f.write(observed_trait.to_tsv(nodes="all"))               # OBSERVED, decoupled trait
    with open(os.path.join(args.out, "trait_ground_truth.tsv"), "w") as f:
        f.write(hidden_trait.to_tsv(nodes="all"))                 # the HIDDEN driving trait
    _write_coupling_manifest(args.out, coupling)
    _write_null_manifest(
        args.out, "traits:genes", "cid", cut="trait -> panel loss/gain",
        preserved="a hidden trait (trait_ground_truth.tsv) shaped the panel; traits.tsv is a second, "
                  "independent neutral trait, decoupled from the panel",
        extra={"observed_trait": args.model, "responsive_families": coupling.n_responsive})
    return (f"wrote traits:genes [cid null] to {args.out}/ (observed trait decoupled from a "
            f"{coupling.n_families}-family panel shaped by a hidden trait, "
            f"{len(tree.extant_leaves())} tips) in {dt:.3g} s")

def _run_coevolve_mode(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """Run the ``coevolve`` umbrella over the six directed edges (``--couple``): ``traits:species``
    (SSE), ``species:traits`` (cladogenetic) and their combination = **ClaSSE**; ``genes:species``
    (key innovations); ``species:genes`` (cladogenetic genome); ``genes:traits`` (gene-conditioned
    trait); and ``traits:genes`` (trait-conditioned genes). Each node-pair's two edges can also be
    combined: ``genes:species``+``species:genes`` = **co-diversification** and
    ``traits:genes``+``genes:traits`` = **trait-gene feedback** (as ``traits:species``+
    ``species:traits`` = ClaSSE). Whether the tree is grown (an arrow into species) or read from
    ``-t`` follows the arrows-into-S rule."""
    # --couple accepts both repeated flags and space-separated lists (append + nargs); flatten
    raw = args.couple or [["traits:species"]]
    edges = [e.strip().lower() for group in raw for e in group]
    for e in edges:
        if e not in _COEVOLVE_EDGES:
            parser.error(f"unknown --couple edge {e!r}: expected 'driver:target' over "
                         f"{{{', '.join(_COEVOLVE_NODES)}}} (e.g. traits:species); see "
                         "docs/guide/coevolution.md for the full edge set")
    eset = set(edges)

    # --null cuts a single directed arrow; a joint/both-arrows model has two, so decline it there
    # (run the null for one arrow at a time). Each edge below applies its own archetype.
    if args.null != "none" and len(eset) > 1:
        parser.error("--null cuts a single directed edge; a joint (both-arrows) model has two — "
                     "run the null for one arrow at a time (e.g. --couple traits:species --null cid)")

    # ---- joint (both-arrow) models: a node-pair with BOTH its directed edges on ----
    # species<->genes: driver gene content drives diversification AND bursts at each speciation
    # (one arrow into S -> the tree is an OUTPUT), the genomic analogue of ClaSSE.
    if eset == {"genes:species", "species:genes"}:
        return _run_co_diversification(args, parser)

    # traits<->genes: the trait and a coupled panel modulate each other (no arrow into S -> an
    # overlay on a given tree). The closed feedback loop writes a trait<->gene tip association.
    if eset == {"traits:genes", "genes:traits"}:
        return _run_trait_gene_feedback(args, parser)

    # traits:genes — a trait conditions a gene-family panel's loss/gain (formerly the standalone
    # 'coevolve-genetrait' command). An overlay edge (no arrow into S), so the tree is an INPUT.
    if "traits:genes" in eset:
        if eset != {"traits:genes"}:
            parser.error("traits:genes combines only with genes:traits (the trait-gene feedback "
                         "joint model); other combinations are future work — see "
                         "docs/guide/coevolution.md")
        return _run_traits_genes(args, parser)

    # genes:species — gene content drives diversification (a forward joint loop). Combines with
    # species:genes above (co-diversification); other combinations are still on the roadmap.
    if "genes:species" in eset:
        if eset != {"genes:species"}:
            parser.error("genes:species combines only with species:genes (the co-diversification "
                         "joint model); other combinations are future work — see "
                         "docs/guide/coevolution.md")
        return _run_genes_species(args, parser)

    # genes:traits — gene content conditions a trait (a modifier gene switches the trait's OU
    # optimum). An overlay edge (no arrow into S), so the tree is an INPUT; runs on a given -t tree.
    if "genes:traits" in eset:
        if eset != {"genes:traits"}:
            parser.error("genes:traits combines only with traits:genes (the trait-gene feedback "
                         "joint model); other combinations are future work — see "
                         "docs/guide/coevolution.md")
        return _run_genes_traits(args, parser)

    # species:genes — speciation drives gene content (cladogenetic genome). An overlay edge (no
    # arrow into S), so the tree is an INPUT; runs on a given -t tree.
    if "species:genes" in eset:
        if eset != {"species:genes"}:
            parser.error("species:genes combines only with genes:species (the co-diversification "
                         "joint model); other combinations are future work — see "
                         "docs/guide/coevolution.md")
        return _run_species_genes(args, parser)

    traits_species = "traits:species" in eset      # SSE arrow (trait -> diversification), into S
    species_traits = "species:traits" in eset      # cladogenetic arrow (speciation -> trait)
    clado = (Cladogenesis(jump_sigma2=args.clado_jump, shift=args.clado_shift)
             if species_traits else None)

    # species:traits ALONE — no arrow into S, so the tree is an INPUT (nothing to grow): evolve the
    # trait along the given tree with cladogenetic jumps at its speciation nodes.
    if species_traits and not traits_species:
        if not args.tree:
            parser.error("species:traits alone runs on a GIVEN tree — pass -t/--tree (no "
                         "diversification happens on this edge, so there is nothing to grow)")
        if args.age is not None or args.tips is not None:
            parser.error("species:traits alone uses the given -t tree; --age/--tips only apply to "
                         "the into-species edges that grow a tree")
        with open(args.tree) as f:
            tree = read_newick(f.read())
        model, kind = _build_anagenetic_trait(args, parser)
        nman = None
        if args.null == "neutral":                         # drop the at-speciation jump entirely
            clado = None
            nman = dict(cut="speciation -> trait jump", preserved="removed (anagenetic change only)")
        elif args.null == "timing":                        # spread the jump along branches (matched)
            clado = None
            model, extra = _null_species_traits_timing(args, tree, parser)
            nman = dict(cut="speciation -> trait jump",
                        preserved="spread along branches (matched anagenetic rate)", extra=extra)
        elif args.null == "cid":
            parser.error("species:traits has no 'cid' null (its driver is the speciation process, "
                         "not a state); use --null neutral or --null timing")
        t0 = time.perf_counter()
        res = simulate_traits(tree, model, cladogenesis=clado,
                              root_state=args.root_state, seed=args.seed)
        dt = time.perf_counter() - t0
        _write_coevolve_outputs(args.out, tree, res)
        if nman is not None:
            _write_null_manifest(args.out, "species:traits", args.null, **nman)
        tag = f" [{args.null} null]" if args.null != "none" else ""
        mode = "anagenetic" if args.null != "none" else "cladogenetic"
        return (f"wrote species:traits{tag} ({mode} {kind}) to {args.out}/ "
                f"({len(tree.extant_leaves())} tips{_sse_tip_signal(res)}) in {dt:.3g} s")

    # traits:species (SSE) or traits:species + species:traits (ClaSSE): an arrow INTO S, so the
    # tree is an OUTPUT — grow it forward with exactly one stopping condition (no input -t tree).
    if args.tree:
        parser.error("traits:species grows the tree (it is an OUTPUT); give --age/--tips, not an "
                     "input -t tree (that is the species:traits-alone edge)")
    if (args.age is None) == (args.tips is None):
        parser.error("traits:species grows the tree — give exactly one of --age or --tips")

    model = _build_sse_model(args, parser)
    nman = None
    if args.null != "none":                                # cut trait -> diversification
        model, nman = _null_sse_model(model, args, parser)
    t0 = time.perf_counter()
    res = simulate_sse(model, age=args.age, n_tips=args.tips, root_state=args.root_state,
                       cladogenesis=clado, seed=args.seed)
    dt = time.perf_counter() - t0
    _write_coevolve_outputs(args.out, res.tree, res)
    if nman is not None:
        _write_null_manifest(args.out, "traits:species", args.null, **nman)
    n_extant = len(res.tree.extant_leaves())
    tag = f" [{args.null} null]" if args.null != "none" else ""
    edge_label = ("traits:species+species:traits" if clado is not None else "traits:species") + tag
    model_label = f"ClaSSE {args.sse_model}" if clado is not None else f"SSE {args.sse_model}"
    return (f"wrote {edge_label} ({model_label}) to {args.out}/ "
            f"({n_extant} extant tips{_sse_tip_signal(res)}) in {dt:.3g} s")

def _write_drivers_manifest(out: str, model: GeneDiversification) -> None:
    """Write ``drivers_manifest.tsv`` — the per-driver effect sizes and rates behind the tree, so
    the gene↔diversification linkage that shaped the profiles is on record for inference."""
    root = ",".join(f"D{i}" for i in sorted(model.root_set)) or "-"
    lines = [f"# lambda0\t{model.lambda0:g}", f"# mu0\t{model.mu0:g}",
             f"# loss\t{model.loss:g}", f"# origination\t{model.origination:g}",
             f"# transfer\t{model.transfer:g}", f"# root_drivers\t{root}",
             "driver\tbeta_speciation\tbeta_extinction"]
    for i in range(model.n_drivers):
        lines.append(f"D{i}\t{model.beta_lambda[i]:.6g}\t{model.beta_mu[i]:.6g}")
    with open(os.path.join(out, "drivers_manifest.tsv"), "w") as f:
        f.write("\n".join(lines) + "\n")

def _run_genes_species(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """Grow a tree whose diversification is driven by a panel of binary key-innovation gene
    families (``genes:species``); the neutral genome is overlaid afterward with ``zombi2 genomes``
    on the resulting tree (exact under independent families)."""
    if args.tree:
        parser.error("genes:species grows the tree (it is an OUTPUT); give --age/--tips, not an "
                     "input -t tree")
    if (args.age is None) == (args.tips is None):
        parser.error("genes:species grows the tree — give exactly one of --age or --tips")

    model = GeneDiversification(
        args.drivers, lambda0=args.lambda0, mu0=args.mu0,
        driver_speciation=args.driver_speciation, driver_extinction=args.driver_extinction,
        loss=args.driver_loss, origination=args.driver_origination,
        transfer=args.driver_transfer, root_drivers=args.root_drivers)
    if args.null == "cid":
        return _run_genes_species_cid_null(args, parser)
    if args.null == "timing":
        parser.error("genes:species has no 'timing' null; use --null neutral or --null cid")
    if args.null == "neutral":
        model = model.null("neutral")
    t0 = time.perf_counter()
    res = simulate_gene_diversification(model, age=args.age, n_tips=args.tips, seed=args.seed)
    dt = time.perf_counter() - t0

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
        f.write(res.tree.to_newick() + "\n")          # the tree the drivers shaped
    with open(os.path.join(args.out, "drivers.tsv"), "w") as f:
        f.write(res.to_tsv(nodes="all"))              # per-node driver presence (0/1 columns)
    _write_drivers_manifest(args.out, model)
    n_extant = len(res.tree.extant_leaves())
    prev = " ".join(f"D{i}:{100 * p:.0f}%" for i, p in enumerate(res.tip_prevalence()))
    print(f"  overlay the neutral genome with: zombi2 genomes -t {args.out}/species_tree.nwk "
          f"--trans 1 --loss 0.5 --write profiles trees -o {args.out}")
    if args.null == "neutral":
        _write_null_manifest(args.out, "genes:species", "neutral", cut="gene content -> (lambda, mu)",
                             preserved="removed (drivers no longer set the rates)")
    tag = " [neutral null]" if args.null == "neutral" else ""
    return (f"wrote genes:species{tag} (key innovations) to {args.out}/ "
            f"({n_extant} extant tips, {model.n_drivers} drivers, tip prevalence {prev}) "
            f"in {dt:.3g} s")

def _run_co_diversification(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """The species<->genes JOINT model (``genes:species`` + ``species:genes``): the same driver
    families both drive diversification AND are reshuffled by a cladogenetic burst at every
    speciation. One arrow points into S, so the tree is an OUTPUT (grown jointly)."""
    if args.tree:
        parser.error("genes:species+species:genes grows the tree (it is an OUTPUT); give "
                     "--age/--tips, not an input -t tree")
    if (args.age is None) == (args.tips is None):
        parser.error("the species<->genes joint model grows the tree — give exactly one of "
                     "--age or --tips")
    if args.driver_clado_loss <= 0.0 and args.driver_clado_gain <= 0.0:
        parser.error("the species:genes arrow needs a cladogenetic burst on the drivers: set "
                     "--driver-clado-loss and/or --driver-clado-gain > 0 (with both 0 there is no "
                     "species:genes coupling — that is plain genes:species)")

    model = GeneDiversification(
        args.drivers, lambda0=args.lambda0, mu0=args.mu0,
        driver_speciation=args.driver_speciation, driver_extinction=args.driver_extinction,
        loss=args.driver_loss, origination=args.driver_origination,
        transfer=args.driver_transfer, root_drivers=args.root_drivers,
        cladogenetic_loss=args.driver_clado_loss, cladogenetic_gain=args.driver_clado_gain)
    t0 = time.perf_counter()
    res = simulate_co_diversification(model, age=args.age, n_tips=args.tips, seed=args.seed)
    dt = time.perf_counter() - t0

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
        f.write(res.tree.to_newick() + "\n")          # the tree the drivers shaped and burst
    with open(os.path.join(args.out, "drivers.tsv"), "w") as f:
        f.write(res.to_tsv(nodes="all"))              # per-node driver presence (0/1 columns)
    _write_drivers_manifest(args.out, model)
    n_extant = len(res.tree.extant_leaves())
    prev = " ".join(f"D{i}:{100 * p:.0f}%" for i, p in enumerate(res.tip_prevalence()))
    print(f"  overlay the neutral genome with: zombi2 genomes -t {args.out}/species_tree.nwk "
          f"--trans 1 --loss 0.5 --write profiles trees -o {args.out}")
    return (f"wrote genes:species+species:genes (co-diversification) to {args.out}/ "
            f"({n_extant} extant tips, {model.n_drivers} drivers, tip prevalence {prev}) "
            f"in {dt:.3g} s")

def _run_species_genes(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """Evolve a genome down a GIVEN tree with a cladogenetic ('punctuational') burst of gene loss
    and gain at every speciation (the ``species:genes`` edge — speciation drives gene content)."""
    if not args.tree:
        parser.error("species:genes runs on a GIVEN tree — pass -t/--tree (speciation drives the "
                     "genome; there is no diversification to grow here)")
    if args.age is not None or args.tips is not None:
        parser.error("species:genes uses the given -t tree; --age/--tips only apply to the "
                     "into-species edges that grow a tree")
    with open(args.tree) as f:
        tree = read_newick(f.read())
    model = CladogeneticGenome(
        initial_families=args.genome_size, loss=args.gene_loss, origination=args.gene_origination,
        cladogenetic_loss=args.clado_gene_loss, cladogenetic_gain=args.clado_gene_gain)
    nman = None
    if args.null == "cid":
        parser.error("species:genes has no 'cid' null (its driver is the speciation process, not a "
                     "state); use --null neutral or --null timing")
    if args.null == "neutral":
        model = model.null("neutral")
        nman = dict(cut="speciation -> gene burst", preserved="removed (anagenetic turnover only)")
    elif args.null == "timing":
        model = model.null("timing", tree=tree)
        nman = dict(cut="speciation -> gene burst",
                    preserved="spread along branches (matched anagenetic rate)",
                    extra={"loss": f"{model.loss:.6g}", "origination": f"{model.origination:.6g}"})
    t0 = time.perf_counter()
    res = simulate_cladogenetic_genome(tree, model, seed=args.seed)
    dt = time.perf_counter() - t0

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick() + "\n")              # the given tree, for provenance
    pm = res.profile_matrix()
    with open(os.path.join(args.out, "Profiles.tsv"), "w") as f:
        f.write(pm.to_tsv())
    with open(os.path.join(args.out, "Presence.tsv"), "w") as f:
        f.write(pm.to_tsv(presence=True))
    sizes = res.genome_sizes()
    with open(os.path.join(args.out, "genome_sizes.tsv"), "w") as f:
        f.write("node\tgenome_size\n")
        for node in tree.nodes():
            f.write(f"{node.name}\t{sizes[node]}\n")
    tips = tree.extant_leaves()
    mean_size = sum(sizes[t] for t in tips) / len(tips) if tips else 0
    if nman is not None:
        _write_null_manifest(args.out, "species:genes", args.null, **nman)
    tag = f" [{args.null} null]" if args.null != "none" else ""
    return (f"wrote species:genes{tag} (cladogenetic genome) to {args.out}/ "
            f"({len(tips)} tips, {len(pm.families)} families, mean genome {mean_size:.0f}) "
            f"in {dt:.3g} s")

def _run_genes_traits(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """Evolve a trait down a GIVEN tree whose OU optimum is switched by a modifier gene's presence
    (the ``genes:traits`` edge — gene content conditions the trait)."""
    if not args.tree:
        parser.error("genes:traits runs on a GIVEN tree — pass -t/--tree (gene content conditions "
                     "the trait; there is no diversification to grow here)")
    if args.age is not None or args.tips is not None:
        parser.error("genes:traits uses the given -t tree; --age/--tips only apply to the "
                     "into-species edges that grow a tree")
    with open(args.tree) as f:
        tree = read_newick(f.read())
    model = GeneConditionedTrait(
        gene_gain=args.modifier_gain, gene_loss=args.modifier_loss, root_gene=args.root_modifier,
        theta_absent=args.theta_absent, theta_present=args.theta_present,
        alpha=args.trait_alpha, sigma2=args.trait_sigma2, x0=args.trait_x0)
    if args.null == "cid":
        return _run_genes_traits_cid_null(args, parser)
    if args.null == "timing":
        parser.error("genes:traits has no 'timing' null; use --null neutral")
    if args.null == "neutral":
        model = model.null("neutral")
    t0 = time.perf_counter()
    res = simulate_gene_conditioned_trait(tree, model, seed=args.seed)
    dt = time.perf_counter() - t0

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick() + "\n")              # the given tree, for provenance
    with open(os.path.join(args.out, "traits.tsv"), "w") as f:
        f.write(res.to_tsv(nodes="all"))              # per-node modifier presence + trait value
    with open(os.path.join(args.out, "trait_tree.nwk"), "w") as f:
        f.write(res.to_newick() + "\n")               # trait annotated on every node

    tips = tree.extant_leaves()
    tv, gp = res.trait_values(), res.gene_presence()
    car = [tv[t] for t in tips if gp[t]]
    non = [tv[t] for t in tips if not gp[t]]
    car_m = f"{sum(car) / len(car):.2g}" if car else "-"
    non_m = f"{sum(non) / len(non):.2g}" if non else "-"
    if args.null == "neutral":
        _write_null_manifest(args.out, "genes:traits", "neutral", cut="gene presence -> trait optimum",
                             preserved="removed (theta_present = theta_absent; plain OU)")
    tag = " [neutral null]" if args.null == "neutral" else ""
    return (f"wrote genes:traits{tag} (gene-conditioned trait) to {args.out}/ "
            f"({len(tips)} tips; carrier trait mean {car_m} vs non-carrier {non_m}) in {dt:.3g} s")

def _run_trait_gene_feedback(args: argparse.Namespace, parser: argparse.ArgumentParser) -> str:
    """The traits<->genes JOINT model (``traits:genes`` + ``genes:traits``): the trait and a coupled
    panel evolve together, each modulating the other. An overlay (no arrow into S), so it needs a
    given ``-t`` tree. Reuses --panel/--effect-loss/--loss/--trans for the panel and
    --theta-absent/--theta-present/--trait-alpha/--trait-sigma2 for the trait (see --panel-root-fraction)."""
    if not args.tree:
        parser.error("traits:genes+genes:traits runs on a GIVEN tree — pass -t/--tree (neither "
                     "arrow points into S, so there is nothing to grow)")
    if args.age is not None or args.tips is not None:
        parser.error("the traits<->genes joint model uses the given -t tree; --age/--tips only "
                     "apply to the into-species edges that grow a tree")
    with open(args.tree) as f:
        tree = read_newick(f.read())
    model = TraitGeneFeedback(
        n_families=args.panel, effect_loss=args.effect_loss, base_loss=args.loss, gain=args.trans,
        theta_low=args.theta_absent, theta_high=args.theta_present,
        alpha=args.trait_alpha, sigma2=args.trait_sigma2, x0=args.trait_x0,
        root_fraction=args.panel_root_fraction, steps=args.trait_steps)
    t0 = time.perf_counter()
    res = simulate_trait_gene_feedback(tree, model, seed=args.seed)
    dt = time.perf_counter() - t0

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick() + "\n")                  # the given tree, for provenance
    with open(os.path.join(args.out, "Profiles.tsv"), "w") as f:
        f.write(res.profiles.to_tsv(presence=True))       # panel presence at the extant tips
    with open(os.path.join(args.out, "traits.tsv"), "w") as f:
        f.write("node\ttrait\tpanel_occupancy\n")         # the coupled trait + panel at every node
        for n in tree.nodes():
            f.write(f"{n.name}\t{res.node_trait[n]:.6g}\t{res.node_presence[n].mean():.6g}\n")

    corr = res.trait_gene_correlation()
    corr_s = f"{corr:.2f}" if corr == corr else "n/a"
    return (f"wrote traits:genes+genes:traits (trait-gene feedback) to {args.out}/ "
            f"({len(tree.extant_leaves())} tips; tip trait-panel corr {corr_s}) in {dt:.3g} s")


def run(args, parser):
    summary = _run_coevolve_mode(args, parser)
    print(summary)
    _write_params_log(os.path.join(args.out, "coevolve.log"), args, summary)
    return 0
