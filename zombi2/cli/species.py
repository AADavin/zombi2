"""zombi2 species command."""
from __future__ import annotations

import argparse
import os
import time



from zombi2.species.ghosts import add_ghost_lineages
from zombi2.species.model import (
    BirthDeath, CladeShiftBirthDeath, ClaDS, DiversityDependent, EpisodicBirthDeath,
)
from zombi2.species.sim import simulate_species_tree
from zombi2.tree import prune

from zombi2.cli.framework import _add_params_arg, _write_params_log

def _add_species_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("general")
    _add_params_arg(g)
    g.add_argument("--mode", dest="model", choices=("backward", "forward"), default="backward",
                   metavar="MODE",
                   help="backward: reconstructed tree conditioned on --tips extant species "
                        "(default); forward: complete tree grown in time, keeping extinct "
                        "lineages (and fossils)")
    g.add_argument("--tips", type=int, default=None, metavar="N",
                   help="number of extant species (backward default 50; forward: --tips OR --age)")
    g.add_argument("--age", type=float, default=None, metavar="T",
                   help="tree age / timescale, in the same time units as the rates "
                        "(backward default 1.0; forward: --tips OR --age)")
    g.add_argument("--age-type", choices=("crown", "stem"), default="crown", metavar="KIND",
                   help="[backward] interpret --age as crown (default) or stem age")
    g.add_argument("--seed", type=int, default=None, metavar="N",
                   help="RNG seed for reproducibility")
    g.add_argument("-o", "--out", required=True, metavar="DIR", help="output directory")

    g = p.add_argument_group("diversification model",
                             "the rate process, chosen by --diversification (forward only)")
    g.add_argument("--diversification", choices=("constant", "clads", "diversity-dependent"),
                   default="constant", metavar="PROCESS",
                   help="constant-rate birth-death (default); clads = per-lineage rates that "
                        "shift at each speciation (ClaDS); diversity-dependent = rates decline "
                        "toward a carrying capacity")
    g.add_argument("--birth", type=float, nargs="+", default=[1.0], metavar="RATE",
                   help="speciation rate (default 1.0); several values with --shifts give an "
                        "episodic (skyline) model. For clads/diversity-dependent it is the "
                        "root/intrinsic rate λ₀ (a single value)")
    g.add_argument("--death", type=float, nargs="+", default=[0.3], metavar="RATE",
                   help="extinction rate (default 0.3); several values with --shifts give an "
                        "episodic (skyline) model. The constant μ for --diversification "
                        "diversity-dependent (clads uses --turnover instead)")
    g.add_argument("--shifts", type=float, nargs="+", default=None, metavar="AGE",
                   help="[episodic] rate-shift ages, present -> past (K-1 ages for K rate values)")

    g = p.add_argument_group("clads model", "with --diversification clads")
    g.add_argument("--clads-alpha", type=float, default=0.9, metavar="ALPHA",
                   help="speciation-rate trend per branch; α<1 = rates slow toward the present "
                        "(default 0.9)")
    g.add_argument("--clads-sigma", type=float, default=0.1, metavar="SIGMA",
                   help="lognormal spread of the per-branch rate shift (default 0.1)")
    g.add_argument("--turnover", type=float, default=0.0, metavar="EPS",
                   help="extinction/speciation ratio ε=μ/λ, in [0,1) (0 = pure birth; default 0.0)")

    g = p.add_argument_group("diversity-dependent model",
                             "with --diversification diversity-dependent")
    g.add_argument("--carrying-capacity", "-K", type=float, default=None, metavar="K",
                   help="carrying capacity K; the speciation rate is λ₀·(1−n/K) (required for "
                        "this model)")

    g = p.add_argument_group("clade-specific shifts", "forward only")
    g.add_argument("--clade-shift", action="append", nargs=3, type=float,
                   metavar=("AGE", "BIRTH", "DEATH"), default=None, dest="clade_shift",
                   help="a clade-specific rate shift: at AGE before the present, one random "
                        "lineage then alive (and its descendants) switches to speciation BIRTH / "
                        "extinction DEATH. Repeat for several shifting clades, e.g. "
                        "--clade-shift 3.0 2.5 0.1")

    g = p.add_argument_group("forward sampling & fossils", "only with --mode forward")
    g.add_argument("--sampling-fraction", type=float, default=1.0, metavar="RHO",
                   help="fraction of extant species sampled, 0<rho<=1 (default 1.0)")
    g.add_argument("--fossilization", type=float, default=0.0, metavar="PSI",
                   help="fossil (serial) sampling rate psi — fossilized birth-death "
                        "(default 0 = no fossils)")
    g.add_argument("--removal", type=float, default=1.0, metavar="R",
                   help="removal probability on sampling, 0<=r<=1 (r<1 keeps sampled ancestors; "
                        "default 1.0)")
    g.add_argument("--mass-extinction", action="append", nargs=2, type=float,
                   metavar=("AGE", "FRACTION"), default=None, dest="mass_extinction",
                   help="a mass extinction: at AGE before the present, each lineage dies with "
                        "probability FRACTION (0<FRACTION<=1). Repeat for several pulses, e.g. "
                        "--mass-extinction 1.0 0.75 --mass-extinction 2.5 0.5")

    g = p.add_argument_group("ghost lineages", "backward only")
    g.add_argument("--ghosts", action="store_true",
                   help="graft the extinct/unsampled 'ghost' lineages back onto the tree")
    g.add_argument("--ghost-method", choices=("rejection", "htransform"), default="rejection",
                   metavar="METHOD",
                   help="ghost-subtree sampler used with --ghosts (default rejection)")

    g = p.add_argument_group("run limits", "forward only")
    g.add_argument("--max-attempts", type=int, default=10000, metavar="N",
                   help="retries before giving up when the process goes extinct (default 10000)")
    g.add_argument("--max-lineages", type=int, default=1_000_000, metavar="N",
                   help="abort a run exceeding this many live lineages (default 1000000)")

def _build_species_model(args: argparse.Namespace, parser: argparse.ArgumentParser):
    """Construct a species-tree model (BirthDeath / EpisodicBirthDeath / ClaDS /
    DiversityDependent) from the CLI args (validated)."""
    if args.model == "backward" and (args.fossilization or args.removal != 1.0
                                     or args.sampling_fraction != 1.0):
        parser.error("--fossilization / --removal / --sampling-fraction require --mode forward "
                     "(the backward reconstructed sampler assumes complete sampling)")
    if args.model == "backward" and args.mass_extinction:
        parser.error("--mass-extinction requires --mode forward (mass extinctions kill real "
                     "lineages forward in time; the backward reconstructed sampler never sees them)")
    # [(age, fraction), ...] pulses, or None; carried by whichever model is built
    mass_ext = args.mass_extinction

    if args.clade_shift and args.diversification != "constant":
        parser.error("--clade-shift is its own constant-background model; it does not combine "
                     "with --diversification clads/diversity-dependent")
    if args.diversification != "constant":
        return _build_heterogeneous_model(args, parser, mass_ext)
    if args.clade_shift:
        return _build_clade_shift_model(args, parser, mass_ext)

    episodic = args.shifts is not None or len(args.birth) > 1 or len(args.death) > 1
    if not episodic:
        return BirthDeath(args.birth[0], args.death[0], fossilization=args.fossilization,
                          sampling_fraction=args.sampling_fraction, removal=args.removal,
                          mass_extinctions=mass_ext)
    shifts = args.shifts or []
    if len(args.birth) != len(args.death) or len(shifts) != len(args.birth) - 1:
        parser.error("episodic model needs len(--birth) == len(--death) == len(--shifts)+1 "
                     f"(got {len(args.birth)} birth, {len(args.death)} death, {len(shifts)} shifts)")
    return EpisodicBirthDeath(birth=args.birth, death=args.death, shifts=shifts,
                              fossilization=(args.fossilization or None),
                              sampling_fraction=args.sampling_fraction, removal=args.removal,
                              mass_extinctions=mass_ext)

def _build_heterogeneous_model(args: argparse.Namespace, parser: argparse.ArgumentParser,
                               mass_ext):
    """Build a ClaDS or DiversityDependent model — both forward-only, per-lineage/diversity-
    dependent rate processes selected by ``--diversification``."""
    if args.model != "forward":
        parser.error(f"--diversification {args.diversification} is a forward-in-time process; "
                     "add --mode forward")
    if args.shifts is not None or len(args.birth) > 1 or len(args.death) > 1:
        parser.error(f"--diversification {args.diversification} takes a single --birth/--death "
                     "(no --shifts / multiple rates — those are the episodic model)")
    if args.fossilization or args.removal != 1.0:
        parser.error(f"--fossilization / --removal are not supported by --diversification "
                     f"{args.diversification}")
    if args.diversification == "clads":
        return ClaDS(args.birth[0], alpha=args.clads_alpha, sigma=args.clads_sigma,
                     turnover=args.turnover, sampling_fraction=args.sampling_fraction,
                     mass_extinctions=mass_ext)
    # diversity-dependent
    if args.carrying_capacity is None:
        parser.error("--diversification diversity-dependent needs --carrying-capacity/-K")
    return DiversityDependent(args.birth[0], args.death[0],
                              carrying_capacity=args.carrying_capacity,
                              sampling_fraction=args.sampling_fraction,
                              mass_extinctions=mass_ext)

def _build_clade_shift_model(args: argparse.Namespace, parser: argparse.ArgumentParser,
                             mass_ext):
    """Build a CladeShiftBirthDeath — constant background plus scheduled clade-specific rate
    shifts (forward-only, age mode)."""
    if args.model != "forward":
        parser.error("--clade-shift requires --mode forward (the shifts play out forward in time)")
    if args.shifts is not None or len(args.birth) > 1 or len(args.death) > 1:
        parser.error("--clade-shift takes a single background --birth/--death (no --shifts; those "
                     "are the episodic model)")
    if args.fossilization or args.removal != 1.0:
        parser.error("--fossilization / --removal are not supported with --clade-shift")
    shifts = [(a, b, d) for a, b, d in args.clade_shift]
    return CladeShiftBirthDeath(args.birth[0], args.death[0], clade_shifts=shifts,
                                sampling_fraction=args.sampling_fraction,
                                mass_extinctions=mass_ext)


def run(args, parser):
    if args.model == "forward" and args.ghosts:
        parser.error("--ghosts un-prunes a reconstructed (backward) tree; forward trees "
                     "already include extinct lineages")
    model = _build_species_model(args, parser)
    common = dict(age_type=args.age_type, max_attempts=args.max_attempts,
                  max_lineages=args.max_lineages, seed=args.seed)

    t0 = time.perf_counter()
    if args.model == "backward":
        n_tips = args.tips if args.tips is not None else 50
        age = args.age if args.age is not None else 1.0
        tree = simulate_species_tree(model, n_tips=n_tips, age=age,
                                     direction="backward", **common)
        if args.ghosts:
            add_ghost_lineages(tree, model, method=args.ghost_method, seed=args.seed)
    else:  # forward
        if (args.tips is None) == (args.age is None):
            parser.error("forward model needs exactly one of --tips or --age "
                         "(--tips to stop at that many extant species; "
                         "--age to grow for that long)")
        if args.mass_extinction and args.age is None:
            parser.error("--mass-extinction needs --age: its times are ages before a fixed "
                         "present, which --tips (random age) leaves undefined")
        if args.clade_shift and args.age is None:
            parser.error("--clade-shift needs --age: its times are ages before a fixed "
                         "present, which --tips (random age) leaves undefined")
        try:
            tree = simulate_species_tree(model, n_tips=args.tips, age=args.age,
                                         direction="forward", **common)
        except RuntimeError:
            raise RuntimeError(
                f"forward simulation kept going extinct in {args.max_attempts} attempts. "
                f"With --death {args.death} vs --birth {args.birth}, most runs die out — "
                f"lower --death, raise --max-attempts, or use --mode backward.") from None
    dt = time.perf_counter() - t0

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick() + "\n")               # the complete tree (extinct/ghost tips kept)
    leaves = tree.leaves()
    n_extant = sum(1 for n in leaves if n.is_extant)
    n_unsampled = sum(1 for n in leaves if n.name.startswith("u"))   # ghost tips (u*), from ρ<1
    n_extinct = len(leaves) - n_extant - n_unsampled
    wrote = "species_tree.nwk"
    if n_extant and n_extant < len(leaves):            # dead tips present: also the pruned tree
        with open(os.path.join(args.out, "species_tree_extant.nwk"), "w") as f:
            f.write(prune(tree, keep="extant").to_newick() + "\n")
        wrote += " + species_tree_extant.nwk"
    with open(os.path.join(args.out, "species_nodes.tsv"), "w") as f:   # node metadata table
        f.write("name\ttime\tis_leaf\tis_extant\n")
        for node in tree.nodes():
            is_leaf = not node.children
            f.write(f"{node.name}\t{node.time:.10g}\t{is_leaf}\t{bool(node.is_extant)}\n")
    wrote += " + species_nodes.tsv"
    parts = [f"{n_extant} extant"]
    if n_extinct:
        parts.append(f"{n_extinct} extinct")
    if n_unsampled:
        parts.append(f"{n_unsampled} unsampled")
    summary = " + ".join(parts) + " tips"
    print(f"wrote {args.out}/{wrote} ({summary}) in {dt:.3g} s")
    _write_params_log(os.path.join(args.out, "species_tree.log"), args, summary)
    return 0
