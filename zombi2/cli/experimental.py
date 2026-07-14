"""zombi2 experimental command."""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np


from zombi2.tree import read_newick

from zombi2.cli.framework import ZombiHelpFormatter, _examples, _write_params_log

def _add_experimental_args(p: argparse.ArgumentParser) -> None:
    """The ``experimental`` command groups unstable, not-yet-validated models (the
    ``zombi2.experimental`` layer). Each is a sub-subcommand: ``selection`` (ESM2 codon dN/dS) and
    ``ils`` (multispecies coalescent)."""
    esub = p.add_subparsers(dest="experimental_command", metavar="<model>", required=True)
    sp = esub.add_parser(
        "selection",
        help="language-model (ESM2) codon selection on a real annotated genome (emergent dN/dS)",
        description=(
            "Evolve a real annotated genome down a species tree with protein-language-model "
            "selection on its coding genes. The nucleotide genome model runs the structural "
            "simulation (inversion/duplication/loss/transfer/...); each gene evolves as coding DNA "
            "along its own gene tree under a codon mutation-selection process whose selection comes "
            "from an ESM2 critic (mutation on DNA, selection on the encoded protein -> emergent "
            "dN/dS), while intergenic DNA drifts neutrally. Genomes are reconstructed at every node.\n\n"
            "EXPERIMENTAL: APIs and outputs may change; needs the optional deps "
            "(pip install 'zombi2[selection]': torch, fair-esm, scipy)."
        ),
        usage="zombi2 experimental selection -t FILE --gff FILE --genome-fasta FILE -o DIR [options]",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # evolve a real genome with ESM2 purifying selection on its genes",
            "  zombi2 experimental selection -t species_tree.nwk --gff genome.gff --genome-fasta genome.fna \\",
            "      --beta 1.0 --dup 0.01 --loss 0.01 --seed 1 -o out/",
            "",
            "  # ...or calibrate the selection strength to a target genome-wide dN/dS, with the big ESM2",
            "  zombi2 experimental selection -t species_tree.nwk --gff genome.gff --genome-fasta genome.fna \\",
            "      --target-dnds 0.2 --esm-model esm2_t33_650M_UR50D -o out/",
        ),
    )
    _add_experimental_selection_args(sp)

    sp_ils = esub.add_parser(
        "ils",
        help="incomplete lineage sorting: gene trees under the multispecies coalescent",
        description=(
            "Simulate gene trees under the multispecies coalescent, so gene lineages need not "
            "coalesce at the nodes they pass through; deep coalescence makes gene trees disagree with "
            "the containing tree -- incomplete lineage sorting (ILS). The amount of ILS is set by "
            "--population-size N, in the tree's own time units: it grows with branch_length / N.\n\n"
            "Two modes. Plain: the coalescent inside the species tree -t (single-copy orthologs). "
            "DTL + ILS: add --events-trace from a 'zombi2 genomes' run and the coalescent runs inside "
            "each gene family's locus tree (duplications/transfers/losses), one gene tree per family; "
            "a duplication's new copy, a transferred copy and the family origination are single-copy "
            "foundings (bounded coalescent), while speciations allow deep coalescence.\n\n"
            "EXPERIMENTAL: APIs and outputs may change. Pure numpy (no optional dependencies)."
        ),
        usage="zombi2 experimental ils -t FILE -N POP [--events-trace FILE] [-n R] [-k C] -o DIR",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # plain ILS: 1000 gene trees under the MSC on a species tree, one allele per species",
            "  zombi2 experimental ils -t species_tree.nwk -N 0.5 -n 1000 --seed 1 -o out/",
            "",
            "  # DTL + ILS: a coalescent gene tree per family from a 'genomes' run (write it with --write trace)",
            "  zombi2 experimental ils -t species_tree.nwk --events-trace run/Events_trace.tsv -N 0.5 -o out/",
        ),
    )
    _add_experimental_ils_args(sp_ils)

def _add_experimental_selection_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input / output")
    g.add_argument("-t", "--tree", required=True, metavar="FILE",
                   help="dated species-tree Newick (as written by 'zombi2 species')")
    g.add_argument("--gff", required=True, metavar="FILE",
                   help="GFF3 with CDS features (keeps strand + reading frame); defines the coding "
                        "genes. Single-exon CDS only")
    g.add_argument("--gff-seqid", default=None, metavar="ID",
                   help="which sequence/contig of the GFF (and FASTA) to use (default: the GFF's "
                        "sole sequence; required if it has several)")
    g.add_argument("--genome-fasta", required=True, metavar="FILE",
                   help="the root genome FASTA (optionally .gz); its length sets the chromosome "
                        "length and it seeds the root that then evolves")
    g.add_argument("-o", "--out", required=True, metavar="DIR", help="output directory")
    g.add_argument("--seed", type=int, default=None, metavar="N", help="RNG seed for reproducibility")

    g = p.add_argument_group("selection (the language-model critic)")
    g.add_argument("--critic", choices=["esm2"], default="esm2", metavar="NAME",
                   help="the protein-language-model critic (default: esm2; the Critic API is "
                        "pluggable from Python for other models)")
    g.add_argument("--esm-model", default="esm2_t6_8M_UR50D", metavar="NAME",
                   help="ESM2 model: small default esm2_t6_8M_UR50D (8M params); go big with e.g. "
                        "esm2_t33_650M_UR50D (650M, GPU recommended)")
    sel = g.add_mutually_exclusive_group()
    sel.add_argument("--beta", type=float, default=None, metavar="B",
                     help="selection strength (>= 0; 0 = neutral). Default 1.0 unless --target-dnds "
                          "is given. Larger = stronger purifying selection (lower dN/dS)")
    sel.add_argument("--target-dnds", type=float, default=None, metavar="W",
                     help="instead of --beta, calibrate beta so the genome-wide expected dN/dS is "
                          "about W (in (0, 1)); measured on the root proteins")

    g = p.add_argument_group("mutation model (nucleotide; codon backbone + intergene)")
    g.add_argument("--subst-model", default="hky85", metavar="MODEL",
                   choices=["jc69", "k80", "hky85", "gtr"],
                   help="nucleotide substitution model: jc69 | k80 | hky85 | gtr (default hky85)")
    g.add_argument("--kappa", type=float, default=2.0, metavar="K",
                   help="[k80/hky85] transition/transversion ratio (default 2.0)")
    g.add_argument("--base-freqs", type=float, nargs=4, default=None, metavar=("A", "C", "G", "T"),
                   help="[hky85/gtr] equilibrium base frequencies (default equal)")
    g.add_argument("--gtr-rates", type=float, nargs=6, default=None,
                   metavar=("AC", "AG", "AT", "CG", "CT", "GT"),
                   help="[gtr] the 6 exchangeabilities (default all 1)")
    g.add_argument("--subst-rate", type=float, default=1.0, metavar="R",
                   help="overall divergence scale: neutral substitutions/site at the root (default 1.0)")
    g.add_argument("--gamma-shape", type=float, default=None, metavar="ALPHA",
                   help="discrete-Gamma across-site rate heterogeneity for INTERGENE blocks "
                        "(coding-block heterogeneity is emergent from selection; default: none)")

    g = p.add_argument_group("genome structural events (per-nucleotide rates)")
    g.add_argument("--inversion", type=float, default=0.001, metavar="R",
                   help="inversion rate (default 0.001)")
    g.add_argument("--duplication", "--dup", type=float, default=0.0, dest="duplication",
                   metavar="R", help="segmental duplication rate (default 0)")
    g.add_argument("--loss", type=float, default=0.0, metavar="R",
                   help="loss / deletion rate (default 0)")
    g.add_argument("--transfer", "--trans", type=float, default=0.0, dest="transfer",
                   metavar="R", help="transfer rate (default 0)")
    g.add_argument("--transposition", type=float, default=0.0, metavar="R",
                   help="transposition rate (default 0)")
    g.add_argument("--origination", "--orig", type=float, default=0.0, dest="origination",
                   metavar="R", help="per-branch novel-gene origination rate (default 0)")
    g.add_argument("--pseudogenization", type=float, default=0.0, metavar="P",
                   help="probability a loss demotes a gene to intergene rather than deleting it "
                        "(default 0). Note: pseudogenized lineages currently stay under selection")

def _add_experimental_ils_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input / output")
    g.add_argument("-t", "--tree", required=True, metavar="FILE",
                   help="species-tree Newick (as written by 'zombi2 species') -- the coalescent "
                        "container, and (with --events-trace) the frame the locus trees live in")
    g.add_argument("--events-trace", default=None, metavar="FILE", dest="events_trace",
                   help="a 'genomes' run's Events_trace.tsv (write it with 'zombi2 genomes ... "
                        "--write trace'). When given, run DTL + ILS: the coalescent within each gene "
                        "family's locus tree, one gene tree per family. Without it, plain species-tree ILS")
    g.add_argument("-o", "--out", required=True, metavar="DIR", help="output directory")
    g.add_argument("--seed", type=int, default=None, metavar="N", help="RNG seed for reproducibility")

    g = p.add_argument_group("coalescent")
    g.add_argument("-N", "--population-size", type=float, required=True, metavar="POP",
                   dest="population_size",
                   help="effective population size in the tree's time units: pairwise coalescence "
                        "rate 1/POP per unit time. Larger POP => more ILS (governed by branch / POP)")
    g.add_argument("-n", "--replicates", type=int, default=1, metavar="R",
                   help="independent gene trees to draw (default 1); with --events-trace, per family")
    g.add_argument("-k", "--samples", type=int, default=1, metavar="C",
                   help="alleles sampled per species tip / per extant gene copy (default 1)")

def _selection_cds_protein(genome: str, c, translate, reverse_complement):
    """The clean 5'->3' protein of one CDS, or ``None`` if it is out-of-frame / has an internal stop /
    is not ACGT (so calibration only sees translatable coding sequence)."""
    sub = genome[c.start:c.end]
    coding = sub if c.strand == 1 else reverse_complement(sub)
    if c.phase != 0 or len(coding) % 3 or any(ch not in "ACGT" for ch in coding):
        return None
    if len(coding) >= 3 and translate(coding[-3:]) == "*":
        coding = coding[:-3]
    if not coding:
        return None
    prot = translate(coding)
    return None if "*" in prot else prot

def _gff_contigs(path: str) -> set:
    """The set of sequence ids carrying a CDS feature in a GFF3 (to pair the GFF with the FASTA)."""
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    contigs: set = set()
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if len(f) >= 3 and f[2] == "CDS":
                contigs.add(f[0])
    return contigs

def _calibrate_beta_genomewide(critic, proteins, target_dnds, nuc_model, *,
                               hi: float = 64.0, tol: float = 1e-3, max_iter: int = 48) -> float:
    """Find one ``beta`` whose LENGTH-WEIGHTED-MEAN expected dN/dS over ``proteins`` is ~ ``target_dnds``.

    Each protein is profiled by the critic **once** (bounded by the CDS length, so tractable), then dN/dS
    is analytic per beta — unlike a single genome-length concatenation, which would blow up the critic's
    O(L^2) attention. dN/dS is monotone decreasing in beta, so a bisection on ``[0, hi]`` converges.
    """
    from zombi2.experimental.codon_selection import CodonSelection
    from zombi2.experimental.selection import FixedProfileCritic
    if not 0.0 < target_dnds < 1.0:
        raise ValueError(f"--target-dnds must be in (0, 1), got {target_dnds}")
    # one critic call per CDS, then a reusable analytic model per CDS (beta is set live in the loop)
    models = [(CodonSelection(FixedProfileCritic(critic.profile(p)), beta=1.0, nuc_model=nuc_model), p)
              for p in proteins]
    total = float(sum(len(p) for _, p in models))

    def omega(b: float) -> float:
        acc = 0.0
        for sel, p in models:
            sel.beta = b
            acc += len(p) * sel.dnds(p)
        return acc / total

    if omega(hi) > target_dnds:
        raise ValueError(f"--target-dnds {target_dnds} needs beta > {hi}; choose a less extreme target")
    lo, high = 0.0, float(hi)
    for _ in range(max_iter):
        mid = 0.5 * (lo + high)
        w = omega(mid)
        if w > 0.0 and abs(w - target_dnds) <= tol:
            return mid
        lo, high = (mid, high) if w > target_dnds else (lo, mid)
    mid = 0.5 * (lo + high)
    if abs(omega(mid) - target_dnds) > tol:
        raise ValueError(f"--target-dnds calibration did not converge within {max_iter} iterations; "
                         "try a less extreme target")
    return mid

def _run_experimental_selection(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    print("zombi2: 'experimental selection' is unstable — APIs and outputs may change "
          "(zombi2.experimental).", file=sys.stderr)
    import importlib.util
    missing = [m for m in ("torch", "esm", "scipy") if importlib.util.find_spec(m) is None]
    if missing:                                        # probe the ACTUAL optional deps (all lazy inside)
        raise RuntimeError(f"experimental selection needs optional dependencies {missing}; install them "
                           "with:  pip install 'zombi2[selection]'  (torch, fair-esm, scipy)")
    from zombi2.experimental import read_cds_gff, simulate_nucleotide_selection
    from zombi2.experimental.codon_selection import translate
    from zombi2.experimental.selection import ESM2Critic
    from zombi2.sequences.models import GammaRates, make_model, read_fasta, reverse_complement

    with open(args.tree) as f:
        tree = read_newick(f.read())
    fa = read_fasta(args.genome_fasta)
    # pick ONE sequence and use the same id for both the GFF and the FASTA, so CDS coordinates can
    # never be silently applied to the wrong contig
    if args.gff_seqid is not None:
        seqid = args.gff_seqid
    else:
        contigs = _gff_contigs(args.gff)
        if len(contigs) != 1:
            raise ValueError(f"the GFF spans {len(contigs)} sequences {sorted(contigs)}; "
                             "pass --gff-seqid to pick one")
        seqid = next(iter(contigs))
    if seqid not in fa:
        raise ValueError(f"the GFF is annotated on {seqid!r} but --genome-fasta has no such sequence "
                         f"(have: {', '.join(fa)}); supply the matching FASTA or --gff-seqid")
    genome = fa[seqid].upper()
    cds = read_cds_gff(args.gff, seqid=seqid)
    if not cds:
        raise ValueError(f"no CDS features found for {seqid!r} in {args.gff!r}")

    model = make_model(args.subst_model, kappa=args.kappa, freqs=args.base_freqs, rates=args.gtr_rates)
    gamma = GammaRates(args.gamma_shape) if args.gamma_shape else None
    critic = ESM2Critic(args.esm_model)                # args.critic == "esm2" (the only choice)

    if args.target_dnds is not None:
        proteins = [p for p in (_selection_cds_protein(genome, c, translate, reverse_complement)
                                for c in cds) if p]
        if not proteins:
            raise ValueError("no cleanly-translatable CDS to calibrate --target-dnds on")
        beta = _calibrate_beta_genomewide(critic, proteins, args.target_dnds, model)
        print(f"zombi2: calibrated beta = {beta:.4g} for a genome-wide dN/dS ~ {args.target_dnds} "
              f"(length-weighted over {len(proteins)} CDS)", file=sys.stderr)
    else:
        beta = 1.0 if args.beta is None else args.beta

    result, report = simulate_nucleotide_selection(
        tree, genome, cds, critic=critic, beta=beta, nuc_model=model, gamma=gamma,
        subst_rate=args.subst_rate, seed=args.seed,
        inversion=args.inversion, duplication=args.duplication, loss=args.loss,
        transfer=args.transfer, transposition=args.transposition, origination=args.origination,
        pseudogenization=args.pseudogenization)

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick())
    _write_selection_outputs(args.out, result, tree, report, beta)

    n_nodes = sum(1 for _ in tree.nodes_preorder())
    summary = (f"experimental selection: {report.n_selected}/{report.n_gene_blocks} gene blocks under "
               f"selection ({report.n_neutral_fallback} fell back to neutral), "
               f"{report.n_intergene} intergene blocks; beta={beta:.4g}. "
               f"Genomes for {n_nodes} nodes -> {args.out}/")
    print(summary)
    _write_params_log(os.path.join(args.out, "selection.log"), args, summary)
    return 0

def _run_experimental_ils(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    print("zombi2: 'experimental ils' is unstable — APIs and outputs may change "
          "(zombi2.experimental).", file=sys.stderr)
    from zombi2.experimental.ils import MultispeciesCoalescent, is_concordant

    if args.replicates < 1:
        raise ValueError("--replicates must be >= 1")
    if args.samples < 1:
        raise ValueError("--samples must be >= 1")
    with open(args.tree) as f:
        tree = read_newick(f.read())
    msc = MultispeciesCoalescent(population_size=args.population_size)
    rng = np.random.default_rng(args.seed)
    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick())

    if args.events_trace:                          # DTL + ILS: the coalescent within each locus tree
        from zombi2.genomes.reconciliation import extant_species_from_records
        from zombi2.genomes.simulation import read_events_trace
        with open(args.events_trace) as f:
            families = read_events_trace(f.read(), tree)
        if not families:
            raise ValueError(f"no gene-family events in {args.events_trace!r}")
        gid2species = extant_species_from_records(families, tree)
        fam_trees = msc._family_trees(families, gid2species, tree.total_age,
                                      args.samples, args.replicates, rng)
        gdir = os.path.join(args.out, "gene_trees")
        os.makedirs(gdir, exist_ok=True)
        for family, trees in sorted(fam_trees.items()):
            with open(os.path.join(gdir, f"{family}.nwk"), "w") as fh:
                for t in trees:
                    fh.write(t.to_newick(include_internal_names=False) + "\n")
        reps = "" if args.replicates == 1 else f", x{args.replicates} reps"
        summary = (f"experimental ils (DTL + ILS): coalescent gene trees for {len(fam_trees)} of "
                   f"{len(families)} families -> {args.out}/gene_trees/ (N={args.population_size:g}{reps})")
        print(summary)
        _write_params_log(os.path.join(args.out, "ils.log"), args, summary)
        return 0

    genes = msc.sample_gene_trees(tree, args.replicates, samples=args.samples, rng=rng)
    with open(os.path.join(args.out, "gene_trees.nwk"), "w") as f:
        for g in genes:
            f.write(g.to_newick(include_internal_names=False) + "\n")
    copies = "1 copy/species" if args.samples == 1 else f"{args.samples} copies/species"
    summary = (f"experimental ils: {len(genes)} gene tree(s) under the multispecies coalescent "
               f"(N={args.population_size:g}, {copies}) -> {args.out}/gene_trees.nwk")
    if args.samples == 1:
        conc = sum(is_concordant(g, tree) for g in genes) / len(genes)
        summary += f"; {conc:.1%} match the species-tree topology"
    print(summary)
    _write_params_log(os.path.join(args.out, "ils.log"), args, summary)
    return 0

def _write_selection_outputs(out: str, result, tree, report, beta: float) -> None:
    """Write the per-node genomes + architecture, the extant gene alignments, and the selection report."""
    from zombi2.sequences.models import write_fasta
    adir = os.path.join(out, "Architecture")
    gdir = os.path.join(out, "Genomes")
    os.makedirs(adir, exist_ok=True)
    os.makedirs(gdir, exist_ok=True)
    for node in tree.nodes_preorder():
        name = node.name
        lines = ["order\tblock\tkind\tgene_id\tstrand\tlength"]
        for i, (aid, strand) in enumerate(result.node_mosaic(node)):
            a = result._block_by_id[aid]
            lines.append(f"{i}\tblock{aid}\t{a.kind}\t{a.gene_id or '-'}\t"
                         f"{'+' if strand > 0 else '-'}\t{a.length}")
        with open(os.path.join(adir, f"{name}.tsv"), "w") as f:
            f.write("\n".join(lines) + "\n")
        write_fasta(os.path.join(gdir, f"{name}.fasta.gz"), {name: result.node_sequence(node)},
                    gzip_out=True)

    aln_dir = os.path.join(out, "Gene_alignments")
    os.makedirs(aln_dir, exist_ok=True)
    for gene, aln in result.gene_alignments().items():
        write_fasta(os.path.join(aln_dir, f"{gene}.fasta"), aln)

    with open(os.path.join(out, "Selection_report.tsv"), "w") as f:
        f.write("metric\tvalue\n")
        for k in ("n_blocks", "n_gene_blocks", "n_selected", "n_neutral_fallback",
                  "n_intergene", "n_empty"):
            f.write(f"{k}\t{getattr(report, k)}\n")
        f.write(f"beta\t{beta:.6g}\n")
    if report.fallbacks:
        with open(os.path.join(out, "Selection_fallbacks.tsv"), "w") as f:
            f.write("block_id\tgene_id\treason\n")
            for block_id, gene_id, reason in report.fallbacks:
                f.write(f"{block_id}\t{gene_id or '-'}\t{reason}\n")


def run(args, parser):
    if args.experimental_command == "selection":
        return _run_experimental_selection(args, parser)
    if args.experimental_command == "ils":
        return _run_experimental_ils(args, parser)
    parser.error(f"unknown experimental model {args.experimental_command!r}")   # unreachable
