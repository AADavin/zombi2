"""zombi2 genomes command."""
from __future__ import annotations

import argparse
import os
import sys
import time

from zombi2.genomes.nucleotide_sim import simulate_nucleotide_genomes
from zombi2.genomes.profiles import ProfileMatrix
from zombi2.genomes.genome import OrderedGenome
from zombi2.genomes.rates import LineageRates, FamilySampledRates, Rates
from zombi2.genomes.conversion import ConversionModel
from zombi2.genomes.read_rates import read_lineage_rates, read_family_rates
from zombi2.genomes.simulation import Genomes, simulate_genomes
from zombi2.genomes.transfers import TransferModel
from zombi2.tree import Tree, read_newick
from zombi2.cli.framework import _add_params_arg, _int_or_float, _write_params_log
from zombi2.cli.tools import _write_reconciliation_likelihoods


def _add_rate_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("general")
    _add_params_arg(g)
    g.add_argument("-t", "--tree", required=True, metavar="FILE",
                   help="input species tree in Newick format (e.g. species_tree.nwk)")
    g.add_argument("--genome-resolution", "--genome-model", dest="genome_model",
                   choices=("unordered", "ordered", "nucleotide"), default="unordered",
                   metavar="RESOLUTION",
                   help="genome resolution: unordered (default) evolves gene families with no "
                        "positional structure; ordered places genes on a chromosome where order "
                        "matters (adds inversion/transposition on gene segments; distance counted "
                        "in genes, not nucleotides); nucleotide evolves nucleotide-resolution "
                        "genomes by variable-length structural events, genes emerge as 'blocks' "
                        "(see the nucleotide sections). --genome-model is a deprecated alias")
    g.add_argument("--rate-per", "--per", choices=("copy", "lineage", "shared", "genome"), default=None,
                   dest="rate_per", metavar="UNIT",
                   help="what each rate is counted per — the opportunity that scales it "
                        "(unordered/ordered resolutions): copy = per gene copy, so total rates grow with "
                        "genome size (default; Rust for unordered); lineage = a constant rate per "
                        "lineage (the whole genome as one unit), giving linear rather than "
                        "exponential growth (Python); shared = one tree-wide clock per family "
                        "(constant TOTAL duplication/loss rate however many lineages carry it — "
                        "unordered only, Python). Per-family rates come from --family-rates; "
                        "nucleotide genomes are always per nucleotide. Rearrangements "
                        "(--inversion/--transposition) need --rate-per copy. (genome = deprecated "
                        "alias of lineage)")
    g.add_argument("--rate-model", choices=("shared", "per-genome", "family"), default=None,
                   metavar="MODEL",
                   help=argparse.SUPPRESS)  # deprecated -> --rate-per (still accepted; warns on use)
    g.add_argument("--seed", type=int, default=None, metavar="N",
                   help="RNG seed for reproducibility")
    g.add_argument("-o", "--out", required=True, metavar="DIR", help="output directory")

    g = p.add_argument_group(
        "gene-family rates",
        "per gene copy (--rate-per copy) or per lineage (--rate-per lineage); "
        "per nucleotide for --genome-resolution nucleotide")
    g.add_argument("--dup", type=float, default=0.0, metavar="RATE", help="duplication rate")
    g.add_argument("--trans", type=float, default=0.0, metavar="RATE", help="transfer (HGT) rate")
    g.add_argument("--loss", type=float, default=0.0, metavar="RATE", help="loss/deletion rate")
    g.add_argument("--orig", type=float, default=0.0, metavar="RATE",
                   help="origination rate (per lineage)")
    g.add_argument("--initial-families", type=int, default=None, metavar="N",
                   dest="initial_families",
                   help="number of gene families seeded at the root, for the unordered and ordered "
                        "genome resolutions (--genome-resolution unordered/ordered) (default: 20)")
    g.add_argument("--max-family-size", type=_int_or_float, default=None, metavar="CAP",
                   help="bound family growth: integer = absolute cap, decimal = fraction of the "
                        "number of species (e.g. 0.5) [not used by --genome-resolution nucleotide]")

    g = p.add_argument_group(
        "gene conversion",
        "intra-genome (ectopic) gene conversion; per-copy rates on unordered genomes")
    g.add_argument("--conversion", type=float, default=0.0, metavar="RATE",
                   help="per-copy intra-genome gene-conversion rate: one copy of a family overwrites "
                        "another copy of the SAME family (concerted evolution), pulling their "
                        "coalescence toward the present. Fires only on families with >= 2 copies")
    g.add_argument("--conversion-bias", type=float, default=0.0, metavar="B",
                   dest="conversion_bias",
                   help="donor directionality in [0, 1]: 0 = donor drawn uniformly (default); "
                        "1 = always the family's oldest copy (homogenise toward the founder)")

    g = p.add_argument_group(
        "custom rate tables",
        "user-supplied per-family and per-lineage rates (unordered genomes; Python engine "
        "except a receptivity-only --lineage-rates, which stays on Rust)")
    g.add_argument("--family-rates", metavar="FILE", dest="family_rates",
                   help="TSV of explicit per-family duplication/transfer/loss rates (per copy; "
                        "columns: family duplication transfer loss) — the per-family rate source; "
                        "families not listed fall back to --dup/--trans/--loss")
    g.add_argument("--lineage-rates", "--branch-rates", metavar="FILE", dest="lineage_rates",
                   help="TSV of per-lineage transfer emission (donation-rate factor) and/or "
                        "receptivity (absorption weight) (columns: lineage emission receptivity). "
                        "Emission scales that lineage's transfer rate; receptivity biases which "
                        "lineage receives. (--branch-rates is a deprecated alias)")

    g = p.add_argument_group("output")
    g.add_argument("--write", dest="write", nargs="+", metavar="PART",
                   choices=(*Genomes.WRITE_PARTS, "ancestral", "bed", "geneorder", "all"),
                   default=["profiles", "trees"],
                   help="which output files to write — any of {profiles, trace, trees, events, "
                        "transfers, summary, branch_events, reconciliations, layout, karyotype} or "
                        "'all' (default: profiles trees). 'reconciliations' writes "
                        "Reconciled_complete/extant.nwk (tips <species>|<gid>) — the truth input for "
                        "'tools recon-accuracy' and 'tools reconcile'. "
                        "species_tree.nwk is always written; 'profiles' alone takes the fast Rust "
                        "counts-only path; 'trace' (optionally with 'profiles') writes the compact "
                        "single-file event log events_trace.tsv near counts-only speed, from which "
                        "gene trees can be reconstructed later on demand; 'branch_events' writes "
                        "branch_events.tsv, the per-species-branch event counts (with an is_extant "
                        "flag). [ordered] 'layout' writes gene_order.tsv (which chromosome each gene "
                        "sits on) and 'karyotype' writes karyotype_trace.tsv (fission/fusion/"
                        "origination/loss) — both added automatically for a multi-chromosome or "
                        "fission/fusion run. [nucleotide] 'ancestral' simulates DNA and reconstructs the genome "
                        "(architecture + gzipped FASTA) at every node; 'bed' writes BED gene "
                        "annotations — genes.bed for the root genome and BED/<node>.bed per node "
                        "(needs --genes/--gff); 'geneorder' writes geneorder_events.tsv, the "
                        "structural-event log with physical breakpoints (chrom/start/length/strand) "
                        "per branch — the input for gene-order / breakpoint export)")
    g.add_argument("--sparse", action="store_true",
                   help="write the profile as a sparse long table (profiles_sparse.tsv: "
                        "family/species/copies, present cells only) instead of the dense matrix — "
                        "the scalable output for huge trees (needs 'profiles' in --write)")
    g.add_argument("--threads", type=int, default=1, metavar="N",
                   help="parallelise the counts-only profile simulation across N cores (only with "
                        "--write profiles; Poisson-thins the gene families into N independent "
                        "copies and sums them — a different but statistically identical "
                        "realisation, whose output depends on N). Default 1 (serial)")
    g.add_argument("--annotate-species", action="store_true",
                   help="label internal gene-tree nodes <gid>|<species-branch> (e.g. g570|i5)")

    g = p.add_argument_group(
        "reconciliation likelihoods (ALE)",
        "score every extant gene tree (forces the full gene-family path)")
    g.add_argument("--score-likelihoods", action="store_true",
                   help="also write reconciliation_likelihoods.tsv: the ALE reconciliation "
                        "log-likelihood of every extant family's gene tree, under each "
                        "--score-model, at the --dup/--trans/--loss rates (zombi2.tools ALElite)")
    g.add_argument("--score-model", nargs="+", metavar="MODEL",
                   choices=("dated", "undated", "reldated"), default=["dated", "undated"],
                   help="ALE model(s) to score with (default: dated undated). dated = faithful "
                        "time-sliced likelihood (rates are per-unit-time); undated/reldated use "
                        "per-branch odds")
    g.add_argument("--score-nsteps", type=int, default=100, metavar="N",
                   help="dated model time-grid resolution (sub-steps per slice; default 100)")
    g.add_argument("--score-origination", choices=("root", "uniform"), default="root",
                   metavar="WHERE",
                   help="where the family enters the tree: 'root' (default; exact for root-seeded "
                        "families) or 'uniform' over branches")

    g = p.add_argument_group("structural events (rearrangements)",
                             "--genome-resolution ordered/nucleotide")
    g.add_argument("--inversion", type=float, default=None, metavar="RATE",
                   help="inversion rate — per gene copy for --genome-resolution ordered (default 0), "
                        "per nucleotide for nucleotide (default 0.001)")
    g.add_argument("--transposition", type=float, default=None, metavar="RATE",
                   help="transposition rate — per gene copy for --genome-resolution ordered, per "
                        "nucleotide for nucleotide (default 0)")
    g.add_argument("--mean-length", type=float, default=None, metavar="L", dest="mean_length",
                   help="mean length of an inversion/transposition segment (geometric): in genes "
                        "for --genome-resolution ordered (default 1 = single-gene events), in "
                        "nucleotides for nucleotide (default 100)")
    g.add_argument("--transposition-flip", type=float, default=0.0, metavar="P",
                   dest="transposition_flip",
                   help="probability a transposed segment reinserts reverse-complemented "
                        "(gene order reversed and strands flipped), for --genome-resolution ordered "
                        "(default 0 = always keep orientation)")
    g.add_argument("--translocation", type=float, default=0.0, metavar="RATE",
                   help="[ordered/nucleotide] translocation rate: a segment (ordered, per gene) or "
                        "arc (nucleotide, per nucleotide) moves to a different chromosome of the same "
                        "genome (needs >1 chromosome; distinct from transposition, which stays on one "
                        "chromosome) (default 0 = off)")
    g.add_argument("--n-chromosomes", type=int, default=1, metavar="N", dest="n_chromosomes",
                   help="number of chromosomes seeded at the root, for --genome-resolution "
                        "ordered/nucleotide (default 1). [ordered] the root's initial families are "
                        "spread across them; [nucleotide] each is an independent full-length copy of "
                        "the root chromosome. Rearrangements stay within a chromosome (see "
                        "--fission/--fusion for chromosome-level changes)")
    g.add_argument("--linear-chromosomes", action="store_true", dest="linear_chromosomes",
                   help="ordered chromosomes are linear (segments never wrap the origin), for "
                        "--genome-resolution ordered (default: circular, as for bacteria). Nucleotide "
                        "chromosomes are always circular")
    # chromosome-tier events (ordered + nucleotide genomes; off by default). When any is set — or
    # with more than one chromosome — the run also writes the karyotype (gene_order.tsv /
    # chromosomes.tsv layout) + karyotype_trace.tsv.
    g.add_argument("--fission", type=float, default=0.0, metavar="RATE",
                   help="[ordered/nucleotide] chromosome fission rate, per chromosome: a chromosome "
                        "splits in two (linear: one breakpoint; circular: two) (default 0 = off)")
    g.add_argument("--fusion", type=float, default=0.0, metavar="RATE",
                   help="[ordered/nucleotide] chromosome fusion rate, per chromosome: two "
                        "chromosomes merge into one (default 0 = off)")
    g.add_argument("--chromosome-origination", type=float, default=0.0, metavar="RATE",
                   dest="chromosome_origination",
                   help="[ordered/nucleotide] chromosome origination rate, per genome: a de-novo "
                        "replicon (a plasmid) appears (default 0 = off)")
    g.add_argument("--chromosome-loss", type=float, default=0.0, metavar="RATE",
                   dest="chromosome_loss",
                   help="[ordered/nucleotide] chromosome loss rate, per chromosome: a whole "
                        "chromosome and its genes are lost (default 0 = off)")

    g = p.add_argument_group("nucleotide model", "with --genome-resolution nucleotide")
    g.add_argument("--initial-chromosomes", type=int, default=None, metavar="N",
                   dest="initial_chromosomes",
                   help=argparse.SUPPRESS)  # deprecated -> --n-chromosomes (still accepted; warns on use)
    g.add_argument("--root-length", type=int, default=1000, metavar="BP",
                   help="length of the root chromosome, in nucleotides (default 1000)")
    g.add_argument("--insertion", type=float, default=0.0, metavar="RATE",
                   help="per-nucleotide intergenic insertion rate: lay down a run of novel "
                        "nucleotides (a fresh block) inside an intergene (default 0)")
    g.add_argument("--deletion", type=float, default=0.0, metavar="RATE",
                   help="per-nucleotide intergenic deletion rate: remove a run from within a single "
                        "intergene, never touching a gene, never below the min-genome floor "
                        "(default 0)")
    g.add_argument("--indel-mean-length", type=float, default=10.0, metavar="L",
                   dest="indel_mean_length",
                   help="mean length (in nucleotides) of an insertion/deletion run — geometric, a "
                        "separate knob from --mean-length (default 10)")

    g = p.add_argument_group("genes & intergenes",
                             "--genome-resolution nucleotide; declare genes to enable genic mode")
    g.add_argument("--gff", metavar="FILE", default=None,
                   help="a GFF3 annotation (optionally .gz) to start from a real genome: sets the "
                        "chromosome length and the gene coordinates (intergenes are the gaps). "
                        "Overlapping genes are trimmed to be disjoint. Supersedes --genes/--root-length")
    g.add_argument("--gff-seqid", metavar="ID", default=None,
                   help="which GFF sequence to read (default: the most-annotated one — the "
                        "chromosome of a single-chromosome bacterium)")
    g.add_argument("--genes", metavar="FILE", default=None,
                   help="BED/TSV of gene intervals on the root chromosome (columns: start end "
                        "[name], 0-based half-open) — an alternative to --gff. Event breakpoints "
                        "fall only in intergene positions so genes are never split; genes and "
                        "intergenes are recovered as separate tree sets")
    g.add_argument("--pseudogenization", type=float, default=0.0, metavar="P",
                   help="[genic] probability a loss hitting a gene demotes it to intergene "
                        "(sequence retained, a state change) instead of deleting it (default 0)")
    g.add_argument("--replacement", type=float, default=0.0, metavar="P",
                   help="[genic] probability a transfer is a homologous replacement (the copy "
                        "replaces the recipient's syntenic locus, located via flanking genes; "
                        "additive when no homolog) (default 0)")

    g = p.add_argument_group("sequences & ancestral genomes",
                             "--genome-resolution nucleotide, with --write ancestral")
    g.add_argument("--subst-model", choices=("jc69", "k80", "hky85", "gtr"), default="hky85",
                   metavar="MODEL",
                   help="nucleotide substitution model for the sequences (default hky85)")
    g.add_argument("--kappa", type=float, default=2.0, metavar="K",
                   help="transition/transversion ratio for k80/hky85 (default 2.0)")
    g.add_argument("--base-freqs", type=float, nargs=4, default=None, metavar=("A", "C", "G", "T"),
                   help="equilibrium base frequencies for hky85/gtr (default equal)")
    g.add_argument("--gtr-rates", type=float, nargs=6, default=None,
                   metavar=("AC", "AG", "AT", "CG", "CT", "GT"),
                   help="the 6 GTR exchangeabilities (default all 1)")
    g.add_argument("--gamma-shape", type=float, default=None, metavar="ALPHA",
                   help="discrete-Gamma across-site rate heterogeneity shape (default: none)")
    g.add_argument("--subst-rate", type=float, default=1.0, metavar="RATE",
                   help="overall substitutions/site per unit time — scales sequence divergence "
                        "(default 1.0)")
    g.add_argument("--genome-fasta", metavar="FILE", default=None,
                   help="the input genome's DNA (FASTA, optionally .gz) to seed the root sequence; "
                        "without it the root is drawn at random")


def _write_profiles_only(out: str, tree: Tree, profiles, sparse: bool = False) -> None:
    """Emit the reduced profiles-only output: tree + copy-number/presence matrices.

    With ``sparse=True`` the profile is written as a single COO long table
    (``profiles_sparse.tsv``) that is O(present cells), so the output scales to trees
    where the dense families x species matrix would be astronomically large.
    """
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick() + "\n")
    if sparse:
        with open(os.path.join(out, "profiles_sparse.tsv"), "w") as f:
            f.write(profiles.to_coo_tsv())
        return
    with open(os.path.join(out, "profiles.tsv"), "w") as f:
        f.write(profiles.to_tsv())
    with open(os.path.join(out, "presence.tsv"), "w") as f:
        f.write(profiles.to_tsv(presence=True))


def _extension_from_mean_length(mean_length: float | None) -> float | None:
    """User-facing knob → engine parameter. The user gives the *mean* segment length L (genes or
    nucleotides); the engine wants the geometric continuation probability. ``None`` keeps the
    per-level default; otherwise ``extension = 1 - 1/L`` (L=1 → single-element events)."""
    if mean_length is None:
        return None
    if mean_length < 1.0:
        raise ValueError(f"--mean-length must be >= 1 (a segment spans at least one unit), got {mean_length}")
    return 1.0 - 1.0 / mean_length


def _run_genomes(tree: Tree, args: argparse.Namespace,
                 parser: argparse.ArgumentParser) -> str:
    """Simulate gene families along ``tree``, write output, and return a one-line summary.

    The default per-copy rates run on the Rust engine automatically (``simulate_genomes``
    raises a build hint if the extension is missing); ``--rate-per lineage`` runs on Python.
    """
    args.extension = _extension_from_mean_length(args.mean_length)   # mean-length knob → engine p
    parts = set(Genomes.WRITE_PARTS) if "all" in args.write else set(args.write)
    if args.sparse and "profiles" not in parts:
        raise ValueError("--sparse affects the profile output; add 'profiles' to --write")
    if args.threads > 1:
        # --threads parallelises ONLY the counts-only Rust fast path: built-in per-copy rates on an
        # unordered genome, profiles-only, no conversion / lineage-rates / scoring. Reject every other
        # combination up front with a flag-level message — otherwise per-lineage/family rates crash
        # deep in the engine and --conversion / ordered silently run serial (the threads ignored).
        reason = None
        if parts != {"profiles"}:
            reason = "use it with exactly --write profiles"
        elif args.genome_model != "unordered":
            reason = f"--genome-resolution {args.genome_model} runs serially; use --genome-resolution unordered"
        elif (args.rate_per in ("lineage", "genome") or args.rate_model in ("per-genome", "family")
              or args.family_rates):
            reason = ("the built-in per-copy rates are required "
                      "(--rate-per lineage / --family-rates run on Python)")
        elif args.conversion:
            reason = "--conversion runs on the full (serial) path"
        elif args.lineage_rates:
            reason = "--lineage-rates runs on the full (serial) path"
        elif getattr(args, "score_likelihoods", False):
            reason = "--score-likelihoods forces the full (serial) gene-tree path"
        if reason is not None:
            parser.error(f"--threads > 1 parallelises only the counts-only path: {reason}")

    # ---- opportunity axis (--rate-per), folding in deprecated spellings (--rate-model, "genome") ----
    rate_per = args.rate_per  # None (default per-copy) | "copy" | "lineage" | "genome" (deprecated)
    if rate_per == "genome":  # "genome" is the deprecated spelling of "lineage" (one genome / lineage)
        print("warning: --rate-per genome is deprecated; use --rate-per lineage "
              "(one genome per lineage — the same measure).", file=sys.stderr)
        rate_per = "lineage"
    if args.rate_model is not None:
        print("warning: --rate-model is deprecated; use --rate-per {copy,lineage} "
              "(and --family-rates for per-family rates).", file=sys.stderr)
        if args.rate_model == "family":
            if args.family_rates is None:
                parser.error("--rate-model family is deprecated; supply per-family rates with "
                             "--family-rates FILE")
        else:
            mapped = "copy" if args.rate_model == "shared" else "lineage"
            if rate_per is not None and rate_per != mapped:
                parser.error("pass --rate-per or the deprecated --rate-model, not both")
            rate_per = mapped
    delattr(args, "rate_model")  # log the single canonical field (rate_per), not both spellings

    if args.genome_model == "nucleotide":
        if rate_per is not None:
            parser.error("--rate-per (and the deprecated --rate-model) apply to the unordered/"
                         "ordered genome resolutions; the nucleotide model is per nucleotide by "
                         "construction")
        if args.initial_families is not None:
            parser.error("--initial-families is for the unordered genome resolution "
                         "(--genome-resolution unordered); the nucleotide model uses --n-chromosomes")
        if getattr(args, "score_likelihoods", False):
            parser.error("--score-likelihoods scores reconstructed gene-family trees, which the "
                         "nucleotide genome model does not produce; use --genome-resolution "
                         "unordered/ordered to score reconciliation likelihoods")
        return _run_nucleotides(tree, args, parts)

    if args.initial_chromosomes is not None:
        parser.error("--initial-chromosomes is only for --genome-resolution nucleotide; the "
                     "unordered and ordered genome resolutions use --initial-families")
    if args.translocation and args.genome_model != "ordered":
        parser.error("--translocation needs chromosomes: use --genome-resolution ordered or nucleotide "
                     "(it moves a segment/arc between chromosomes of a multi-chromosome genome)")

    ordered = args.genome_model == "ordered"
    initial_families = 20 if args.initial_families is None else args.initial_families
    args.initial_families = initial_families  # record the effective value in the params log
    # user-supplied custom rate tables are unordered-only (Python engine, except a
    # receptivity-only --lineage-rates on a plain PerCopyRates, which stays on Rust)
    if (args.family_rates or args.lineage_rates) and args.genome_model != "unordered":
        parser.error("--family-rates / --lineage-rates are only for --genome-resolution unordered")
    per_lineage = rate_per == "lineage"
    per_shared = rate_per == "shared"
    if per_shared and args.genome_model != "unordered":
        parser.error("--rate-per shared is currently supported for --genome-resolution unordered only")
    args.rate_per = rate_per if rate_per in ("lineage", "shared") else "copy"  # effective value in the log
    if not (0.0 <= args.transposition_flip <= 1.0):
        parser.error("--transposition-flip must be a probability in [0, 1]")
    if args.transposition_flip and not ordered:
        parser.error("--transposition-flip applies to transpositions on an ordered chromosome; "
                     "use --genome-resolution ordered")
    if args.n_chromosomes < 1:
        parser.error("--n-chromosomes must be >= 1")
    if args.n_chromosomes != 1 and not ordered:  # nucleotide is handled in its own branch above
        parser.error("--n-chromosomes applies to --genome-resolution ordered or nucleotide")
    if args.linear_chromosomes and not ordered:
        parser.error("--linear-chromosomes applies to --genome-resolution ordered")
    chrom_tier = bool(args.fission or args.fusion or args.chromosome_origination
                      or args.chromosome_loss)
    if chrom_tier and not ordered:
        parser.error("--fission / --fusion / --chromosome-origination / --chromosome-loss apply to "
                     "--genome-resolution ordered or nucleotide")
    if ordered and (args.n_chromosomes > 1 or chrom_tier):
        # auto-surface the karyotype outputs when non-trivial, so a multi-chromosome or fission /
        # fusion run captures its layout (and genealogy) without the user asking; single-chromosome
        # runs are untouched.
        parts.add("layout")
        if chrom_tier:
            parts.add("karyotype")
    family_mode = args.family_rates is not None  # per-family rates come from the table
    if family_mode and per_lineage:
        parser.error("--family-rates supplies per-family per-copy rates; it does not combine with "
                     "--rate-per lineage")
    if args.conversion and family_mode:
        parser.error("gene conversion (--conversion) needs per-copy rates; the per-family "
                     "table carries no conversion rate")

    rates = None  # None => D/T/L/O shorthand (plain per-copy, unordered — the Rust fast path)
    if per_lineage:
        if args.conversion:
            parser.error("gene conversion (--conversion) needs --rate-per copy; "
                         "per-lineage rates do not carry it")
        if ordered and (args.inversion is not None or args.transposition is not None):
            parser.error("rearrangements (--inversion/--transposition) need --rate-per copy; "
                         "per-lineage rates do not carry them")
        rates = Rates(args.dup, args.trans, args.loss, args.orig, per="lineage")
    elif per_shared:
        if args.conversion:
            parser.error("gene conversion (--conversion) needs --rate-per copy")
        if family_mode:
            parser.error("--family-rates does not combine with --rate-per shared")
        if args.trans:
            parser.error("--rate-per shared does not yet support transfer; use --dup/--loss "
                         "(and --orig, which stays per lineage)")
        if args.lineage_rates:
            parser.error("--rate-per shared does not support --lineage-rates (per-lineage transfer "
                         "emission/receptivity); the shared clock has no transfer channel")
        rates = Rates(args.dup, 0.0, args.loss, args.orig, per="shared")
    elif ordered:  # per-copy rates + rearrangements on an ordered chromosome
        if args.conversion:
            parser.error("gene conversion (--conversion) is only supported on unordered genomes "
                         "(--genome-resolution unordered)")
        inv = 0.0 if args.inversion is None else args.inversion
        tps = 0.0 if args.transposition is None else args.transposition
        args.inversion, args.transposition = inv, tps  # record effective values in the params log
        rates = Rates(args.dup, args.trans, args.loss, args.orig,
                      inversion=inv, transposition=tps, translocation=args.translocation,
                      chromosome_origination=args.chromosome_origination,
                      chromosome_loss=args.chromosome_loss,
                      fission=args.fission, fusion=args.fusion)
    elif family_mode:  # each family its own rates, from the table (unlisted -> --dup/--trans/--loss)
        rates = FamilySampledRates(duplication=args.dup, transfer=args.trans, loss=args.loss,
                                   origination=args.orig,
                                   rates=read_family_rates(args.family_rates))

    # --lineage-rates overlay: per-lineage transfer emission (LineageRates) + receptivity (TransferModel)
    transfers = None
    if args.lineage_rates is not None:
        emission, receptivity = read_lineage_rates(args.lineage_rates)
        if rates is None:  # plain per-copy base — carry the conversion rate through the overlay
            rates = Rates(args.dup, args.trans, args.loss, args.orig,
                          conversion=args.conversion)
        if emission:
            rates = LineageRates(rates, factors=emission, events=("transfer",))
        if receptivity:
            transfers = TransferModel(receptivity=receptivity)

    if rates is None:  # plain shared (unordered): D/T/L/O shorthand + optional gene conversion
        model_kw = dict(duplication=args.dup, transfer=args.trans, loss=args.loss,
                        origination=args.orig, conversion=args.conversion)
    else:
        model_kw = dict(rates=rates)
    if args.conversion:
        model_kw["conversions"] = ConversionModel(bias=args.conversion_bias)
    if transfers is not None:
        model_kw["transfers"] = transfers
    rate_kw = dict(**model_kw, initial_families=initial_families,
                   max_family_size=args.max_family_size, seed=args.seed)
    if ordered:
        ext = args.extension  # ordered event length is counted in genes; None -> single-gene events
        flip = args.transposition_flip  # P(a transposed segment reinserts reverse-complemented)
        n_chrom = args.n_chromosomes     # number of chromosomes to seed at the root
        circular = not args.linear_chromosomes  # linear = eukaryotic ends (no wrap); circular default
        rate_kw["genome_factory"] = (
            lambda ids, _e=ext, _f=flip, _n=n_chrom, _c=circular: OrderedGenome(
                ids, extension=_e, transposition_flip=_f, n_chromosomes=_n, circular=_c)
        )

    # scoring reconciliation likelihoods needs the full gene-family genealogy, so it forces the
    # full path (the fast counts-only / trace paths don't reconstruct gene trees).
    score = getattr(args, "score_likelihoods", False)

    t0 = time.perf_counter()
    if parts == {"profiles"} and not score and not ordered and not args.conversion:
        # counts-only Rust fast path: no genealogy reconstructed (parallel when --threads > 1).
        # Conversion only reshapes gene trees (copy numbers are unchanged), so it has no bearing on
        # a profiles-only run and stays on the full path below.
        profiles = simulate_genomes(tree, output="profiles", threads=args.threads, **rate_kw)
        dt = time.perf_counter() - t0
        _write_profiles_only(args.out, tree, profiles, sparse=args.sparse)
        n_families = len(profiles.families)
    elif ("trace" in parts and parts <= {"trace", "profiles"} and not score and not ordered
          and not args.conversion):
        # event-trace fast path: compact events_trace.tsv (+ profile), no per-event objects,
        # no gene-tree reconstruction — near counts-only speed, trees reconstructable later
        trace = simulate_genomes(tree, output="trace", **rate_kw)
        dt = time.perf_counter() - t0
        trace.write(args.out, include=parts, sparse=args.sparse)
        n_families = len(trace.profiles.families)
    else:
        genomes = simulate_genomes(tree, **rate_kw)
        dt = time.perf_counter() - t0
        genomes.write(args.out, include=parts, sparse=args.sparse,
                      annotate_species=args.annotate_species)
        n_families = len(genomes.profiles.families)
        if score:
            _write_reconciliation_likelihoods(genomes, args)
    suffix = " + reconciliation_likelihoods.tsv" if score else ""
    return (f"wrote [{' '.join(sorted(parts))}]{suffix} to {args.out}/ "
            f"({len(tree.leaves())} tips, {n_families} gene families) in {dt:.3g} s")


def _run_nucleotides(tree: Tree, args: argparse.Namespace, parts: set) -> str:
    """Simulate nucleotide-resolution genomes (variable-length structural events) along ``tree``.

    Genes are not atomic here — they emerge as **blocks** (maximal intervals with one shared
    history). ``profiles`` writes the emergent block-by-species profile (plus ``blocks.tsv`` and
    the per-leaf ``mosaics.tsv``); ``trees`` writes the per-block gene trees and their
    reconciliations. Only ``profiles``/``trees`` apply here (the family-model ``events`` /
    ``transfers`` / ``summary`` do not). ``profiles`` alone takes the fast Rust path.
    """
    want = parts & {"profiles", "trees", "reconciliations", "ancestral", "bed", "geneorder"}
    if not want:
        raise ValueError("the nucleotide model writes 'profiles', 'trees', 'reconciliations', "
                         "'ancestral', 'bed' and/or 'geneorder'; --write events/transfers/summary/"
                         "branch_events do not apply to it")
    ancestral = "ancestral" in want
    bed = "bed" in want
    # --n-chromosomes is the unified flag (both models); --initial-chromosomes is a deprecated alias.
    # The canonical flag wins; the deprecated one is honoured only when --n-chromosomes is at its
    # default, and a genuine conflict is an error rather than a silent override.
    initial_chromosomes = args.n_chromosomes
    if args.initial_chromosomes is not None:
        print("warning: --initial-chromosomes is deprecated; use --n-chromosomes.", file=sys.stderr)
        if args.n_chromosomes != 1 and args.n_chromosomes != args.initial_chromosomes:
            raise ValueError("pass --n-chromosomes or the deprecated --initial-chromosomes, not both")
        initial_chromosomes = args.initial_chromosomes
    if initial_chromosomes < 1:
        raise ValueError("--n-chromosomes must be >= 1")
    args.n_chromosomes = args.initial_chromosomes = initial_chromosomes  # effective value in the log
    # the structural knobs are shared with the ordered level, so their defaults are resolved here
    args.inversion = 0.001 if args.inversion is None else args.inversion
    args.transposition = 0.0 if args.transposition is None else args.transposition
    args.extension = 0.99 if args.extension is None else args.extension
    if args.gff and args.genes:
        raise ValueError("give either --gff or --genes (not both) to set the gene coordinates")
    gff_info = None
    gff_all = None
    root_chromosomes = None
    if args.gff:                              # start from a real genome: length + gene coordinates
        from zombi2.genomes.gff import read_gff, read_gff_all
        if args.gff_seqid:                    # one named sequence -> a single chromosome
            gff_info = read_gff(args.gff, seqid=args.gff_seqid)
            genes = gff_info.genes
            args.root_length = gff_info.length
        else:
            gff_all = read_gff_all(args.gff)  # every sequence becomes its own chromosome
            if len(gff_all) == 1:
                gff_info = gff_all[0]
                genes = gff_info.genes
                args.root_length = gff_info.length
            else:                             # a multi-replicon genome -> heterogeneous chromosomes
                if initial_chromosomes != 1:
                    raise ValueError("--gff with several sequences already sets the chromosomes; "
                                     "drop --n-chromosomes (or --gff-seqid to pick one sequence)")
                # each replicon keeps its own topology from the GFF `Is_circular` flag, so a real
                # genome's linear chromosome + circular plasmids seed as a mixed-topology genome;
                # --linear-chromosomes forces them all linear.
                root_chromosomes = [(g.length, g.genes, False if args.linear_chromosomes else g.circular)
                                    for g in gff_all]
                genes = None
    elif args.genes:
        genes = _read_gene_intervals(args.genes)
    else:
        genes = None
    genic = bool(genes) or root_chromosomes is not None
    if bed and not genic:
        raise ValueError("--write bed annotates genes on each genome, so it needs gene "
                         "coordinates: supply --genes or --gff")
    transfers = TransferModel(replacement=0.0) if genic else None  # homologous repl. is genome-side
    indels = bool(args.insertion or args.deletion)
    chrom_tier = bool(args.fission or args.fusion or args.chromosome_origination
                      or args.chromosome_loss)
    sim_kw = dict(inversion=args.inversion, loss=args.loss, duplication=args.dup,
                  transfer=args.trans, transposition=args.transposition,
                  translocation=args.translocation,
                  origination=args.orig, insertion=args.insertion, deletion=args.deletion,
                  indel_mean_length=args.indel_mean_length, root_length=args.root_length,
                  extension=args.extension, initial_chromosomes=initial_chromosomes, seed=args.seed,
                  circular=not args.linear_chromosomes,
                  fission=args.fission, fusion=args.fusion,
                  chromosome_origination=args.chromosome_origination,
                  chromosome_loss=args.chromosome_loss,
                  gene_intervals=genes, root_chromosomes=root_chromosomes,
                  pseudogenization=args.pseudogenization,
                  replacement=args.replacement, transfers=transfers,
                  retain_internal=ancestral or bed)  # bed annotates every node's genome

    t0 = time.perf_counter()
    # the Python engine is needed for the event log, the genic model, indels, the chromosome tier,
    # translocation, and explicit heterogeneous replicons (the Rust profiles path is single-chromosome).
    if ("trees" in want or "reconciliations" in want or "geneorder" in want or genic or ancestral
            or bed or indels or chrom_tier or args.translocation or args.linear_chromosomes
            or root_chromosomes is not None):
        result = simulate_nucleotide_genomes(tree, output="genomes", **sim_kw)
    else:                                     # profiles only -> Rust fast path (Python fallback)
        try:
            result = simulate_nucleotide_genomes(tree, output="profiles", **sim_kw)
        except (ImportError, RuntimeError):
            result = simulate_nucleotide_genomes(tree, output="genomes", **sim_kw)
    dt = time.perf_counter() - t0

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "species_tree.nwk"), "w") as f:
        f.write(tree.to_newick() + "\n")
    if genic:
        _write_genes_table(args.out, result.registry)

    if "profiles" in want:
        _write_blocks_table(args.out, result.blocks)
        block_ids, species, matrix = result.profile_matrix()
        pm = ProfileMatrix([f"block{a}" for a in block_ids], species, matrix)
        if args.sparse:
            with open(os.path.join(args.out, "profiles_sparse.tsv"), "w") as f:
                f.write(pm.to_coo_tsv())
        else:
            with open(os.path.join(args.out, "profiles.tsv"), "w") as f:
                f.write(pm.to_tsv())
            with open(os.path.join(args.out, "presence.tsv"), "w") as f:
                f.write(pm.to_tsv(presence=True))
        _write_mosaics(args.out, result)
    if "trees" in want or "reconciliations" in want:
        result.write_reconciliations(args.out)   # Reconciled_complete/extant.nwk + events.tsv
    if "trees" in want:
        _write_block_gene_trees(args.out, result, genic=genic)
        if genic:
            _write_pseudogenizations(args.out, result)
    if ancestral:
        _write_ancestral(args.out, result, tree, args, gff_info, gff_all)
    if bed:
        _write_bed(args.out, result, tree, gff_info, gff_all)
    if "geneorder" in want:
        from zombi2.genomes.simulation import geneorder_events_from_log
        with open(os.path.join(args.out, "geneorder_events.tsv"), "w") as f:
            f.write(geneorder_events_from_log(result.event_log))
    # karyotype: when the run is multi-chromosome or uses the chromosome tier, surface the layout
    # (chromosomes.tsv) and the fission/fusion/origination/loss genealogy (karyotype_trace.tsv).
    # Needs the Python engine (the Rust profiles path is single-chromosome and event-log-free).
    py_engine = any(isinstance(getattr(g, "chromosomes", None), dict)
                    for g in result.leaf_genomes.values())
    multi_chrom = initial_chromosomes > 1 or (root_chromosomes is not None
                                              and len(root_chromosomes) > 1)
    if py_engine and (chrom_tier or multi_chrom):
        _write_nucleotide_karyotype(args.out, result)

    if gff_all is not None and len(gff_all) > 1:
        print(f"  GFF: {len(gff_all)} sequences -> {len(gff_all)} chromosomes "
              f"({', '.join(f'{g.seqid}:{g.length}bp/{len(g.genes)}g' for g in gff_all)})")
    elif gff_info is not None:
        print(f"  GFF {gff_info.seqid}: {gff_info.length} bp, {gff_info.n_features} genes "
              f"-> {len(gff_info.genes)} after trimming ({gff_info.n_trimmed} trimmed, "
              f"{gff_info.n_dropped} dropped as overlapping)")
    extra = f", {len(result.gene_blocks())} genes" if genic else ""
    return (f"wrote [{' '.join(sorted(want))}] (nucleotide{'/genic' if genic else ''}) to "
            f"{args.out}/ ({len(result.leaf_genomes)} tips, {len(result.blocks)} blocks{extra}) "
            f"in {dt:.3g} s")


def _write_nucleotide_karyotype(out: str, result) -> None:
    """Write the karyotype of a multi-chromosome / chromosome-tier nucleotide run.

    ``chromosomes.tsv`` — per extant leaf, which chromosome each segment sits on and in what order,
    with the chromosome's topology (``species chromosome topology position source start end strand``);
    ``karyotype_trace.tsv`` — the fission / fusion / origination / loss genealogy
    (``time event branch parents children``), one row per chromosome-tier event (header-only if the
    karyotype never changed).
    """
    lay = ["species\tchromosome\ttopology\tposition\tsource\tstart\tend\tstrand"]
    for leaf, genome in sorted(result.leaf_genomes.items(), key=lambda kv: kv[0].name):
        chroms = getattr(genome, "chromosomes", None)
        if not isinstance(chroms, dict):
            continue
        for chrom in chroms.values():
            topology = "circular" if chrom.circular else "linear"
            for pos, s in enumerate(chrom.elements):
                strand = "+" if s.strand >= 0 else "-"
                lay.append(f"{leaf.name}\t{chrom.chrom_id}\t{topology}\t{pos}\t{s.source}\t"
                           f"{s.src_start}\t{s.src_end}\t{strand}")
    with open(os.path.join(out, "chromosomes.tsv"), "w") as f:
        f.write("\n".join(lay) + "\n")

    kar = ["time\tevent\tbranch\tparents\tchildren"]
    for r in result.event_log.chromosome_records:
        parents = ";".join(str(p) for p in r.parents)
        children = ";".join(str(c) for c in r.children)
        kar.append(f"{r.time:.10g}\t{r.event.value}\t{r.branch}\t{parents}\t{children}")
    with open(os.path.join(out, "karyotype_trace.tsv"), "w") as f:
        f.write("\n".join(kar) + "\n")


def _write_ancestral(out: str, result, tree, args, gff_info, gff_all=None) -> None:
    """Simulate sequences and write the genome (architecture + gzipped DNA) at every node.

    ``architecture/<node>.tsv`` — the ordered, oriented gene/intergene mosaic of the node's genome
    (a ``chromosome`` column keeps replicons apart); ``Genomes/<node>.fasta.gz`` — its assembled DNA,
    one FASTA record per chromosome for a multi-chromosome genome (else one record, the whole
    genome); ``gene_alignments/<gene>.fasta`` — the extant per-gene alignments. The root sequence is
    seeded from ``--genome-fasta`` when given (a multi-record FASTA for a multi-sequence GFF, matched
    to each replicon by sequence name), else drawn at random.
    """
    from zombi2.sequences.models import make_model, GammaRates, read_fasta, write_fasta
    model = make_model(args.subst_model, kappa=args.kappa,
                       freqs=args.base_freqs, rates=args.gtr_rates)
    gamma = GammaRates(args.gamma_shape) if args.gamma_shape else None
    multi = gff_all is not None and len(gff_all) > 1
    root_fasta = None
    if args.genome_fasta:
        fa = read_fasta(args.genome_fasta)
        if multi:
            # a multi-record FASTA: match each replicon's record (by GFF sequence name) to its source
            root_chroms = list(result.node_genomes[tree.root].chromosomes.values())
            root_fasta = {}
            for g, chrom in zip(gff_all, root_chroms):
                if g.seqid not in fa:
                    raise ValueError(f"--genome-fasta has no record {g.seqid!r} "
                                     f"(have: {', '.join(sorted(fa))})")
                dna = fa[g.seqid]
                if len(dna) != g.length:
                    raise ValueError(f"--genome-fasta record {g.seqid} is {len(dna)} bp but the GFF "
                                     f"says {g.length} bp")
                root_fasta[chrom.elements[0].source] = dna
        else:
            seqid = gff_info.seqid if gff_info is not None else None
            root_fasta = fa[seqid] if seqid in fa else next(iter(fa.values()))
            if len(root_fasta) != args.root_length:
                raise ValueError(f"--genome-fasta sequence is {len(root_fasta)} bp but the genome is "
                                 f"{args.root_length} bp; supply the matching chromosome FASTA")
    result.simulate_sequences(model, gamma=gamma, root_fasta=root_fasta,
                              subst_rate=args.subst_rate, seed=args.seed)

    adir = os.path.join(out, "architecture")
    gdir = os.path.join(out, "genomes")
    os.makedirs(adir, exist_ok=True)
    os.makedirs(gdir, exist_ok=True)
    for node in tree.nodes_preorder():
        name = node.name
        mosaics = result.node_mosaics(node)      # {chrom_id: [(block_id, strand), ...]}
        seqs = result.node_sequences(node)        # {chrom_id: dna}
        lines = ["chromosome\torder\tblock\tkind\tgene_id\tstrand\tlength"]
        for cid, mosaic in mosaics.items():
            for i, (aid, strand) in enumerate(mosaic):
                a = result._block_by_id[aid]
                lines.append(f"{cid}\t{i}\tblock{aid}\t{a.kind}\t{a.gene_id or '-'}\t"
                             f"{'+' if strand > 0 else '-'}\t{a.length}")
        with open(os.path.join(adir, f"{name}.tsv"), "w") as f:
            f.write("\n".join(lines) + "\n")
        # one FASTA record per chromosome (multi-chromosome), else one record for the whole genome
        records = ({name: next(iter(seqs.values()))} if len(seqs) == 1
                   else {f"{name}_chr{cid}": dna for cid, dna in seqs.items()})
        write_fasta(os.path.join(gdir, f"{name}.fasta.gz"), records, gzip_out=True)

    aln_dir = os.path.join(out, "gene_alignments")
    os.makedirs(aln_dir, exist_ok=True)
    for gene, aln in result.gene_alignments().items():
        write_fasta(os.path.join(aln_dir, f"{gene}.fasta"), aln)


def _read_gene_intervals(path: str) -> list[tuple]:
    """Read a BED/TSV of gene intervals: ``start end [name]`` per line (0-based half-open).

    Blank lines and ``#`` comments are skipped; a leading ``track``/``chrom``-style header is
    tolerated (any line whose first field is not an integer).
    """
    out: list[tuple] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            try:
                a, b = int(fields[0]), int(fields[1])
            except (ValueError, IndexError):
                continue  # header / non-numeric row
            name = fields[2] if len(fields) > 2 else None
            out.append((a, b, name))
    if not out:
        raise ValueError(f"no gene intervals found in {path!r} (expected 'start end [name]' lines)")
    return out


def _write_bed(out: str, result, tree, gff_info, gff_all=None) -> None:
    """Write BED gene annotations — one BED6 feature per gene on each genome.

    ``genes.bed`` is the root (seed) genome's annotation, using the input sequence name as the
    chromosome (the GFF/FASTA seqid when a real genome was supplied) so it loads against the
    original genome. ``BED/<node>.bed`` is the annotation of every node's genome *after*
    rearrangements — genes at their coordinates on that node's chromosome, whose chromosome name
    matches the corresponding ``Genomes/<node>.fasta.gz`` record (written by ``--write ancestral``).

    Each chromosome is annotated separately, with coordinates that restart at 0 per chromosome — so
    a multi-chromosome genome (a chromosome plus its plasmids) gets one BED contig per replicon,
    named to line up with its FASTA record: ``<seqid>`` at the root (the input names) and
    ``<node>_chr<id>`` at every node (single-chromosome runs keep the plain ``<seqid>`` / ``<node>``).

    Columns are standard BED6: ``chrom  chromStart  chromEnd  name  score  strand`` — 0-based,
    half-open, the same coordinate convention ZOMBI2 uses internally, so no conversion is needed.
    Only gene blocks are emitted (intergenes are the gaps); the score field is a constant 0.

    ``strand`` is the gene's orientation *relative to the root genome* — every gene is ``+`` at the
    root, and an inversion during evolution flips it. It is not a GFF-annotated coding strand: the
    genic model does not track that (``read_gff`` reads only coordinates), so ``genes.bed`` is
    always all ``+``.
    """
    multi = gff_all is not None and len(gff_all) > 1
    # root chromosome id -> input sequence name (for genes.bed, so it loads on the real genome)
    seqid_by_cid = {}
    if multi:
        root_chroms = list(result.node_genomes[tree.root].chromosomes.values())
        seqid_by_cid = {c.chrom_id: g.seqid for g, c in zip(gff_all, root_chroms)}

    def bed_rows(node, chrom_name) -> list[str]:
        """One BED6 row per gene, coordinates restarting at 0 on each chromosome. ``chrom_name(cid)``
        gives the contig name to match this node's FASTA record(s)."""
        rows: list[str] = []
        for cid, mosaic in result.node_mosaics(node).items():
            chrom, offset = chrom_name(cid), 0
            for block_id, strand in mosaic:
                block = result._block_by_id[block_id]
                if block.kind == "gene":
                    rows.append(f"{chrom}\t{offset}\t{offset + block.length}\t"
                                f"{block.gene_id}\t0\t{'+' if strand > 0 else '-'}")
                offset += block.length
        return rows

    def write_bed(path: str, rows: list[str]) -> None:
        with open(path, "w") as f:
            f.write("\n".join(rows) + ("\n" if rows else ""))

    # root (seed) annotation — chromosome named after the input sequence(s)
    root_seqid = gff_info.seqid if gff_info is not None else "root_chromosome"
    root_name = (lambda cid: seqid_by_cid.get(cid, f"chr{cid}")) if multi else (lambda cid: root_seqid)
    write_bed(os.path.join(out, "genes.bed"), bed_rows(tree.root, root_name))

    # every node's genome (ancestral + extant), each contig keyed to its FASTA record id
    bdir = os.path.join(out, "bed")
    os.makedirs(bdir, exist_ok=True)
    for node in tree.nodes_preorder():
        node_name = ((lambda cid, n=node.name: f"{n}_chr{cid}") if multi
                     else (lambda cid, n=node.name: n))
        write_bed(os.path.join(bdir, f"{node.name}.bed"), bed_rows(node, node_name))


def _write_genes_table(out: str, registry) -> None:
    """Write ``genes.tsv`` — the gene annotation (seed genes + any originated novel genes)."""
    lines = ["gene\tsource\tstart\tend\tlength"]
    for source in sorted(registry.genes):
        for gi in registry.genes[source]:
            lines.append(f"{gi.gene_id}\t{gi.source}\t{gi.start}\t{gi.end}\t{gi.length}")
    with open(os.path.join(out, "genes.tsv"), "w") as f:
        f.write("\n".join(lines) + "\n")


def _write_pseudogenizations(out: str, result) -> None:
    """Write ``pseudogenizations.tsv`` — every gene->intergene state flip (branch, time, lineage)."""
    lines = ["block\tgene\tspecies_branch\ttime\tgene_lineage"]
    for block_id, gene_id, species, t, gid in result.pseudogenizations():
        lines.append(f"block{block_id}\t{gene_id}\t{species}\t{t:.10g}\t{gid}")
    with open(os.path.join(out, "pseudogenizations.tsv"), "w") as f:
        f.write("\n".join(lines) + "\n")


def _write_blocks_table(out: str, blocks) -> None:
    """Write ``blocks.tsv`` — the emergent gene families (uncut ancestral intervals).

    Carries the ``kind`` (gene/intergene) and ``gene_id`` classification (``-`` for intergene).
    """
    lines = ["block\tsource\tstart\tend\tlength\tkind\tgene_id"]
    for a in blocks:
        lines.append(f"block{a.block_id}\t{a.source}\t{a.start}\t{a.end}\t{a.length}\t"
                     f"{a.kind}\t{a.gene_id or '-'}")
    with open(os.path.join(out, "blocks.tsv"), "w") as f:
        f.write("\n".join(lines) + "\n")


def _write_mosaics(out: str, result) -> None:
    """Write ``mosaics.tsv`` — each extant genome as an ordered, signed sequence of blocks."""
    lines = ["leaf\tmosaic"]
    for leaf in sorted(result.leaf_genomes, key=lambda n: n.name):
        seq = " ".join(("+" if s > 0 else "-") + f"block{aid}"
                       for aid, s in result.leaf_mosaic(leaf))
        lines.append(f"{leaf.name}\t{seq}")
    with open(os.path.join(out, "mosaics.tsv"), "w") as f:
        f.write("\n".join(lines) + "\n")


def _write_block_gene_trees(out: str, result, genic: bool = False) -> None:
    """Write per-block trees to ``block<id>_complete.nwk`` / ``_extant.nwk``.

    Plain nucleotide model: everything under ``gene_trees/``. Genic mode: gene blocks under
    ``gene_trees/`` and intergene blocks under ``intergene_trees/`` (both tree sets recovered).
    """
    def dump(tdir: str, trees: dict) -> None:
        os.makedirs(tdir, exist_ok=True)
        for block_id, (complete, extant) in trees.items():
            if complete:
                with open(os.path.join(tdir, f"block{block_id}_complete.nwk"), "w") as f:
                    f.write(complete + "\n")
            if extant:
                with open(os.path.join(tdir, f"block{block_id}_extant.nwk"), "w") as f:
                    f.write(extant + "\n")

    if genic:
        dump(os.path.join(out, "gene_trees"), result.gene_trees())
        dump(os.path.join(out, "intergene_trees"), result.intergene_trees())
    else:
        dump(os.path.join(out, "gene_trees"), result.block_gene_trees())


def run(args, parser):
    with open(args.tree) as f:
        tree = read_newick(f.read())
    if len(tree.leaves()) < 2:
        parser.error(f"{args.tree} is not a usable species tree — fewer than 2 tips "
                     "(is it a valid Newick file?)")
    summary = _run_genomes(tree, args, parser)
    print(summary)
    _write_params_log(os.path.join(args.out, "genomes.log"), args, summary)
    return 0
