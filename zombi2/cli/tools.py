"""zombi2 tools command."""
from __future__ import annotations

import argparse
import os



from zombi2.tree import read_newick

from zombi2.cli.framework import ZombiHelpFormatter, _examples

def _write_reconciliation_likelihoods(genomes, args: argparse.Namespace) -> None:
    """Score every extant family's gene tree (ALElite) and write Reconciliation_likelihoods.tsv."""
    from zombi2.tools.reconciliation import write_scores_tsv

    models = list(dict.fromkeys(args.score_model))  # de-dupe, keep order
    rows = genomes.reconciliation_likelihoods(
        args.dup, args.trans, args.loss, models=models,
        origination=args.score_origination, n_steps=args.score_nsteps,
    )
    write_scores_tsv(rows, os.path.join(args.out, "Reconciliation_likelihoods.tsv"), models=models)

def _run_tools_reconcile(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    """``zombi2 tools reconcile`` — the ALE reconciliation log-likelihood (ALElite) of one or
    more given gene trees, *evaluated* at fixed DTL rates (no rate fitting)."""
    from zombi2.tools import (GeneTree, SpeciesTree, FamilyScore, reconciliation_likelihood,
                        write_scores_tsv)

    with open(args.species_tree) as f:
        tree = read_newick(f.read())
    if len(tree.leaves()) < 2:
        parser.error(f"{args.species_tree} is not a usable species tree — fewer than 2 tips "
                     "(is it a valid Newick file?)")
    species_names = {n.name for n in tree.leaves()}
    sp = SpeciesTree.from_tree(tree)               # build the dated species index once

    with open(args.gene_tree) as f:                # one Newick per non-blank, non-comment line
        newicks = [ln.strip() for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    if not newicks:
        raise ValueError(f"no gene trees found in {args.gene_tree}")

    models = list(dict.fromkeys(args.model))       # de-dupe, keep order
    rows = []
    for i, nwk in enumerate(newicks, 1):
        gt = GeneTree.from_newick(nwk)
        unknown = gt.species_set() - species_names
        if unknown:
            raise ValueError(
                f"gene tree {i} references species absent from the species tree: "
                f"{', '.join(sorted(unknown))} — tip labels must be '<species>|<gid>' with "
                "<species> a species-tree leaf.")
        tips = sum(g.is_leaf for g in gt.nodes)
        logliks = {m: reconciliation_likelihood(
                        gene_tree=gt, species_tree=sp,
                        duplication=args.dup, transfer=args.trans, loss=args.loss,
                        model=m, origination=args.origination, n_steps=args.n_steps)
                   for m in models}
        rows.append(FamilyScore(family=str(i), extant_tips=tips, logliks=logliks))

    if args.out:
        os.makedirs(args.out, exist_ok=True)
        path = os.path.join(args.out, "Reconciliation_likelihoods.tsv")
        write_scores_tsv(rows, path, models=tuple(models))
        print(f"wrote {path} ({len(rows)} gene tree(s) x {len(models)} model(s))")
    elif len(rows) == 1 and len(models) == 1:
        print(f"{rows[0].logliks[models[0]]:.6f}")     # bare number — scripting-friendly
    else:                                              # same columns as write_scores_tsv
        print("family\textant_copies\t" + "\t".join(f"{m}_loglik" for m in models))
        for r in rows:
            print(f"{r.family}\t{r.extant_tips}\t"
                  + "\t".join(f"{r.logliks[m]:.6f}" for m in models))
    return 0

def _run_tools_simulate(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    """``zombi2 tools simulate`` — forward-simulate gene families under the ALE undated/reldated
    model (the generative twin of 'reconcile') and write per-family ground-truth reconciliations."""
    from zombi2.tools.reconciliation.undated import UndatedDTL
    from zombi2.tools.reconciliation.undated_sim import simulate_undated

    with open(args.species_tree) as f:
        tree = read_newick(f.read())
    if len(tree.leaves()) < 2:
        parser.error(f"{args.species_tree} is not a usable species tree — fewer than 2 tips "
                     "(is it a valid Newick file?)")
    transfers = "dated" if args.model == "reldated" else "global"
    model = UndatedDTL(args.dup, args.trans, args.loss)
    try:
        res = simulate_undated(tree, model, n_families=args.families,
                               origination=args.origination, transfers=transfers,
                               seed=args.seed, max_events=args.max_events)
    except (ValueError, RuntimeError) as e:
        parser.error(str(e))

    c = res.event_counts
    print(f"simulated {res.n_families} families under {args.model} "
          f"(d={args.dup}, t={args.trans}, l={args.loss}): "
          f"{res.n_surviving} survived, {res.n_extinct} went extinct")
    print(f"events: {c.get('D', 0)} D, {c.get('T', 0)} T, {c.get('L', 0)} L, {c.get('S', 0)} S")

    if args.score:
        from zombi2.tools.reconciliation.undated import undated_joint_loglik
        ll = undated_joint_loglik(res.gene_trees(), res.species_tree, model,
                                  origination=args.origination, transfers=transfers,
                                  n_extinct=res.n_extinct)
        print(f"joint undated log-likelihood of the {res.n_surviving} survivors "
              f"(+{res.n_extinct} extinct) under the generating odds: {ll:.6f}")

    if args.out:
        os.makedirs(args.out, exist_ok=True)
        n_extant = _write_undated_sim(res, args.out)
        print(f"wrote Reconciled_complete.nwk, Reconciled_extant.nwk ({n_extant} survivors), "
              f"Reconciliation_events.tsv and Gene_family_profiles.tsv into {args.out}/")
    return 0

def _write_undated_sim(res, out: str) -> int:
    """Write a simulated result's ground-truth reconciliations: bare annotated Newicks (one family
    per line — the format 'recon-accuracy' reads) plus a flat S/D/T/L event table. Returns the
    number of surviving (extant) families written."""
    complete_lines, extant_lines = [], []
    ev_lines = ["family\tevent\tspecies\trecipient\ttime\tgene"]
    for i, recon in enumerate(res.reconciliations, 1):
        if recon.complete is not None:
            complete_lines.append(recon.complete)
        if recon.extant is not None:
            extant_lines.append(recon.extant)
        for e in recon.events:
            ev_lines.append(f"{i}\t{e.event}\t{e.species}\t{e.recipient or ''}\t"
                            f"{e.time:.10g}\t{e.gene or ''}")
    prof_lines = ["family\t" + "\t".join(res.leaf_names)]
    for fam, counts in res.profile_rows():
        prof_lines.append(fam + "\t" + "\t".join(str(c) for c in counts))
    for name, lines in (("Reconciled_complete.nwk", complete_lines),
                        ("Reconciled_extant.nwk", extant_lines),
                        ("Reconciliation_events.tsv", ev_lines),
                        ("Gene_family_profiles.tsv", prof_lines)):
        with open(os.path.join(out, name), "w") as f:
            f.write("\n".join(lines) + ("\n" if lines else ""))
    return len(extant_lines)

_TREEDIST_COLS = ("tree", "n_leaves", "rf", "rf_norm", "rf_unrooted",
                  "branch_score", "quartet", "quartet_norm", "matching", "matching_norm")

def _treedist_row(label: str, c) -> str:
    """One TSV line for a :class:`~zombi2.tools.treedist.TreeComparison` (blank when a metric was
    skipped: quartet over max-leaves, or matching over max-leaves / SciPy missing)."""
    quartet = "" if c.quartet is None else str(c.quartet)
    quartet_norm = "" if c.quartet_normalized is None else f"{c.quartet_normalized:.6f}"
    matching = "" if c.matching is None else str(c.matching)
    matching_norm = "" if c.matching_normalized is None else f"{c.matching_normalized:.6f}"
    return "\t".join((
        label, str(c.n_leaves), str(c.rf), f"{c.rf_normalized:.6f}", str(c.rf_unrooted),
        f"{c.branch_score:.6f}", quartet, quartet_norm, matching, matching_norm,
    ))

def _run_tools_treedist(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    """``zombi2 tools treedist`` — RF, branch-score, and quartet distances between a reference
    tree and one or more comparison trees (e.g. a simulated truth vs. inferred trees)."""
    from zombi2.tools.treedist import compare_trees

    def _read_trees(path):
        with open(path) as f:
            return [ln.strip() for ln in f if ln.strip() and not ln.lstrip().startswith("#")]

    ref = _read_trees(args.reference)
    if len(ref) != 1:
        parser.error(f"{args.reference} must contain exactly one reference tree (found {len(ref)})")
    reference = read_newick(ref[0])

    estimates = _read_trees(args.estimate)
    if not estimates:
        raise ValueError(f"no trees found in {args.estimate}")

    rows = []
    for i, nwk in enumerate(estimates, 1):
        try:
            c = compare_trees(reference, read_newick(nwk), quartet=not args.no_quartet,
                              max_leaves=args.max_leaves, branch_score_order=args.branch_order)
        except ValueError as e:
            parser.error(f"tree {i} in {args.estimate}: {e}")
        rows.append(_treedist_row(str(i) if len(estimates) > 1 else "1", c))

    header = "\t".join(_TREEDIST_COLS)
    if args.out:
        os.makedirs(args.out, exist_ok=True)
        path = os.path.join(args.out, "Tree_distances.tsv")
        with open(path, "w") as f:
            f.write(header + "\n" + "\n".join(rows) + "\n")
        print(f"wrote {path} ({len(rows)} comparison(s))")
    else:
        print(header)
        for r in rows:
            print(r)
    return 0

_RECONACC_COLS = ("family", "n_nodes", "event_acc", "mapping_acc", "joint_acc",
                  "transfers", "transfers_recovered")

def _run_tools_recon_accuracy(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    """``zombi2 tools recon-accuracy`` — node-by-node accuracy of an inferred reconciliation
    against a true one, per family (paired by line) and pooled over all families."""
    from zombi2.tools.recon_accuracy import reconciliation_accuracy

    def _read(path):
        with open(path) as f:
            return [ln.strip() for ln in f if ln.strip() and not ln.lstrip().startswith("#")]

    truth, inferred = _read(args.truth), _read(args.inferred)
    if len(truth) != len(inferred):
        parser.error(f"--truth has {len(truth)} tree(s) but --inferred has {len(inferred)}; "
                     "they are paired line by line and must match")
    if not truth:
        raise ValueError(f"no reconciled trees found in {args.truth}")

    accs = []
    for i, (t, e) in enumerate(zip(truth, inferred), 1):
        try:
            accs.append(reconciliation_accuracy(t, e))
        except ValueError as err:
            parser.error(f"family {i}: {err}")

    rows = []
    for i, a in enumerate(accs, 1):
        rows.append("\t".join((
            str(i), str(a.n_nodes), f"{a.event_accuracy:.6f}", f"{a.mapping_accuracy:.6f}",
            f"{a.joint_accuracy:.6f}", str(a.transfer.n_true), str(a.transfer.both_correct),
        )))

    # pooled (micro-averaged over all nodes) summary
    N = sum(a.n_nodes for a in accs)
    ev = sum(round(a.event_accuracy * a.n_nodes) for a in accs)
    mp = sum(round(a.mapping_accuracy * a.n_nodes) for a in accs)
    jt = sum(round(a.joint_accuracy * a.n_nodes) for a in accs)
    nT = sum(a.transfer.n_true for a in accs)
    det = sum(a.transfer.detected for a in accs)
    both = sum(a.transfer.both_correct for a in accs)
    pooled = (
        f"# pooled over {len(accs)} family(ies), {N} node(s): "
        f"event_acc={ev / N:.4f} mapping_acc={mp / N:.4f} joint_acc={jt / N:.4f}"
        if N else "# pooled: no internal nodes to score"
    )
    if nT:
        pooled += f" | transfers: {det}/{nT} detected, {both}/{nT} donor+recipient recovered"

    header = "\t".join(_RECONACC_COLS)
    if args.out:
        os.makedirs(args.out, exist_ok=True)
        path = os.path.join(args.out, "Reconciliation_accuracy.tsv")
        with open(path, "w") as f:
            f.write(header + "\n" + "\n".join(rows) + "\n" + pooled + "\n")
        print(f"wrote {path} ({len(rows)} family(ies))")
        print(pooled)
    else:
        print(header)
        for r in rows:
            print(r)
        print(pooled)
    return 0

def _run_tools_parse(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    """``zombi2 tools parse`` — read an external reconciliation run (ALE or AleRax) and print a
    summary (rates, log-likelihood, top transfers); with -o, also write the tables as TSV."""
    try:                                       # optional 'reconparser' extra (pandas)
        from zombi2.tools.reconparser import ALEParser, AleRaxRun
    except ImportError as e:
        raise RuntimeError(str(e)) from e

    tool = args.tool
    if tool == "auto":
        tool = "alerax" if os.path.isdir(args.path) else "ale"

    if args.out:
        os.makedirs(args.out, exist_ok=True)

    if tool == "ale":
        p = ALEParser(args.path)
        present = [k for k, ok in p.files_exist().items() if ok]
        if not present:
            raise FileNotFoundError(
                f"no ALE output files found for base path {args.path!r} "
                "(expected .ucons_tree / .uTs / .uml_rec) — is --tool right?")
        print(f"ALE reconciliation: {p.base_path}")
        print(f"  files present: {', '.join(present)}")
        try:
            print(f"  log-likelihood: {p.get_log_likelihood():.6f}")
        except (FileNotFoundError, ValueError):
            pass
        try:
            r = p.get_ml_rates()
            print(f"  ML rates:  D={r['duplications']:.4g}  "
                  f"T={r['transfers']:.4g}  L={r['losses']:.4g}")
        except (FileNotFoundError, ValueError):
            pass
        try:
            s = p.get_summary_statistics()
            print(f"  total events:  D={s['total_duplications']:g}  "
                  f"T={s['total_transfers']:g}  L={s['total_losses']:g}  "
                  f"S={s['total_speciations']:g}")
        except (FileNotFoundError, ValueError):
            pass
        try:
            tr = p.get_transfers()
            print(f"  transfers: {len(tr)} edge(s)"
                  + (f" (top {args.top} by frequency)" if len(tr) else ""))
            for _, row in tr.nlargest(args.top, "freq").iterrows():
                print(f"     {row['from']} -> {row['to']}   {row['freq']:.3f}")
            if args.out:
                tp = os.path.join(args.out, "ale_transfers.tsv")
                tr.to_csv(tp, sep="\t", index=False)
                print(f"  wrote {tp}")
        except (FileNotFoundError, ValueError):
            pass
        if args.out:
            try:
                bs = p.get_branch_statistics()
                bp = os.path.join(args.out, "ale_branch_statistics.tsv")
                bs.to_csv(bp, sep="\t", index=False)
                print(f"  wrote {bp}")
            except (FileNotFoundError, ValueError):
                pass
        return 0

    # tool == "alerax"
    run = AleRaxRun(args.path)                  # raises NotADirectoryError on a non-dir path
    print(f"AleRax run: {run.output_dir}")
    try:
        info = run.get_run_info()
        bits = [f"version: {info.get('version', '?')}"]
        if "num_families" in info:
            bits.append(f"families: {info['num_families']}")
        if "num_species" in info:
            bits.append(f"species: {info['num_species']}")
        print("  " + "   ".join(bits))
    except (FileNotFoundError, ValueError):
        pass
    try:
        print(f"  total log-likelihood: {run.get_total_log_likelihood():.6f}")
    except (FileNotFoundError, ValueError):
        pass
    try:
        tr = run.get_transfers()
        score = "score" if "score" in tr.columns else tr.columns[-1]
        print(f"  global transfers: {len(tr)} edge(s)"
              + (f" (top {args.top} by {score})" if len(tr) else ""))
        for _, row in tr.nlargest(args.top, score).iterrows():
            print(f"     {row['from']} -> {row['to']}   {row[score]:.3f}")
        if args.out:
            tp = os.path.join(args.out, "alerax_transfers.tsv")
            tr.to_csv(tp, sep="\t", index=False)
            print(f"  wrote {tp}")
    except (FileNotFoundError, ValueError):
        pass
    if args.out:
        try:
            lk = run.get_per_family_likelihoods()
            lp = os.path.join(args.out, "alerax_per_family_likelihoods.tsv")
            lk.to_csv(lp, sep="\t", index=False)
            print(f"  wrote {lp}")
        except (FileNotFoundError, ValueError):
            pass
    return 0

def _run_tools_red(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    """``zombi2 tools red`` — Relative Evolutionary Divergence of every node of a tree."""
    from zombi2.tools import relative_evolutionary_divergence

    with open(args.tree) as f:
        tree = read_newick(f.read())
    if len(tree.leaves()) < 2:
        parser.error(f"{args.tree} is not a usable tree — fewer than 2 tips "
                     "(is it a valid Newick file?)")
    red = relative_evolutionary_divergence(tree)
    rows = [(n.name, n.is_leaf(), red[n]) for n in tree.nodes_preorder()]

    if args.out:
        os.makedirs(args.out, exist_ok=True)
        path = os.path.join(args.out, "RED.tsv")
        with open(path, "w") as f:
            f.write("node\tis_leaf\tred\n")
            for name, leaf, r in rows:
                f.write(f"{name}\t{leaf}\t{r:.6f}\n")
        print(f"wrote {path} ({len(rows)} node(s))")
    else:
        print("node\tis_leaf\tred")
        for name, leaf, r in rows:
            print(f"{name}\t{leaf}\t{r:.6f}")
    return 0

def _run_tools_export(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    """``zombi2 tools export`` — gene-order study formats from a nucleotide genomes run."""
    from zombi2.tools.geneorder_export import breakpoints_tsv, gff_text, posortho_tsv

    if not os.path.isdir(args.genomes_dir):
        parser.error(f"{args.genomes_dir} is not a directory")
    builders = {"breakpoints": ("Breakpoints.tsv", breakpoints_tsv),
                "gff": ("Genes.gff", gff_text),
                "posortho": ("Positional_orthologs.tsv", posortho_tsv)}
    if args.out:
        os.makedirs(args.out, exist_ok=True)
    for fmt in args.formats:
        filename, build = builders[fmt]
        try:
            text = build(args.genomes_dir)
        except FileNotFoundError as e:
            parser.error(str(e))
        if args.out:
            path = os.path.join(args.out, filename)
            with open(path, "w") as f:
                f.write(text)
            print(f"wrote {path} ({max(0, len(text.splitlines()) - 1)} row(s))")
        else:
            print(text, end="")
    return 0

def _add_tools_args(p: argparse.ArgumentParser) -> None:
    """The ``tools`` command groups analyses that compute on ZOMBI2 outputs (the ``zombi2.tools``
    layer). Each tool is its own sub-subcommand: ``reconcile`` (ALElite likelihood),
    ``simulate`` (its generative twin — sample gene families under the undated model),
    ``treedist`` (tree distances), ``recon-accuracy`` (reconciliation accuracy), ``red`` (RED)
    and ``parse`` (read external ALE / AleRax reconciliation output)."""
    tsub = p.add_subparsers(dest="tools_command", metavar="<tool>", required=True)
    rp = tsub.add_parser(
        "reconcile",
        help="ALE reconciliation log-likelihood of a gene tree given a species tree",
        description=(
            "Compute the ALE reconciliation log-likelihood P(gene tree | species tree, DTL "
            "rates) of one or more gene trees, EVALUATED at the given --dup/--trans/--loss "
            "(ALElite). This is not inference: it scores fixed rates, it does not fit them."
        ),
        usage="zombi2 tools reconcile -g FILE -t FILE --dup D --trans T --loss L [options]",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # log-likelihood of a reconciled gene tree under the faithful dated model",
            "  zombi2 tools reconcile -g gene_trees.nwk -t species_tree.nwk --dup 0.1 --trans 0.05 --loss 0.15",
            "",
            "  # compare all three ALE models and save the table into out/",
            "  zombi2 tools reconcile -g gene_trees.nwk -t species_tree.nwk --dup 0.1 --trans 0.05 --loss 0.15 --model dated undated reldated -o out/",
        ),
    )
    _add_tools_reconcile_args(rp)

    smp = tsub.add_parser(
        "simulate",
        help="simulate gene families under the ALE undated/reldated model (generative twin of 'reconcile')",
        description=(
            "Forward-simulate gene families under the ALEml_undated / GeneRax UndatedDTL model — "
            "the exact generative twin of 'zombi2 tools reconcile' (simulate here, then score there, "
            "and the rates round-trip). --dup/--trans/--loss are per-branch ODDS (dimensionless, "
            "relative to a speciation), NOT per-unit-time rates, and the species tree needs no dates "
            "(a cladogram is fine; unit branches are assumed). Writes a ground-truth reconciliation "
            "per family (complete + extant annotated Newicks and an S/D/T/L event table) — the same "
            "format 'zombi2 tools recon-accuracy' scores. For a dated, contemporaneous-transfer "
            "forward simulation, use 'zombi2 genomes' instead."
        ),
        usage="zombi2 tools simulate -t FILE --dup D --trans T --loss L [-n N] [-o DIR]",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # sample 200 families on a cladogram under undated odds, write ground-truth reconciliations",
            "  zombi2 tools simulate -t species_tree.nwk --dup 0.2 --trans 0.1 --loss 0.3 -n 200 -o truth/",
            "",
            "  # then score an inferred reconciliation against that truth",
            "  zombi2 tools recon-accuracy -t truth/Reconciled_extant.nwk -i inferred_recon.nwk",
        ),
    )
    _add_tools_simulate_args(smp)

    tp = tsub.add_parser(
        "treedist",
        help="RF, branch-score, and quartet distances between two trees",
        description=(
            "Tree distances between a REFERENCE tree (e.g. a simulated truth) and one or more "
            "comparison trees (e.g. inferred estimates), over their shared leaf set: rooted and "
            "unrooted Robinson-Foulds, the Kuhner-Felsenstein branch score, and the quartet "
            "distance. One output row per comparison tree."
        ),
        usage="zombi2 tools treedist -r FILE -e FILE [options]",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # distances between a true species tree and an inferred one",
            "  zombi2 tools treedist -r true_tree.nwk -e inferred_tree.nwk",
            "",
            "  # score many bootstrap/replicate trees against one reference, saved to out/",
            "  zombi2 tools treedist -r true_tree.nwk -e replicates.nwk -o out/",
        ),
    )
    _add_tools_treedist_args(tp)

    ap = tsub.add_parser(
        "recon-accuracy",
        help="accuracy of an inferred reconciliation against a known (simulated) one",
        description=(
            "Node-by-node accuracy of an INFERRED reconciliation against the TRUE one for the "
            "same gene tree: event-type accuracy and per-class precision/recall, species "
            "(MRCA) mapping accuracy, and transfer donor/recipient recovery. Inputs are ZOMBI2 "
            "annotated reconciled Newicks (as written by 'zombi2 tools simulate'), "
            "one family per line, --truth and --inferred paired by line."
        ),
        usage="zombi2 tools recon-accuracy -t FILE -i FILE [-o DIR]",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # score inferred reconciliations against the simulated truth",
            "  zombi2 tools recon-accuracy -t true_recon.nwk -i inferred_recon.nwk",
        ),
    )
    _add_tools_recon_accuracy_args(ap)

    rp = tsub.add_parser(
        "red",
        help="Relative Evolutionary Divergence (RED) of every node of a tree",
        description=(
            "Compute the Relative Evolutionary Divergence (RED, Parks et al. 2018) of every node "
            "of a rooted tree: the root is 0, every leaf is 1, and each internal node sits at its "
            "relative position along the root-to-tip path. RED is invariant to a global rate "
            "rescaling, so on a phylogram it approximates each node's relative age without a "
            "clock — GTDB's rank-normalisation quantity."
        ),
        usage="zombi2 tools red -t FILE [-o DIR]",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # RED of every node (a phylogram recovers relative ages; a dated tree gives them exactly)",
            "  zombi2 tools red -t species_tree.nwk -o out/",
        ),
    )
    _add_tools_red_args(rp)

    pp = tsub.add_parser(
        "parse",
        help="parse external reconciliation output (ALE, AleRax) and summarize it",
        description=(
            "Read the output of an established reconciliation program and print a summary — the "
            "ML DTL rates, the log-likelihood, and the top transfers. Understands classic ALE "
            "(.ucons_tree / .uTs / .uml_rec, v0.4 and v1.0) and AleRax run directories (v1.2+); "
            "the tool is auto-detected from the path (a directory is an AleRax run). With -o it "
            "also writes the transfer / per-branch tables as TSV. This is the reconparser interop "
            "bridge — needs the optional extra:  pip install 'zombi2[reconparser]'."
        ),
        usage="zombi2 tools parse PATH [--tool auto|ale|alerax] [--top N] [-o DIR]",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # summarize a classic ALE result (base path, without the .uml_rec extension)",
            "  zombi2 tools parse results.ale",
            "",
            "  # summarize an AleRax run directory and save its transfer tables into out/",
            "  zombi2 tools parse alerax_output/ --top 20 -o out/",
        ),
    )
    _add_tools_parse_args(pp)

    xp = tsub.add_parser(
        "export",
        help="export gene-order study formats from a nucleotide 'zombi2 genomes' run",
        description=(
            "Derive gene-order study formats from a nucleotide genomes output directory (the "
            "complement of the fork's zombiExporter). 'breakpoints' (adjacencies broken per tree "
            "edge), 'gff' (every node's genes as one GFF3) and 'posortho' (positional ortholog "
            "sets over the leaves) come from the per-node gene orders in BED/, so the run needs "
            "'bed' in --write. breakpoints / posortho are exact for content-conserving "
            "rearrangements (inversion / transposition); under duplication / loss gene content "
            "changes, so interpret those accordingly. ('dupinfo' and 'ffgc' are planned.)"
        ),
        usage="zombi2 tools export GENOMES_DIR --format {breakpoints,gff,posortho} [-o DIR]",
        formatter_class=ZombiHelpFormatter,
        epilog=_examples(
            "  # simulate with the gene-order outputs, then export the broken adjacencies + GFF",
            "  zombi2 genomes -t species_tree.nwk --genome-model nucleotide --genes genes.tsv \\",
            "      --root-length 3000 --inversion 0.01 --transposition 0.005 --write bed geneorder -o run/",
            "  zombi2 tools export run/ --format breakpoints gff posortho -o export/",
        ),
    )
    _add_tools_export_args(xp)

def _add_tools_export_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input / output")
    g.add_argument("genomes_dir", metavar="GENOMES_DIR",
                   help="a 'zombi2 genomes --genome-model nucleotide' output directory")
    g.add_argument("--format", dest="formats", nargs="+", required=True,
                   choices=("breakpoints", "gff", "posortho"), metavar="FORMAT",
                   help="which format(s) to export: breakpoints / gff / posortho (all need "
                        "--write bed)")
    g.add_argument("-o", "--out", metavar="DIR", default=None,
                   help="write the export file(s) into DIR (default: print to stdout)")

def _add_tools_reconcile_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input / output")
    g.add_argument("-g", "--gene-tree", required=True, metavar="FILE",
                   help="Newick file of one or more reconciled gene trees (one per line); tip "
                        "labels '<species>|<gid>', <species> matching a species-tree leaf")
    g.add_argument("-t", "--species-tree", required=True, metavar="FILE",
                   help="dated species-tree Newick (as written by 'zombi2 species')")
    g.add_argument("-o", "--out", metavar="DIR", default=None,
                   help="write Reconciliation_likelihoods.tsv into DIR (default: print to "
                        "stdout — a bare number for one tree and one model, else a table)")

    g = p.add_argument_group("DTL rates")
    g.add_argument("--dup", type=float, default=0.0, metavar="RATE",
                   help="duplication rate (per-unit-time for dated; per-branch odds for undated/reldated)")
    g.add_argument("--trans", type=float, default=0.0, metavar="RATE", help="transfer rate")
    g.add_argument("--loss", type=float, default=0.0, metavar="RATE", help="loss rate")

    g = p.add_argument_group("model")
    g.add_argument("--model", nargs="+", default=["dated"], metavar="MODEL",
                   choices=("dated", "undated", "reldated"),
                   help="ALE model(s) to score with (default: dated). dated = faithful "
                        "time-sliced likelihood (rates per-unit-time); undated = GeneRax "
                        "UndatedDTL (per-branch odds); reldated = time-overlap-constrained undated")
    g.add_argument("--n-steps", type=int, default=100, metavar="N",
                   help="dated model time-grid resolution (sub-steps per slice; default 100)")
    g.add_argument("--origination", choices=("root", "uniform"), default="root", metavar="WHERE",
                   help="where the family originates: 'root' (default; exact for root-seeded "
                        "families) or 'uniform' over branches")

def _add_tools_simulate_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input / output")
    g.add_argument("-t", "--species-tree", required=True, metavar="FILE",
                   help="species-tree Newick; a cladogram with no branch lengths is fine for the "
                        "undated model (unit branches are assumed) — reldated needs real dates")
    g.add_argument("-o", "--out", metavar="DIR", default=None,
                   help="write Reconciled_complete.nwk, Reconciled_extant.nwk and "
                        "Reconciliation_events.tsv into DIR (default: print a summary to stdout)")

    g = p.add_argument_group("DTL odds (per-branch, relative to a speciation — NOT per-unit-time)")
    g.add_argument("--dup", type=float, default=0.0, metavar="ODDS", help="duplication odds d")
    g.add_argument("--trans", type=float, default=0.0, metavar="ODDS", help="transfer odds t")
    g.add_argument("--loss", type=float, default=0.0, metavar="ODDS", help="loss odds l")

    g = p.add_argument_group("model")
    g.add_argument("--model", choices=("undated", "reldated"), default="undated", metavar="MODEL",
                   help="undated = a transfer may land on any branch (default); reldated = only on "
                        "a branch that overlaps the donor in time (needs a dated tree)")
    g.add_argument("--origination", choices=("root", "uniform"), default="root", metavar="WHERE",
                   help="where each family originates: 'root' (default) or 'uniform' over branches")
    g.add_argument("-n", "--families", type=int, default=100, metavar="N",
                   help="number of families to simulate (default 100)")
    g.add_argument("--seed", type=int, default=None, metavar="INT",
                   help="random seed for a reproducible draw")
    g.add_argument("--max-events", type=int, default=1_000_000, metavar="N",
                   help="per-family event cap; guards against runaway families at supercritical "
                        "odds (default 1000000)")
    g.add_argument("--score", action="store_true",
                   help="also report the joint undated log-likelihood of the simulated survivors "
                        "(with the true extinct count) under the generating odds — a round-trip check")

def _add_tools_treedist_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input / output")
    g.add_argument("-r", "--reference", required=True, metavar="FILE",
                   help="Newick file with exactly one reference tree (e.g. a simulated truth)")
    g.add_argument("-e", "--estimate", required=True, metavar="FILE",
                   help="Newick file of one or more comparison trees (one per line); each is "
                        "compared to the reference and must share its leaf label set")
    g.add_argument("-o", "--out", metavar="DIR", default=None,
                   help="write Tree_distances.tsv into DIR (default: print the table to stdout)")

    g = p.add_argument_group("metrics")
    g.add_argument("--no-quartet", action="store_true",
                   help="skip the quartet distance (it is O(n^4) in the number of leaves)")
    g.add_argument("--max-leaves", type=int, default=100, metavar="N",
                   help="quartet-distance guard: skip it above N leaves (default 100); raise to "
                        "force it on larger trees")
    g.add_argument("--branch-order", type=int, choices=(1, 2), default=2, metavar="P",
                   help="branch-score norm: 2 = L2 / Kuhner-Felsenstein (default), 1 = L1")

def _add_tools_recon_accuracy_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input / output")
    g.add_argument("-t", "--truth", required=True, metavar="FILE",
                   help="annotated reconciled Newick(s) of the TRUE reconciliation, one family "
                        "per line (labels '<species>|<EVENT>', '<donor>|T>recipient', tips "
                        "'<species>|<gid>')")
    g.add_argument("-i", "--inferred", required=True, metavar="FILE",
                   help="annotated reconciled Newick(s) of the INFERRED reconciliation, paired "
                        "with --truth line by line (same gene-tree topology and tip labels)")
    g.add_argument("-o", "--out", metavar="DIR", default=None,
                   help="write Reconciliation_accuracy.tsv into DIR (default: print to stdout)")

def _add_tools_red_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input / output")
    g.add_argument("-t", "--tree", required=True, metavar="FILE",
                   help="Newick tree (one tree). Branch lengths are read as-is: pass a phylogram "
                        "(substitutions) to recover relative ages, or a dated tree for exact "
                        "relative ages. Works with the trees 'zombi2 species'/'sequence' write.")
    g.add_argument("-o", "--out", metavar="DIR", default=None,
                   help="write RED.tsv (node, is_leaf, red) into DIR (default: print the table to stdout)")

def _add_tools_parse_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input / output")
    g.add_argument("path", metavar="PATH",
                   help="the reconciliation output: an ALE base path (e.g. results.ale, or any "
                        "of its .ucons_tree/.uTs/.uml_rec files) or an AleRax run directory")
    g.add_argument("--tool", choices=("auto", "ale", "alerax"), default="auto", metavar="NAME",
                   help="which reconciliation tool produced PATH (default: auto — a directory "
                        "is treated as an AleRax run, anything else as classic ALE)")
    g.add_argument("--top", type=int, default=10, metavar="N",
                   help="how many top transfers (by frequency/score) to print (default: 10)")
    g.add_argument("-o", "--out", metavar="DIR", default=None,
                   help="also write the transfer and per-branch/per-family tables as TSV into DIR")


def run(args, parser):
    if args.tools_command == "reconcile":
        return _run_tools_reconcile(args, parser)
    if args.tools_command == "simulate":
        return _run_tools_simulate(args, parser)
    if args.tools_command == "treedist":
        return _run_tools_treedist(args, parser)
    if args.tools_command == "recon-accuracy":
        return _run_tools_recon_accuracy(args, parser)
    if args.tools_command == "red":
        return _run_tools_red(args, parser)
    if args.tools_command == "export":
        return _run_tools_export(args, parser)
    if args.tools_command == "parse":
        return _run_tools_parse(args, parser)
    parser.error(f"unknown tool {args.tools_command!r}")   # unreachable: subparsers validate
