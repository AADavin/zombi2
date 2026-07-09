"""Rule: aggregate every per-run metrics.tsv into summary.tsv + the two figures."""


rule summarize:
    input:
        # keep the full file list as the dependency (so the DAG is correct), but the script
        # discovers them under --runs-dir rather than receiving them all on argv (avoids ARG_MAX
        # on large sweeps — see summarize.py).
        metrics=all_metrics,
    output:
        summary="results/summary.tsv",
        curve="results/figures/red_curve.png",
        scatter="results/figures/red_scatter.png",
    params:
        showcase=lambda w: json.dumps(CFG.get("showcase", {})),
    shell:
        "{PYTHON} {SCRIPTS}/summarize.py "
        "--out-summary {output.summary} --fig-dir results/figures "
        "--showcase {params.showcase:q} --runs-dir results/runs"
