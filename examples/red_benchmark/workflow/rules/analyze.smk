"""Rule: perturb a time tree's branch rates and recover node ages with RED (one realization)."""


rule perturb_and_red:
    input:
        tree="results/trees/{treeset}/s{tseed}/time_tree.nwk",
    output:
        metrics="results/runs/{treeset}/s{tseed}/{pert}/r{cseed}/metrics.tsv",
        points="results/runs/{treeset}/s{tseed}/{pert}/r{cseed}/points.csv",
    params:
        spec=lambda w: json.dumps(PERTS[w.pert]),
        model=lambda w: TREES[w.treeset]["model"],
        n_tips=lambda w: TREES[w.treeset]["n_tips"],
    shell:
        "{PYTHON} {SCRIPTS}/perturb_and_red.py "
        "--tree {input.tree} --spec {params.spec:q} "
        "--treeset {wildcards.treeset} --model {params.model} --n-tips {params.n_tips} "
        "--pert {wildcards.pert} --tseed {wildcards.tseed} --cseed {wildcards.cseed} "
        "--out-metrics {output.metrics} --out-points {output.points}"
