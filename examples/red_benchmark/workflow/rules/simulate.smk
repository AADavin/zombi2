"""Rule: simulate one ultrametric time tree (known node ages) per (tree config, seed)."""


rule simulate_timetree:
    output:
        tree="results/trees/{treeset}/s{tseed}/time_tree.nwk",
    params:
        model=lambda w: TREES[w.treeset]["model"],
        n_tips=lambda w: TREES[w.treeset]["n_tips"],
        birth=lambda w: TREES[w.treeset].get("birth", 1.0),
        death=lambda w: TREES[w.treeset].get("death", 0.0),
    shell:
        "{PYTHON} {SCRIPTS}/simulate_timetree.py "
        "--model {params.model} --n-tips {params.n_tips} "
        "--birth {params.birth} --death {params.death} "
        "--treeset {wildcards.treeset} --tseed {wildcards.tseed} --out {output.tree}"
