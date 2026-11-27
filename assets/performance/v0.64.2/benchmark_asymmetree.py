"""Recheck the shared DTL task, including materialized gene trees in both timers.

Run with benchmark/.venv/bin/python and --source-root pointing at a git archive.
Tree generation, output traversal/counting, and serialization are outside both timers.
The ZOMBI2 simulation-only time is also retained to quantify the old endpoint mismatch.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime, timezone
import gc
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import random
import resource
import statistics
import subprocess
import sys
import time

def rss_mb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024**2 if sys.platform == "darwin" else 1024)

def tree_pair(tips, seed):
    random.seed(seed)
    np.random.seed(seed)
    s = te.species_tree_n_age(tips, 1.0, model="yule")
    planted = s.root
    assert len(planted.children) == 1
    age = planted.tstamp
    root = next(iter(planted.children))
    ordered = [n for n in s.preorder() if n is not planted]
    ids = {id(n): i for i, n in enumerate(ordered)}
    nodes = {}
    for n in ordered:
        kids = tuple(ids[id(c)] for c in n.children)
        i = ids[id(n)]
        nodes[i] = Node(id=i, parent=None if n is root else ids[id(n.parent)],
                        birth_time=0.0 if n is root else float(age-n.parent.tstamp),
                        end_time=float(age-n.tstamp), children=kids or None,
                        fate="speciation" if kids else "extant")
    t = Tree(nodes=nodes, root=ids[id(root)])
    at_length = sum(n.dist for n in s.preorder() if n.parent is not None)
    z_length = sum(n.end_time-n.birth_time for n in t.nodes.values())
    assert len(t.leaves()) == tips
    assert abs(at_length-z_length) < 1e-9
    serial = [(n.id,n.parent,n.birth_time,n.end_time,n.children) for n in t.nodes.values()]
    return s, t, {"tips":tips,"seed":seed,"age":1.0,"total_branch_length":z_length,
                  "sha256":hashlib.sha256(json.dumps(serial).encode()).hexdigest()}

def run(tool, s, t, families, seed):
    random.seed(seed)
    np.random.seed(seed)
    gc.collect()
    start = time.perf_counter()
    if tool == "zombi2":
        g = simulate_genomes_family(t,initial_families=families,duplication=.5,
                                   transfer=.3,loss=.6,origination=0.,seed=seed)
        simulation_s = time.perf_counter()-start
        trees = g.gene_trees
        complete_s = time.perf_counter()-start
        leaves = 0
        for gt in trees.values():
            stack = [gt.complete]
            while stack:
                n = stack.pop()
                if n.children:
                    stack.extend(n.children)
                else:
                    leaves += 1
        events = Counter(e.kind for e in g.events)
        return {"simulation_s":simulation_s,"complete_s":complete_s,
                "gene_tree_leaves":leaves,"families":len(trees),"events":dict(events)}
    sim = te.GeneTreeSimulator(s)
    trees = [sim.simulate(dupl_rate=.5,hgt_rate=.3,loss_rate=.6,
                          prohibit_extinction=False) for _ in range(families)]
    complete_s = time.perf_counter()-start
    leaves = sum(1 for gt in trees for n in gt.preorder() if not n.children)
    events = Counter(n.event for gt in trees for n in gt.preorder() if n.event)
    return {"complete_s":complete_s,"gene_tree_leaves":leaves,"families":len(trees),
            "events":dict(events)}

def env():
    return {"timestamp_utc":datetime.now(timezone.utc).isoformat(),
            "python":sys.version,"numpy":np.__version__,"zombi2":zombi2.__version__,
            "zombi2_source":str(Path(zombi2.__file__).resolve()),
            "source_ref":"v0.64.2","source_commit":"047ce05fd96362439a9b9f633b13ce0cfd45f94e",
            "asymmetree":importlib.metadata.version("asymmetree"),
            "platform":platform.platform(),"processor":platform.processor(),
            "cpu_count":os.cpu_count(),"loadavg":os.getloadavg(),
            "gc_enabled":gc.isenabled(),
            "harness_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}

def main(a):
    if a.only:
        s,t,info = tree_pair(a.tips,20261023+a.tips)
        baseline = rss_mb()
        r = run(a.only,s,t,a.families,1000)
        print(json.dumps({"tool":a.only,"tree":info,"families":a.families,
                          "baseline_rss_mb":baseline,"peak_rss_mb":rss_mb(),"result":r}))
        return
    result = {"env":env(),"task":"DTL simulation including complete gene-tree materialization",
              "rates":{"duplication":.5,"transfer":.3,"loss":.6,"origination":0.},
              "replicates":a.reps,"rows":[],"summary":[],"memory":[]}
    cases = [(200,f) for f in (50,100,200,400)]+[(n,100) for n in (100,400,800)]
    for tips,fam in cases:
        s,t,info = tree_pair(tips,20261023+tips)
        # Warm both implementations on the same smaller workload, outside the timer.
        for tool in ("zombi2","asymmetree"):
            run(tool,s,t,5,999)
        pairs=[]
        for rep in range(a.reps):
            row={"tips":tips,"families":fam,"replicate":rep,"tree":info,"seed":1000+rep,
                 "loadavg":os.getloadavg(),"order":["zombi2","asymmetree"] if rep%2==0
                 else ["asymmetree","zombi2"]}
            for tool in row["order"]:
                row[tool] = run(tool,s,t,fam,1000+rep)
            pairs.append(row)
            result["rows"].append(row)
        z=statistics.median(r["zombi2"]["complete_s"] for r in pairs)
        at=statistics.median(r["asymmetree"]["complete_s"] for r in pairs)
        old=statistics.median(r["zombi2"]["simulation_s"] for r in pairs)
        zl=statistics.median(r["zombi2"]["gene_tree_leaves"] for r in pairs)
        al=statistics.median(r["asymmetree"]["gene_tree_leaves"] for r in pairs)
        summary={"tips":tips,"families":fam,"zombi2_complete_median_s":z,
                 "asymmetree_complete_median_s":at,"speedup_complete":at/z,
                 "zombi2_simulation_only_median_s":old,"speedup_old_endpoint":at/old,
                 "zombi2_gene_tree_leaves_median":zl,"asymmetree_gene_tree_leaves_median":al,
                 "work_ratio_A_over_Z":al/zl}
        result["summary"].append(summary)
        print(json.dumps(summary),flush=True)
        Path(a.out).write_text(json.dumps(result,indent=2))
    if not a.no_memory:
        for tips in (400,800):
            for rep in range(3):
                for tool in ("zombi2","asymmetree"):
                    cmd=[sys.executable,str(Path(__file__).resolve()),"--source-root",a.source_root,
                         "--only",tool,"--tips",str(tips),"--families","200"]
                    p=subprocess.run(cmd,capture_output=True,text=True,check=True)
                    result["memory"].append({"replicate":rep,**json.loads(p.stdout)})
                    Path(a.out).write_text(json.dumps(result,indent=2))
        for tips in (400,800):
            rows=[r for r in result["memory"] if r["tree"]["tips"]==tips]
            assert len({r["tree"]["sha256"] for r in rows})==1
    print("Saved "+a.out,flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-root",required=True)
    p.add_argument("--out",default="asymmetree-v0.64.2.json")
    p.add_argument("--reps",type=int,default=7)
    p.add_argument("--only",choices=["zombi2","asymmetree"])
    p.add_argument("--tips",type=int,default=100)
    p.add_argument("--families",type=int,default=100)
    p.add_argument("--no-memory",action="store_true")
    args=p.parse_args()
    sys.path.insert(0,str(Path(args.source_root).resolve()))
    import numpy as np
    import asymmetree.treeevolve as te
    import zombi2
    from zombi2.genomes import simulate_genomes_family
    from zombi2.tree import Node,Tree
    assert Path(zombi2.__file__).resolve().is_relative_to(Path(args.source_root).resolve())
    assert zombi2.__version__=="0.64.2"
    main(args)
