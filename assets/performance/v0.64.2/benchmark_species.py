"""Species-tree timing and RSS in fresh processes on an archived release."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import gc
import json
import os
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import sys
import time

p=argparse.ArgumentParser(description=__doc__)
p.add_argument("--source-root",required=True)
p.add_argument("--out")
p.add_argument("--tips",type=int)
p.add_argument("--seed",type=int,default=1001)
a=p.parse_args()
sys.path.insert(0,str(Path(a.source_root).resolve()))
import numpy  # noqa: E402 - import only after selecting the archived source
import zombi2  # noqa: E402
from zombi2.species import simulate_species_tree  # noqa: E402
assert Path(zombi2.__file__).resolve().is_relative_to(Path(a.source_root).resolve())
assert zombi2.__version__ == "0.64.2"
if a.tips:
    gc.collect()
    baseline=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2
    gc.disable()
    t=time.perf_counter()
    try:
        r=simulate_species_tree(birth=1.,death=.3,n_extant=a.tips,seed=a.seed)
        dt=time.perf_counter()-t
    finally:
        gc.enable()
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2
    print(json.dumps({"tips":a.tips,"seed":a.seed,"seconds":dt,"peak_rss_mib":rss,
                      "baseline_rss_mib":baseline,"n_extant":r.n_extant,
                      "nodes":len(r.complete_tree.nodes)}))
else:
    out={"env":{"timestamp_utc":datetime.now(timezone.utc).isoformat(),
                "python":sys.version,"numpy":numpy.__version__,"zombi2":zombi2.__version__,
                "source_commit":"047ce05fd96362439a9b9f633b13ce0cfd45f94e",
                "zombi2_source":zombi2.__file__,"platform":platform.platform(),
                "cpu_count":os.cpu_count(),"loadavg":os.getloadavg()},
         "task":"birth=1, death=.3, stop at N extant; GC disabled inside timed call",
         "rows":[],"summary":[]}
    for n in (1000,10000,100000,500000,1000000):
        rows=[]
        for rep in range(3):
            cmd=[sys.executable,str(Path(__file__).resolve()),"--source-root",a.source_root,
                 "--tips",str(n),"--seed",str(1001+rep)]
            r=subprocess.run(cmd,capture_output=True,text=True,check=True)
            row=json.loads(r.stdout)
            assert row["n_extant"]==n
            rows.append(row)
            out["rows"].append(row)
        summary={"tips":n,"median_s":statistics.median(r["seconds"] for r in rows),
                 "peak_rss_mib":max(r["peak_rss_mib"] for r in rows)}
        out["summary"].append(summary)
        Path(a.out).write_text(json.dumps(out,indent=2))
        print(json.dumps(summary),flush=True)
