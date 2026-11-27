"""Genome scaling at three resolutions; simulation endpoints and raw times."""
from pathlib import Path
import argparse,json,gc,sys,time,statistics,platform,os
p=argparse.ArgumentParser()
p.add_argument("--source-root",required=True);p.add_argument("--out",required=True)
a=p.parse_args();sys.path.insert(0,a.source_root)
import numpy as np
from zombi2 import species,genomes,__version__
assert __version__=="0.64.2"
rates=dict(duplication=.2,transfer=.1,loss=.25,origination=.5)
out={"version":__version__,"commit":"047ce05fd96362439a9b9f633b13ce0cfd45f94e",
     "python":sys.version,"numpy":np.__version__,"platform":platform.platform(),
     "rates":rates,"tree_rates":{"birth":1.,"death":.3},"rows":[],"summary":[],
     "endpoint":"simulate_genomes return; gene-tree reconstruction and writing excluded",
     "nucleotide_input":"synthetic circular genome: 4411532 nt, 4000 genes of 1000 nt"}
for mode,sizes in (("family",[100,300,1000,3000,10000]),
                   ("ordered",[100,300,1000,3000]),
                   ("nucleotide",[100,300,1000,3000])):
 for n in sizes:
  tree=species.simulate_species_tree(birth=1.,death=.3,n_extant=n,seed=100+n).extant_tree
  times=[]
  for rep in range(3):
   gc.collect();gc.disable();start=time.perf_counter()
   try:
    if mode=="family": g=genomes.simulate_genomes_family(tree,initial_families=20,**rates,seed=7001+rep)
    elif mode=="ordered": g=genomes.simulate_genomes_ordered(tree,initial_families=20,inversion=.1,**rates,seed=7001+rep)
    else: g=genomes.simulate_genomes_nucleotide(tree,genes=4000,root_length=4411532,gene_length=1000,chromosomes=1,topology="circular",inversion=.1,**rates,seed=7001+rep)
    dt=time.perf_counter()-start
   finally: gc.enable()
   times.append(dt)
   out["rows"].append({"resolution":mode,"tips":n,"seed":7001+rep,"seconds":dt,"loadavg":os.getloadavg(),"gene_edges":len(g.edges) if hasattr(g,"edges") else None})
   del g
  s={"resolution":mode,"tips":n,"median_s":statistics.median(times),"min_s":min(times),"max_s":max(times)}
  out["summary"].append(s);Path(a.out).write_text(json.dumps(out,indent=2));print(s,flush=True)
