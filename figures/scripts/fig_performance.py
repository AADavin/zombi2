"""Figure 8 from raw release-pinned benchmark records."""
from pathlib import Path
import argparse
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--version',default='0.64.2')
parser.add_argument('--data-dir',type=Path,default=ROOT/'assets/performance/v0.64.2')
parser.add_argument('--output-dir',type=Path,default=ROOT/'assets')
args=parser.parse_args();D=args.data_dir;F=args.output_dir;version=args.version
F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,axs=plt.subplots(1,3,figsize=(12,3.6),layout='constrained')
s=json.loads((D/f'species-v{version}.json').read_text());a=axs[0]
x=[r['tips'] for r in s['summary']];y=[r['median_s'] for r in s['summary']]
a.loglog(x,y,'o-',color='#187f83');a.set(xlabel='Extant species tips',ylabel='Runtime (s)',title='Species trees')
g=json.loads((D/f'genomes-v{version}.json').read_text());a=axs[1]
for mode,col,marker in [('family','#187f83','o'),('ordered','#c38220','s'),('nucleotide','#765da1','^')]:
 rr=[r for r in g['summary'] if r['resolution']==mode]
 xx=np.array([r['tips'] for r in rr]);yy=np.array([r['median_s'] for r in rr])
 a.errorbar(xx,yy,yerr=[yy-np.array([r['min_s'] for r in rr]),np.array([r['max_s'] for r in rr])-yy],fmt=marker+'-',color=col,label='gene order' if mode=='ordered' else mode,capsize=2)
a.set(xscale='log',yscale='log',xlabel='Extant species tips',ylabel='Runtime (s)',title='Genome resolutions');a.legend(frameon=False,fontsize=11)
d=json.loads((D/f'asymmetree-v{version}.json').read_text());a=axs[2];rr=d['summary'];xx=np.arange(len(rr))
bars=[]
for off,tool,col,label in [(-.18,'zombi2','#187f83',f'ZOMBI2 {version}'),(.18,'asymmetree','#8a8a8a','AsymmeTree 2.3.1')]:
 yy=[r[tool+'_complete_median_s'] for r in rr]
 bars.extend(a.bar(xx+off,yy,width=.34,color=col,label=label))
a.set(yscale='log',ylabel='Runtime (s)',title='Complete gene trees',xlabel='Species tips / initial families')
a.set_xticks(xx,[str(r['tips'])+'/'+str(r['families']) for r in rr],rotation=50,ha='right',fontsize=11)
a.set_ylim(top=max(b.get_height() for b in bars)*4)
legend=a.legend(frameon=False,fontsize=11,loc='upper left')
for label,a in zip('ABC',axs):
 a.text(-.18,1.07,label,transform=a.transAxes,fontweight='bold',fontsize=14)
fig.canvas.draw()
renderer=fig.canvas.get_renderer()
assert legend.get_window_extent(renderer).y0 > max(b.get_window_extent(renderer).y1 for b in bars), 'Legend overlaps a bar'
for suffix in ['png','pdf','svg']:
 fig.savefig(F/f'performance-overview.{suffix}',dpi=300)

svg=F/'performance-overview.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
