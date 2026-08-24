# ALE on Euler for the transfers example

Cluster home: `/cluster/work/sunagawa/aadria/ZOMBI2/transfers/`

```
transfers/
├── ale.sif                  ALE container (copied from TAYLOR_WDJ)
├── data/                    uploaded by rsync from analyses/transfers/data/
│   ├── species_extant.nwk       internal nodes labelled n<id> (ours; ALE refuses it)
│   ├── species_extant_ale.nwk   the same tree, internal labels stripped (ALE's input)
│   └── gene_trees/fam_*.nwk     718 extant gene trees, leaves n<species>_g<copy>
├── ale/                     .ale files from step 01
├── recs/                    .uml_rec and .uTs files from step 02
├── logs/
├── scripts/                 this folder, mirrored
└── smoke/                   fam_0 run end-to-end by hand; validated the pipeline
```

Order, from the `transfers/` folder on Euler:

```bash
sbatch scripts/01_observe.slurm     # ~minutes: builds 718 .ale files
sbatch scripts/02_undated.slurm     # array 0-71, 10 families per task
```

Pull the reconciliations back to the laptop:

```bash
rsync -az euler:/cluster/work/sunagawa/aadria/ZOMBI2/transfers/recs/ analyses/transfers/recs/
```

Reading the output: a `.uTs` line is `from  to  freq` where a bare number is ALE's own
internal-node id and `n<id>` is a tip. The `S:` line of the `.uml_rec` carries ALE's species
tree with those ids as internal labels, so ALE nodes map to `data/branches.tsv` branches by
clade content. Scoring happens on the laptop, against the truth tables `experiment.py` wrote.

Upload after regenerating the data locally:

```bash
rsync -az data/species_extant.nwk data/species_extant_ale.nwk data/manifest.json \
    data/gene_trees euler:/cluster/work/sunagawa/aadria/ZOMBI2/transfers/data/
rsync -az cluster/ euler:/cluster/work/sunagawa/aadria/ZOMBI2/transfers/scripts/
```
