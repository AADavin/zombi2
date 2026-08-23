#!/bin/bash
#SBATCH --job-name=z2t_observe
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/observe_%j.out
#SBATCH --error=logs/observe_%j.err

# ALEobserve over every gene tree: one short job, the step is near-instant per family.
# Each family is a single true tree, so the .ale holds one observed topology.

BASE=/cluster/work/sunagawa/aadria/ZOMBI2/transfers
cd "$BASE"

for f in data/gene_trees/fam_*.nwk; do
    id=$(basename "$f" .nwk)
    cp "$f" "ale/${id}.nwk"
    apptainer exec -B /cluster "$BASE/ale.sif" ALEobserve "ale/${id}.nwk" > /dev/null
    rm "ale/${id}.nwk"
done

echo "done: $(ls ale/*.ale | wc -l) ale files"
