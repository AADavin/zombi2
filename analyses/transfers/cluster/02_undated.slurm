#!/bin/bash
#SBATCH --job-name=z2t_undated
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-71
#SBATCH --output=logs/undated_%A_%a.out
#SBATCH --error=logs/undated_%A_%a.err

# ALEml_undated on the unlabeled extant species tree, 10 families per array task
# (718 .ale files -> tasks 0-71). Outputs land in recs/: per family a .uml_rec and
# a .uTs (the transfer list the scoring reads). separators=_ splits n343_g203846
# at the first underscore, so the species is n343.

BASE=/cluster/work/sunagawa/aadria/ZOMBI2/transfers
cd "$BASE/recs"

FILES=($(find "$BASE/ale" -maxdepth 1 -name "*.ale" | sort))
PER=10
START=$((SLURM_ARRAY_TASK_ID * PER))

for i in $(seq "$START" $((START + PER - 1))); do
    F=${FILES[$i]}
    [ -z "$F" ] && break
    if [ -f "$(basename "$F").uml_rec" ]; then
        echo "skip $(basename "$F"): already reconciled"
        continue
    fi
    echo "reconciling $(basename "$F")"
    apptainer exec -B /cluster "$BASE/ale.sif" \
        ALEml_undated "$BASE/data/species_extant_ale.nwk" "$F" separators=_
done
