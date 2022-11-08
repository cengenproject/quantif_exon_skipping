#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=suppa_genEvents
#SBATCH -c 1
#SBATCH --mem=5G
#SBATCH --time=5:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


echo "Start $(date)"


module load miniconda
conda activate SUPPA2


WS="WS281"
ref_dir="/gpfs/ycga/project/ysm/hammarlund/aw853/references/$WS"
ref_gtf=$ref_dir"/c_elegans.PRJNA13758.${WS}.canonical_geneset.gtf"

out_dir="data/suppa2_data/221107_events"

mkdir -p $out_dir



suppa.py generateEvents \
    --mode DEBUG \
    --event-type SE \
    -f ioe \
    -i $ref_gtf \
    -o $out_dir/221107_SE





echo "End $(date)"

