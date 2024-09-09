#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=merge_index_bam
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH --time=00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


# Transfer the bam files merging the technical replicates. Index the results (bai)
# Also merge the SJ technical replicates.

echo "Prepare alignments bsn12 starting $(date)"


input_sj="/gpfs/gibbs/pi/hammarlund/CeNGEN/bulk/bulk_alignments/bsn12_junctions"
out_dir="/vast/palmer/scratch/hammarlund/aw853/2024-09-09_bsn12/"


### Prepare metadata
mkdir -p $out_dir
mkdir -p $out_dir/junctions




echo "-----   Merging SJ   -----"
module load R//4.2.0-foss-2020b

Rscript src/combine_sj_files.R \
  --sj_dir $input_sj \
  --out_dir $out_dir/junctions

echo
echo "All done $(date)"


echo "TO DO manually: if needed, remove IL1r107 (no bam, only sj); remove outliers (AVKr113, PVMr122, RICr133)"
