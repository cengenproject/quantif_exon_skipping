#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=schafer_ind
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


echo "Make index for Schafer quantification method. Starting $(date)"


WS="WS281"

ref_dir="/gpfs/ycga/project/ysm/hammarlund/aw853/references/$WS"
ref_gff=$ref_dir"/c_elegans.PRJNA13758.${WS}.canonical_geneset.gtf"
tmp="/home/aw853/scratch60/"

module load HTSeq/0.13.5-foss-2020b-Python-3.8.6
python ~/.utilities/DEXSeq/dexseq_prepare_annotation.py $ref_gff $tmp.flattened.gff

# select "exonic_part" entries, rename them gene_id:exonic_part_number (e.g. WBGene00000001:001)
awk '{OFS="\t"}{if ($3 == "exonic_part") print $1,$2,$3,$4,$5,$6,$7,$8,$14":"$12}' $tmp.flattened.gff \
    | sed 's=[";]==g' \
    > $tmp.exonic_parts.gff

# Filter exons smaller than 5 bp
awk '$5-$4>4{print}' $tmp.exonic_parts.gff \
    > $ref_dir"/${WS}_schafer_index.gff"

rm $tmp.exonic_parts.gff $tmp.flattened.gff


echo "Done $(date)"
