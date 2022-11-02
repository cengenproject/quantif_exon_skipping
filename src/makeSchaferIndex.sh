#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=schafer_ind
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH --time=1-00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


ref_dir="/gpfs/ycga/project/ysm/hammarlund/aw853/references/WS277"

schafer_index=$ref_dir"/schafer_index"
ref_gff=$ref_dir"/c_elegans.PRJNA13758.WS277.canonical_geneset.gtf"

module load HTSeq/0.6.1p1-foss-2016b-Python-2.7.13
python ~/.utilities/DEXSeq/dexseq_prepare_annotation.py $ref_gff WS277_flattened.gff

# select "exonic_part" entries, rename them gene_id:exonic_part_number (e.g. WBGene00000001:001)
awk '{OFS="\t"}{if ($3 == "exonic_part") print $1,$2,$3,$4,$5,$6,$7,$8,$14":"$12}' WS277_flattened.gff | sed 's=[";]==g' > WS277_exonic_parts.gff

# Filter exons smaller than 5 bp
awk '$5-$4>4{print}' WS277_exonic_parts.gff > $ref_dir"/schafer_index.gff"

rm WS277_exonic_parts.gff WS277_flattened.gff


echo "Done."
