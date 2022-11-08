#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=psi
#SBATCH -c 1
#SBATCH --mem=80G
#SBATCH --time=1-00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


# Use the SUPPA2 annotation and an approach inspired from Schafer (but very different)
# to compute PSI for exon skipping, for 1 sample.
# 
# Idea:
# Inclusion reads: make ioe into bed, use bedtools coverage vs bam to count read coverage
# Exclusion reads: takethe SJ.tab file from STAR, filter it to keep only the introns
# PSI: divide the two


echo "Starting PSI quantification $(date)"

set -ue


## Inits ----


module load BEDTools


WS="WS281"
readLength=101


SE_annot="data/suppa2_data/221108_events/${WS}_"
SE_annot="/home/aw853/project/exon_skipping/data/suppa2_data/221108_events/${WS}_"

data_dir="/home/aw853/scratch60/2022-11-03_bsn9"
out_dir="data/2022-11-08_SE_PSI"
tmp_dir="/home/aw853/scratch60/221108_SE_PSI"

mkdir -p $tmp_dir
mkdir -p $out_dir


rm $tmp_dir/*



# Single sample for testing
sample="ADLr173"
my_bam=$data_dir/bams/$sample.bam
my_SJ_tab=$data_dir/junctions/$sample.SJ.tab

s=$tmp_dir/$sample



## Inclusion reads ----

coverageBed -split -s \
            -a ${SE_annot}SE_central_exon.bed \
            -b $my_bam \
  > $s.inclusion


## Exclusion reads ----

grep -F \
     -f ${SE_annot}SE_spanning_intron.tab \
     $my_SJ_tab \
  > $s.exclusion


## PSI ----

module swap BEDTools R

Rscript src/assemble_psi.R \
        --se ${SE_annot}SE_coords.tab \
        --read_len 101 \
        --inputs $tmp_dir \
        --out $out_dir
        



echo "Ending $(date)"
