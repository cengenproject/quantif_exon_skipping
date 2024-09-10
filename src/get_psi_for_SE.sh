#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=psi
#SBATCH -c 1
#SBATCH --mem=310G
#SBATCH --time=23:50:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


# Use the SUPPA2 annotation and an approach inspired from Schafer (but very different)
# to compute PSI for exon skipping, for all samples.
# 
# Idea:
# Inclusion reads: make ioe into bed, use bedtools coverage vs bam to count read coverage
# Exclusion reads: takethe SJ.tab file from STAR, filter it to keep only the introns
# PSI: divide the two


echo "Starting PSI quantification $(date)"

set -ue


## Inits ----


module load BEDTools


WS="WS289"
readLength=101


SE_annot="data/suppa2_data/240910_events/${WS}_"


bam_dir="/gpfs/gibbs/pi/hammarlund/CeNGEN/bulk/bulk_alignments/bsn12_bams"
sj_dir="/vast/palmer/scratch/hammarlund/aw853/2024-09-09_bsn12/junctions"
out_dir="data/2024-09-10_SE_PSI"
tmp_dir="/vast/palmer/scratch/hammarlund/aw853/240910_SE_PSI_all_samples"

mkdir -p $out_dir

mkdir -p $tmp_dir
rm -r $tmp_dir
mkdir -p $tmp_dir


mapfile -t bamList < <(ls $bam_dir/*.bam \
                          | xargs basename -a -s .bam \
                          | cut -f1 -d't'\
                          | uniq)

mapfile -t sjList < <(ls $sj_dir/*.SJ.tab \
                          | xargs basename -a -s .SJ.tab \
                          | cut -f1 -d't'\
                          | uniq)

if [ "${bamList[*]}" != "${sjList[*]}" ]
then
  echo "Error: the BAM and SJ files differ."
  echo
  echo "Number of BAM files: ${#bamList[@]}"
  echo "Number of SJ files: ${#sjList[@]}"
  exit 1
else
  echo "Treating ${#sjList[@]} files."
  echo "Estimate 2.5 min per sample."
fi

#bamList=(${sjList[@]:171})

# Loop on samples ----

for sample in "${bamList[@]}"
do
  
  echo "------   Sample $sample   ------"
  
  my_bam=$bam_dir/$sample.bam
  my_SJ_tab=$sj_dir/$sample.SJ.tab
  
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
  
done


## PSI ----


echo
echo "#######################################"
echo "Assembling the files and computing PSI."
echo "#######################################"
echo

module swap BEDTools R/4.2.0-foss-2020b

Rscript src/assemble_psi.R \
        --se ${SE_annot}SE_coords.tab \
        --read_len $readLength \
        --inputs $tmp_dir \
        --out $out_dir



echo "Ending $(date)"

