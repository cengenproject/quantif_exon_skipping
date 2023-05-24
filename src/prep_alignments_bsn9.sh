#!/bin/bash
#SBATCH --partition=week
#SBATCH --job-name=merge_index_bam
#SBATCH -c 1
#SBATCH --mem=25G
#SBATCH --time=5-00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


# Transfer the bam files merging the technical replicates. Index the results (bai)
# Also merge the SJ technical replicates.

echo "Prepare alignments bsn9 starting $(date)"


input_bams="/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/bsn9_bams/"
input_sj="/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/bsn9_junctions/"
out_dir="/home/aw853/scratch60/2023-05-24_bsn9/"


### Prepare metadata
mkdir -p $out_dir
mkdir -p $out_dir/bams
mkdir -p $out_dir/junctions


mapfile -t sampleList < <(ls $input_bams/*.bam \
                          | xargs basename -a -s .bam \
                          | cut -f1 -d't'\
                          | uniq)

echo ${#sampleList[@]}" samples"



echo "-----   Merging BAMs   -----"

module load SAMtools

for sample in ${sampleList[@]}
do
  samtools merge -@ $SLURM_CPUS_PER_TASK $out_dir/bams/$sample".bam" $(echo $input_bams/$sample"*.bam")
done

echo "-----   BAM files merged. Indexing.   -----"

for sample in ${sampleList[@]}
do
  if [ -f "$out_dir/bams/$sample.bam.bai" ]; then
    echo "  $sample already indexed."
  else 
    echo "  Indexing $sample"
    samtools index -@ $SLURM_CPUS_PER_TASK $out_dir/bams/$sample".bam" $out_dir/bams/$sample".bam.bai"
  fi
done


echo "-----   Merging SJ   -----"
module swap SAMtools R

Rscript src/combine_sj_files.R \
  --sj_dir /SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/bsn9_junctions/ \
  --out_dir $out_dir/junctions

echo
echo "All done $(date)"


echo "TO DO manually: if needed, remove IL2r107 (no bam, only sj); remove outliers (AVKr113, PVMr122, RICr133)"
