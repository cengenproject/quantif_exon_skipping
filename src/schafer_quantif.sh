#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=schafer_quantif
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH --time=1-00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


echo "Schafer quantification. Starting $(date)"



WS="WS281"
readLength=101


ref_dir="/gpfs/ycga/project/ysm/hammarlund/aw853/references/$WS"                                                                                                                                                                                                                                                                                                                                                                   
gff=$ref_dir"/${WS}_schafer_index.gff"

data_dir="/home/aw853/scratch60/2022-11-03_bsn9/"
out_dir="2022-11-03_schafer_psi/"
tmp_dir="/home/aw853/scratch60/221103_schafer/"





module load BEDTools/2.23.0-foss-2016a

# For each sample, ex ADLr173
sample="ADLr173"
my_bam=$data_dir/bams/$sample.bam
my_SJ_tab=$data_dir/junctions/$sample.SJ.tab



$s=$tmp_dir/${sample}


coverageBed -split -abam $my_bam -b $gff \
	| awk 'BEGIN{OFS="\t"}{print $1,$4,$5,$5-$4+1,$9,$10}' \
	| sort -k 5 \
	> $tmp_dir/$sample.inclusion


# Convert STAR SJ output to emulate TopHat
awk 'BEGIN{OFS="\t"}{                                           \
      print $1, $2-20-1, $3+20, "JUNCBJ"NR, $7,                 \
            ($4 == 1)? "+":"-",$2-20-1, $3+20, "255,0,0", 2,    \
            "20,20", "0,300"                                    \
    }' $my_SJ_tab \
  > $s.junctions.bed


# "left" covers the 20bp at the exon start, "right" covers the last 20bp of the exon. Here junctions are found by STAR.
sed 's/,/\t/g' $s.junctions.bed \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+$13,$4,$5,$6}' \
  > $s.left.bed
sed 's/,/\t/g' $s.junctions.bed \
  | awk 'BEGIN{OFS="\t"}{print $1,$3-$14,$3,$4,$5,$6}' \
  > $s.right.bed

# Filter to keep only introns with both ends annotated
intersectBed -u -s -a $s.left.bed -b $gff \
  > $s.left.overlap
intersectBed -u -s -a $s.right.bed -b $gff \
  > $s.right.overlap

cat $s.left.overlap $s.right.overlap \
  | cut -f4 \
  | sort \
  | uniq -c \
  | awk '{ if($1 == 2) print $2 }' \
  > $s.filtered_junctions.txt

grep -F -f $s.filtered_junctions.txt $s.junctions.bed \
  > $s.filtered_junctions.bed



# remove the overhang
sed 's/,/\t/g' $s.filtered_junctions.bed \
  | grep -v description \
  | awk '{OFS="\t"}{print $1,$2+$13,$3-$14,$4,$5,$6}' \
  > $s.intron.bed


# Count number of junction reads spanning the exonic part (i.e. excluded reads)
intersectBed -wao -f 1.0 -s -a $gff -b $s.intron.bed \
  |	awk 'BEGIN{OFS="\t"}                             \
              {$16 == 0? s[$9] += 0:s[$9] += $14}    \
           END{for (i in s) {print i,s[i]}}'         \
	| sort -k 1 \
	> $s.exonic_parts.exclusion




paste $s.exonic_parts.inclusion $s.exonic_parts.exclusion \
  | awk -v "len=$readLength" \
      'BEGIN{OFS="\t";                                                         \
             print "exon_ID" , "length" , "inclusion" , "exclusion" , "PSI"}   \
       {NIR=$6/($4+len-1) ; NER=$8/(len-1)}                                    \
       {print $5,$4,$6,$8,(NIR+NER<=0)? "NA":NIR / (NIR + NER)}' \
  > $out_dir/$sample.psi




echo "Done $(date)"

