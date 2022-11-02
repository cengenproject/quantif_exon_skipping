module load BEDTools/2.23.0-foss-2016a

ref_dir="/gpfs/ycga/project/ysm/hammarlund/aw853/references/WS277"                                                                                                                                                                                                                                                                                                                                                                    schafer_index=$ref_dir"/schafer_index"                                                                                                                                                                                ref_gff=$ref_dir"/c_elegans.PRJNA13758.WS277.canonical_geneset.gtf" 


# Flatten exons (i.e. create counting bins)

module load HTSeq/0.6.1p1-foss-2016b-Python-2.7.13
python ~/.utilities/DEXSeq/dexseq_prepare_annotation.py $ref_gff WS277_flattened.gff

# select "exonic_part" entries, rename them gene_id:exonic_part_number (e.g. WBGene00000001:001)
awk '{OFS="\t"}{if ($3 == "exonic_part") print $1,$2,$3,$4,$5,$6,$7,$8,$14":"$12}' WS277_flattened.gff | sed 's=[";]==g' > WS277_exonic_parts.gff

# Filter exons smaller than 5 bp
awk '$5-$4>4{print}' WS277_exonic_parts.gff > $ref_dir"/schafer_index.gff"

rm WS277_exonic_parts.gff WS277_flattened.gff




module load BEDTools/2.23.0-foss-2016a

# For each sample, ex AIYr65
my_bam="/home/aw853/scratch60/counts/work/e5/2fedf2aab21f9bd6832581c138a474/Aligned.sortedByCoord.out.bam"
my_SJ_tab="/home/aw853/scratch60/counts/work/e5/2fedf2aab21f9bd6832581c138a474/SJ.out.tab"

coverageBed -split -abam $my_bam -b WS277_exonic_parts.gff | \
	awk 'BEGIN{OFS="\t"}{print $1,$4,$5,$5-$4+1,$9,$10}' | \
	sort -k 5 > AIYr65.inclusion

# Convert STAR output to emulate TopHat
awk 'BEGIN{OFS="\t"}{print $1, $2-20-1, $3+20, "JUNCBJ"NR, $7, ($4 == 1)? "+":"-",$2-20-1, $3+20, "255,0,0", 2, "20,20", "0,300" }' $my_SJ_tab > AIYr65_junctions.bed

# "left" covers the 20bp at the exon start, "right" covers the last 20bp of the exon. Here junctions are found by STAR.
sed 's/,/\t/g' AIYr65_junctions.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+$13,$4,$5,$6}' > left.bed
sed 's/,/\t/g' AIYr65_junctions.bed | awk 'BEGIN{OFS="\t"}{print $1,$3-$14,$3,$4,$5,$6}' > right.bed
# Filter to keep only introns with both ends annotated
intersectBed -u -s -a left.bed -b WS277_exonic_parts.gff > left.overlap
intersectBed -u -s -a right.bed -b WS277_exonic_parts.gff > right.overlap
cat left.overlap right.overlap | cut -f4 | sort | uniq -c | awk '{ if($1 == 2) print $2 }' > filtered_junctions.txt
grep -F -f filtered_junctions.txt AIYr65_junctions.bed > AIYr65_filtered_junctions.bed



# remove the overhang
sed 's/,/\t/g' AIYr65_filtered_junctions.bed | grep -v description | awk '{OFS="\t"}{print $1,$2+$13,$3-$14,$4,$5,$6}' > AIYr65_intron.bed


# Count number of junction reads spanning the exonic part (i.e. excluded reads)
intersectBed -wao -f 1.0 -s -a WS277_exonic_parts.gff -b AIYr65_intron.bed | \
	awk 'BEGIN{OFS="\t"}{$16 == 0? s[$9] += 0:s[$9] += $14}END{for (i in s) {print i,s[i]}}' | \
	sort -k 1 > exonic_parts.exclusion



readLength=101
paste exonic_parts.inclusion exonic_parts.exclusion | awk -v "len=$readLength" 'BEGIN{OFS="\t"; print "exon_ID" , "length" , "inclusion" , "exclusion" , "PSI"}{NIR=$6/($4+len-1) ; NER=$8/(len-1)}{print $5,$4,$6,$8,(NIR+NER<=0)? "NA":NIR / (NIR + NER)}' > AIY_exonic_parts.psi


# Clean up
rm left.bed right.bed left.overlap right.overlap filtered_junctions.txt 
rm filtered_junctions.bed 
rm intron.bed
rm exonic_parts.inclusion
rm exonic_parts.exclusion 




