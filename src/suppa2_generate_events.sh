#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=suppa_genEvents
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --time=10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


echo "Start $(date)"


module load miniconda
conda activate SUPPA2


WS="WS281"
ref_dir="/gpfs/ycga/project/ysm/hammarlund/aw853/references/$WS"
ref_gtf=$ref_dir"/c_elegans.PRJNA13758.${WS}.canonical_geneset.gtf"

out_dir="data/suppa2_data/221108_events"

mkdir -p $out_dir


# SUPPA2 generate events

suppa.py generateEvents \
    --mode DEBUG \
    --event-type SE \
    -f ioe \
    -i $ref_gtf \
    -o $out_dir/${WS}

### Reformat result for easy conversion to BED
# From event name, extract coords of the skipped exon and the skipping intron, to tsv

# Initial version used an intermediate file, not necessary in the end, but useful 
# for diagnostics
# echo -e "name\tgene\tchr\tstrand\tstartLong\tendLong\tstartExon\tendExon" > SE_coords.tab
# cat $out_dir/${WS}_SE_strict.ioe \
#     | sed -nr 's/.*\t(WBGene[0-9]{8});SE:([IVX]{1,3}):([0-9]+)-([0-9]+):([0-9]+)-([0-9]+):([+\-])\t.*/\1\t\2\t\7\t\3\t\6\t\4\t\5/p' \
#     | awk -v'OFS=\t' '{print "SE_"NR, $0}' \
#      >> $out_dir/${WS}_SE_coords.tab


# BED file for inclusion (exonic) reads, contains coords of central exon (the one that can be skipped)
cat $out_dir/${WS}_SE_strict.ioe \
    | sed -nr 's/.*\t(WBGene[0-9]{8});SE:([IVX]{1,3}):([0-9]+)-([0-9]+):([0-9]+)-([0-9]+):([+\-])\t.*/\1\t\2\t\7\t\3\t\6\t\4\t\5/p' \
    | awk -v'OFS=\t' '{print "SE_"NR, $0}' \
    | awk -v'OFS="\t"' '{print($3, $7, $8, $1, 0, $4)}' \
    > $out_dir/${WS}_SE_central_exon.bed

# junctions tab file to match the ones created by STAR, with the coords of the intron
cat $out_dir/${WS}_SE_strict.ioe \
    | sed -nr 's/.*\t(WBGene[0-9]{8});SE:([IVX]{1,3}):([0-9]+)-([0-9]+):([0-9]+)-([0-9]+):([+\-])\t.*/\1\t\2\t\7\t\3\t\6\t\4\t\5/p' \
    | awk -v'OFS=\t' '{print "SE_"NR, $0}' \
    | awk -v 'OFS=\t' '{print($3,$5+1, $6-1, $4=="+"? 1: 2)}' \
    > $out_dir/${WS}_SE_spanning_intron.tab



echo "End $(date)"

