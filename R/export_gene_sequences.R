# Extract the sequence of each gene, and the coordinates of the exons within that tx
# Copied from same script in `stringtie_quantif` Project.


# Inits ----
library(GenomicRanges)
library(tidyverse)
library(wbData)



#~ Load reference data ----
fasta <- Rsamtools::FaFile(str_remove(wb_get_genome_path(289), "\\.gz$"))
if(! file.exists(fasta$path)){
  R.utils::gunzip(wb_get_genome_path(289))
}
gene_coords <- wb_load_gene_coords(289)


# which genes to keep
gene_expression <- read.csv("data/thresholded_gene_expression/bsn9_subtracted_integrated_binarized_expression_withVDDD_FDR_0.1_092022.tsv",
                            sep = "\t")
gene_expression2 <- apply(gene_expression, 1, sum)
genes_expressed_in_neurons <- names(gene_expression2)[gene_expression2 > 0]

rm(gene_expression)
rm(gene_expression2)

gene_coords_filt <- gene_coords |>
  filter(gene_id %in% genes_expressed_in_neurons)


# convert ----
gr <- GRanges(seqnames = gene_coords_filt$chr,
              strand = gene_coords_filt$strand,
              IRanges(gene_coords_filt$start, gene_coords_filt$end),
              names = gene_coords_filt$gene_id)

gseq <- Biostrings::getSeq(fasta, gr)
names(gseq) <- gene_coords_filt$gene_id

# # Save to disk
# Biostrings::writeXStringSet(gseq, "data/intermediates_for_DL/220920_gene_sequences.fa")
# 
# R.utils::gzip("data/intermediates_for_DL/220920_gene_sequences.fa")


# # Save to disk unfiltered
gr <- GRanges(seqnames = gene_coords$chr,
              strand = gene_coords$strand,
              IRanges(gene_coords$start, gene_coords$end),
              names = gene_coords$gene_id)

gseq <- Biostrings::getSeq(fasta, gr)
names(gseq) <- gene_coords$gene_id

Biostrings::writeXStringSet(gseq, "data/intermediates_for_DL/240906_all_gene_sequences.fa")

R.utils::gzip("data/intermediates_for_DL/240906_all_gene_sequences.fa")

