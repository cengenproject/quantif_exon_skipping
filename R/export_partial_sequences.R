# Export sequence up to the exon, for structure

# Partially copied from `export_gene_sequences.R`


# Inits ----
library(GenomicRanges)
library(tidyverse)
library(wbData)



#~ Load reference data ----
fasta <- Rsamtools::FaFile(str_remove(wb_get_genome_path(281), "\\.gz$"))
if(! file.exists(fasta$path)){
  R.utils::gunzip(wb_get_genome_path(281))
}
gene_coords <- wb_load_gene_coords(281)




# genes
gr <- GRanges(seqnames = gene_coords$chr,
              strand = gene_coords$strand,
              IRanges(gene_coords$start, gene_coords$end),
              names = gene_coords$gene_id)

gseq <- Biostrings::getSeq(fasta, gr)
names(gseq) <- gene_coords$gene_id



# Events ----
event_coords <- read_tsv("data/export_for_arman/230912_events_coordinates.tsv") |>
  left_join(gene_coords,
                           by = "gene_id") |>
# to genomic coordinates
  mutate(start_g = if_else(strand == "+",
                                start,
                           end - exon_end + 1),
         end_g = if_else(strand == "+",
                              start + exon_end - 1,
                              end))

event_coords |>
  select(event_id, name, strand, start_g, end_g) |>
  mutate(width = end_g - start_g)

egr <- GRanges(seqnames = event_coords$chr,
               strand = event_coords$strand,
               IRanges(event_coords$start_g, event_coords$end_g),
               names = event_coords$event_id)

eseq <- Biostrings::getSeq(fasta, egr)
names(eseq) <- event_coords$event_id




Biostrings::writeXStringSet(gseq, "data/intermediates_for_DL/240221_partial_sequences.fa")

R.utils::gzip("data/intermediates_for_DL/240221_partial_sequences.fa")









