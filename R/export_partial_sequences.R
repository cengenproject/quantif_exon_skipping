# Export sequence up to the exon, for structure

# Partially copied from `export_gene_sequences.R`


# Inits ----
library(GenomicRanges)
library(tidyverse)
library(wbData)

WS <- 289
path_ev_coords <- "data/export_for_arman/240308_events_coordinates.tsv"
path_export <- "data/intermediates_for_DL/240718_partial_sequences_match_240308_events.fa"

#~ Load reference data ----
fasta <- Rsamtools::FaFile(str_remove(wb_get_genome_path(WS), "\\.gz$"))
if(! file.exists(fasta$path)){
  R.utils::gunzip(wb_get_genome_path(WS))
}
gene_coords <- wb_load_gene_coords(WS)





# Events ----
event_coords <- read_tsv(path_ev_coords) |>
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




Biostrings::writeXStringSet(eseq, path_export)

R.utils::gzip(path_export)









