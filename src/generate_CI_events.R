
# From the WS gtf, find constitutive exons for each gene of interest (that contains a SE)

message("Starting generate CI events, ", date())

# Inits
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
})
library(wbData)



# read arguments
library(optparse)

option_list <- list( 
  make_option(c("-t", "--tab_file"), type = "character",
              help = "Path to input SE_coords.tab file containing skipped events records, CI records will be appended to that file"),
  make_option(c("-w", "--ws"), type = "character",
              help = "Wormbase WS version")
)


# ~ validate inputs ----
if(interactive()){
  opt <- list(tab_file = "data/suppa2_data/221108_events/WS281_SE_coords.tab",
              ws = "WS281")
} else{
  opt <- parse_args(OptionParser(option_list=option_list))
}

stopifnot(file.exists(opt$tab_file))

options(wb_dir_cache = paste0("/gpfs/ycga/project/ysm/hammarlund/aw853/references/", opt$ws))



# List genes that need an exon ----
list_SE <- read_tsv(opt$tab_file, col_select = "gene", show_col_types = FALSE)



# Prepare reference data ----
txdb <- wb_load_TxDb(opt$ws)
tx2g <- wb_load_tx2gene(opt$ws)


all_introns <- GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)
all_exons <- GenomicFeatures::exonsBy(txdb, "tx", use.names = TRUE)



# main function
select_ci_exon <- function(my_gene, all_introns, all_exons){
  # get all introns in that gene
  introns <- wb_g2tx(my_gene, tx2g, warn_missing = TRUE) |>
    unlist() |>
    {\(.x) all_introns[.x]}() |>
    unlist() |>
    unique()
  
  
  # get all exons in that gene
  exons <- wb_g2tx(my_gene, tx2g, warn_missing = TRUE) |>
    unlist() |>
    {\(.x) all_exons[.x]}() |>
    unlist() |>
    unique()
  
  
  #~ Find internal, constitutive, exons ----
  
  # exons which overlap with a different exon (e.g. alt 5/3' splice sites)
  exons_with_alt_ss <- findOverlaps(exons, exons) |>
    {\(.x) .x[from(.x) != to(.x)]}() |>
    from()
  
  # exons that are included in intron (e.g. skipped or mutually exclusive exons)
  exons_skippable <- which(exons %within% introns)
  
  
  # exons that can have two different flanking introns OR that are terminal
  # in practice, find the flanking introns for each exon, ensure there are exactly 1 on each side
  ex_intr_overlap <- findOverlaps(exons, introns, minoverlap = 1L) |> as.data.frame()
  
  left_flanking <- findOverlaps(exons, shift(introns, 1L)) |>
    as.data.frame() |>
    anti_join(ex_intr_overlap, by = c("queryHits", "subjectHits")) |>
    dplyr::count(queryHits)
  
  right_flanking <- findOverlaps(exons, shift(introns, -1L)) |>
    as.data.frame() |>
    anti_join(ex_intr_overlap, by = c("queryHits", "subjectHits")) |>
    dplyr::count(queryHits)
  
  # gather all this info
  internal_constitutive_exons <- left_flanking$queryHits[left_flanking$n == 1L] |>
    intersect(right_flanking$queryHits[right_flanking$n == 1L]) |>
    unique() |>
    setdiff(exons_with_alt_ss) |>
    setdiff(exons_skippable)
  
  if(length(internal_constitutive_exons) == 0L){
    return(tibble())
  }
  
  # select one exon randomly
  selected_exon <- sample(internal_constitutive_exons, 1)
  
  
  left_intron <- findOverlaps(exons[selected_exon,], shift(introns, 1L)) |>
    as.data.frame() |>
    anti_join(ex_intr_overlap, by = c("queryHits", "subjectHits"))
  right_intron <- findOverlaps(exons[selected_exon,], shift(introns, -1L)) |>
    as.data.frame() |>
    anti_join(ex_intr_overlap, by = c("queryHits", "subjectHits"))
  
  flanking_introns <- introns[c(left_intron$subjectHits, right_intron$subjectHits),]
  
  # export
  tibble(gene = my_gene,
         chr = seqnames(exons) |> runValue() |> as.character(),
         strand = strand(exons) |> runValue() |> as.character(),
         startLong = min(start(flanking_introns)) - 1L,
         endLong = max(end(flanking_introns)) + 1L,
         startExon = start(exons[selected_exon,]),
         endExon = end(exons[selected_exon,]))
}



# run function
set.seed(123)
CI_coords <- map_dfr(list_SE$gene, select_ci_exon, all_introns, all_exons)

CI_coords <- mutate(CI_coords,
                    name = paste0("CI_", seq_len(nrow(CI_coords))),
                    .before = gene)

write_tsv(CI_coords,
          opt$tab_file,
          col_names = FALSE,
          append = TRUE)

message("Done generating CI events, ", date())
