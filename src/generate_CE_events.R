
# From the WS gtf, find constitutive non-exons for each gene of interest (that contains a SE)
#
# More precisely, find splice donor/acceptor motifs inside an intron

message("Starting generate CE events, ", date())

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
  opt <- list(tab_file = "data/from_ruddle/WS281_SE_coords.tab",
    # tab_file = "data/suppa2_data/221108_events/WS281_SE_coords.tab",
              ws = "WS281")
} else{
  opt <- parse_args(OptionParser(option_list=option_list))
}

stopifnot(file.exists(opt$tab_file))

options(wb_dir_cache = paste0("/gpfs/ycga/project/ysm/hammarlund/aw853/references/", opt$ws))



# List genes that need an exon ----
list_SE_genes <- read_tsv(opt$tab_file, col_select = "gene", show_col_types = FALSE)$gene



# Prepare reference data ----
txdb <- wb_load_TxDb(opt$ws)
tx2g <- wb_load_tx2gene(opt$ws)

dna_seq <- Biostrings::readDNAStringSet(wb_get_genome_path(opt$ws))


# Make sets for our genes of interest
all_genes_coords <- genes(txdb)[list_SE_genes]
all_genes_seqs <- dna_seq[all_genes_coords]


all_introns <- map(wb_g2tx(list_SE_genes[1:3], tx2g, warn_missing = TRUE),
                   ~ GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)[.x])

all_exons <- GenomicFeatures::exonsBy(txdb, "gene")[list_SE_genes]


# Look at real data to define scores
unlist(all_exons) |>
  width() |>
  log10() |>
  hist(breaks = 100)
abline(v = log10(c(90, 180)), col = 'green4')
abline(v = log10(c(60, 300)), col = 'red3')
abline(v = log10(c(35, 800)), col = 'purple2')
#> ideal exons between 90 and 180 bp,
#> typical between 60 and 300 bp,
#> acceptable between 35 and 800 bp
#> not acceptable outside 20-4000 bp.


dna_seq[unlist(all_exons)] |>
  Biostrings::letterFrequency("GC", as.prob = TRUE) |>
  hist(breaks = 100)
abline(v = c(.42, .48), col = 'green4')
abline(v = c(.38, .52), col = 'red3')
abline(v = c(.3, .6), col = 'purple2')
#> ideal exons between .42 and .48 bp,
#> typical between .38 and .52 bp,
#> acceptable between .3 and .6 bp
#> not acceptable outside .2-.7 bp.







# introns that don't intersect with exon
introns <- all_introns[[1]] |>
  unlist() |>
  unique()
intersecting <- GenomicRanges::findOverlaps(introns,
                            all_exons[[1]] |> unique()) |>
  from()
introns <- introns[-intersecting]



find_internal_exon <- function(intron_seq){
  matches <- Biostrings::matchLRPatterns(Lpattern = "GT", Rpattern = "AG",
                              max.gaplength = length(intron_seq),
                              subject = intron_seq)
  matches[start(matches) != 1 & end(matches) != length(intron_seq)]
}

potential_exons <- lapply(dna_seq[introns],
       find_internal_exon) |>
  set_names(NULL)

sapply(potential_exons, length) |> sum()
# we have a list of 16,176 views, now we should give them scores (e.g. length typical of real exons) and pick one

imap_dfr(potential_exons,
        ~ tibble(intron = .y,
                 gc_content = Biostrings::letterFrequency(.x, "GC",as.prob = TRUE)[,1],
                 width = width(.x)) |>
          mutate(view_nb_in_intron = row_number())) |>
  mutate(score_width = case_when(width > 90 & width < 180 ~ 10,
                                 width > 60 & width < 300 ~ 7,
                                 width > 35 & width < 800 ~ 2,
                                 width > 20 & width < 4000 ~ 1,
                                 .default = 0),
         score_gc = case_when(gc_content > .42 & gc_content < .48 ~ 10,
                              gc_content > .38 & gc_content < .52 ~ 7,
                              gc_content > .3 & gc_content < .6 ~ 2,
                              gc_content > .2 & gc_content < .7 ~ 1,
                              .default = 0),
         score_total = score_width + score_gc) |>
  filter(score_total > 0) |>
  slice_sample(n = 1, weight_by = exp(score_total)) |>
  (\(.x) potential_exons[[.x$intron]][.x$view_nb_in_intron])()


# table(xx$score_total)
# 
# 0    1    2    3    4    7    8    9   10   11   12   14   17   20 
# 19   81  194 7946 3449   85  605 1641   46  316  848  414  386  146

# run it 100 times
# table(unlist(xx))
# 14 17 20 
# 2 10 88 


Biostrings::vmatchPattern("GT", dna_seq[introns])
Biostrings::matchLRPatterns("GT", "AG", max.gaplength = max(width(introns)), subject = dna_seq[introns][[1]])

Biostrings::matchPDict(Biostrings::PDict(c("GT","AG")), subject = dna_seq[introns][[1]])

donors <- map(all_genes_seqs,
              ~ Biostrings::matchPattern("GT", .x) |>
                start())
acceptor <- map(all_genes_seqs,
              ~ Biostrings::matchPattern("AG", .x) |>
                start())


grep("GT", dna_seq[all_genes][[1]])

Biostrings::matchPattern("GT", dna_seq[all_genes][[1]])
exons_starts <- lapply(all_exons[[my_gene]] |>
  start() |>
  unique())


## TODO!!! Check on a gene on reverse strand! -----




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
         startLong = min(start(flanking_introns)),
         endLong = max(end(flanking_introns)),
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
