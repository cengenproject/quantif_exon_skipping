
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
list_SE_genes <- read_tsv(opt$tab_file, show_col_types = FALSE) |>
  filter(startsWith(name, "SE_")) |>
  pull(gene)



# Prepare reference data ----
txdb <- wb_load_TxDb(opt$ws)
tx2g <- wb_load_tx2gene(opt$ws)

dna_seq <- Biostrings::readDNAStringSet(wb_get_genome_path(opt$ws))


# Make sets for our genes of interest
all_genes_coords <- genes(txdb)[list_SE_genes]
all_genes_seqs <- dna_seq[all_genes_coords]

all_introns_genome <- GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)

all_introns <- purrr::map(wb_g2tx(list_SE_genes, tx2g, warn_missing = TRUE),
                   ~ all_introns_genome[.x])


all_exons <- GenomicFeatures::exonsBy(txdb, "gene")[list_SE_genes]


# # Look at real data to define scores
# unlist(all_exons) |>
#   width() |>
#   log10() |>
#   hist(breaks = 100)
# abline(v = log10(c(90, 180)), col = 'green4')
# abline(v = log10(c(60, 300)), col = 'red3')
# abline(v = log10(c(35, 800)), col = 'purple2')
# abline(v = log10(c(18, 3000)), col = 'chocolate1')
# #> ideal exons between 90 and 180 bp,
# #> typical between 60 and 300 bp,
# #> acceptable between 35 and 800 bp
# #> not acceptable outside 18-3000 bp.
# 
# 
# dna_seq[unlist(all_exons)] |>
#   Biostrings::letterFrequency("GC", as.prob = TRUE) |>
#   hist(breaks = 100)
# abline(v = c(.42, .48), col = 'green4')
# abline(v = c(.38, .52), col = 'red3')
# abline(v = c(.3, .6), col = 'purple2')
# abline(v = c(.25, .65), col = 'chocolate1')
# #> ideal exons between .42 and .48 bp,
# #> typical between .38 and .52 bp,
# #> acceptable between .3 and .6 bp
# #> not acceptable outside .25-.65 bp.





# Functions ----

find_internal_ss_pair <- function(intron_seq, strand){
  if(strand == "+"){
    matches <- Biostrings::matchLRPatterns(Lpattern = "GT", Rpattern = "AG",
                                           max.gaplength = length(intron_seq),
                                           subject = intron_seq)
  } else if(strand == "-"){
    matches <- Biostrings::matchLRPatterns(Lpattern = "CT", Rpattern = "AC",
                                           max.gaplength = length(intron_seq),
                                           subject = intron_seq)
  } else{
    stop("Error determining strand")
  }
  
  matches[start(matches) != 1 & end(matches) != length(intron_seq)]
}



select_ce_exon <- function(my_gene, all_introns, all_exons){
  
  # introns that don't intersect with exon
  introns <- all_introns[[my_gene]] |>
    unlist() |>
    unique()
  intersecting <- GenomicRanges::findOverlaps(introns,
                                              all_exons[[my_gene]] |> unique()) |>
    from()
  introns <- introns[-intersecting]
  
  
  potential_exons <- map2(dna_seq[introns],
                          strand(introns) |> as.factor(),
                          find_internal_ss_pair) |>
    setNames(NULL)
  
  
  
  # we have a list of many potential exons, now we should give them scores and pick one
  
  imap_dfr(potential_exons,
           ~ tibble(intron = .y,
                    gc_content = Biostrings::letterFrequency(.x, "GC",as.prob = TRUE)[,1],
                    width = width(.x)) |>
             mutate(view_nb_in_intron = row_number())
  ) |>
    mutate(score_width = case_when(width > 90 & width < 180 ~ 10,
                                   width > 60 & width < 300 ~ 7,
                                   width > 35 & width < 800 ~ 2,
                                   width > 18 & width < 3000 ~ 1,
                                   .default = 0),
           score_gc = case_when(gc_content > .42 & gc_content < .48 ~ 10,
                                gc_content > .38 & gc_content < .52 ~ 7,
                                gc_content > .3 & gc_content < .6 ~ 2,
                                gc_content > .25 & gc_content < .65 ~ 1,
                                .default = 0),
           score_total = score_width + score_gc) |>
    filter(score_total > 0) |>
    slice_sample(n = 1, weight_by = exp(score_total) ) |>
    (\(.x) tibble(gene = my_gene,
                  chr = seqnames(introns) |> runValue() |> as.character(),
                  strand = strand(introns) |> runValue() |> as.character(),
                  startLong = start(introns[.x$intron]),
                  endLong = end(introns[.x$intron]),
                  startExon = startLong + start(potential_exons[[.x$intron]])[[.x$view_nb_in_intron]] - 1,
                  endExon = startLong + end(potential_exons[[.x$intron]])[[.x$view_nb_in_intron]])
    )()
}








# run function
set.seed(123)
CE_coords <- map_dfr(list_SE_genes[1:10], select_ce_exon, all_introns, all_exons,
                     .progress = TRUE)

CE_coords <- mutate(CE_coords,
                    name = paste0("CE_", row_number()),
                    .before = gene)

write_tsv(CE_coords,
          opt$tab_file,
          col_names = FALSE,
          append = TRUE)

message("Done generating CE events, ", date())
