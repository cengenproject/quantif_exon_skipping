# Compute PSI per sample
# 
# Started with SUPPA2 generateEvents for SE, made it into a tab file.
# For each sample, we already created files with inclusion and exclusion counts,
# now for each sample we merge the incl and excl to the ref, and also compute PSI.


# Initializations ----


cat("Starting: ", date(),"\n")


#~ packages ----
suppressPackageStartupMessages({
  library(tidyverse)
})

library(optparse)


# Inputs
if(interactive()){
  opt <- list(
    se = "/home/aw853/project/quantif_exon_skipping/data/suppa2_data/230531_events/WS281_SE_coords.tab",
    read_len = 101,
    inputs = "/home/aw853/scratch60/230531_SE_PSI_all_samples",
    out_dir = "/home/aw853/project/quantif_exon_skipping/data/2023-05-31_SE_PSI"
  )
} else{
  
  
  option_list <- list( 
    make_option(c("-s", "--se"), type = "character",
                help = "Path to SE_coords.tab file adapted from SUPPA2 ioe."),
    make_option(c("-l", "--read_len"), type = "integer", default = 101,
                help = "Read length."),
    make_option(c("-i", "--inputs"), type = "character",
                help = "Path to dir with the sample.inclusion and sample.exclusion files. Note the sample names will be deduced from these file names."),
    make_option(c("-o", "--out_dir"), type="character",
                help="Output directory.")
  )
  opt <- parse_args(OptionParser(option_list=option_list))
  
}


# ~ validate inputs ----
stopifnot(file.exists(opt$se))
stopifnot(dir.exists(opt$inputs))
stopifnot(dir.exists(opt$out_dir))




#~ Functions ----

read_incl <- function(path){
  incl <- read_tsv(path,
                   col_names = c("chr","startExon","endExon", "name","bedScore", "strand",
                                 "depth", "length_at_depth", "total_length", "fraction_covered"),
                   col_types = cols(
                     chr = col_factor(levels = c("I","II","III","IV","V","X")),
                     startExon = col_integer(),
                     endExon = col_integer(),
                     name = col_character(),
                     bedScore = col_integer(),
                     strand = col_factor(levels = c("+","-")),
                     depth = col_integer(),
                     length_at_depth = col_integer(),
                     total_length = col_integer(),
                     fraction_covered = col_double()
                   ),
                   col_select = -c(name, bedScore, length_at_depth, fraction_covered))
  stop_for_problems(incl)
  incl
}


read_excl <- function(path){
  # note: the multimapper and overhang columns of the SJ.tab file were removed during merging
  excl <- read_tsv(path,
                   col_names = c("chr", "startLong", "endLong", "strand",
                                 "intron_motif", "annotated","unimappers"),
                   col_types = cols(
                     chr = col_factor(levels = c("I","II","III","IV","V","X")),
                     startLong = col_integer(),
                     endLong = col_integer(),
                     strand = col_factor(levels = c("1", "2")),
                     intron_motif = col_factor(levels = as.character(0:6)),
                     annotated = col_factor(levels = c("0", "1")),
                     unimappers = col_integer()
                   )) |>
    mutate(strand = recode_factor(strand,
                                  "1" = "+",
                                  "2" = "-"),
           startLong = startLong - 1,
           endLong = endLong + 1) |>
    filter(annotated == "1")
  stop_for_problems(excl)
  
  
  excl
}



# note: the same exon can be in several events (with diff intron),
# so we have duplicates in `included` df.
# On the other hand, there are fewer lines in `excluded` df as the SJ is not reported by STAR
# when there are no skipping reads; so we get NA after merging.
compute_psi <- function(incl, excl){
  se_coords |>
    full_join(unique(incl),
              by = c("chr", "strand", "startExon", "endExon")) |>
    full_join(excl,
              by = c("chr", "strand", "startLong", "endLong")) |>
    select(name:endExon,
           exon_length = total_length,
           IR_raw = depth,
           ER_raw = unimappers) |>
    mutate(ER_raw = replace_na(ER_raw, 0L),
           IRn = IR_raw / (exon_length + opt$read_len - 1),
           ERn = ER_raw / (opt$read_len - 1),
           PSI = IRn / (IRn + ERn)) |>
    select(- IRn, -ERn)
}



# Main code ----

se_coords <- read_tsv(opt$se,
                      col_types = cols(
                        name = col_character(),
                        gene = col_character(),
                        chr = col_factor(levels = c("I","II","III","IV","V","X")),
                        strand = col_factor(levels = c("+","-")),
                        startLong = col_integer(),
                        endLong = col_integer(),
                        startExon = col_integer(),
                        endExon = col_integer()
                      ))
stop_for_problems(se_coords)




paths_incl <- list.files(opt$input, pattern = "\\.inclusion$")
paths_excl <- list.files(opt$input, pattern = "\\.exclusion$")

samples <- str_remove(paths_incl, pattern = "\\.inclusion$")

stopifnot(all.equal(samples,
                    str_remove(paths_excl, pattern = "\\.exclusion$")))

message("Treating ",length(samples), " samples: ", head(samples), "...")


all_res <- tibble(sample = samples,
                  included_reads_df = map(file.path(opt$input, paths_incl),
                                          read_incl),
                  excluded_reads_df = map(file.path(opt$input, paths_excl),
                                          read_excl)) |>
  mutate(res = map2(included_reads_df,
                    excluded_reads_df,
                    compute_psi)) |>
  select(sample, res) |>
  unnest(res) |>
  pivot_wider(values_from = c(IR_raw, ER_raw, PSI),
              names_from = sample)


write_tsv(all_res,
          file.path(opt$out_dir, "assembled_psi.tsv"))







