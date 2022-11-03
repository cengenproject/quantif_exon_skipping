# For SJ.out.tab files output by STAR, combine technical replicates


# Initializations ----


cat("Starting: ", date(),"\n")


#~ packages ----
suppressPackageStartupMessages({
  library(tidyverse)
})

library(optparse)

option_list <- list( 
  make_option(c("-s", "--sj_dir"), type = "character",
              default = "/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/bsn9_junctions/",
              help = "Path to input directory containing JS files created by STAR"),
  make_option(c("-o", "--out_dir"), type="character",
              help="Output directory")
)


# ~ validate inputs ----
opt <- parse_args(OptionParser(option_list=option_list))

stopifnot(dir.exists(opt$sj_dir))
stopifnot(dir.exists(opt$out_dir))



#~ functions ----

read_sj_file <- function(cur_path){
  read_table(cur_path,
             col_names=c("chr", "start","end","strand","motif", "annotated",
                         "count_unique","count_multimap", "max_overhang"),
             col_types = cols(
               chr = col_factor(levels = c("I","II","III","IV","MtDNA", "V", "X")),
               start = col_integer(),
               end = col_integer(),
               strand = col_factor(levels = as.character(0:2)),
               motif = col_factor(levels = as.character(0:6)),
               annotated = col_factor(levels = as.character(0:1)),
               count_unique = col_integer(),
               count_multimap = col_integer(),
               max_overhang = col_integer()
             )) |>
    select(-max_overhang, -count_multimap)
}



# fn must be a function operating on the rows of a matrix (e.g. rowSums, rowMaxs, etc...)
combine_sj <- function(sj_file, fn){
  reduce(sj_file, full_join, by = c("chr", "start", "end", "strand", "motif", "annotated")) |>
    mutate(across(starts_with("count"), replace_na, 0)) %>%
    mutate(count_unique = fn(as.matrix(select(., starts_with("count_unique"))))) |>
    select(-starts_with("count_unique.")) |>
    arrange(chr, start, end)
}



#~ process samples ----
cat("Processing.\n")

tibble(path = list.files(opt$sj_dir, full.names = TRUE),
       replicate = stringr::str_split_fixed(basename(path), "\\.", 2)[,1],
       sample = stringr::str_split_fixed(replicate, "t", 2)[,1]) |>
  mutate(sj_file = map(path, read_sj_file)) |>
  group_by(sample) |>
  summarize(sj_file_combined = list(combine_sj(sj_file, rowSums))) |>
  mutate(out_path = file.path(opt$out_dir, sample) |>
           paste0(".SJ.tab")) |>
  select(x = sj_file_combined,
         file = out_path) |>
  pwalk(write_tsv)


