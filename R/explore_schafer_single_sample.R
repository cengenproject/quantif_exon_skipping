# Explore data, for a single sample

library(tidyverse)


# load reference
ref <- rtracklayer::readGFF("data/from_ruddle/WS281_schafer_index.gff") |>
  select(exon_ID = group,
         chrom = seqid,
         start, end, strand)

# Annotate the cases where exonic_part is at start or end of tx
tx_coords <- wbData::wb_get_gtf_path(281) |>
  rtracklayer::readGFF(filter = list(type = "transcript")) |>
  select(transcript_id, seqid,
         tx_start = start, tx_end = end,
         strand)

# Expand all tx for each exonic_part
exons_with_tx_coords <- ref |>
  mutate(tx_id = str_split_fixed(ref$exon_ID, ":", n=2)[,2]) |>
  separate_rows(tx_id, sep = "\\+") |>
  left_join(tx_coords,
            by = c(tx_id = "transcript_id",
                   chrom = "seqid",
                   strand = "strand"))
  
exons_terminals <- exons_with_tx_coords |>
  mutate(is_terminal = (start == tx_start | end == tx_end)) |>
  group_by(exon_ID) |>
  summarize(is_terminal = any(is_terminal))

ref <- ref |>
  left_join(exons_terminals, by = "exon_ID")


psi <- read_tsv("data/from_ruddle/ADLr173.psi",
                col_types = cols(
                  exon_ID = col_character(),
                  length = col_integer(),
                  inclusion = col_integer(),
                  exclusion = col_integer(),
                  PSI = col_double()
                )) |>
  left_join(ref, by = "exon_ID") |>
  arrange(chrom, start, end)

psi |>
  filter(PSI > .1 & PSI < .9,
         inclusion > 15 & exclusion > 15,
         ! is_terminal) |>
  mutate(coords = paste0(chrom,":",start,"..",end),
         .after = exon_ID) |>
  View()


xx <- psi$PSI
xx[xx>0 & xx<1] <- 0.5
table(xx, useNA = 'ifany')
hist(psi$PSI)

xx <- psi$PSI
xx[xx == 0 | xx == 1] <- NA
hist(xx, breaks = 100)


