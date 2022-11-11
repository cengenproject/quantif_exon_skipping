# Explore data, for a single sample

library(tidyverse)




psi <- read_tsv("data/from_ruddle/221108_ADLr173_psi.tsv",
                col_types = cols(
                  name = col_character(),
                  gene = col_character(),
                  chr = col_factor(levels = c("I","II","III","IV","V","X")),
                  strand = col_factor(levels = c("+","-")),
                  startLong = col_integer(),
                  endLong = col_integer(),
                  startExon = col_integer(),
                  endExon = col_integer(),
                  exon_length = col_integer(),
                  IR_raw_ADLr173 = col_integer(),
                  ER_raw_ADLr173 = col_integer(),
                  PSI_ADLr173 = col_double()
                )) |>
  arrange(chr, startExon, endExon)



xx <- psi$PSI_ADLr173
xx[xx>0 & xx<1] <- 0.5
table(xx, useNA = 'ifany')
hist(psi$PSI_ADLr173)

xx <- psi$PSI_ADLr173
xx[xx == 0 | xx == 1] <- NA
hist(xx, breaks = 100)




psi |>
  # filter(PSI_ADLr173 > .1 & PSI_ADLr173 < .9,
  #        IR_raw_ADLr173 > 15 & ER_raw_ADLr173 > 15) |>
  mutate(coords = paste0(chr,":",startLong,"..",endLong),
         .after = name) |>
  relocate(IR_raw_ADLr173, ER_raw_ADLr173, PSI_ADLr173, .after = coords) |>
  View()



