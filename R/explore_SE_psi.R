# Load, explore, format, and export the SE PSI data


library(tidyverse)
library(wbData)


# Load and format data ----

psi <- read_tsv("data/from_ruddle/2022-11-10_SE_PSI/assembled_psi.tsv",
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
                  .default = col_double()
                )) |>
  rename(event_id = name,
         gene_id = gene) |>
  mutate(across(starts_with(c("ER", "IR")), as.integer)) |>
  arrange(chr, startExon, endExon)

psi_long <- pivot_longer(psi,
                         starts_with(c("IR","ER","PSI")),
                         names_to = c(".value", "sample_id"),
                         names_pattern = "(ER_raw|IR_raw|PSI)_([A-Z0-9ef]{2,4}r[0-9]{1,5})") |>
  rename(IR = IR_raw,
         ER = ER_raw) |>
  mutate(nb_reads = IR + ER,
         .before = PSI) |>
  mutate(neuron_id = str_match(sample_id, "([A-Z0-9ef]{2,4})r[0-9]{1,5}")[,2],
         .before = sample_id)


psi_long

# Explorations ----
xx <- psi_long$PSI
xx[xx>0 & xx<1] <- 0.5
table(xx, useNA = 'ifany')
hist(psi_long$PSI)

xx <- psi_long$PSI
xx[xx == 0 | xx == 1] <- NA
hist(xx, breaks = 100)


psi_long |>
  filter(IR > 15 | ER > 15) |>
  filter(event_id == "SE_1") |>
  pull(PSI) |> range() |> diff()

psi_long |>
  filter(IR > 15 | ER > 15) |>
  mutate(coords = paste0(chr,":",startLong,"..",endLong),
         .after = event_id) |>
  View()


# Look for differentially skipped exons
potential_DSE <- psi_long |>
  filter(IR > 10 | ER > 10) |>
  group_by(event_id) |>
  summarize(Dpsi = diff(range(PSI)))

hist(potential_DSE$Dpsi, breaks = 50)

potential_DSE |>
  filter(Dpsi > .99)



# Look for genes with multiple events
psi_long |>
  filter(nb_reads > 10) |>
  select(gene_id, event_id) |>
  distinct() |>
  count(gene_id) |>
  arrange(desc(n))

# ex: mel-11  WBGene00003196   II:9,357,531..9,369,596



