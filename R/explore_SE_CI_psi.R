# Very quick exploration of the PSI values (note, deeper analyses in other projects)

library(tidyverse)

psi <- read_tsv("data/export_for_arman/230228_PSI_quantifications.tsv") |>
  mutate(type = str_match(event_id, "^(CI|SE)_[0-9]+$")[,2])

table(psi$type, useNA = 'ifany')

# Look at distribution

ggplot(psi) +
  theme_classic() +
  geom_histogram(aes(x = PSI), color = 'white') +
  facet_wrap(~type)

table(psi$PSI[psi$type == "CI"])


ggplot(psi) +
  theme_classic() +
  geom_histogram(aes(x = nb_reads), color = 'white') +
  facet_wrap(~type) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Number of supporting reads")


# Compare with previous version (without CI), to make sure rerunning everything didn't introduce new problems

psi_old <- read_tsv("data/export_for_arman/221110_PSI_quantifications.tsv")
psi_new <- psi |> filter(type == "SE") |> select(-type)

psi_old
psi_new

waldo::compare(psi_old, psi_new)





