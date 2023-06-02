# Very quick exploration of the PSI values (note, deeper analyses in other projects)

library(tidyverse)

psi <- read_tsv("data/export_for_arman/230531_PSI_quantifications.tsv") |>
  mutate(type = str_match(event_id, "^(CI|SE|CE)_[0-9]+$")[,2],
         type = factor(type, levels = c("SE","CI","CE")))

table(psi$type, useNA = 'ifany')

# Look at distribution

ggplot(psi) +
  theme_classic() +
  geom_histogram(aes(x = PSI), color = 'white') +
  facet_wrap(~type)# +scale_y_log10()

table(psi$PSI[psi$type == "CI"] == 1)


ggplot(psi) +
  theme_classic() +
  geom_histogram(aes(x = nb_reads), color = 'white') +
  facet_wrap(~type) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Number of supporting reads")



psi |>
  filter(nb_reads > 4) |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = PSI), color = 'white') +
  facet_wrap(~type)



# Compare with previous version (without CE), to make sure rerunning everything didn't introduce new problems

psi_old <- read_tsv("data/export_for_arman/230228_PSI_quantifications.tsv") |>
  as_tibble() # necessray to remove {readr}-specific spec_tbl_df
psi_new <- psi |> filter(type == "SE"| type == "CI") |> select(-type)

psi_old
psi_new

psi_old |> tail()
psi_new |> tail()

all.equal(psi_old[startsWith(psi_old$event_id, "SE_"),],
          psi_new[startsWith(psi_new$event_id, "SE_"),])


all.equal(psi_old, psi_new)
waldo::compare(psi_old, psi_new)

table(psi_old$nb_reads == psi_new$nb_reads)
head(which(psi_old$nb_reads != psi_new$nb_reads))

psi_old[310,]
psi_new[310,]

plot(psi_old$nb_reads[startsWith(psi_old$event_id, "CI_")],
     psi_new$nb_reads[startsWith(psi_new$event_id, "CI_")])
plot(jitter(psi_old$PSI[startsWith(psi_old$event_id, "CI_")]),
     psi_new$PSI[startsWith(psi_new$event_id, "CI_")], col = alpha("black", 0.1))

table(psi_new$PSI[startsWith(psi_new$event_id, "CI_")] < 1)


psi_new


# Check coordinates ----

coords <- read_tsv("data/export_for_arman/230531_events_coordinates.tsv")

# frm-5.2
coords |> filter(gene_id == "WBGene00021406")


# WBGene00016022
coords |> filter(gene_id == "WBGene00016022")


gene_coords <- Biostrings::readDNAStringSet("../stringtie_quantif/data/intermediates_for_DL/220920_gene_sequences.fa.gz")

gene_coords[["WBGene00016022"]]

# SE_936 intron
gene_coords[["WBGene00016022"]][2331:4384]
# exon
gene_coords[["WBGene00016022"]][3552:3722]


# CI_844 intron
gene_coords[["WBGene00016022"]][6975:7267]
# exon
gene_coords[["WBGene00016022"]][7029:7159]


# CE_922 intron
gene_coords[["WBGene00016022"]][2331:3551]
# exon
gene_coords[["WBGene00016022"]][2375:2488]



