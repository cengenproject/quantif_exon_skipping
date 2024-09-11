# Very quick exploration of the PSI values (note, deeper analyses in other projects)

# Inits ----
library(tidyverse)

psi <- read_tsv("data/export_for_arman/240910_PSI_quantifications.tsv") |>
  mutate(type = str_match(event_id, "^(CI|SE|CE)_[0-9]+$")[,2],
         type = factor(type, levels = c("SE","CI","CE")))

coords <- read_tsv("data/export_for_arman/240910_events_coordinates.tsv")





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


# only high-quality events
psi |>
  filter(nb_reads > 10) |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = PSI), color = 'white') +
  facet_wrap(~type)




#~ Check sequence at boundaries ----

gseq <- Biostrings::readDNAStringSet("data/intermediates_for_DL/240906_all_gene_sequences.fa.gz")

# 0         intron_start       exon_start     exon_end         intron_end   gene_length
# |              |                    |         |               |            |
# ---------------GT-----------------AG-----------GT------------AG-------------

correct_boundary_seq <- logical(length = nrow(coords))
for(i in 1:nrow(coords)){
  seq <- gseq[[ coords$gene_id[[i]] ]]
  
  irg <- IRanges::IRanges(
    start = c(
      coords$intron_start[[i]],
      coords$intron_end[[i]] - 1L,
      coords$exon_start[[i]] - 2L,
      coords$exon_end[[i]] + 1L
    ),
    end = c(
      coords$intron_start[[i]] + 1L,
      coords$intron_end[[i]],
      coords$exon_start[[i]] - 1L,
      coords$exon_end[[i]] + 2L
    )
  )
  
  correct_boundary_seq[[i]] <- identical(
    Biostrings::extractAt(seq, irg) |> as.character(),
    c("GT","AG","AG","GT")
  )
}

table(correct_boundary_seq)
# correct_boundary_seq
# FALSE  TRUE 
#    21  2279 

#> manual checking SE_640, CE_692 (and CI_664 same place), SE_975: non-canonical sites in the annotation.






# Compare with older versions ----

#~ 240308



psi_old <- read_tsv("data/export_for_arman/240308_PSI_quantifications.tsv") |>
  mutate(type = str_match(event_id, "^(CI|SE|CE)_[0-9]+$")[,2],
         type = factor(type, levels = c("SE","CI","CE")))
psi_new <- psi

old_ev_coords <- read_tsv("data/export_for_arman/240308_events_coordinates.tsv")

new_old_correspondance <- right_join(coords, old_ev_coords |> rename(old_id = event_id),
                                     by = c("intron_start", "intron_end", "exon_start", "exon_end",
                                            "gene_length", "gene_id"),
                                     relationship = "many-to-many") |>
  filter(!is.na(event_id)) |>
  select(new_id = event_id,
         old_id)

psi_old <- psi_old |>
  left_join(new_old_correspondance,
            by = c(event_id = "old_id"),
            relationship = "many-to-many") |>
  select(-event_id) |>
  rename(event_id = new_id)




#~ 230912 ----

psi_old <- read_tsv("data/export_for_arman/230912_PSI_quantifications.tsv") |>
  mutate(type = str_match(event_id, "^(CI|SE|CE)_[0-9]+$")[,2],
         type = factor(type, levels = c("SE","CI","CE")))
psi_new <- psi

old_ev_coords <- read_tsv("data/export_for_arman/230912_events_coordinates.tsv")

new_old_correspondance <- right_join(coords, old_ev_coords |> rename(old_id = event_id),
          by = c("intron_start", "intron_end", "exon_start", "exon_end",
                 "gene_length", "gene_id"),
          relationship = "many-to-many") |>
  filter(!is.na(event_id)) |>
  select(new_id = event_id,
         old_id)

psi_old <- psi_old |>
  left_join(new_old_correspondance,
            by = c(event_id = "old_id"),
            relationship = "many-to-many") |>
  select(-event_id) |>
  rename(event_id = new_id)


left_join(
  psi_old |> filter(type == "CI"),
  psi_new |> filter(type == "CI"),
  by = c("event_id", "sample_id", "type")
) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = PSI.x, y = PSI.y), alpha = .1)

left_join(
  psi_old |> filter(type == "CE"),
  psi_new |> filter(type == "CE"),
  by = c("event_id", "sample_id", "type")
) |> filter(PSI.x <.01, PSI.y == 1)






#~ 230907

psi_old <- read_tsv("data/export_for_arman/230907_PSI_quantifications.tsv") |>
  mutate(type = str_match(event_id, "^(CI|SE|CE)_[0-9]+$")[,2],
         type = factor(type, levels = c("SE","CI","CE")))
psi_new <- psi

# only CE changed
table(psi_old$type, useNA = 'ifany')
table(psi_new$type, useNA = 'ifany')

all.equal(
  psi_old |> filter(type != "CE"),
  psi_new |> filter(type != "CE")
)


ggplot(psi_old) +
  theme_classic() +
  geom_histogram(aes(x = PSI), color = 'white') +
  facet_wrap(~type)

ggplot(psi_new) +
  theme_classic() +
  geom_histogram(aes(x = PSI), color = 'white') +
  facet_wrap(~type)

bind_rows(
  psi_old |> filter(type == "CE") |> add_column(version = "old"),
  psi_new |> filter(type == "CE") |> add_column(version = "new")
) |>
  ggplot() +
  theme_classic() +
  geom_density(aes(x = PSI, fill = version), alpha = .2)



# the cases that shouldn't happen too much
psi_new |> filter(type == "CE", PSI > .9, nb_reads > 10)
psi_old |> filter(type == "CE", PSI > .9, nb_reads > 10)

psi_old |> filter(type == "CE", PSI > .9, nb_reads > 10) |>
  count(event_id)

psi_new |> filter(type == "CE", PSI > .9, nb_reads > 10) |>
  count(event_id)




#~ 230531 ----

psi_old <- read_tsv("data/export_for_arman/230531_PSI_quantifications.tsv") |>
  mutate(type = str_match(event_id, "^(CI|SE|CE)_[0-9]+$")[,2],
         type = factor(type, levels = c("SE","CI","CE")))
psi_new <- psi

# only CE changed
table(psi_old$type, useNA = 'ifany')
table(psi_new$type, useNA = 'ifany')

all.equal(
  psi_old |> filter(type != "CE"),
  psi_new |> filter(type != "CE")
)


ggplot(psi_old) +
  theme_classic() +
  geom_histogram(aes(x = PSI), color = 'white') +
  facet_wrap(~type)

ggplot(psi_new) +
  theme_classic() +
  geom_histogram(aes(x = PSI), color = 'white') +
  facet_wrap(~type)

bind_rows(
  psi_old |> filter(type == "CE") |> add_column(version = "old"),
  psi_new |> filter(type == "CE") |> add_column(version = "new")
) |>
  ggplot() +
  theme_classic() +
  geom_density(aes(x = PSI, fill = version), alpha = .2)

# the cases that shouldn't happen too much
psi_new |> filter(type == "CE", PSI > .9, nb_reads > 10)
psi_old |> filter(type == "CE", PSI > .9, nb_reads > 10)

psi_old |> filter(type == "CE", PSI > .9, nb_reads > 10) |>
  count(event_id)

psi_new |> filter(type == "CE", PSI > .9, nb_reads > 10) |>
  count(event_id)




#~ 230228 ----

#  Compare with previous version (without CE), to make sure rerunning everything didn't introduce new problems

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



# frm-5.2
coords |> filter(gene_id == "WBGene00021406")


# WBGene00016022
coords |> filter(gene_id == "WBGene00016022")


gene_coords <- Biostrings::readDNAStringSet("../stringtie_quantif/data/intermediates_for_DL/220920_gene_sequences.fa.gz")

gene_coords[["WBGene00016022"]]

# SE_946 intron
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


# Check weird cases
psi |> filter(type == "CE", PSI > .9, nb_reads > 10) |>
  filter(event_id != "CE_657")


# e.g. CE_657 in CEPr278
coords |> filter(event_id == "CE_657")
#> due to gene on other strand inside the intron

# e.g. CE_677 === CE_678 in I5r210
coords |> filter(event_id == "CE_677")
#> due to noise (there is a weird bump around there,
#> though gene likely not expressed, however there is an expressed ncRNA in its 3'UTR
#> that may appear like expression in scRNA-Seq (though doesn't seem polyA, so unsure how))


# e.g. CE_695 in ADLr172
coords |> filter(event_id == "CE_695")
#> some noise in an unexpressed region of the gene (so there is no exclusion read)
#> it's an unexpressed region because of alt first exon: these neurons start the gene downstream








# older versions below (230912)

# e.g. CE_980 in ADLr171
coords |> filter(event_id == "CE_980")
#> due to alt 5' ss unannotated but highly used in ADL

# e.g. CE_671 in LUAr
coords |> filter(event_id == "CE_671")
#> gene not expressed, so no excl reads,
#> but noisy for some reason, so some incl reads

# e.g. CE_55
coords |> filter(event_id == "CE_55")
#> the fake exon ends up on top of a highly expressed snoRNA


# CE_1010
coords |> filter(event_id == "CE_1010")

# CE_649, 597, 669
coords |> filter(event_id == "CE_649")
coords |> filter(event_id == "CE_597")
coords |> filter(event_id == "CE_669")



# Exons properties ----
eulerr::euler(list(coords = unique(coords$event_id), psi = unique(psi$event_id))) |> plot()


#~ length ----




coords_annot <- coords |>
  mutate(event_type = str_match(event_id, "^(CI|SE|CE)_[0-9]+$")[,2],
         exon_length = exon_end - exon_start + 1,
         upstream_intron_length = intron_end - exon_end + 1,
         downstream_intron_length = exon_start - intron_start + 1) |>
  left_join(psi |>
              filter(nb_reads > 10) |>
              group_by(event_id) |>
              summarize(dPSI = diff(range(PSI, na.rm = TRUE))),
            by = "event_id")

ggplot(coords_annot) +
  theme_classic() +
  facet_wrap(~event_type) +
  geom_histogram(aes(x = dPSI))

ggplot(coords_annot) +
  theme_classic() +
  geom_freqpoly(aes(x = exon_length, color = event_type)) +
  scale_x_log10()

