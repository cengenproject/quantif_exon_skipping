# Explore the exported data


library(tidyverse)
library(wbData)


# Load and format data ----

psi <- read_tsv("data/export_for_arman/230531_PSI_quantifications.tsv",
                col_types = cols(
                  event_id = col_character(),
                  sample_id = col_character(),
                  nb_reads = col_integer(),
                  PSI = col_double()
                ))

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

psi_long$event_id |> unique() |> str_split_i("_", 1) |> table()

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
potential_DSE |>
  filter(Dpsi > .99) |> pull(event_id) |> unique() |> str_split_i("_",1) |> table()


# Look for genes with multiple events
psi_long |>
  filter(nb_reads > 10) |>
  select(gene_id, event_id) |>
  distinct() |>
  count(gene_id) |>
  arrange(desc(n))

# ex: mel-11  WBGene00003196   II:9,357,531..9,369,596



# Compare exons properties ----
coords <- read_tsv("data/export_for_arman/230531_events_coordinates.tsv")

eulerr::euler(list(coords = unique(coords$event_id), psi = unique(psi$event_id))) |> plot()
table(unique(coords$event_id) %in% unique(psi_long$event_id))
table(unique(psi_long$event_id) %in% unique(coords$event_id))
which(! unique(psi_long$event_id) %in% unique(coords$event_id)) |> head()

unique(psi_long$event_id)[4]

psi$event_id |> unique() |> str_split_i("_", 1) |> table()




psi |> filter(event_id == "CE_659") |> pull(IR_raw_AVLr232)
[1] 3
psi |> filter(event_id == "CE_659") |> pull(ER_raw_AVLr232)
[1] 0