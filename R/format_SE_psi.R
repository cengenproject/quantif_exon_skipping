# Load, explore, format, and export the SE PSI data


library(tidyverse)
library(wbData)


# Load and format data ----

psi <- read_tsv("data/from_ruddle/2024-09-10_SE_PSI/assembled_psi.tsv",
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

stop_for_problems(psi)

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



# Filter out the gene-neuron pairs that have no expression ----

# load Alec's sc-bulk integration
gene_expression <- read.csv("data/thresholded_gene_expression/bsn12_subtracted_integrated_binarized_expression_withVDDD_FDR0.05_030424.tsv",
                            sep = "\t") |>
  rownames_to_column("gene_id") |>
  pivot_longer(-gene_id,
               names_to = "neuron_id",
               values_to = "expression")



# Clean up: if not present in gene_expression, almost surely not expressed in any neuron
psi_filt <- psi_long |>
  filter(neuron_id %in% gene_expression$neuron_id,
         gene_id %in% gene_expression$gene_id)


psi_joint <- left_join(psi_filt,
                      gene_expression,
                      by = c("gene_id", "neuron_id"))



ggplot(psi_joint) +
  # ggbeeswarm::geom_quasirandom(aes(x = as.character(expression), y = nb_reads)) +
  geom_violin(aes(x = as.character(expression), y = nb_reads)) +
  geom_boxplot(aes(x = as.character(expression), y = nb_reads), width = .2) +
  scale_y_log10()

psi_export <- psi_joint |>
  filter(expression == 1L)

# Remove outlier samples if not done previously
outlier_samples <- c("AVKr113","RICr133","PVMr122", "ASIr154")

psi_export <- psi_export |>
  filter(! sample_id %in% outlier_samples)

# Convert into gene coordinates ----

gene_coords <- wb_load_gene_coords(289)


psi_export <- psi_export |>
  left_join(gene_coords |>
              select(gene_id,
                     startGene = start, endGene = end),
            by = c("gene_id")) |>
  mutate(intron_start = if_else(strand == "+",
                                as.integer(startLong - startGene + 2),
                                as.integer(endGene - endLong + 2)),
         intron_end   = if_else(strand == "+",
                                as.integer(endLong - startGene),
                                as.integer(endGene - startLong)),
         exon_start = if_else(strand == "+",
                              as.integer(startExon - startGene + 1),
                              as.integer(endGene - endExon + 1)),
         exon_end   = if_else(strand == "+",
                              as.integer(endExon - startGene + 1),
                              as.integer(endGene - startExon + 1)),
         gene_length = as.integer(endGene - startGene +1))


psi_export |>
  select(event_id, intron_start, intron_end, exon_start, exon_end, gene_length, gene_id) |>
  distinct() |>
  write_tsv("data/export_for_arman/240910_events_coordinates.tsv")

psi_export |>
  select(event_id, sample_id, nb_reads, PSI) |>
  write_tsv("data/export_for_arman/240910_PSI_quantifications.tsv")


# check ----
evc <- read_tsv("data/export_for_arman/240910_events_coordinates.tsv")



stopifnot(all(evc$exon_end <= evc$gene_length))
stopifnot(all(evc$intron_end <= evc$gene_length))
stopifnot(all(evc$exon_start <= evc$gene_length))
stopifnot(all(evc$intron_start <= evc$gene_length))
stopifnot(all(evc$exon_end <= evc$intron_end))
stopifnot(all(evc$exon_start >= evc$intron_start))

evc |>
  filter(exon_end > intron_end)
psi_export |>
  filter(exon_end > intron_end) |> slice_head(n = 1) |> as.data.frame()




# Compare previous ----
library(tidyverse)

evc_old <- read_tsv("data/export_for_arman/240308_events_coordinates.tsv")
evc_new <- read_tsv("data/export_for_arman/240910_events_coordinates.tsv")

# non-CI identical
waldo::compare(evc_old |> filter(!startsWith(event_id, "CI")),
               evc_new |> filter(!startsWith(event_id, "CI")))

evc_old |> filter(startsWith(event_id, "CI"))
evc_new |> filter(startsWith(event_id, "CI"))

# event names identical
all.equal(evc_old$event_id |> sort(), evc_new$event_id |> sort())

# gene IDs identical (i.e. CI_n corresponds to the same gene)
full_join(
  evc_old |> filter(startsWith(event_id, "CI")) |> select(event_id, gene_id),
  evc_new |> filter(startsWith(event_id, "CI")) |> select(event_id, gene_id),
  by = "event_id"
  ) |>
  filter(gene_id.x != gene_id.y)


xx <- full_join(
  evc_old |> filter(startsWith(event_id, "CI")),
  evc_new |> filter(startsWith(event_id, "CI")),
  by = "event_id"
)


stable_CIs <- xx$event_id[xx$intron_start.x == xx$intron_start.y &
                            xx$intron_end.x == xx$intron_end.y &
                            xx$exon_start.x == xx$exon_start.y &
                            xx$exon_end.x == xx$exon_end.y &
                            xx$gene_length.x == xx$gene_length.y]
length(stable_CIs)

unstable_CIs <- setdiff(xx$event_id, stable_CIs)
length(unstable_CIs)

xx |>
  filter(event_id %in% unstable_CIs) |>
  mutate(stab_gene = gene_id.x == gene_id.y) |>
  pull(stab_gene) |>
  table()


psi_old <- read_tsv("data/export_for_arman/240308_PSI_quantifications.tsv")
psi_new <- read_tsv("data/export_for_arman/240910_PSI_quantifications.tsv")

waldo::compare(psi_old |> filter(!startsWith(event_id, "CI")),
               psi_new |> filter(!startsWith(event_id, "CI")))

psi_old |> filter(startsWith(event_id, "CI"))
psi_new |> filter(startsWith(event_id, "CI"))

all.equal(psi_old$event_id |> sort(), psi_new$event_id |> sort())


waldo::compare(
  psi_old |> filter(startsWith(event_id, "CI"), event_id %in% stable_CIs),
  psi_new |> filter(startsWith(event_id, "CI"), event_id %in% stable_CIs)
)


yy <- full_join(
  psi_old |> filter(startsWith(event_id, "CI"), event_id %in% unstable_CIs),
  psi_new |> filter(startsWith(event_id, "CI"), event_id %in% unstable_CIs),
  by = c("event_id", "sample_id")
)

ggplot(yy) +
  theme_classic() +
  scale_x_log10() + scale_y_log10() +
  coord_equal() +
  geom_point(aes(x = nb_reads.x, y = nb_reads.y),
             alpha = .02)


ggplot(yy) +
  theme_classic() +
  coord_equal() +
  geom_point(aes(x = PSI.x, y = PSI.y),
             alpha = .02)
ggplot(yy) +
  theme_classic() +
  coord_equal() +
  geom_hex(aes(x = PSI.x, y = PSI.y))


table(yy$PSI.x < .9, yy$PSI.y < .9)





