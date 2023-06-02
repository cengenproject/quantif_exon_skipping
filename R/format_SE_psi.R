# Load, explore, format, and export the SE PSI data


library(tidyverse)
library(wbData)


# Load and format data ----

psi <- read_tsv("data/from_ruddle/2023-05-31_SE_PSI/assembled_psi.tsv",
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
gene_expression <- read.csv("data/thresholded_gene_expression/bsn9_subtracted_integrated_binarized_expression_withVDDD_FDR_0.1_092022.tsv",
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
outlier_samples <- c("AVKr113","RICr133","PVMr122")

psi_export <- psi_export |>
  filter(! sample_id %in% outlier_samples)

# Convert into gene coordinates ----

gene_coords <- wb_load_gene_coords(281)


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
  write_tsv("data/export_for_arman/230531_events_coordinates.tsv")

psi_export |>
  select(event_id, sample_id, nb_reads, PSI) |>
  write_tsv("data/export_for_arman/230531_PSI_quantifications.tsv")

  








