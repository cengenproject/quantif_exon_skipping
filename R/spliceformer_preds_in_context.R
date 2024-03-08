# To put model prediction into context, check PSI vs expression of high attention weight SFs

# Inits ----
library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(289)

psi <- read_tsv("data/export_for_arman/230912_PSI_quantifications.tsv") |>
  mutate(type = str_match(event_id, "^(CI|SE|CE)_[0-9]+$")[,2],
         type = factor(type, levels = c("SE","CI","CE")),
         neuron_id = str_match(sample_id, "^([A-Z0-9]{2,4})r[0-9]{1,3}$")[,2])

table(psi$type, useNA = 'ifany')
table(psi$neuron_id, useNA = 'ifany')

sf <- read_tsv("data/export_for_arman/tx_expression.tsv.gz")



# cdgs-1 ----
vals_cdgs <- bind_rows(
  psi |>
    filter(event_id == "SE_173") |>
    group_by(event_id, neuron_id) |>
    summarize(PSI = mean(PSI, na.rm = TRUE),
              .groups = 'drop') |>
    select(neuron = neuron_id,
           id = event_id,
           value = PSI) |>
    arrange(value),
  sf |>
    filter(gene_id %in% c("WBGene00022060","WBGene00005105","WBGene00021332","WBGene00000862",
                          "WBGene00000796","WBGene00017004","WBGene00015233")) |>
    group_by(neuron_id, gene_id, sample_id) |>
    summarize(TPM_gene = sum(TPM, na.rm = TRUE),
              .groups = 'drop') |>
    group_by(neuron_id, gene_id) |>
    summarize(TPM = log10(mean(TPM_gene, na.rm = TRUE)),
              .groups = 'drop') |>
    mutate(gene_id = i2s(gene_id, gids)) |>
    select(neuron = neuron_id,
           id = gene_id,
           value = TPM)
)

vals_cdgs |>
  mutate(neuron = fct_inorder(neuron)) |>
  ggplot() +
  theme_classic() +
  geom_tile(aes(y = id, x = neuron, fill = value)) +
  viridis::scale_fill_viridis() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))






# unc-68 ----
vals_unc68 <- bind_rows(
  psi |>
    filter(event_id == "SE_850") |>
    group_by(event_id, neuron_id) |>
    summarize(PSI = mean(PSI, na.rm = TRUE),
              .groups = 'drop') |>
    select(neuron = neuron_id,
           id = event_id,
           value = PSI) |>
    arrange(value),
  sf |>
    filter(gene_id %in% c("WBGene00015233","WBGene00021332","WBGene00003082","WBGene00013964","WBGene00005105")) |>
    group_by(neuron_id, gene_id, sample_id) |>
    summarize(TPM_gene = sum(TPM, na.rm = TRUE),
              .groups = 'drop') |>
    group_by(neuron_id, gene_id) |>
    summarize(TPM = log10(mean(TPM_gene, na.rm = TRUE)),
              .groups = 'drop') |>
    mutate(gene_id = i2s(gene_id, gids)) |>
    select(neuron = neuron_id,
           id = gene_id,
           value = TPM)
)

vals_unc68 |>
  mutate(neuron = fct_inorder(neuron)) |>
  ggplot() +
  theme_classic() +
  geom_tile(aes(y = id, x = neuron, fill = value)) +
  viridis::scale_fill_viridis() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




# madd-2 ----
vals_madd2 <- bind_rows(
  psi |>
    filter(event_id == "SE_819") |>
    group_by(event_id, neuron_id) |>
    summarize(PSI = mean(PSI, na.rm = TRUE),
              .groups = 'drop') |>
    select(neuron = neuron_id,
           id = event_id,
           value = PSI) |>
    arrange(value),
  sf |>
    filter(gene_id %in% c("WBGene00015233","WBGene00021332","WBGene00013964","WBGene00003082",
                          "WBGene00044513","WBGene00005105","WBGene00000890")) |>
    group_by(neuron_id, gene_id, sample_id) |>
    summarize(TPM_gene = sum(TPM, na.rm = TRUE),
              .groups = 'drop') |>
    group_by(neuron_id, gene_id) |>
    summarize(TPM = log10(mean(TPM_gene, na.rm = TRUE)),
              .groups = 'drop') |>
    mutate(gene_id = i2s(gene_id, gids)) |>
    select(neuron = neuron_id,
           id = gene_id,
           value = TPM)
)

vals_madd2 |>
  mutate(neuron = fct_inorder(neuron)) |>
  ggplot() +
  theme_classic() +
  geom_tile(aes(y = id, x = neuron, fill = value)) +
  viridis::scale_fill_viridis() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# nrx-1 ----
vals_nrx <- bind_rows(
  psi |>
    filter(event_id == "SE_773") |>
    group_by(event_id, neuron_id) |>
    summarize(PSI = mean(PSI, na.rm = TRUE),
              .groups = 'drop') |>
    select(neuron = neuron_id,
           id = event_id,
           value = PSI) |>
    arrange(value),
  sf |>
    filter(gene_id %in% c("WBGene00013964","WBGene00015233","WBGene00021332","WBGene00003082",
                          "WBGene00044513","WBGene00000890","WBGene00005105","WBGene00018156")) |>
    group_by(neuron_id, gene_id, sample_id) |>
    summarize(TPM_gene = sum(TPM, na.rm = TRUE),
              .groups = 'drop') |>
    group_by(neuron_id, gene_id) |>
    summarize(TPM = log10(mean(TPM_gene, na.rm = TRUE)),
              .groups = 'drop') |>
    mutate(gene_id = i2s(gene_id, gids)) |>
    select(neuron = neuron_id,
           id = gene_id,
           value = TPM)
)

vals_nrx |>
  mutate(neuron = fct_inorder(neuron)) |>
  ggplot() +
  theme_classic() +
  geom_tile(aes(y = id, x = neuron, fill = value)) +
  viridis::scale_fill_viridis() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))






