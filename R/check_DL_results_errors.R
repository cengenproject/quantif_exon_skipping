# Check quality of prediction on different events/neurons

# inits ----
library(tidyverse)
library(wbData)
gids <- wb_load_gene_ids(281)


#~ load data ----
res_by_event <- read_csv("data/res_from_egbert/error_by_eid_2023-10-16-20-24-30_8E027.csv")

res_by_neuron <- read_csv("data/res_from_egbert/error_by_nid_2023-10-16-20-24-30_8E027.csv")

res_all <- read_csv("data/res_from_egbert/error_all_2023-10-16-20-24-30_8E027.csv") |>
  select(-1)

event_coords <- read_tsv("data/export_for_arman/230912_events_coordinates.tsv") |>
  mutate(gene_name = i2s(gene_id, gids),
         coordinates = paste0(intron_start, "-" , intron_end)) |>
  filter(startsWith(event_id, "SE_"))

psi <- read_tsv("data/export_for_arman/230912_PSI_quantifications.tsv") |>
  filter(startsWith(event_id, "SE_"))


# reproduce Egbert's figures ----
res_by_event |>
  arrange(desc(error)) |>
  mutate(eid = fct_inorder(eid)) |>
  head(15) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = eid, y = error)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("event_id") +
  ggtitle("")

res_by_event |>
  arrange(desc(error)) |>
  mutate(eid = fct_inorder(eid)) |>
  tail(15) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = eid, y = error)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



# Check best/worst events

 |>
  head(5)


worst <- res_by_event |>
  select(-1) |>
  left_join(event_coords,
            by = c("eid" = "event_id")) |>
  arrange(desc(error)) |>
  mutate(eid = fct_inorder(eid)) |>
  head(10)

best <- res_by_event |>
  select(-1) |>
  left_join(event_coords,
            by = c("eid" = "event_id")) |>
  arrange(desc(error)) |>
  mutate(eid = fct_inorder(eid)) |>
  tail(10)

View(worst)

psi |>
  filter(event_id %in% worst$eid) |>
  left_join(res_by_event, by = c(event_id = "eid")) |>
  mutate(event_id2 = paste(event_id, ": ", round(error, 2))) |>
  ggplot() +
  theme_classic() +
  facet_wrap(~event_id2, ncol = 5) +
  geom_point(aes(x = nb_reads, y = PSI)) +
  ggtitle("Worst performers") +
  scale_x_log10()


psi |>
  filter(event_id %in% best$eid) |>
  left_join(res_by_event, by = c(event_id = "eid")) |>
  mutate(event_id2 = paste(event_id, ": ", round(error, 2))) |>
  ggplot() +
  theme_classic() +
  facet_wrap(~event_id2, ncol = 5) +
  geom_point(aes(x = nb_reads, y = PSI)) +
  ggtitle("Best performers") +
  scale_x_log10()

psi |>
  filter(event_id %in% sample(unique(psi$event_id), 10)) |>
  left_join(res_by_event, by = c(event_id = "eid")) |>
  mutate(event_id2 = paste(event_id, ": ", round(error, 2))) |>
  ggplot() +
  theme_classic() +
  facet_wrap(~event_id2, ncol = 5) +
  geom_point(aes(x = nb_reads, y = PSI)) +
  ggtitle("Random 10")

# check amount of variability impact on error
psi |>
  group_by(event_id) |>
  summarize(var = var(PSI, na.rm = TRUE),
            nb_points = sum(!is.na(PSI))) |>
  right_join(res_by_event,
             by = c(event_id = "eid")) |>
  mutate(`error vigintile` = cut(error, breaks = quantile(error, probs = c(0,.05,.95,1)), labels = c("bottom 5%", "mid", "top 5%"))) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = var, y = error, color = `error vigintile`)) +
  xlab("Variance of PSI")

psi |>
  group_by(event_id) |>
  summarize(var = var(PSI, na.rm = TRUE),
            nb_points = sum(!is.na(PSI))) |>
  right_join(res_by_event,
             by = c(event_id = "eid")) |>
  mutate(decile = cut(error, breaks = quantile(error, probs = c(0,.1,.9,1)), labels = c("bottom_decile", "mid", "top_decile"))) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = nb_points, y = error, color = decile))


psi |>
  group_by(event_id) |>
  summarize(var = var(PSI, na.rm = TRUE),
            nb_points = sum(!is.na(PSI))) |>
  filter(event_id %in% worst$eid | event_id %in% best$eid) |>
  mutate(event_id = factor(event_id, levels = levels(worst$eid))) |>
  arrange(event_id)



# by event ----
psi_event <- psi |>
  right_join(res_by_event |> select(-c(`...1`,eid_indx)),
             by = c(event_id = "eid"))
  

# Check 
psi_event |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = nb_reads, y = error), alpha = .2) +
  scale_x_log10()


psi_event |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = PSI, y = error), alpha = .2)



# Mean variance
psi_event |>
  group_by(event_id) |>
  summarize(var = var(PSI, na.rm = TRUE),
            error = mean(error),
            nb_points = sum(!is.na(PSI))) |>
  mutate(decile = cut(error, breaks = quantile(error, probs = c(0,.1,.9,1)), labels = c("bottom_decile", "mid", "top_decile"))) |>
  ggplot() +
  theme_classic() +
  geom_point(aes(x = nb_points, y = error, color = decile))



