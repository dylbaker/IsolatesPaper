### This script uses collection_tax_data.csv and logistic_mod_stats_normCoefs.csv to build figure 2

##### Libraries ####
library(tidyverse)
library(ggpubr)
library(tidytext)
library(patchwork)

#### Read in Data from Master ####
colors <- read.csv("./csv_files/colors.csv")|>
  mutate(host_species = str_to_lower(hostLong))|>
  select(host_species, color)

stats_pure_meanCoefs <- read.csv("./csv_files/collection_tax_data_pureOnly.csv") |>
  mutate(asv = str_replace(string = asv, pattern = "0+", replacement = ''),
         order = factor(parse_number(asv)),
         isolation_day = ifelse(str_detect(Isolate, "D31|DF"), "D31", "D3")) |>
  arrange(order, desc = TRUE) 

asv_presence <- stats_pure_meanCoefs |>
  group_by(asv) |>
  mutate(distinct_hosts = n_distinct(host_species),
         distinct_day = n_distinct(isolation_day))
#### Faceted Plots ####
grPlot_facet <- stats_pure_meanCoefs |>
  mutate(significant = ifelse(pGR_greater <= 0.05 |
                                pGR_less <= 0.05, T, F)) |>
  distinct(growthrate, .keep_all =TRUE) |>
  group_by(asv,host_species) |>
  mutate(distinct_day = as.character(n_distinct(isolation_day))) |>
  ggplot(aes(
    x = logNormGR,
    y = fct_reorder(Genus,desc(order)),
    shape = significant,
    color = distinct_day
  )) +
  geom_vline(aes(xintercept = 0), color = "grey", linetype = "longdash") +
  geom_point(size = 2, alpha = 0.7,stroke = 1) +
  scale_shape_manual(values = c(1, 16), "Significant Effect")+
  facet_grid(host_species ~ isolation_day, scales = "free_y") +
  #log of coculture growthrate/axenic growthrate
  labs(x = "Log Normalized Growth Rate",
       y = "Genus",
       title = "Host Growth Rate") +
  xlim(-0.9,0.9) +
  scale_color_manual(values = c("2" = "royalblue"), "ASVs Isolated at both Timepoints from the Same Host") +
  theme_pubclean() +
  theme(
    strip.background.y = element_blank(),
    strip.text.y = element_blank()
  )

grPlot_facet


ccPlot_facet <- stats_pure_meanCoefs |>
  mutate(significant = ifelse(pGR_greater <= 0.05 |
                                pGR_less <= 0.05, T, F)) |>
  distinct(growthrate, .keep_all =TRUE) |>
  group_by(asv,host_species) |>
  mutate(distinct_hosts = as.character(n_distinct(host_species)),
         distinct_day = as.character(n_distinct(isolation_day))) |>
  ggplot(aes(
    x = logNormCC,
    y = fct_reorder(Genus,desc(order)),
    shape = significant,
    color = distinct_day
  )) +
  geom_vline(aes(xintercept = 0), color = "grey", linetype = "longdash") +
  geom_point(size = 2, alpha = 0.7, stroke = 1) +
  scale_shape_manual(values = c(1, 16), "Significant Effect")+
  facet_grid(host_species ~ isolation_day, scales = "free_y") +
  #log of coculture growthrate/axenic growthrate
  labs(x = "Log Normalized Carrying Capacity",
       y = "Genus",
       title = "Host Carrying Capacity") +
  xlim(-3, 3) +
  scale_color_manual(values = c("2" = "royalblue"), "ASVs Isolated at both Timepoints from the Same Host") +
  # scale_color_manual(values = c("5" = "#D2042D",
  #                               "4" = "#E1AD01",
  #                               "3" = "green",
  #                               "2" = "royalblue",
  #                               "1" = "black"), "ASVs Shared by Hosts") +
  theme_pubclean()
ccPlot_facet


final_plot <-
  (grPlot_facet | (ccPlot_facet + lapply(c("y.text","ylab"),rremove))) +
  plot_annotation(tag_levels = "A") + 
  plot_layout(guides = 'collect') &
  theme(legend.position = "bottom",
        legend.justification = "center")

final_plot
png(filename = "./fig2/GR,CC Isolate Impacts on Host.png",
    res = 450,
    type = "cairo",
    units = "in",
    width = 14,
    height = 12)
final_plot
dev.off()