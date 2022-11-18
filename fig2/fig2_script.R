### This script uses collection_tax_data.csv and logistic_mod_stats_normCoefs.csv to build figure 2

##### Libraries ####
library(tidyverse)
library(ggpubr)
library(patchwork)

#### Read in Data from Master ####
colors <- read.csv("./csv_files/colors.csv")|>
  mutate(host_species = str_to_lower(hostLong))|>
  select(host_species, color)

contaminates <- read.csv("./fig1/Isolate_contamination_report.csv") |>
  mutate(asv = str_replace(string = asv, pattern = "0+", replacement = '')) |>
  filter(phylum == "Bacteroidota" | phylum == "Firmicutes")
 
stats_meanCoefs <- read.csv("./csv_files/collection_tax_data.csv") |>
  mutate(asv = str_replace(string = asv, pattern = "0+", replacement = ''),
         order = factor(parse_number(asv)),
         isolation_day = ifelse(str_detect(Isolate, "D31|DF|S"), "D31", "D3")) |>
  distinct(exact_isolate, .keep_all = T) |>
  distinct(Isolate, .keep_all = T)|>
  anti_join(contaminates, by = "asv") |>
  arrange(order, desc = TRUE) 

stats_pure_meanCoefs <- filter(stats_meanCoefs, mixed == F)

asv_presence <- stats_pure_meanCoefs |>
  group_by(asv) |>
  mutate(distinct_hosts = n_distinct(host_species),
         distinct_day = n_distinct(isolation_day))

summary_data <- stats_meanCoefs |>
  group_by(host_species,isolation_day) |>
  summarise(n_pos_gr = sum(grEffect == "Increased GR")/n(),
            n_0_gr = sum(grEffect == "No Significant GR Change")/n(),
            n_neg_gr = sum(grEffect == "Decreased GR")/n(),
            n_pos_cc = sum(ccEffect == "Increased CC")/n(),
            n_0_cc = sum(ccEffect == "No Significant CC Change")/n(),
            n_neg_cc = sum(ccEffect == "Decreased CC")/n(),
            n_pos = sum(c(n_pos_gr,n_pos_cc)),
            n_neg = sum(c(n_neg_gr, n_neg_cc)),
            n_0 = sum(c(n_0_gr,n_0_cc)),
            total = n()
            ) |>
  mutate(host_species = str_to_title(host_species)) |>
  pivot_longer(!host_species &! isolation_day & !total, names_to = "growth_outcome", values_to = "percent") |>
  mutate(number = round(percent * total),
         label = paste("n =", number),
         metric = case_when(str_detect(growth_outcome, "gr") ~ "Growth Rate",
                            str_detect(growth_outcome, "cc") ~ "Carrying Capacity",
                            T ~ "Total"),
         growth_outcome = case_when(str_detect(growth_outcome, "n_pos") ~ "Positive",
                                    str_detect(growth_outcome, "n_neg") ~ "Negative",
                                    str_detect(growth_outcome, "n_0") ~ "Neutral"))

impact_barplot <-  summary_data |>
  filter(metric == "Total") |>
  ggplot(aes(x = growth_outcome, y = percent, fill = growth_outcome)) +
  geom_col(position = "dodge", show.legend = F) +
  facet_grid(host_species ~ isolation_day) +
  scale_x_discrete("Growth Outcome", 
                   limits = c("Negative", "Neutral" , "Positive")) +
  geom_text(aes(label = label),  vjust = "inward") +
  ylab("Percent of Isolates") +
  scale_fill_manual(values = c("Neutral" = "gray", "Positive" = "springgreen", "Negative" = "#D2042D")) +
  theme_pubclean() 

impact_barplot_sep <-  summary_data |>
  filter(metric %in% c("Growth Rate", "Carrying Capacity")) |>
  ggplot(aes(x = growth_outcome, y = percent, fill = growth_outcome)) +
  geom_col(position = "dodge", show.legend = F) +
  facet_grid(host_species ~ metric + isolation_day) +
  scale_x_discrete("Growth Outcome", 
                   limits = c("Negative","Neutral","Positive")) +
  geom_text(aes(label = label),  vjust = "inward") +
  ylab("Percent of Isolates") +
  scale_fill_manual(values = c("Neutral" = "gray", "Positive" = "springgreen", "Negative" = "#D2042D")) +
  theme_pubclean() 

png(filename = "./fig2/isolate_impact_barplot_sep.png",
    res = 450,
    type = "cairo",
    units = "in",
    width = 8,
    height = 8)
impact_barplot_sep
dev.off()

png(filename = "./fig2/isolate_impact_barplot.png",
    res = 450,
    type = "cairo",
    units = "in",
    width = 6,
    height = 8)
impact_barplot
dev.off()

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