### This script uses collection_tax_data.csv and logistic_mod_stats_normCoefs.csv to build figure 2

##### Libraries ####
library(tidyverse,ggpubr)

#### Read in Data from Master ####
colors <- read.csv("./csv_files/colors.csv")|>
  mutate(host_species = str_to_lower(hostLong))|>
  select(host_species, color)

tax <- read.csv("./csv_files/collection_tax_data.csv") |>
  arrange( "Genus", desc = T)

stats_normCoefs <- read.csv("./csv_files/logistic_mod_stats_normCoefs.csv")|>
  mutate(Isolate = ifelse(str_detect(Isolate, "\\.") == T,
                          str_replace(Isolate, "\\.", ","), Isolate))|>
  # grubbs()|> (Unneeded because grubbs run in master)
  filter(outlier == F)

#### Faceted Plots ####
grPlot_facet <- stats_normCoefs |>
  mutate(significant = ifelse(pGR_greater <= 0.05| pGR_less <= 0.05, T, F))|>
  drop_na()|>
  left_join(colors)|>
  left_join(tax)|>
  ggplot(aes(y = Genus, x = logNormGR, color = host_species, shape = significant))+ # genus
  # ggplot(aes(y = Family, x = logNormGR, color = host_species, shape = significant))+ # fam
  geom_vline(aes(xintercept = 0), color = "grey", linetype = "longdash")+
  geom_point(size = 2, alpha = 0.5, position = "jitter")+
  scale_shape_manual(values = c(1, 16))+
  xlim(-4, 4)+
  labs(x = "log(Coculture Growth Rate / Axenic Growth Rate)",
       y = "Bacterial Genus",
       title = "Genus Level Impact on Algal Growth Rate")+ ## Genus
       # y = "Bacterial Family",
       # title = "Family Level Impact on Algal Growth Rate")+ ## Family
  facet_grid(vars(host_species), scales = "free_y")+
  theme_pubclean() +
  scale_color_manual(values = colors$color)
grPlot_facet

ggsave("./fig2/growthRate_plot_facet.pdf", grPlot_facet, device = "pdf", height = 10, width = 8)

ccPlot_facet <- stats_normCoefs |>
  mutate(significant = ifelse(pCC_greater <= 0.05 | pCC_less <= 0.05, T, F))|>
  drop_na()|>
  left_join( tax)|>
  ggplot( aes(y = Genus, x = logNormCC, color = host_species, shape = significant))+ # genus
  # ggplot( aes(y = Family, x = logNormCC, color = host_species, shape = significant))+ # fam
  geom_vline(aes(xintercept = 0), color = "grey", linetype = "longdash")+
  geom_point(size = 2, alpha = 0.5, position = "jitter")+
  scale_shape_manual(values = c(1, 16))+
  xlim(-4, 4)+
  labs(x = "log(Coculture Carrying Capacity / Axenic Carrying Capacity)",
       y = "Bacterial Genus",
       title = "Genus Level Impact on Algal Carrying Capacity")+ ## genus
       # y = "Bacterial Family",
       # title = "Family Level Impact on Algal Carrying Capacity")+ ## family
  facet_grid(vars(host_species), scales = "free_y")+
  theme_pubclean() +
  scale_color_manual(values = colors$color)
ccPlot_facet

ggsave("./fig2/carryingCapacity_plot_facet.pdf", ccPlot_facet, device = "pdf", height = 10, width = 8)