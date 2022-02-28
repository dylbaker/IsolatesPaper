#### Libraries ####
library(tidyverse)
library(data.table)
library(patchwork)
library(ggpubr)

#### Directory Set ####
setwd("C:/Users/jlaue/OneDrive/Documents/Denef_Lab_Stuff/IsolatesPaper/")
getwd()

#### Colors ####
colors <- read.csv("./csv_files/colors.csv")%>%
  mutate(host_species = str_to_lower(hostLong))%>%
  select(host_species, color)

#### GRBC Data ####
collectionTax <- read.csv("./csv_files/collection_tax_data_pureOnly.csv")%>%
  mutate(Isolate = str_replace_all(string=Isolate,
                                   pattern="D",
                                   replacement="_D"),
         Isolate = str_replace_all(string = Isolate,
                                   pattern = "\\,",
                                   replacement = "\\."))%>%
  select(otu, readType, Isolate, host_species, count,
         Domain, Phylum, Class, Order, Family, Genus)%>%
  mutate(Isolate = gsub("DF", "D31", Isolate))%>%
  left_join(., colors)%>%
  group_by(Isolate, readType)%>%
  mutate(time = ifelse(str_detect(Isolate, "_D31"), "D31", "D3"),
         Genus=str_replace_all(string=Genus,
                               pattern="Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
                               replacement="A-N-P-R"),
         rel_abund = count/sum(count))%>%
  ungroup()

#### Phycosphere Community Data ####
phycosphereTax <- read.csv("./csv_files/phycosphere_tax_data.csv")%>%
  mutate(Isolate = str_replace_all(string=Isolate,
                                   pattern="D",
                                   replacement="_D"),
         Isolate = str_replace_all(string = Isolate,
                                   pattern = "\\,",
                                   replacement = "\\."))%>%
  select(otu, Isolate, host_species, count, Domain, Phylum, Class, Order, Family, Genus)%>%
  mutate(Isolate = gsub("DF", "D31", Isolate))%>%
  group_by(Isolate, Genus)%>%
  mutate(genus_count = sum(count))%>%
  ungroup()%>%
  group_by(Isolate)%>%
  mutate(time = ifelse(str_detect(Isolate, "_D31"), "D31", "D3"),
         Genus=str_replace_all(string=Genus,
                               pattern="Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
                               replacement="A-N-P-R"),
         rel_abund = count/sum(count),
         genus_rel_abund = genus_count/sum(count))%>%
  ungroup()

#### Pond Data (Natural Community) ####
pondTax <- read.csv("./csv_files/natCom_tax_data.csv")%>%
  mutate(Isolate = str_replace_all(string=Isolate,
                                   pattern="D",
                                   replacement="_D"),
         Isolate = str_replace_all(string = Isolate,
                                   pattern = "\\,",
                                   replacement = "\\."))%>%
  select(otu, Isolate, host_species, count, Domain, Phylum, Class, Order, Family, Genus)%>%
  mutate(Isolate = gsub("DF", "D31", Isolate))%>%
  group_by(Isolate)%>%
  mutate(time = ifelse(str_detect(Isolate, "_D31"), "D31", "D3"),
         Genus=str_replace_all(string=Genus,
                               pattern="Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
                               replacement="A-N-P-R"),
         rel_abund = count/sum(count))%>%
  ungroup()

#### Genus Level Data Organization ####
isolatedOTUs <- select(collectionTax, otu, host_species, color, time)%>%
  mutate(isolated = T)%>%
  distinct()

genus_data_D3 <- phycosphereTax %>%
  filter(time == "D3")%>%
  group_by(otu, host_species)%>%
  mutate(agg_rel_abund = sum(rel_abund),
         avg_rel_abund = mean(rel_abund))%>%
  ungroup()%>%
  left_join(., isolatedOTUs)%>%
  mutate(color = if_else(is.na(color), "black", color),
         isolated = if_else(is.na(isolated), F, isolated))

genus_data_D31 <- phycosphereTax %>%
  filter(time == "D31")%>%
  group_by(otu, host_species)%>%
  mutate(agg_rel_abund = sum(rel_abund),
         avg_rel_abund = mean(rel_abund))%>%
  ungroup()%>%
  left_join(., isolatedOTUs)%>%
  mutate(color = if_else(is.na(color), "black", color),
         isolated = if_else(is.na(isolated), F, isolated))

#### Plotting Fxn ####
abundance_plot_isolate <- function(species_name, data) { #, color_data, Genus_species
  df <- filter(data, host_species == species_name)%>%
    select(otu, avg_rel_abund, isolated, color)%>%
    # select(Genus, avg_rel_abund, isolated, color)%>%
    distinct()%>%
    arrange(desc(avg_rel_abund)) %>% slice(1:30) # limit curve to 30 genera
  plot <- df %>%
    #Reorder Genus by average relative abundance to get a smooth curve
    ggplot(aes(x = reorder(otu, -avg_rel_abund), y = avg_rel_abund, fill = isolated)) +
    # ggplot(aes(x = reorder(Genus, -avg_rel_abund), y = avg_rel_abund, fill = isolated)) +
    geom_bar(stat = "Identity") +
    theme_minimal() +
    scale_fill_manual(values = unique(df$color),
                      name = '', labels = c("Not Isolated", "Cultured Isolate"))+
    # scale_fill_manual(values = c(colors[(colors$host_species==species_name),]$color, "black"),
    #                   name =  '', labels = c("Cultured Isolate", "Not Isolated")) +
    # guides(fill = guide_legend(
    #   title = substitute(paste(italic(Genus_species))),
    #   title.position = "top",
    #   nrow = 2, byrow = TRUE
    # )) +
    # labs(y = "Relative Abundance (Genus)",
    #      x = "",
    #      fill = '') +
    # # , title = "") +
    # scale_y_continuous("Average Relative Abundance", 
    #                    limits = c(0, 0.32), 
    #                    breaks = c(0.0, 0.1, 0.2, 0.3)) +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(angle = 70,
                                   hjust = 1,
                                   size = 14),
                                   # color = color_data$color),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom",
        legend.box.just = "center",
        legend.direction = "vertical",
        legend.title = element_text(vjust = .5, size = 12),
        legend.text = element_text(vjust = .5, size = 12))
  # , plot.title = element_text(hjust = .5, size = 16)) +
  # scale_x_reordered()
  
  return(plot)
}
#### D3 Plots ####
chl_D3 <- abundance_plot_isolate(species_name = "chlorella", data = genus_data_D3)
coe_D3 <- abundance_plot_isolate(species_name = "coelastrum", data = genus_data_D3)
sce_D3 <- abundance_plot_isolate(species_name = "scenedesmus", data = genus_data_D3)
mon_D3 <- abundance_plot_isolate(species_name = "monoraphidium", data = genus_data_D3)
sel_D3 <- abundance_plot_isolate(species_name = "selenastrum", data = genus_data_D3)

chl_D3
coe_D3
sce_D3
mon_D3
sel_D3

#### D31 Plots ####
chl_D31 <- abundance_plot_isolate(species_name = "chlorella", data = genus_data_D31)
coe_D31 <- abundance_plot_isolate(species_name = "coelastrum", data = genus_data_D31)
sce_D31 <- abundance_plot_isolate(species_name = "scenedesmus", data = genus_data_D31)
mon_D31 <- abundance_plot_isolate(species_name = "monoraphidium", data = genus_data_D31)
sel_D31 <- abundance_plot_isolate(species_name = "selenastrum", data = genus_data_D31)

chl_D31
coe_D31
sce_D31
mon_D31
sel_D31

#### Natural Community Plot ####
#Special case can't be handled by plot abundance function
nat_comm_D0 <- pondTax %>%
  # filter(time == "D3")%>%
  group_by(otu, host_species)%>%
  mutate(agg_rel_abund = sum(rel_abund),
         avg_rel_abund = mean(rel_abund))%>%
  ungroup()%>%
  filter(!(avg_rel_abund < 10e-4))%>%
  mutate(host_species = ifelse(host_species == "naturalCommunity", "Natural Community", NA))%>%
  filter(host_species == "Natural Community")%>%
  slice(1:60)%>%
  #Reorder Genus by average relative abundance to get a smooth curve
  ggplot(aes(x = reorder(otu, -avg_rel_abund), 
  # ggplot(aes(x = reorder(Genus, -avg_rel_abund), 
             y = avg_rel_abund, 
             fill = host_species)) +
  geom_bar(stat = "Identity") +
  theme_minimal() +
  scale_fill_manual(values = c("#D4D5B0"), name = 'Pond Community Day 0', labels = c("Natural Inhabitant")) +
  labs( y = "Relative Abundance (Genus)", x = "Genus", fill = 'Pond Community Day 0',
        title = "Rank Abundance Curves Based on 16S v4 Tag Sequencing") +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(angle = 70, hjust = 1, size = 14,
                                   color = "black"),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom",
        legend.box.just = "center",
        legend.direction = "vertical",
        legend.title = element_text(vjust = .5, size = 12),
        legend.text = element_text(vjust = .5, size = 12),
        plot.title = element_text(hjust = .5, size = 16)) +
  labs( y = "Relative Abundance (Genus)", x = "Genus",  fill = "" ) +
  scale_y_continuous("Average Relative Abundance", limits = c(0, 0.2), breaks = c(0.0, 0.1, 0.2))
# +scale_x_reordered()

#### Combine Plots with Patchwork ####
patchwork <- (nat_comm_D0 / (chl_D3 | chl_D31 + rremove("ylab")) / (coe_D3 | coe_D31 + rremove("ylab")) / (sce_D3 | sce_D31 + rremove("ylab")) / (mon_D3 | mon_D31 + rremove("ylab")) / (sel_D3 | sel_D31+ rremove("ylab"))) 

patchwork + 
  plot_annotation(
    tag_levels = "A"
  ) + plot_layout(
    guides = 'collect'
  ) & theme(
    legend.spacing = unit(0, 'cm'),
    legend.box.just = 'left',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.text.x = element_text(size = 10)
  )

ggsave(filename = "fig1_rewrite.pdf",  dpi = 300, path = "./fig1", width = 14, height = 24)

#### Percentage of Phycosphere OTUs Isolated ####

grbc_genera <- collectionTax %>%
  select(Genus, host_species)%>% # Genus,
  distinct()%>%
  mutate(isolated = T)

phyc_genera <- phycosphereTax %>% select(otu, Genus) %>% distinct()


## Calculates the percentage of genera isolated for each algal host
### excludes genera that with relative abundance less that 10^-5
percentage_isolated <- phycosphereTax%>%
  filter(rel_abund >= 1e-5)%>%
  select(Genus, host_species)%>%
  distinct()%>%
  left_join(., grbc_genera)%>%
  mutate(isolated = ifelse(is.na(isolated), F, isolated))%>%
  group_by(host_species)%>%
  summarize(isolated_count = sum(isolated),
            total_genera = n(),
            per_isolated = sum(isolated)/n()*100)
