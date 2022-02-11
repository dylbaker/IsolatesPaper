#### Libraries ####
library(tidyverse)
library(data.table)
library(patchwork)
library(ggpubr)

#### Directory Set ####
setwd("C:/Users/jlaue/OneDrive/Documents/Denef_Lab_Stuff/IsolatesPaper/")
getwd()

#### Jar Info ####
jar_no_grbc <- read_csv(file = "./raw_data/Jar_Number_Correction.csv",
                        na = "empty", show_col_types = F) %>%
  mutate(jar = paste0(jar_no, "-"),
         point1s = str_replace(Isolate_Number, "(?<=\\d)_", ".1_"),
         point2s = str_replace(Isolate_Number, "(?<=\\d)_", ".2_"),
         point3s = str_replace(Isolate_Number, "(?<=\\d)_", ".3_"),
         point4s = str_replace(Isolate_Number, "(?<=\\d)_", ".4_"))%>%
  gather(., key = "point_added", value = "Isolate_Number", c("Isolate_Number", "point1s",
                                                             "point2s", "point3s", "point4s"))%>%
  select(Isolate_Number, jar)

jar_num <- tibble("jar_no" = paste0("J", 1:45,"-"))%>%
  mutate(host_species=ifelse(str_detect(jar_no,
                                        paste0("J", 1:9, "-")), "Chlorella",
                      ifelse(str_detect(jar_no,
                                        paste0("J", 10:18, "-")), "Coelastrum",
                      ifelse(str_detect(jar_no,
                                        paste0("J", 19:27, "-")), "Scenedesmus",
                      ifelse(str_detect(jar_no,
                                        paste0("J", 28:36, "-")), "Monoraphidium",
                      ifelse(str_detect(jar_no,
                                        paste0("J", 37:45, "-")), "Selenastrum",
                                                         ""))))))

#### GRBC Phyloseq Data ####
colors <- read.csv("./csv_files/colors.csv")%>%
  mutate(host_species = hostLong)%>%
  select(host_species, color)

GRBC_tax <- read.csv("./csv_files/collection_tax_data_pureOnly.csv")%>%
  mutate(Isolate = str_replace_all(string=Isolate,
                               pattern="D",
                               replacement="_D"),
         Isolate = str_replace_all(string = Isolate,
                               pattern = "\\,",
                               replacement = "\\."))%>%
  select(-mixed, -count, -otuMatches, -numOTUs, -totalReads, -contamFlag)%>%
  left_join(., jar_no_grbc, by = c("Isolate" = "Isolate_Number"))%>%
  # inner_join(jar_no_grbc, by = c("Isolate" = "Isolate_Number"))%>%
  #Filter out genera that were not in the database, not helpful for our analysis
  # filter(!str_detect(Genus, "unclassified")) %>%
  # filter(!str_detect(Genus, "uncultured")) %>%
  mutate(Isolate = gsub("DF", "D31", Isolate))%>%  # replace DF with D31, quicker
  left_join(., jar_num, by = c("jar" = "jar_no"))%>%
  left_join(., colors) %>%
  mutate(time = ifelse(str_detect(Isolate, "_D31"), "D31", "D3"),
         # color = ifelse(is.na(color), "black", color),
         replicate = ifelse(str_detect(jar, paste(jar_num$jar_no[seq(1, 45, 3)],
                                                  collapse = "|")), "1",
                     ifelse(str_detect(jar, paste(jar_num$jar_no[seq(2, 45, 3)],
                                                  collapse = "|")), "2",
                     ifelse(str_detect(jar, paste(jar_num$jar_no[seq(3, 45, 3)],
                                                  collapse = "|")), "3", "NA"))),
         pond = ifelse(str_detect(jar, paste(jar_num$jar_no[c(1:3,10:12,19:21,28:30,37:39)],
                                             collapse = "|")), "1",
                ifelse(str_detect(jar, paste(jar_num$jar_no[c(4:6,13:15,22:24,31:33,40:42)],
                                             collapse = "|")), "2",
                ifelse(str_detect(jar, paste(jar_num$jar_no[c(7:9,16:18,25:27,34:36,43:45)],
                                             collapse = "|")), "3", "NA"))),
         Genus=str_replace_all(string=Genus,
                               pattern="Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
                               replacement="A-N-P-R"))

#### Phycosphere Community Data  (FREAD HERE) ####
# Phycosphere_tax <- fread("Phycosphere_stability.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy") %>%
Phycosphere_tax <- fread("./raw_data/mothur_outputs_phycosphere/Phycosphere_stability.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy")%>%
  rename_all(tolower) %>% 
  mutate(taxonomy=str_replace_all(string=taxonomy,
                                  pattern="\\(\\d*\\)",
                                  replacement=""))%>%
  mutate(taxonomy=str_replace_all(string=taxonomy,
                                  pattern=";$",
                                  replacement=""))%>%
  separate(taxonomy,
           into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep=";")%>%
  filter(!str_detect(Genus, "unclassified"))%>%
  filter(!str_detect(Genus, "uncultured"))%>%
  filter(!str_detect(Genus, "uncultured_fa"))%>%
  filter(!str_detect(Genus, "Unknown")) 


#### Phycosphere Community Data  (FREAD HERE) ####
Phycosphere_otu_count <- fread("./raw_data/mothur_outputs_phycosphere/Phycosphere_stability.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared.txt")%>%
  select(-label, -numOtus) %>%
  pivot_longer(cols=-Group, names_to="otu", values_to="count")%>%
  group_by(Group)%>%
  mutate("total_count_group" = sum(count))%>%
  ungroup()%>%
  mutate(rel_abund = count/total_count_group,
         #Add Host Species Information
         host_species = ifelse(str_detect(Group,
                                          paste(jar_num$jar_no[1:9],
                                                collapse = "|")), "Chlorella",
                        ifelse(str_detect(Group, paste(jar_num$jar_no[10:18],
                                                       collapse = "|")), "Coelastrum",
                        ifelse(str_detect(Group, paste(jar_num$jar_no[19:27],
                                                       collapse = "|")), "Scenedesmus",
                        ifelse(str_detect(Group, paste(jar_num$jar_no[28:36],
                                                       collapse = "|")), "Monoraphidium",
                        ifelse(str_detect(Group, paste(jar_num$jar_no[37:45],
                                                       collapse = "|")), "Selenastrum",
                               "Natural Community"))))),
         #Add Replicate info
         replicate = ifelse(str_detect(Group, paste(jar_num$jar_no[seq(1, 45, 3)],
                                                    collapse = "|")), "1",
                     ifelse(str_detect(Group, paste(jar_num$jar_no[seq(2, 45, 3)],
                                                    collapse = "|")), "2",
                     ifelse(str_detect(Group, paste(jar_num$jar_no[seq(3, 45, 3)],
                                                    collapse = "|")), "3",
                     ifelse(str_detect(Group, "R1") | str_detect(Group, "-1"), "1", 
                     ifelse(str_detect(Group, "R2") | str_detect(Group, "-2"), "2",
                     ifelse(str_detect(Group, "R3"), "3", "NA")))))),
         #Add Pond info
         pond = ifelse(str_detect(Group,
                                  paste(jar_num$jar_no[c(1:3,10:12,19:21,28:30,37:39)],
                                        collapse = "|")), "1",
                       ifelse(str_detect(Group,
                                         paste(jar_num$jar_no[c(4:6,13:15,22:24,31:33,40:42)],
                                               collapse = "|")), "2",
                       ifelse(str_detect(Group,
                                         paste(jar_num$jar_no[c(7:9,16:18,25:27,34:36,43:45)],
                                               collapse = "|")), "3",
                       ifelse(str_detect(Group, "P1"), "1",
                       ifelse(str_detect(Group, "P2"), "2",
                       ifelse(str_detect(Group, "P3"), "3", "NA")))))))



#### Genus Level Data Organization (D3) ####
genus_data_D3 <- Phycosphere_otu_count %>% 
  filter((str_detect(Group,"-D3-") | str_detect(Group, "P")) & str_detect(Group,"022um"))%>%
  inner_join(., Phycosphere_tax) %>% #by OTU
  group_by(Genus, Group, host_species, pond) %>% #Triplicate 1 is J1,2,3 and so on for each host species
  summarize(agg_rel_abund=sum(rel_abund), .groups = "drop") %>% # Sum
  group_by(Genus, host_species) %>%
  summarize(avg_rel_abund = mean(agg_rel_abund), .groups = "drop") %>% # Mean
  mutate(time = "D3",
         Genus=str_replace_all(string = Genus,
                               pattern = "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
                               replacement="A-N-P-R"))%>%
  left_join(., GRBC_tax)%>%
  # inner_join(GRBC_tax, .)%>%
  # filter(!str_detect(otu, "D31"))%>%
  distinct()%>%
  mutate(color = if_else(is.na(color), "black", color),
         isolated = if_else(is.na(otu), F, T))

#### Genus Level Data Organization (D31) ####
genus_data_D31 <- Phycosphere_otu_count %>% 
  filter(str_detect(Group,"-D31-") & str_detect(Group,"022um"))%>%
  left_join(., Phycosphere_tax)%>%
  # inner_join(., Phycosphere_tax) %>% #by OTU
  group_by(Genus, Group, host_species, pond) %>% #Triplicate 1 is J1,2,3 and so on for each host species
  summarize(agg_rel_abund=sum(rel_abund), .groups = "drop") %>% # Sum
  group_by(Genus, host_species) %>%
  summarize(avg_rel_abund = mean(agg_rel_abund), .groups = "drop") %>% # Mean
  mutate(time = "D31",
         Genus=str_replace_all(string = Genus,
                               pattern = "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
                               replacement="A-N-P-R"))%>%
  left_join(., GRBC_tax)%>%
  # inner_join(GRBC_tax, .)%>%
  # filter(str_detect(otu, "D31"))%>%
  distinct()%>%
  mutate(color = if_else(is.na(color), "black", color),
         isolated = if_else(is.na(otu), F, T))

#### Plotting Fxn (Old) ####
abundance_plot_isolate <- function(species_name, data) { #, color_data, Genus_species
  df <- filter(data, host_species == species_name)%>%
    select(Genus, avg_rel_abund, isolated, color)%>%
    distinct()%>%
    arrange(desc(avg_rel_abund)) %>% slice(1:30) # limit curve to 30 genera
  plot <- df %>%
    #Reorder Genus by average relative abundance to get a smooth curve
    ggplot(aes(x = reorder(Genus, -avg_rel_abund), y = avg_rel_abund, fill = isolated)) +
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
chl_D3 <- abundance_plot_isolate(species_name = "Chlorella", data = genus_data_D3)
coe_D3 <- abundance_plot_isolate(species_name = "Coelastrum", data = genus_data_D3)
sce_D3 <- abundance_plot_isolate(species_name = "Scenedesmus", data = genus_data_D3)
mon_D3 <- abundance_plot_isolate(species_name = "Monoraphidium", data = genus_data_D3)
sel_D3 <- abundance_plot_isolate(species_name = "Selenastrum", data = genus_data_D3)
# 
# chl_D3
# coe_D3
# sce_D3
# mon_D3
# sel_D3

#### D31 Plots ####
chl_D31 <- abundance_plot_isolate(species_name = "Chlorella", data = genus_data_D31)
coe_D31 <- abundance_plot_isolate(species_name = "Coelastrum", data = genus_data_D31)
sce_D31 <- abundance_plot_isolate(species_name = "Scenedesmus", data = genus_data_D31)
mon_D31 <- abundance_plot_isolate(species_name = "Monoraphidium", data = genus_data_D31)
sel_D31 <- abundance_plot_isolate(species_name = "Selenastrum", data = genus_data_D31)
# 
# chl_D31
# coe_D31
# sce_D31
# mon_D31
# sel_D31

#### Natural Community Plot ####
#Special case can't be handled by plot abundance function
nat_comm_D0 <- genus_data_D3 %>%
  # agg_Genus_data_D3 %>% 
  #Filter out reads below the detection limit that are not part of our isolate collection
  # group_by(Genus, host_species) %>%
  # summarize(avg_rel_abund = mean(agg_rel_abund)) %>% #Take average of all the groups for one host species
  # ungroup() %>%
  filter(!(avg_rel_abund < 10e-4))%>%
  # filter(!(avg_rel_abund < 10e-4) | Genus %in% Phycosphere_GRBC_D3$Genus) %>%
  # mutate(time = "D3")  %>% 
  # mutate(host_species = factor(host_species, levels = c("Chlorella", "Coelastrum", "Scenedesmus", "Monoraphidium", "Selenastrum", "Natural Community")), 
  # Genus = reorder_within(Genus, -avg_rel_abund, host_species)) %>%
  filter(host_species == "Natural Community") %>%
  #Reorder Genus by average relative abundance to get a smooth curve
  ggplot(aes(x = reorder(Genus, -avg_rel_abund), 
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

ggsave(filename = "fig1_rewrite.png",  dpi = 300, path = "./fig1", width = 14, height = 24)

# ggsave(filename = "Isolates D3-D31 Phycosphere Comparison.pdf",  dpi = 300, path = "./fig1", width = 14, height = 18)

#### Percentage of Phycosphere OTUs Isolated ####

grbc_genera <- GRBC_tax %>%
  select(Genus, host_species)%>% # Genus,
  distinct()%>%
  mutate(isolated = T)

phyc_genera <- Phycosphere_tax %>% select(otu, Genus) %>% distinct()

percentage_isolated <- Phycosphere_otu_count%>%
  # filter(host_species != "Natural Community")%>%
  left_join(., phyc_genera)%>%
  filter(rel_abund >= 1e-5)%>%
  select(Genus, host_species)%>%
  distinct()%>%
  left_join(., grbc_genera)%>%
  mutate(isolated = ifelse(is.na(isolated), F, isolated))%>%
  group_by(host_species)%>%
  summarize(isolated_count = sum(isolated),
            total_genera = n(),
            per_isolated = sum(isolated)/n()*100)
