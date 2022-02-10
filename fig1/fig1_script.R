#### Libraries ####
library(tidyverse)
library(data.table)
library(patchwork)
library(ggpubr)

#### Directory Set ####
setwd("C:/Users/jlaue/OneDrive/Documents/Denef_Lab_Stuff/R_scripts_git/Github Reorganization/")
getwd()

#### Plot Function ####
abundance_plot_isolate <- function(species_name, data, color_data, Genus_species) {
  plot <- data %>%
    filter(host_species == species_name) %>%
    # arrange(desc(avg_rel_abund)) %>% slice(1:30) %>% # limit curve to 30 genera
    #Reorder Genus by average relative abundance to get a smooth curve
    ggplot(aes(x = reorder(Genus, -avg_rel_abund), y = avg_rel_abund, fill = color)) +
    geom_bar(stat = "Identity") +
    theme_minimal() +
    scale_fill_manual(values = c(colors[(colors$host_species==species_name),]$color, "black"),
                      name =  '', labels = c("Cultured Isolate", "Not Isolated")) +
    guides(fill = guide_legend(
      title = substitute(paste(italic(Genus_species))),
      title.position = "top",
      nrow = 2, byrow = TRUE
    )) +
    labs(y = "Relative Abundance (Genus)",
         x = "",
         fill = '') +
    # , title = "") +
    scale_y_continuous("Average Relative Abundance", 
                       limits = c(0, 0.32), 
                       breaks = c(0.0, 0.1, 0.2, 0.3)) +
    theme(axis.title.x = element_text(size = 12),
          axis.text.x = element_text(angle = 70, 
                                     hjust = 1,
                                     size = 14,
                                     color = color_data$color),
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

family_abundance_plot <- function(species_name, data, color_data, Genus_species){
  plot <- data %>%
    filter(host_species == species_name) %>%
    ggplot(aes(x = reorder(Family, -avg_rel_abund), y = avg_rel_abund, fill = color)) +
    geom_bar(stat = "Identity") +
    theme_minimal() +
    scale_fill_manual(values = c(colors[(colors$host_species==species_name),]$color, "black"),
                      name =  '', labels = c("Cultured Isolate", "Not Isolated")) +
    guides(fill = guide_legend(
      title = substitute(paste(italic(Genus_species))),
      title.position = "top",
      nrow = 2, byrow = TRUE
    )) +
    labs(y = "Relative Abundance (Family)",
         x = "",
         fill = '') +
    # , title = "") +
    scale_y_continuous("Average Relative Abundance", 
                       limits = c(0, 0.32), 
                       breaks = c(0.0, 0.1, 0.2, 0.3)) +
    theme(axis.title.x = element_text(size = 12),
          axis.text.x = element_text(angle = 70, 
                                     hjust = 1,
                                     size = 14,
                                     color = color_data$color),
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
#### Jar Info ####
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

#### READ CSV HERE ####
jar_no_grbc <- read_csv(file = "./raw_data/Jar_Number_Correction.csv",
                        na = "empty", show_col_types = F) %>%
  mutate(jar = paste0(jar_no, "-"),
         point1s = str_replace(Isolate_Number, "(?<=\\d)_", ".1_"),
         point2s = str_replace(Isolate_Number, "(?<=\\d)_", ".2_"),
         point3s = str_replace(Isolate_Number, "(?<=\\d)_", ".3_"))%>%
  gather(., key = "point_added", value = "Isolate_Number", c("Isolate_Number", "point1s",
                                                             "point2s", "point3s"))%>%
  select(Isolate_Number, jar)

#### GRBC Phyloseq Data ####
colors <- read.csv("./csv_files/colors.csv")%>%
  mutate(host_species = hostLong)%>%
  select(host_species, color)

GRBC_tax <- read.csv("./csv_files/collection_tax_data_pureOnly.csv")%>%
  mutate(otu = str_replace_all(string=Isolate,
                               pattern="D",
                               replacement="_D"),
         otu = str_replace_all(string = otu,
                               pattern = "\\,",
                               replacement = "\\."))%>%
  select(-Isolate, -mixed, -count, -otuMatches,
         -numOTUs, -totalReads, -contamFlag)%>%
  left_join(jar_no_grbc, by = c("otu" = "Isolate_Number"))%>%
  # inner_join(jar_no_grbc, by = c("otu" = "Isolate_Number"))%>%
  #Filter out genera that were not in the database, not helpful for our analysis
  filter(!str_detect(Genus, "unclassified")) %>%
  filter(!str_detect(Genus, "uncultured")) %>%
  select("otu", "jar", everything()) %>%
  mutate(otu = gsub("DF", "D31", otu),  # replace DF with D31, quicker
         host_species = ifelse(str_detect(jar, 
                                          paste(jar_num$jar_no[1:9], 
                                                collapse = "|")), "Chlorella",
                        ifelse(str_detect(jar,
                                          paste(jar_num$jar_no[10:18],
                                                collapse = "|")), "Coelastrum",
                        ifelse(str_detect(jar,
                                          paste(jar_num$jar_no[19:27],
                                                collapse = "|")), "Scenedesmus",
                        ifelse(str_detect(jar,
                                          paste(jar_num$jar_no[28:36],
                                                collapse = "|")), "Monoraphidium",
                        ifelse(str_detect(jar,
                                          paste(jar_num$jar_no[37:45],
                                                collapse = "|")), "Selenastrum",
                               "NA")))))) %>%
  left_join(., colors) %>%
  mutate(time = ifelse(str_detect(otu, "_D31"), "D31", "D3"),
         color = ifelse(is.na(color), "black", color),
         # color = ifelse(host_species == "Chlorella", "#7FC97F", 
         #         ifelse(host_species == "Coelastrum", "#BEAED4",
         #         ifelse(host_species == "Scenedesmus", "#E0115F",
         #         ifelse(host_species == "Monoraphidium", "#FDC086",
         #         ifelse(host_species == "Selenastrum", "#386CB0", "black"))))),
         replicate = ifelse(str_detect(jar, 
                                       paste(jar_num$jar_no[seq(1, 45, 3)], 
                                             collapse = "|")), "1",
                     ifelse(str_detect(jar,
                                       paste(jar_num$jar_no[seq(2, 45, 3)],
                                             collapse = "|")), "2",
                     ifelse(str_detect(jar,
                                       paste(jar_num$jar_no[seq(3, 45, 3)],
                                             collapse = "|")), "3", "NA"))),
         pond = ifelse(str_detect(jar,
                                  paste(jar_num$jar_no[c(1:3,10:12,19:21,28:30,37:39)],
                                        collapse = "|")), "1",
                ifelse(str_detect(jar,
                                  paste(jar_num$jar_no[c(4:6,13:15,22:24,31:33,40:42)],
                                        collapse = "|")), "2",
                ifelse(str_detect(jar,
                                  paste(jar_num$jar_no[c(7:9,16:18,25:27,34:36,43:45)],
                                        collapse = "|")), "3", "NA"))),
         Genus=str_replace_all(string=Genus,
                               pattern="Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
                               replacement="A-N-P-R")) %>%
  select(otu, jar, replicate, time, pond, color, everything())

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

# Phycosphere_otu_count <- fread("Phycosphere_stability.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared.txt")%>%
Phycosphere_otu_count <- fread("./raw_data/mothur_outputs_phycosphere/Phycosphere_stability.file.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared.txt")%>%
  select(-label, -numOtus) %>%
  pivot_longer(cols=-Group, names_to="otu", values_to="count")%>%
  group_by(Group)%>%
  mutate("total_count_group" = sum(count))%>%
  ungroup()%>%
  mutate(rel_abund = count/total_count_group,
         host_species = ifelse(str_detect(Group,  #Add Host Species Information
                                          paste(jar_num$jar_no[1:9],
                                                collapse = "|")), "Chlorella",
                        ifelse(str_detect(Group,
                                          paste(jar_num$jar_no[10:18],
                                                collapse = "|")), "Coelastrum",
                        ifelse(str_detect(Group,
                                          paste(jar_num$jar_no[19:27],
                                                collapse = "|")), "Scenedesmus",
                        ifelse(str_detect(Group,
                                          paste(jar_num$jar_no[28:36],
                                                collapse = "|")), "Monoraphidium",
                        ifelse(str_detect(Group,
                                          paste(jar_num$jar_no[37:45],
                                                collapse = "|")), "Selenastrum",
                               "Natural Community"))))),
         replicate = ifelse(str_detect(Group, paste(jar_num$jar_no[seq(1, 45, 3)],
                                                    collapse = "|")), "1",
                     ifelse(str_detect(Group, paste(jar_num$jar_no[seq(2, 45, 3)],
                                                    collapse = "|")), "2",
                     ifelse(str_detect(Group, paste(jar_num$jar_no[seq(3, 45, 3)],
                                                    collapse = "|")), "3",
                     ifelse(str_detect(Group, "R1") | str_detect(Group, "-1"), "1", 
                     ifelse(str_detect(Group, "R2") | str_detect(Group, "-2"), "2",
                     ifelse(str_detect(Group, "R3"), "3", "NA")))))),
         pond = ifelse(str_detect(Group,   #Add Pond info
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
agg_Genus_data_D3 <- Phycosphere_otu_count %>% 
  filter((str_detect(Group,"-D3-") | str_detect(Group, "P")) & str_detect(Group,"022um"))%>%
  inner_join(., Phycosphere_tax) %>% #by OTU
  group_by(Genus, Group, host_species, pond) %>% #Triplicate 1 is J1,2,3 and so on for each host species
  summarize(agg_rel_abund=sum(rel_abund), .groups = "drop") %>% # Sum
  # group_by(Genus, triplicate, host_species) %>%
  # summarize(avg_rel_abund = mean(agg_rel_abund), .groups = "drop") %>% # Mean
  mutate(time = "D3",
         Genus=str_replace_all(string = Genus,
                               pattern = "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
                               replacement="A-N-P-R")) 

Phycosphere_GRBC_D3 <- inner_join(GRBC_tax, agg_Genus_data_D3) %>%
  filter(!str_detect(otu, "D31")) %>%
  select(Genus, host_species, color)%>%
  distinct()

agg_Genus_data_D3_color <- agg_Genus_data_D3 %>% 
  #Filter out reads below a threshold (1/500 worked best) that are not part of our isolate collection
  filter(!(agg_rel_abund < 1/500) | Genus %in% Phycosphere_GRBC_D3$Genus) %>%
  group_by(Genus, host_species) %>%
  summarize(avg_rel_abund = mean(agg_rel_abund)) %>% #Take average of all the groups for one host species
  ungroup() %>%
  left_join(., Phycosphere_GRBC_D3, by = c("Genus", "host_species")) %>%
  mutate(color = if_else(condition = is.na(color), true = "black", false = color),
         time = "D3",
         host_species = factor(host_species, levels = c("Chlorella",
                                                        "Coelastrum",
                                                        "Scenedesmus",
                                                        "Monoraphidium",
                                                        "Selenastrum",
                                                        "Natural Community")))

#### Genus Level Data Organization (D31) ####
agg_Genus_data_D31 <- Phycosphere_otu_count %>% 
  filter(str_detect(Group,"-D31-") & str_detect(Group,"022um"))%>%
  left_join(., Phycosphere_tax)%>%
  # inner_join(., Phycosphere_tax) %>% #by OTU
  group_by(Genus, Group, host_species, pond) %>% #Triplicate 1 is J1,2,3 and so on for each host species
  summarize(agg_rel_abund=sum(rel_abund), .groups = "drop") %>% # Sum
  # group_by(Genus, triplicate, host_species) %>%
  # summarize(avg_rel_abund = mean(agg_rel_abund), .groups = drop) %>% # Mean
  mutate(time = "D31",
         Genus=str_replace_all(string = Genus,
                               pattern = "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
                               replacement="A-N-P-R")) 

Phycosphere_GRBC_D31 <- inner_join(GRBC_tax, agg_Genus_data_D31) %>%
  filter(str_detect(otu, "D31"))%>%
  select(Genus, host_species, color) %>%
  distinct()

agg_Genus_data_D31_color <- agg_Genus_data_D31 %>% 
  #Filter out reads below the detection limit that are not part of our isolate collection
  filter(!(agg_rel_abund < 1/500) | Genus %in% Phycosphere_GRBC_D31$Genus) %>%
  group_by(Genus, host_species) %>%
  summarize(avg_rel_abund = mean(agg_rel_abund)) %>% #Take average of all the groups for one host species
  ungroup() %>%
  left_join(., Phycosphere_GRBC_D3, by = c("Genus", "host_species")) %>%
  mutate(color = if_else(condition = is.na(color), true = "black", false = color)) %>%
  mutate(time = "D3")  %>% 
  mutate(host_species = factor(host_species, levels = c("Chlorella",
                                                        "Coelastrum",
                                                        "Scenedesmus",
                                                        "Monoraphidium",
                                                        "Selenastrum",
                                                        "Natural Community")))

#### Needed? ####
# Needed?
agg_Genus_data <- merge(agg_Genus_data_D3, agg_Genus_data_D31, all=TRUE)

Phycosphere_GRBC_cont <- anti_join(GRBC_tax, agg_Genus_data) %>%
  mutate(contaminated = TRUE)
Phycosphere_GRBC_cont %>% group_by(Genus) %>%
  summarize(distinct_Genus = n_distinct(Genus))

#### Genus Level Color Schemes ####
d3_color_schemes_split <- agg_Genus_data_D3_color%>%
  select(host_species, color, Genus, avg_rel_abund)%>%
  group_by(host_species)%>%
  arrange(desc(avg_rel_abund))%>%
  split(., .$host_species)

d31_color_schemes_split <- agg_Genus_data_D31_color%>%
  select(host_species, color, Genus, avg_rel_abund)%>%
  group_by(host_species)%>%
  arrange(desc(avg_rel_abund))%>%
  split(., .$host_species)

#### D3 Plots ####
chl_D3 <- abundance_plot_isolate(species_name = "Chlorella",
                                 data = agg_Genus_data_D3_color,
                                 color_data = d3_color_schemes_split[["Chlorella"]],
                                 # bar_color = "#7FC97F",
                                 Genus_species = "Chlorella sorokiniana" )


coe_D3 <- abundance_plot_isolate(species_name = "Coelastrum",
                                 data = agg_Genus_data_D3_color,
                                 color_data = d3_color_schemes_split[["Coelastrum"]],
                                 # bar_color ="#BEAED4",
                                 Genus_species = "Coelastrum microporum" )

sce_D3 <- abundance_plot_isolate(species_name = "Scenedesmus",
                                 data = agg_Genus_data_D3_color,
                                 color_data = d3_color_schemes_split[["Scenedesmus"]],
                                 # bar_color ="#E0115F",
                                 Genus_species = "Scenedesmus acuminatus")

mon_D3 <- abundance_plot_isolate(species_name = "Monoraphidium",
                                 data = agg_Genus_data_D3_color,
                                 color_data = d3_color_schemes_split[["Monoraphidium"]],
                                 # bar_color = "#FDC086",
                                 Genus_species = "Monoraphidium minutum")


sel_D3 <- abundance_plot_isolate(species_name = "Selenastrum",
                                 data = agg_Genus_data_D3_color,
                                 color_data = d3_color_schemes_split[["Selenastrum"]],
                                 # bar_color = "#386CB0",
                                 Genus_species = "Selenastrum capricornutum")

#### D31 Plots ####

chl_D31 <- abundance_plot_isolate(species_name = "Chlorella",
                                  data = agg_Genus_data_D31_color,
                                  color_data = d31_color_schemes_split[["Chlorella"]],
                                  # bar_color = "#7FC97F",
                                  Genus_species = "Chlorella sorokiniana" )


coe_D31 <- abundance_plot_isolate(species_name = "Coelastrum",
                                  data = agg_Genus_data_D31_color,
                                  color_data = d31_color_schemes_split[["Coelastrum"]],
                                  # bar_color ="#BEAED4",
                                  Genus_species = "Coelastrum microporum" ) 

sce_D31 <- abundance_plot_isolate(species_name = "Scenedesmus",
                                  data = agg_Genus_data_D31_color,
                                  color_data = d31_color_schemes_split[["Scenedesmus"]],
                                  # bar_color ="#E0115F",
                                  Genus_species = "Scenedesmus acuminatus")

mon_D31 <- abundance_plot_isolate(species_name = "Monoraphidium",
                                  data = agg_Genus_data_D31_color,
                                  color_data = d31_color_schemes_split[["Monoraphidium"]],
                                  # bar_color = "#FDC086",
                                  Genus_species = "Monoraphidium minutum") 

sel_D31 <- abundance_plot_isolate(species_name = "Selenastrum",
                                  data = agg_Genus_data_D31_color,
                                  color_data = d31_color_schemes_split[["Selenastrum"]],
                                  # bar_color = "#386CB0",
                                  Genus_species = "Selenastrum capricornutum")

#### Natural Community Plot ####
#Special case can't be handled by plot abundance function
nat_comm_D0 <- agg_Genus_data_D3 %>% 
  #Filter out reads below the detection limit that are not part of our isolate collection
  group_by(Genus, host_species) %>%
  summarize(avg_rel_abund = mean(agg_rel_abund)) %>% #Take average of all the groups for one host species
  ungroup() %>%
  filter(!(avg_rel_abund < 10e-4) | Genus %in% Phycosphere_GRBC_D3$Genus) %>%
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
ggsave(filename = "Isolates D3-D31 Phycosphere Comparison-newTax.png",  dpi = 300, path = "./fig1", width = 14, height = 24)

ggsave(filename = "Isolates D3-D31 Phycosphere Comparison.pdf",  dpi = 300, path = "./fig1", width = 14, height = 18)

#### Family Level Data Organization (D3) ####
agg_family_data_D3 <- Phycosphere_otu_count %>% 
  filter((str_detect(Group,"-D3-") | str_detect(Group, "P")) & str_detect(Group,"022um"))%>%
  inner_join(., Phycosphere_tax) %>% #by OTU
  group_by(Family, Group, host_species, pond) %>% #Triplicate 1 is J1,2,3 and so on for each host species
  summarize(agg_rel_abund=sum(rel_abund), .groups = "drop") %>% # Sum
  mutate(time = "D3")

Phycosphere_GRBC_D3_fam <- inner_join(GRBC_tax, agg_family_data_D3) %>%
  filter(!str_detect(otu, "D31")) %>%
  select(Family, host_species, color)%>%
  distinct()

agg_family_data_D3_color <- agg_family_data_D3 %>% 
  #Filter out reads below a threshold (1/500 worked best) that are not part of our isolate collection
  filter(!(agg_rel_abund < 1/500) | Family %in% Phycosphere_GRBC_D3_fam$Family) %>%
  group_by(Family, host_species) %>%
  summarize(avg_rel_abund = mean(agg_rel_abund)) %>% #Take average of all the groups for one host species
  ungroup() %>%
  left_join(., Phycosphere_GRBC_D3_fam, by = c("Family", "host_species")) %>%
  mutate(color = if_else(condition = is.na(color), true = "black", false = color),
         time = "D3",
         host_species = factor(host_species, levels = c("Chlorella",
                                                        "Coelastrum",
                                                        "Scenedesmus",
                                                        "Monoraphidium",
                                                        "Selenastrum",
                                                        "Natural Community")))

#### Family Level Data Organization (D31) ####
agg_family_data_D31 <- Phycosphere_otu_count %>% 
  filter((str_detect(Group,"-D31-") | str_detect(Group, "P")) & str_detect(Group,"022um"))%>%
  inner_join(., Phycosphere_tax) %>% #by OTU
  group_by(Family, Group, host_species, pond) %>% #Triplicate 1 is J1,2,3 and so on for each host species
  summarize(agg_rel_abund=sum(rel_abund), .groups = "drop") %>% # Sum
  mutate(time = "D31")

Phycosphere_GRBC_D31_fam <- inner_join(GRBC_tax, agg_family_data_D31) %>%
  filter(str_detect(otu, "D31"))%>%
  select(Family, host_species, color) %>%
  distinct()

agg_family_data_D31_color <- agg_family_data_D31 %>% 
  #Filter out reads below a threshold (1/500 worked best) that are not part of our isolate collection
  filter(!(agg_rel_abund < 1/500) | Family %in% Phycosphere_GRBC_D31_fam$Family) %>%
  group_by(Family, host_species) %>%
  summarize(avg_rel_abund = mean(agg_rel_abund)) %>% #Take average of all the groups for one host species
  ungroup() %>%
  left_join(., Phycosphere_GRBC_D31_fam, by = c("Family", "host_species")) %>%
  mutate(color = if_else(condition = is.na(color), true = "black", false = color),
         time = "D3",
         host_species = factor(host_species, levels = c("Chlorella",
                                                        "Coelastrum",
                                                        "Scenedesmus",
                                                        "Monoraphidium",
                                                        "Selenastrum",
                                                        "Natural Community")))


#### Family Level Color Schemes ####
d3_famColors_split <- agg_family_data_D3_color%>%
  select(host_species, color, Family, avg_rel_abund)%>%
  group_by(host_species)%>%
  arrange(desc(avg_rel_abund))%>%
  split(., .$host_species)

d31_famColors_split <- agg_family_data_D31_color%>%
  select(host_species, color, Family, avg_rel_abund)%>%
  group_by(host_species)%>%
  arrange(desc(avg_rel_abund))%>%
  split(., .$host_species)
#### D3 Plots (Family) ####
chl_D3_fam <- family_abundance_plot(species_name = "Chlorella",
                                     data = agg_family_data_D3_color,
                                     color_data = d3_famColors_split[["Chlorella"]],
                                     Genus_species = "Chlorella sorokiniana" )


coe_D3_fam <- family_abundance_plot(species_name = "Coelastrum",
                                     data = agg_family_data_D3_color,
                                     color_data = d3_famColors_split[["Coelastrum"]],
                                     Genus_species = "Coelastrum microporum" )

sce_D3_fam <- family_abundance_plot(species_name = "Scenedesmus",
                                     data = agg_family_data_D3_color,
                                     color_data = d3_famColors_split[["Scenedesmus"]],
                                     Genus_species = "Scenedesmus acuminatus")

mon_D3_fam <- family_abundance_plot(species_name = "Monoraphidium",
                                 data = agg_family_data_D3_color,
                                 color_data = d3_famColors_split[["Monoraphidium"]],
                                 Genus_species = "Monoraphidium minutum")


sel_D3_fam <- family_abundance_plot(species_name = "Selenastrum",
                                 data = agg_family_data_D3_color,
                                 color_data = d3_famColors_split[["Selenastrum"]],
                                 Genus_species = "Selenastrum capricornutum")

#### D31 Plots (Family) ####
chl_D31_fam <- family_abundance_plot(species_name = "Chlorella",
                                  data = agg_family_data_D31_color,
                                  color_data = d31_famColors_split[["Chlorella"]],
                                  Genus_species = "Chlorella sorokiniana" )


coe_D31_fam <- family_abundance_plot(species_name = "Coelastrum",
                                  data = agg_family_data_D31_color,
                                  color_data = d31_famColors_split[["Coelastrum"]],
                                  Genus_species = "Coelastrum microporum" ) 

sce_D31_fam <- family_abundance_plot(species_name = "Scenedesmus",
                                  data = agg_family_data_D31_color,
                                  color_data = d31_famColors_split[["Scenedesmus"]],
                                  Genus_species = "Scenedesmus acuminatus")

mon_D31_fam <- family_abundance_plot(species_name = "Monoraphidium",
                                  data = agg_family_data_D31_color,
                                  color_data = d31_famColors_split[["Monoraphidium"]],
                                  Genus_species = "Monoraphidium minutum") 

sel_D31_fam <- family_abundance_plot(species_name = "Selenastrum",
                                  data = agg_family_data_D31_color,
                                  color_data = d31_famColors_split[["Selenastrum"]],
                                  Genus_species = "Selenastrum capricornutum")

#### Natural Community Plot ####
nat_comm_D0_fam <- agg_family_data_D3 %>%
  group_by(Family, host_species) %>%
  summarize(avg_rel_abund = mean(agg_rel_abund)) %>%
  ungroup() %>%
  filter(!(avg_rel_abund < 10e-4) | Family %in% Phycosphere_GRBC_D3_fam$Family) %>%
  filter(host_species == "Natural Community") %>%
  ggplot(aes(x = reorder(Family, -avg_rel_abund), 
             y = avg_rel_abund, 
             fill = host_species)) +
  geom_bar(stat = "Identity") +
  theme_minimal() +
  scale_fill_manual(values = c("#D4D5B0"),
                    name = 'Pond Community Day 0',
                    labels = c("Natural Inhabitant")) +
  labs( y = "Relative Abundance (Family)",
        x = "Family", fill = 'Pond Community Day 0',
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
  labs( y = "Relative Abundance (Family)", x = "Family",  fill = "" ) +
  scale_y_continuous("Average Relative Abundance", limits = c(0, 0.2), breaks = c(0.0, 0.1, 0.2))

#### Combine Plots with Patchwork ####
patchwork_fam <- (nat_comm_D0_fam / (chl_D3_fam | chl_D31_fam + rremove("ylab")) / (coe_D3_fam | coe_D31_fam + rremove("ylab")) / (sce_D3_fam | sce_D31_fam + rremove("ylab")) / (mon_D3_fam | mon_D31_fam + rremove("ylab")) / (sel_D3_fam | sel_D31_fam + rremove("ylab"))) 

patchwork_fam + 
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
ggsave(filename = "Isolates D3-D31 Phycosphere Comparison-Family Level.png",  dpi = 300, path = "./fig1", width = 14, height = 24)

ggsave(filename = "Isolates D3-D31 Phycosphere Comparison Family Level.pdf",  dpi = 300, path = "./fig1", width = 14, height = 18)

