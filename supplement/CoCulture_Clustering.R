library(tidyverse)
library(ggplot2)
library(lme4)
library(kmlShape)
source("plate_functions.R")

set.seed(09042018)

# # Chlorella
# species_name <- "chlorella"
# optimal_number <- 6
# clusNums <- c(1:6)
# facet_plot_name <- "chlorella_facet_clusters.png"

# Coelastrum
species_name <- "coelastrum"
optimal_number <- 4
clusNums <- c(1:4)
facet_plot_name <- "coelastrum_facet_clusters.png"

# # Scenedesmus
# species_name <- "scenedesmus"
# optimal_number <- 5
# clusNums <- c(1:5)
# facet_plot_name <- "scenedesmus_facet_clusters.png"

# # Monoraphidium
# species_name <- "monoraphidium"
# optimal_number <- 5
# clusNums <- c(1:5)
# facet_plot_name <- "monoraphidium_facet_clusters.png"

# # Selenastrum
# species_name <- "selenastrum"
# optimal_number <- 4
# clusNums <- c(1:4)
# facet_plot_name <- "selenastrum_facet_clusters.png"

#### Contaminant and Effect Data ####
cont_data <- read.csv("./CoCulture_data/GRBC_Match&Cont.csv")%>%
  mutate(Isolate = as.character(str_replace_all(otu, "_", "")))%>%
  group_by(Isolate)%>%
  summarize(contaminated = as.logical(sum(contaminated)))%>%
  ungroup()

all_effects <- read.csv("./CoCulture_Data/all_effects.csv")%>%
  select(Isolate, effect_type, Day_isolated)

#### Set Paths and Prep Data ####

species_path <- paste("CoCulture_Figures/cluster_plots/", species_name, sep = "")


## this could be done without the stats (coculture_data.csv)
species_data <- read.csv("./CoCulture_Data/all_data_with_stats.csv")%>%
  mutate(Isolate = as.character(Isolate))%>%
  filter(host_species == species_name)%>%
  flag_algae()%>%
  mutate(algae_flag = ifelse(is.na(algae_flag), FALSE, algae_flag))

likely_autotrophic_isolates <- species_data %>%
  filter(algae_flag == TRUE)%>%
  select(Isolate, algae_flag)%>%
  distinct()
  
data_4clus <- species_data %>%
  filter(Day_isolated != "Ctrl" & Outlier != TRUE & algae_flag != TRUE)%>%
  group_by(Isolate, read_timeHours)%>%
  summarize(ChlA_100_mean = mean(ChlA_100, na.rm = TRUE))%>%
  ungroup()%>%
  as.data.frame()

# Cool looking plot with final clusters superimposed over all trajectories
myClds <- kmlShape::cldsLong(data_4clus)
traj_plot <- kmlShape::plotTraj(myClds)

# #### Cluster Optimization TAKES FOREVER, Comment once optimal number determined####
# 
# cNums <- c(1:20)
# 
# list_distances <- lapply(cNums, cluster_optimize)%>%
#   bind_rows()
# 
# plot_title <- paste("Cluster Optimization for\n", species_name, " Growth Curves", sep = "")
# 
# in_clus_dist_plot <- ggplot(data = list_distances,
#                             aes(x = Number_of_Clus,
#                                 y = Mean_in_Clus_Dist))+
#   geom_vline(xintercept = optimal_number, color = "red")+
#   geom_point()+
#   geom_line()+
#   labs(x = "Number of Clusters",
#        y = "In-Cluster Distance",
#        title = plot_title)+
#   theme_classic()
#   # theme(legend.position = "none",
#   #       panel.background = element_rect(fill = "transparent"),
#   #       plot.background = element_rect(fill = "transparent", color = NA),
#   #       axis.line = element_line(colour = "black"),
#   #       axis.title.x = element_text(size = 16),
#   #       axis.title.y = element_text(size = 16),
#   #       axis.text.y = element_text(size = 14),
#   #       axis.text.x = element_text(size = 14))
# in_clus_dist_plot
# 
# optimize_filename <- paste(species_name, "cluster", "optimize", "plot", sep = "_")%>%
#   paste(., "png", sep = ".")
# ggsave(optimize_filename,
#        plot = in_clus_dist_plot,
#        device = "png",
#        path = species_path)
# 
#### Using KmlShape to Cluster ####

par(ask = FALSE)
kmlShape::kmlShape(myClds, optimal_number)
par(ask = TRUE)
kmlShape::plot(myClds)

clus_data <- myClds@clusters

Isolate <- names(myClds@clusters)
clusters <- data.frame(Isolate, clus_data)%>%
  mutate(Isolate = as.character(Isolate))

clust_csv_name <- paste(species_name, "clusters", sep = "_")%>%
  paste(., "csv", sep = ".")
write.csv(clusters, clust_csv_name)

#### AC Plot ####
AC <- species_data %>%
  filter(Isolate == "AC")%>%
  group_by(Isolate, read_day)%>%
  summarize(avg_timeHours = mean(read_timeHours, na.rm = TRUE),
            ChlA_100_mean = mean(ChlA_100, na.rm = TRUE),
            ChlA_100_sd = sd(ChlA_100, na.rm = TRUE))

AC_plot <- ggplot(data = AC, aes(x = avg_timeHours, y = ChlA_100_mean)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = ChlA_100_mean - ChlA_100_sd,
                    ymax = ChlA_100_mean + ChlA_100_sd),
                width = .2) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 800)) +
  labs(x = "Time (hours)", y = "Chlorophyll A Fluorescence (RFU)")
plot(AC_plot)

#### Cluster Plots ####

data_4plots <- species_data %>%
  filter(Day_isolated != "Ctrl" & Outlier != TRUE & algae_flag != TRUE)%>%
  group_by(Isolate, read_day, data_sig)%>%
  summarize(avg_timeHours = mean(read_timeHours, na.rm = TRUE),
            ChlA_100_mean = mean(ChlA_100, na.rm = TRUE),
            ChlA_100_sd = sd(ChlA_100, na.rm = TRUE))%>%
  ungroup()%>%
  full_join(., clusters)

cluster_plots <- lapply(clusNums, cluster_plot)
cluster_plot_ggsave(cluster_plots)


#### Faceted Plot ####

# facet_title <- paste()

data_4facet <- left_join(data_4plots, cont_data, by = "Isolate")%>%
  left_join(., all_effects, by = "Isolate")
# %>%
#   filter(contaminated == FALSE | is.na(contaminated))

faceted_plot <- data_4facet %>%
  mutate(clus_data = as.numeric(clus_data))%>%
  # filter(clus_data == clusNum)%>%
  ggplot(data = ., 
         aes(x=avg_timeHours, y=ChlA_100_mean, group=Isolate, color=effect_type), 
         fill = "transparent")+
  facet_wrap(vars(clus_data))+
  geom_line()+
  geom_line(aes(x = avg_timeHours, y = ChlA_100_mean),
            data = AC,
            color = "black",
            size = 1.5)+
  geom_point()+
  geom_errorbar(aes(ymin = ChlA_100_mean - ChlA_100_sd,
                    ymax = ChlA_100_mean + ChlA_100_sd),
                width = .2)+
                # color="green")+
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 800)) +
  labs(x = "Time (hours)", y = "Chlorophyll A Fluorescence (RFU)")+
  scale_color_manual(name = "Effect Observed",
                     values = c("Not Significant" = "gray",
                                "Positive" = "cornflowerblue", 
                                "Negative" = "orange"),
                    guide = guide_legend(reverse = TRUE))+
  theme_classic()
  # theme(legend.position = "none",
  # theme(panel.background = element_rect(fill = "transparent"),
  #       plot.background = element_rect(fill = "transparent", color = NA),
  #       axis.line = element_line(colour = "black"),
  #       axis.title.x = element_text(size = 16),
  #       axis.title.y = element_text(size = 16),
  #       axis.text.y = element_text(size = 14),
  #       axis.text.x = element_text(size = 14))

faceted_plot

ggsave(facet_plot_name, 
       plot = faceted_plot, 
       device = "png", 
       path = species_path, 
       scale = 1, 
       width = 13.5, 
       height = 8, 
       dpi = 600, 
       units = "in", 
       limitsize = TRUE, 
       bg = "transparent")

