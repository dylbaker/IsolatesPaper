library(tidyverse)
library(ggplot2)
library(lme4)
library(kmlShape)

set.seed(09042018)

#### Functions ####
# Attempts Clustering with any number of clusters N
cluster_optimize <- function(n, species, df){
  cNum <- n
  
  species_data <- filter(df, host_species == species)%>%
    select(Isolate, read_timeHours, ChlA_100_mean)
  
  traj <- kmlShape::cldsLong(species_data)
  
  # traj <- df
  
  par(ask = FALSE)
  clds <- kmlShape::kmlShape(traj, cNum, toPlot = "none")
  par(ask = TRUE)
  
  clusters <- data.frame(clds@clusters) %>%
    rownames_to_column("Isolate") %>%
    rename(Cluster = clds.clusters)
  
  traj_with_clus <- left_join(traj@trajLong, clusters, by = "Isolate")%>%
    split(., .$Isolate)
  
  fDist <- list()
  for(i in 1:length(traj_with_clus)){
    cluster <- unique(as.numeric(traj_with_clus[[i]]$Cluster))
    clus_traj_mean <- clds@trajMeans %>%
      filter(iCenters == cluster)
    fDist[[i]] <- data.frame(distFrechet = distFrechet(traj_with_clus[[i]]$read_timeHours,
                                                       traj_with_clus[[i]]$ChlA_100_mean,
                                                       clus_traj_mean$times,
                                                       clus_traj_mean$traj))
  }
  
  mean_in_clus_dist <- data.frame(Mean_in_Clus_Dist = mean(bind_rows(fDist)$distFrechet),
                                  Number_of_Clus = cNum)
  return(mean_in_clus_dist)
}

# Plots clustering attempts with different N for implimentation of elbow method
optimization_plot <- function(species, df, optimal_number){
  plot_title <- paste("Cluster Optimization for\n", str_to_title(species), " Growth Curves", sep = "")
  
  in_clus_dist_plot <- ggplot(data = df,
                              aes(x = Number_of_Clus,
                                  y = Mean_in_Clus_Dist))+
    geom_vline(xintercept = optimal_number, color = "red")+
    geom_point()+
    geom_line()+
    labs(x = "Number of Clusters",
         y = "In-Cluster Distance",
         title = plot_title)+
    theme_classic()
  # theme(legend.position = "none",
  #       panel.background = element_rect(fill = "transparent"),
  #       plot.background = element_rect(fill = "transparent", color = NA),
  #       axis.line = element_line(colour = "black"),
  #       axis.title.x = element_text(size = 16),
  #       axis.title.y = element_text(size = 16),
  #       axis.text.y = element_text(size = 14),
  #       axis.text.x = element_text(size = 14))
  in_clus_dist_plot
  
  optimize_filename <- paste(species, "cluster", "optimize", "plot", sep = "_")%>%
    paste(., "png", sep = ".")
  
  ggsave(optimize_filename,
         plot = in_clus_dist_plot,
         device = "png",
         path = "./supplement/cluster_plots/")
  
  return(in_clus_dist_plot)
}

# Takes optimal number of clusters, species, and plate reader data and generates clusters
cluster <- function(species, optimal_number, df){
  species_data <- filter(df, host_species == species)%>%
    select(Isolate, read_timeHours, ChlA_100_mean)
  
  traj <- kmlShape::cldsLong(species_data)
  
  par(ask = FALSE)
  clds <- kmlShape::kmlShape(traj, optimal_number, toPlot = "none")
  par(ask = TRUE)
  
  return(clds)
}

# Generates Plots for one cluster at a time
cluster_plot <- function(df){
  species <- unique(df$host_species)
  
  cluster <- unique(df$Cluster)
  
  plot_title <- str_to_title(paste(species, "cluster", cluster, sep = " "))
  
  acData <- filter(AC, host_species == species)
  
  plot <- ggplot(data = df, 
                 aes(x=avg_timeHours, y=ChlA_100_mean, color = annotation, grpup=Isolate), 
                 fill = "transparent")+
    geom_line()+
    geom_line(aes(x = avg_timeHours, y = ChlA_100_mean),
              data = acData,
              color = "black",
              size = 1.5)+
    geom_point()+
    geom_errorbar(aes(ymin = ChlA_100_mean - ChlA_100_sd,
                      ymax = ChlA_100_mean + ChlA_100_sd),
                  width = .2)+
    # facet_wrap(~ Cluster) +
    coord_cartesian(xlim = c(0, 250), ylim = c(0, 800)) +
    labs(x = "Time (hours)", y = "Chlorophyll A Fluorescence (RFU)", title = plot_title)+
    theme_classic()+
    theme(panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = NA),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14))
  # legend.position = "none",
  return(plot)
  
}

# Generates Faceted plot of all clusters for a given species
# Uses different labelling for effect observed (Lumps together CC and GR)
# NOTE: can change effect labelling, but must edit data_4facet dataframe as well
facet_plot <- function(df){
  host <- unique(df$host_species)
  
  ac_data <- filter(AC, host_species == host)
  
  ggplot(data = df, 
         aes(x=avg_timeHours, y=ChlA_100_mean, group=Isolate, color=pnnEffect), 
         fill = "transparent")+
    facet_wrap(vars(Cluster))+
    geom_line()+
    geom_line(aes(x = avg_timeHours, y = ChlA_100_mean),
              data = ac_data,
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
                       values = c("No Significant Effect" = "gray",
                                  "Positive" = "cornflowerblue", 
                                  "Negative" = "orange",
                                  "Inconclusive" = "green"),
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
}

#### Data Prep ####
reader_data_raw <- read.csv("./csv_files/coculture_data.csv")

stats_data <- read.csv("./csv_files/logistic_mod_stats_coefs.csv")%>%
  select(Isolate, ccEffect, grEffect, annotation)

AC <- filter(reader_data_raw, Isolate == "AC" & Outlier != TRUE)%>%
  group_by(Isolate, host_species, read_day)%>%
  summarize(avg_timeHours = mean(read_timeHours, na.rm = TRUE),
            ChlA_100_mean = mean(ChlA_100, na.rm = TRUE),
            ChlA_100_sd = sd(ChlA_100, na.rm = TRUE))

reader_data <- reader_data_raw %>%
  filter(Day_isolated != "Ctrl" & Outlier != TRUE)%>%
  group_by(Isolate, host_species, read_timeHours)%>%
  summarize(ChlA_100_mean = mean(ChlA_100, na.rm = TRUE),
            ChlA_100_sd = sd(ChlA_100, na.rm = TRUE))%>%
  ungroup()%>%
  as.data.frame()

algal_flag <- reader_data %>%
  filter(read_timeHours == 0)%>%
  group_by(host_species)%>%
  mutate(day0_mean = mean(ChlA_100_mean),
         day0_sd = sd(ChlA_100_mean),
         flag_algae = ifelse(ChlA_100_mean - 2*day0_sd >= day0_mean, T, F))%>%
  ungroup()%>%
  select(Isolate, flag_algae)

data_4clus <- left_join(reader_data, algal_flag)%>%
  filter(flag_algae != TRUE)

#### Axenic Control Plot (Not Needed for Clustering) ####
AC_plot <- ggplot(data = AC, aes(x = avg_timeHours, y = ChlA_100_mean, color = host_species)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = ChlA_100_mean - ChlA_100_sd,
                    ymax = ChlA_100_mean + ChlA_100_sd),
                width = .2) +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 800)) +
  labs(x = "Time (hours)", y = "Chlorophyll A Fluorescence (RFU)")
#plot(AC_plot)

#### Cluster Optimization ####

cNums <- c(1:15)

cs_distances <- bind_rows(lapply(cNums, cluster_optimize, species = "chlorella", df = data_4clus))
cm_distances <- bind_rows(lapply(cNums, cluster_optimize, species = "coelastrum", df = data_4clus))
sa_distances <- bind_rows(lapply(cNums, cluster_optimize, species = "scenedesmus", df = data_4clus))
mm_distances <- bind_rows(lapply(cNums, cluster_optimize, species = "monoraphidium", df = data_4clus))
sm_distances <- bind_rows(lapply(cNums, cluster_optimize, species = "selenastrum", df = data_4clus))

#### Optimization Plots ####
# NOTE: Optimal Number is determined by looking at the plots and identifying the "Elbow"
# The plotting function requires an optimal number to plot a vertical line at that point
# If running script for the first time, guess any number 1:15 then refine after plot is made
plot_cs <- optimization_plot("chlorella", cs_distances, optimal_number = 6)
plot_cm <- optimization_plot("coelastrum", cm_distances, optimal_number = 4)
plot_sa <- optimization_plot("scenedesmus", sa_distances, optimal_number = 5)
plot_mm <- optimization_plot("monoraphidium", mm_distances, optimal_number = 5)
plot_sm <- optimization_plot("selenastrum", sm_distances, optimal_number = 4)

# These optimal numbers are fed into the clustering portion of the script for each species

#### Clustering and Cluster Plots ####
species_list <- c("chlorella", "coelastrum", "scenedesmus", "monoraphidium", "selenastrum")
optimals <- c(6, 4, 5, 5, 4) # Optimal Numbers from Optimization Portion

# Clustering Happens in this For Loop Using cluster() Function
cluster_data <- list()
avg_traj <- list()
for(i in 1:length(species_list)){
  species <- species_list[i]
  optimal_number <- optimals[i]
  cld <- cluster(species, optimal_number, data_4clus)
  cluster_data[[i]] <- data.frame(cld@clusters) %>% 
    rownames_to_column("Isolate")
  avg_traj[[i]] <- cld@trajMeans %>%
    mutate(host_species = species)
}

cluster_data_all <- bind_rows(cluster_data) %>% rename(Cluster = cld.clusters)
avg_traj_all <- bind_rows(avg_traj)

data_4plots <- left_join(reader_data_raw, algal_flag)%>%
  filter(Day_isolated != "Ctrl" & Outlier != TRUE & flag_algae == FALSE)%>%
  left_join(., stats_data)%>%
  group_by(Isolate, host_species, read_day, ccEffect, grEffect, annotation)%>% # trying this, nix annotation if necessary
  summarize(avg_timeHours = mean(read_timeHours, na.rm = TRUE),
            ChlA_100_mean = mean(ChlA_100, na.rm = TRUE),
            ChlA_100_sd = sd(ChlA_100, na.rm = TRUE))%>%
  ungroup()%>%
  full_join(., cluster_data_all)

data_4plots_split <- data_4plots %>%
  mutate(species_cluster = paste(host_species, Cluster, sep = "_"),
         annotation = ifelse(is.na(annotation), "no_annotation", annotation)) %>%
  split(., .$species_cluster)

cluster_plots <- lapply(data_4plots_split, cluster_plot)

plot_path <- "./supplement/cluster_plots/"
for (i in 1:length(cluster_plots)){
  name <- names(cluster_plots[i])
  file_name <- paste(name, "plot", sep = '_') %>%
    paste(., "png", sep = ".")
  ggsave(file_name, 
         plot = cluster_plots[[i]], 
         device = "png", 
         path = plot_path, 
         scale = 1, 
         width = 4.5, 
         height = 4, 
         dpi = 600, 
         units = "in", 
         limitsize = TRUE,
         bg = "transparent")
}

#### Facet Plots ####
data_4facet <- data_4plots %>%
  mutate(annotation = ifelse(is.na(annotation), "no_annotation", annotation),
         pnnEffect = ifelse((ccEffect == "Positive" & grEffect != "Negative")|(ccEffect != "Negative" & grEffect == "Positive"),
                            "Positive",
                            ifelse((ccEffect == "Negative" & grEffect != "Positive")|(ccEffect != "Positive" & grEffect == "Negative"),
                                   "Negative",
                                   ifelse(ccEffect == "Not Significant" & grEffect == "Not Significant",
                                          "No Significant Effect", "Inconclusive")))) %>%
  split(., .$host_species)

facet_plots <- lapply(data_4facet, facet_plot)

for(i in 1:length(facet_plots)){
  name <- names(facet_plots[i])
  file_name <- paste(name, "facet_plot", sep = '_') %>%
    paste(., "png", sep = ".")
  ggsave(file_name, 
         plot = facet_plots[[i]], 
         device = "png", 
         path = plot_path, 
         scale = 1, 
         width = 13.5, 
         height = 8, 
         dpi = 600, 
         units = "in", 
         limitsize = TRUE,
         bg = "transparent")
}
