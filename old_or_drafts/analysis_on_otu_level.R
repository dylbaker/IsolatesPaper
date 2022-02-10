#### Data Needed for Older Analysis ####


library(tidyverse)

setwd("C:/Users/jlaue/OneDrive/Documents/Denef_Lab_Stuff/R_scripts_git")

isolateTax <- read.csv("./CoCulture_data/collection_all_cultures.csv")%>%
  select(-X)

isolateKey <- read.csv("./CoCulture_data/isolateKey.csv")%>%
  select(-X)%>%
  filter(Isolate != "MC", Isolate != "AC")

readData <- read.csv("./CoCulture_data/coculture_data.csv")%>%
  select(-X)

axenic_data <- filter(readData, Isolate == "AC")%>%
  select(Isolate, host_species, read_day, read_timeHours, ChlA_100, Outlier)%>%
  filter(Outlier == F)%>%
  group_by(Isolate, host_species, read_day)%>%
  mutate(mean_rfu = mean(ChlA_100),
         mean_readTime = mean(read_timeHours))%>%
  ungroup()%>%
  mutate(numOTUs = "Axenic Control")
  
all <- full_join(isolateTax, isolateKey, by = "Isolate")%>%
  left_join(., readData)%>%
  mutate(numOTUs = as_factor(numOTUs))

mixed <- filter(all, mixed == T)%>%
  mutate(numOTUs = as.integer(numOTUs))

otu_list <- split(mixed, mixed$otu)

#### OTU plots ####

plot_otu <- function(df){
  host <- unique(df$host_species)
  otu <-  unique(df$otu)
  genus <- unique(df$genus)
  species <- unique(df$species)
  plotTitle <- paste(host, otu, genus, species, sep = " ")
  filename <- paste("./otu_plots/", host, sep = "")%>%
    paste(., otu, "plot.png", sep = "_")
  
  axenic <- filter(axenic_data, host_species == host)
  
  plot <- ggplot(df, aes(x = read_timeHours, y = ChlA_100, shape = numOTUs, color = Isolate, linetype = numOTUs))+
    geom_point()+
    geom_smooth(alpha = 0)+
    geom_smooth(data = axenic, aes(x = read_timeHours, y = ChlA_100), color = "green", linetype = 1)+
    # geom_point(data = axenic, aes(x = mean_readTime, y = mean_rfu), color = "black", shape = 1)+
    geom_point(data = axenic, aes(x = read_timeHours, y = ChlA_100), color = "black", shape = 1)+
    labs(title = plotTitle)+
    facet_wrap(vars(numOTUs))
  
  # color = Isolate,
  
  ggsave(filename, plot, device = "png")
  
  return(plot)
}

# otu_plots <- lapply(otu_list, plot_otu)

#### Chlorella ####
chlorella_otus <- filter(all, host_species == "chlorella")%>%
  split(., .$otu)

chlorella_otu_plots <- lapply(chlorella_otus, plot_otu)

#### Coelastrum ####
coelastrum_otus <- filter(all, host_species == "coelastrum")%>%
  split(., .$otu)

coelastrum_otu_plots <- lapply(coelastrum_otus, plot_otu)

#### Scenedesmus ####
scenedesmus_otus <- filter(all, host_species == "scenedesmus")%>%
  split(., .$otu)

scenedesmus_otu_plots <- lapply(scenedesmus_otus, plot_otu)

#### Monoraphidium ####
monoraphidium_otus <- filter(all, host_species == "monoraphidium")%>%
  split(., .$otu)

monoraphidium_otu_plots <- lapply(monoraphidium_otus, plot_otu)

#### Selenastrum ####
selenastrum_otus <- filter(all, host_species == "selenastrum")%>%
  split(., .$otu)

selenastrum_otu_plots <- lapply(selenastrum_otus, plot_otu)
