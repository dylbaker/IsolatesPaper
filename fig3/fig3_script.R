## This script pulls plate reader data from the "Away From Home" series of Coculture experiments
## This script ultimately produces figure 3, and uses area under the growth curve to determine statistical differences in algal growth.

#### Directory Set ####

setwd("C:/Users/jlaue/OneDrive/Documents/Denef_Lab_Stuff/IsolatesPaper/")
getwd()

#### Libraries ####

library(tidyverse)
library(lubridate)
library(DescTools)
library(broom)

#### Read in Data From Master ####
colors <- read.csv("./csv_files/colors.csv")%>%
  mutate(speciesCode = hostShort)%>%
  select(speciesCode, color)

tax <- read.csv("./csv_files/collection_tax_data.csv")%>%
  mutate(collectionCode = Isolate,
         facet_title = paste(Genus, sep = " "))%>%
  select(collectionCode, facet_title, mixed, Genus)

#### Plate Data #####
afh_plate_list <- list.files(path = "./raw_data/AFH_data", pattern = "^Plate_AW(\\d\\d|\\d)_Day_(\\d\\d|\\d).csv",
                         full.names = TRUE, ignore.case = TRUE)

read_plate_afh <- function(afh_plate_list) {
  read_csv(afh_plate_list, show_col_types = F)%>%
    mutate(plate_no = str_extract(str_extract(afh_plate_list,
                                              pattern = "Plate_AW(\\d\\d|\\d)"),
                                  pattern = regex("(\\d\\d|\\d)")),
           filename = afh_plate_list,
           read_day = str_extract(str_extract(afh_plate_list,
                                              pattern = "(\\d\\d|\\d).csv"),
                                  pattern = regex("(\\d\\d|\\d)")),
           read_time = as.POSIXct(`Reading Date/Time`, 
                                  tryFormats = c("%m/%d/%Y %H:%M",
                                                 "%m/%d/%Y %H:%M:%S")))
}

afh_plateDataAll <-
  afh_plate_list %>% 
  map_df(~read_plate_afh(.))%>%
  rename(ChlA_100 = "Mean RFU [ChlA_100:460,685]")%>%
  select(Well, ChlA_100, read_time, read_day, plate_no)

#### Map Data ####
afh_map_list <- list.files(path = "./raw_data/AFH_data",
                       pattern = "^plate_AW(\\d\\d|\\d)_map.csv",
                       full.names = TRUE, ignore.case = TRUE)

read_map_afh <- function(afh_map_list) {
  read_csv(afh_map_list, col_types = cols(.default = col_character())) %>%
    mutate(Sample = as.character(Sample),
           plate_no = str_extract(str_extract(afh_map_list, pattern = "plate_AW(\\d\\d|\\d)"),pattern = regex("(\\d\\d|\\d)")), 
           filename = afh_map_list) %>%
    select(-filename)
}

afh_mapDataAll <- 
  afh_map_list %>%
  map_df(~read_map_afh(.))

#### Key Files ####
isolate_key <- read_csv("./raw_data/AFH_data/plate_map_isolate_codes.csv")

species_key <- read_csv("./raw_data/AFH_data/plate_map_species_codes.csv")

#### Assemble Data ####
afhData <- inner_join(afh_plateDataAll, afh_mapDataAll,
                   by = c("Well", "plate_no")) %>%
  mutate(Sample = ifelse(Sample == "MC", "MC_MC", Sample))%>%
  separate(Sample, c("sampleCode", "speciesCode"), remove = F)%>%
  left_join(., isolate_key, by = "sampleCode") %>%
  left_join(., species_key, by = "speciesCode") %>%
  mutate(plateCode = paste(plate_no, speciesCode, sep = "_"),
         collectionCode = ifelse(sampleCode == "AC", "AC", collectionCode),
         hostCode = ifelse(sampleCode == "AC", speciesCode, hostCode),
         sample_exact = paste("plate", plate_no, "well", Well, Sample, sep = "_"),
         afh = ifelse(sampleCode == "AC", "Axenic", 
                      ifelse(sampleCode == "MC", "Media",
                             ifelse(speciesCode == hostCode, "Native", "Non-Native"))))%>%
  # rename(Sample = Isolate)%>%
  # calc time replacement
  group_by(sample_exact)%>%
  mutate(begin = min(read_time),
         read_interval = begin %--% read_time,
         read_timeHours = as.numeric(as.duration((read_interval)/dhours(1))))%>%
  select(-begin)%>%
  group_by(plate_no, Sample, read_day)%>%
  mutate(day_mean = mean(ChlA_100, na.rm=TRUE),
         day_sd = sd(ChlA_100, na.rm=TRUE))%>%
  ungroup()%>%
  # calc time replacement end
  # calc_time()%>%
  mutate(read_day = as.numeric(read_day))%>%
  filter(read_day <= 12)
 
#### Calc AUC ####
afhData_split <- split(afhData, afhData$sample_exact) 

aucData_split <- as.list(1:length(afhData_split))
for(i in 1:length(afhData_split)){
  df <- afhData_split[[i]]
  auc <- AUC(x = df$read_time, y = df$ChlA_100, method = "spline")
  aucData_split[[i]] <- data.frame("sample_exact" = df$sample_exact,
                              "Sample" = df$Sample,
                              "sampleCode" = df$sampleCode,
                              "plateCode" = df$plateCode,
                              "auc" = auc)%>%
    distinct()
}
aucData <- bind_rows(aucData_split)

#### AUC T-Test ####
ctrl_aucData <- filter(aucData, sampleCode == "AC")

sample_aucData_split <- filter(aucData, sampleCode != "AC", sampleCode != "MC")%>%
  split(., .$Sample)

tTestData_split <- as.list(1:length(sample_aucData_split))
for(i in 1:length(sample_aucData_split)){
  df <- sample_aucData_split[[i]]
  plate <- unique(df$plateCode)
  sample <- unique(df$Sample)
  ctrl <- filter(ctrl_aucData, plateCode == plate)
  stat_greater <- tidy(t.test(df$auc, ctrl$auc, alternative = "greater"))%>%
    mutate(Sample = sample,
           p_greater = p.value) %>%
    select(Sample, p_greater)
  stat_less <- tidy(t.test(df$auc, ctrl$auc, alternative = "less"))%>%
    mutate(Sample = sample,
           p_less = p.value) %>%
    select(Sample, p_less)
  tTestData_split[[i]] <- full_join(stat_greater, stat_less)%>%
    mutate(Effect = ifelse(p_greater >= 0.05 & p_less >= 0.05, "Not Significant",
                    ifelse(p_greater <= 0.05, "Positive",
                    ifelse(p_less <= 0.05, "Negative",
                    ifelse(p_greater <= 0.05 & p_less <= 0.05, "Error", NA)))))
}

tTestData <- bind_rows(tTestData_split)

#### Normalization ####
acData_auc <- group_by(ctrl_aucData, plateCode)%>%
  mutate(auc_mean = mean(auc),
         auc_sd = sd(auc))%>%
  ungroup()%>%
  select(Sample, sampleCode, plateCode, auc_mean, auc_sd)%>%
  distinct()

acData <- filter(afhData, sampleCode == "AC")%>%
  select(Sample, sampleCode, plateCode, read_day, read_timeHours, day_mean, day_sd)%>%
  distinct()%>%
  left_join(acData_auc)%>%
  mutate(acDay_mean = day_mean,
         acDay_sd = day_sd,
         acAuc_mean = auc_mean,
         acAuc_sd = auc_sd)%>%
  select(plateCode, read_day, read_timeHours, acDay_mean, acDay_sd, acAuc_mean, acAuc_sd)

afhData_all <- left_join(afhData, aucData)%>%
  left_join(., acData)%>%
  mutate(ChlA_100norm = ChlA_100/acDay_mean,
         aucNorm = auc/acAuc_mean)

#### Plots ####
## Faceted Bar Plot
### Idea for this plot: bars colored by algal host, stat/effect with annotation above/below bar
bar_chart <- afhData_all%>%
  filter(sampleCode != "AC", sampleCode != "MC")%>%
  group_by(Sample)%>%
  mutate(aucLogNorm_mean = mean(log(aucNorm)),
         aucLogNorm_sd = sd(log(aucNorm)))%>%
  left_join(., tax)%>%
  select(Sample, collectionCode, facet_title, algalSpecies, afh, plateCode,
         aucLogNorm_mean, aucLogNorm_sd)%>%
  distinct()%>%
  left_join(., tTestData)%>%
  ggplot(., aes(x = algalSpecies, y = aucLogNorm_mean, fill = Effect, color = afh))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = aucLogNorm_mean-(2*aucLogNorm_sd), ymax = aucLogNorm_mean+(2*aucLogNorm_sd)))+
  facet_wrap(vars(facet_title))
bar_chart

ggsave("./fig3/AFH_experiments/barChart.png", plot = bar_chart, device = "png",
       scale = 1, width = 16, height = 8, dpi = 300, units = "in", limitsize = TRUE)

## Line Plot
line_plot <- afhData_all%>%
  filter(sampleCode != "AC", sampleCode != "MC")%>%
  mutate(logNormRFU = log(ChlA_100norm))%>%
  left_join(., tax)%>%
  left_join(., colors)%>%
  # group_by(Sample, read_day)%>%
  # mutate()
  select(Sample, collectionCode, facet_title, read_day, algalSpecies, color,
         plateCode, afh, read_timeHours, ChlA_100norm, logNormRFU)%>%
  left_join(., tTestData)%>%
  ggplot(., aes(x = read_timeHours, y = logNormRFU, color = algalSpecies, linetype = afh))+
  geom_hline(yintercept = 0)+
  # geom_point()+
  geom_smooth(alpha = 0)+
  facet_wrap(vars(facet_title))+
  scale_color_manual(values = colors$color)
  # facet_wrap(vars(collectionCode))
line_plot

ggsave("./fig3/AFH_experiments/linePlot_geomSmooth.png", plot = line_plot, device = "png",
       scale = 1, width = 16, height = 8, dpi = 300, units = "in", limitsize = TRUE)
