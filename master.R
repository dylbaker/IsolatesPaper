## Master Script for CoCulture Manuscript

## This script pulls plate reader data and mothur outputs into R and combines them to make a variety of data files used to build figures

#### Libraries ####

## Definitely Needed
# Data Processing
library(tidyverse)
library(lubridate)
library(growthrates)
library(rstatix)

# Phyloseq Processing
library(phyloseq)

# Logistic Mod
library(broom)

#Area Under Curve
library(DescTools)

# for Foreach/doParallel
library(foreach)
library(doParallel)
numCores <- detectCores()
registerDoParallel(numCores-1)

# Needed for Grubbs Test (used for normalization)
library(outliers)

# ## Potentially Needed
# # If I make the tree Figs in Phyloseq (potentially used for fig 4)
# library(ape) ## Possibly Unnecessary
# library(patchwork) ## Possibly Unnecessary
# library(ggpubr) ## Possibly Unnecessary
# library(ggtree) ## Possibly Unnecessary

#### Functions ####
## For Data Processing
read_plate <- function(plate_list) {
  read_csv(plate_list, show_col_types = F)|>
    rename(ChlA_100 = `Mean RFU [ChlA_100:460,685]`, 
           read_time = `Reading Date/Time`)|>
    #Extract plate information and species from each file to aid in isolate mapping.
    mutate(plate_no = str_extract(str_extract(plate_list,
                                              pattern = "Plate_(\\d\\d|\\d)"),
                                  pattern = regex("(\\d\\d|\\d)")), 
           host_species = str_extract(plate_list,
                                      pattern = "(chlorella|coelastrum|scenedesmus|monoraphidium|selenastrum)"),
           filename = plate_list,
           read_day = str_extract(str_extract(plate_list,
                                              pattern = "(\\d\\d|\\d).csv"),
                                  pattern = regex("(\\d\\d|\\d)")),
           read_time = as.POSIXct(read_time, 
                                  format = "%Y-%m-%d %H:%M:%S")) |>
    dplyr::select(Well, ChlA_100, read_time, read_day, plate_no, host_species)
}

read_map <- function(map_list) {
  read_csv(map_list, col_types = cols(.default = col_character())) |>
    mutate(Isolate = as.character(Isolate),
           plate_no = str_extract(str_extract(map_list,
                                              pattern = "Plate_(\\d\\d|\\d)"),
                                  pattern = regex("(\\d\\d|\\d)")), 
           host_species = str_extract(map_list,
                                      pattern = "(chlorella|coelastrum|scenedesmus|monoraphidium|selenastrum)"),
           filename = map_list) |>
    dplyr::select(-filename)
}

extract_dfs <- function(x, full_results) { 
  tryCatch(
   expr = {
    num <- x
    name <- names(full_results[num]) 
    df <- as.data.frame(full_results[[num]]$coefficients) %>%
      mutate(sample = name) %>%
      dplyr::select(sample, `Std. Error`, Estimate) %>%
      rename(growthrate = Estimate,
           std_error = `Std. Error`) %>%
      slice(2) }, 
  finally = {
  return(df)
  }
  )
}

## Grubbs Test checks for outliers and gives statistical basis for exclusion
## Used just before normalization
grubbs <- function(df){
  data_grubbs <- split(df, df$host_species)
  for(i in 1:length(data_grubbs)){
    x <- data_grubbs[[i]]
    data_grubbs[[i]] <- x |>
      mutate(p_high = grubbs.test(x$asymptote, type = 10)[[3]],
             p_low = grubbs.test(x$asymptote, type = 10, opposite = T)[[3]],
             maxVal = max(x$asymptote),
             minVal = min(x$asymptote),
             outlier = ifelse(p_high <= 0.05 & maxVal == asymptote, T,
                              ifelse(p_low <= 0.05 & minVal == asymptote, T, F)))
    # |>
    #   filter(outlier == F)|>
    #   mutate(p_high = grubbs.test(x$asymptote, type = 10)[[3]],
    #          p_low = grubbs.test(x$asymptote, type = 10, opposite = T)[[3]],
    #          maxVal = max(x$asymptote),
    #          minVal = min(x$asymptote),
    #          outlier = ifelse(p_high <= 0.05 & maxVal == asymptote, T,
    #                           ifelse(p_low <= 0.05 & minVal == asymptote, T, F)))
  }
  df_final <- bind_rows(data_grubbs)|>
    dplyr::select(-p_high, -p_low, -maxVal, -minVal)
  return(df_final)
}

#### Set Up Colors ####
color_df <- data.frame(hostLong = c("chlorella", "coelastrum", "scenedesmus",
                                    "monoraphidium", "selenastrum"),
                       hostShort = c("CS", "CM", "MM", "SA", "SC"),
                       color = c("#7FC97F", "#BEAED4", "#E0115F", "#FDC086", "#386CB0"))
#### Read in Plate Reader Data and Combine with Plate Maps ####
## Plate Reader Data
folder_list <- list.dirs(path = "./raw_data", full.names = TRUE, recursive = FALSE) %>%
  ifelse(str_detect(. , pattern = regex("plate_(\\d|\\d\\d)_(chlorella|coelastrum|scenedesmus|monoraphidium|selenastrum)$", ignore_case = TRUE)), .,NA) %>%
  na.omit()

plate_list <- list.files(path = folder_list, 
                         pattern = "^Plate_(\\d\\d|\\d)_Day_(\\d\\d|\\d).csv",
                         full.names = TRUE, 
                         ignore.case = TRUE)

plate_data_all <- map_df(plate_list, ~read_plate(.))

## Plate Map Data (for Plate Reader)
map_list <- list.files(path = folder_list,
                       pattern = "^plate_(\\d\\d|\\d)_map.csv",
                       full.names = TRUE,
                       ignore.case = TRUE)

map_data_all <- map_df(map_list, ~read_map(.))


## Combine Plate and Plate Map Data
data <- inner_join(plate_data_all, map_data_all,
                   by = c("Well", "plate_no", "host_species"))|>
  mutate(exact_isolate = paste("Plate", plate_no, Isolate, Well, sep = "_"),
         Day_isolated = ifelse(str_detect(Isolate, "D3")==TRUE, "D3",
                               ifelse(str_detect(Isolate, "DF")==TRUE, "D31", "Ctrl")))|>
  group_by(exact_isolate)|>
  mutate(begin = min(read_time),
         read_interval = begin %--% read_time,
         read_timeHours = as.numeric(as.duration((read_interval)/dhours(1))))|>
  dplyr::select(-begin)|>
  group_by(plate_no, Isolate, read_day)|>
  mutate(day_mean = mean(ChlA_100, na.rm=TRUE),
         day_sd = sd(ChlA_100, na.rm=TRUE))|>
  ungroup()|>
# This chunk flags outliers - defined as any point more than 2 std. dev. from day mean
  mutate(Outlier = ifelse((day_mean-(2*day_sd))>ChlA_100|(day_mean+(2*day_sd))<ChlA_100,
                          TRUE, FALSE))|>
  dplyr::select(Well, host_species, plate_no, Day_isolated, 
         Isolate, exact_isolate, read_day, read_time,
         read_timeHours, ChlA_100, day_mean, day_sd, Outlier)

#Evaporation coefficient for the plates
evap_coef <- 60/200/max(data$read_timeHours)

#Correct for plate evaporation
data <- data |>
  mutate(ChlA_100 = (1-evap_coef*read_timeHours)*ChlA_100) |>
  arrange(read_timeHours)

#### Statistical Analysis of Plate Reader Data ####
algal_species <- unique(data$host_species)

## data_flagged appends a T/F flag to any wells that demonstrate unreasonably high Chlorophyll A fluorescence on day 0 (defined as more than 2* the axenic control wells) as these wells likely contain algal contaminants or cyanobacteria that will disrupt fluorescence measurement throughout the experiment
data_flagged <- as.list(1:length(algal_species))

for(i in 1:length(algal_species)){
  host <- algal_species[i]
  host_data <- data |> filter(host_species == host)
  AC_day0_avg <- filter(host_data, Isolate == "AC" & read_day == 0)|>
    summarise(ChlA_100_mean = mean(ChlA_100))|>
    as.numeric()
  flags <- filter(host_data, read_day == 0)|>
    mutate(algae_flag = ifelse(ChlA_100 > 2*AC_day0_avg, T, F))|>
    dplyr::select(Isolate, algae_flag)|>
    distinct()
  data_flagged[[i]] <- left_join(host_data, flags)
}

gr_chlA <- data |>
  filter(Isolate != "MC") |>
  all_easylinear(ChlA_100 ~ read_timeHours | exact_isolate)

results_full_chlA <- summary(gr_chlA) %>%
  discard(~ length(.x) == 3)

results_chlA<- results(gr_chlA) |>
  rename(growthrate = mumax)

std_err_num <- as.list(1:length(results_full_chlA)) 

std_err_chlA <- lapply(std_err_num, extract_dfs,
                       full_results = results_full_chlA) |>
  bind_rows() |>
  mutate(yplus = growthrate + std_error,
         yminus = growthrate - std_error) |>
  dplyr::select(-growthrate)

#Add standard error to the dataframe
growthrates <- inner_join(results_chlA, std_err_chlA, by = c("exact_isolate" = "sample")) |>
  dplyr::select(exact_isolate, growthrate,
         r2, std_error) |>
  mutate(growthrate = ifelse(is.na(growthrate), 0, growthrate))

####AUC, GR, CC Stats####

data_split <- bind_rows(data_flagged) |>
  mutate(Day_isolated = ifelse(Isolate == "AC", "AC", Day_isolated))|>
  filter(algae_flag != T,
         Day_isolated != "Ctrl") %>%
  split(., .$exact_isolate)

aucData_split <- as.list(1:length(data_split))
for(i in 1:length(data_split)){
  df <- data_split[[i]]
  auc <- AUC(x = df$read_time, y = df$ChlA_100, method = "spline")
  auc <- as.numeric(auc)
  aucData_split[[i]] <- data.frame("exact_isolate" = df$exact_isolate,
                                   "Isolate" = df$Isolate,
                                   "host_species" = df$host_species,
                                   "plate_no" = df$plate_no,
                                   "auc" = auc)|>
    distinct()
}
#Combine area under the curve data into one list.
aucData <- bind_rows(aucData_split)

#### AUC T-Test ####
ctrl_aucData <- filter(aucData, Isolate == "AC")
sample_aucData <- filter(aucData, Isolate != "AC")

acAUC_outliers <- ctrl_aucData |>
  group_by(host_species, plate_no) |>
  identify_outliers(auc)

ctrl_aucData <- anti_join(ctrl_aucData, acAUC_outliers, by = "auc")

sample_aucData_split <- filter(aucData, Isolate != "AC", Isolate != "MC") %>%
  split( .$Isolate)

tTestData_split <- foreach(i = 1:length(sample_aucData_split), .packages = c("tidyverse","broom")) %dopar% {
  tryCatch({
  df <- sample_aucData_split[[i]]
  plate <- unique(df$plate_no)
  sample <- unique(df$Isolate)
  ctrl <- filter(ctrl_aucData, plate_no == plate)
  stat_greater <- tidy(t.test(df$auc, ctrl$auc, alternative = "greater"))|>
    mutate(Isolate = sample,
           p_greater = p.value) |>
    select(Isolate, p_greater)
  stat_less <- tidy(t.test(df$auc, ctrl$auc, alternative = "less"))|>
    mutate(Isolate = sample,
           p_less = p.value) |>
    select(Isolate, p_less)
  tTestData <- full_join(stat_greater, stat_less)|>
    mutate(Effect = case_when((p_greater > 0.05 & p_less > 0.05) ~ "Not Significant",
                              p_greater <= 0.05 ~ "Positive",
                              p_less <= 0.05 ~ "Negative",
                              T ~ "Error"))
  return(tTestData)
  
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

tTestData <- bind_rows(tTestData_split)

#### AUC Normalization ####
AC_auc_mean <- ctrl_aucData |>
  group_by(host_species,plate_no) |>
  mutate(acAuc_mean = mean(auc)) |>
  select(host_species, plate_no, acAuc_mean, acAuc_se) |>
  distinct()

auc_data_norm <- sample_aucData |>
  group_by(Isolate, host_species, plate_no) |>
  left_join(AC_auc_mean)|>
  mutate(auc_norm = auc/acAuc_mean,
         log_auc_norm = log(auc_norm),
         se_lan = sd(log_auc_norm)/sqrt(n())) |>
  distinct(Isolate, .keep_all = T) |>
  filter(!is.na(se_lan)) |>
  left_join(tTestData) 

coefList <- foreach(i = 1:length(data_split), .packages = "tidyverse") %dopar%{
  tryCatch({
    df <- data_split[[i]]
    aPara <- max(df$ChlA_100)+ 0.05*max(df$ChlA_100)
    parameters <- coef(lm(qlogis(ChlA_100/aPara) ~ read_timeHours, data = df))
    iPara <- parameters[[1]]
    xPara <- parameters[[2]]
    
    mod <- nls(ChlA_100 ~ A/(1 + exp(-(I + X * read_timeHours))),
               start = list(A = aPara, I = iPara, X = xPara),
               data = df, trace = T, nls.control(maxiter = 100, warnOnly = T))
    modSummary <- summary(mod)
    
    aMod <- coef(mod)[1]
    iMod <- coef(mod)[2]
    xMod <- coef(mod)[3]
    
    coefDF <- data.frame(exact_isolate = unique(df$exact_isolate),
                         asymptote = aMod,
                         intercept = iMod,
                         growthParam = xMod)
    return(coefDF)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

iKey <- data |>
  dplyr::select(exact_isolate, Isolate, host_species, plate_no)|>
  distinct()

coefs <- bind_rows(coefList)|>
  left_join(iKey)|>
  left_join(growthrates) |>
  mutate(growthrate = ifelse(is.na(growthrate), growthParam, growthrate))

ACcoefs <- coefs |> filter(Isolate == "AC")

coefList_4t <- split(coefs, list(coefs$Isolate, coefs$plate_no), drop = T)#remove empty lists from resulting split dataframe

statsList <- foreach(i = 1:length(coefList_4t), .packages = c("tidyverse", "broom")) %dopar%{
  tryCatch({
    tt_df <- coefList_4t[[i]]
    host <- unique(tt_df$host_species)
    plate <- unique(tt_df$plate_no)
    Isolate <- unique(tt_df$Isolate)
    ac_df <- ACcoefs |> filter(host_species == host & plate_no == plate)
    cc_greater <- tidy(t.test(tt_df$asymptote, 
                              ac_df$asymptote, 
                              alternative = "greater"))|>
      mutate(pCC_greater = p.value, Isolate = Isolate)|>
      dplyr::select(pCC_greater, Isolate)
    cc_less <- tidy(t.test(tt_df$asymptote, 
                           ac_df$asymptote, 
                           alternative = "less"))|>
      mutate(pCC_less = p.value, Isolate = Isolate)|>
      dplyr::select(pCC_less, Isolate)
    gr_greater <- tidy(t.test(tt_df$growthrate, 
                              ac_df$growthrate, 
                              alternative = "greater"))|>
      mutate(pGR_greater = p.value, Isolate = Isolate)|>
      dplyr::select(pGR_greater, Isolate)
    gr_less <- tidy(t.test(tt_df$growthrate, 
                           ac_df$growthrate, 
                           alternative = "less"))|>
      mutate(pGR_less = p.value, Isolate = Isolate)|>
      dplyr::select(pGR_less, Isolate)
    cc_stat <- full_join(cc_greater, cc_less)
    gr_stat <- full_join(gr_greater, gr_less)
    stats_df <- full_join(cc_stat, gr_stat)|>
      mutate(plate_no = plate,
             host_species = host)
    return(stats_df)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

stats <- bind_rows(statsList)|>
  dplyr::select(Isolate, host_species, plate_no, pCC_greater, pCC_less, pGR_greater, pGR_less)|>
  mutate(ccEffect = case_when(pCC_greater <= 0.05 & pCC_less > 0.05 ~ "Positive",
                              pCC_greater > 0.05 & pCC_less <= 0.05 ~ "Negative",
                              pCC_greater > 0.05 & pCC_less > 0.05 ~ "Not Significant",
                              T ~ "Error"),
         grEffect = case_when(pGR_greater <= 0.05 & pGR_less > 0.05 ~ "Positive",
                              pGR_greater > 0.05 & pGR_less <= 0.05 ~ "Negative",
                              pGR_greater > 0.05 & pGR_less > 0.05 ~ "Not Significant",
                              T ~ "Error"),
         annotation = paste(ccEffect,"CC", grEffect,"GR" ,sep = " "))

stats_coefs <- full_join(stats, coefs)|>
  grubbs() # This test added here when normalization was added 01/19/2022 -JL

#### Normalization of Plate Reader Data to Axenic Controls ####
ac_coefs <- stats_coefs |> 
  filter(Isolate == "AC" & outlier != T)|>
  group_by(host_species, plate_no)|>
  mutate(nSamples = n(),
         meanCC = mean(asymptote),
         meanGR = mean(growthrate))|>
  ungroup()|>
  # These columns are empty, as there is no statistics data for axenic controls samples
  dplyr::select(-pCC_greater, -pCC_less, -pGR_greater, -pGR_less,
         -ccEffect, -grEffect, -annotation)
# running grubbs() again doesnt make much sense to me
# |>grubbs()

stats_coefs_split <- stats_coefs |>
  filter(Isolate != "AC" & Isolate != "MC")|>
  mutate(split_var = paste(Isolate, plate_no, sep = "_"))%>%
  split(., .$split_var)

for(i in 1:length(stats_coefs_split)){
  df <- stats_coefs_split[[i]]
  plate <- unique(df$plate_no)
  host <- unique(df$host_species)
  ac <- ac_coefs |>
    filter(host_species == host & plate_no == plate & asymptote >= 0)|>
    summarize(acCC = mean(asymptote),
              acGR = mean(growthrate))
  acCC <- as.numeric(ac$acCC)
  acGR <- as.numeric(ac$acGR)
  stats_coefs_split[[i]] <- df |>
    mutate(normCC = asymptote/acCC,
           normGR = growthrate/acGR,
           logNormCC = log(normCC),
           logNormGR = log(normGR))
}
stats_normCoefs <- bind_rows(stats_coefs_split)

#### Calculate Triplicate Mean of Stats Coefficients ####
stats_meanCoefs <- stats_normCoefs |>
  group_by(Isolate, host_species, plate_no) |>
  mutate(n = n(),
    meanCC = mean(asymptote),
    seCC = sd(asymptote) / sqrt(n),
    meanGR = mean(growthrate),
    seGR = sd(growthrate) / sqrt(n),
    mean_normCC = mean(normCC),
    mean_normGR = mean(normGR),
    mean_logNormCC = mean(logNormCC),
    mean_logNormGR = mean(logNormGR)
  ) |>
  filter(n >= 3) |>
  distinct(Isolate, plate_no, .keep_all = T) |>
  left_join(auc_data_norm) |>
  ungroup()

#### Read in Mothur Outputs and Combine with Coculture Data and Stats ####
# Define mothur outputs
list.file <- "./mothur_outputs/all_seqs/final.merge.asv.list"
constax.file <- "./mothur_outputs/all_seqs/final.merge.asv.ASV.cons.taxonomy"
shared.file <- "./mothur_outputs/all_seqs/final.merge.asv.shared"

# Read mothur outputs into phyloseq/R
mo.data <- import_mothur(mothur_list_file = list.file,
                         mothur_constaxonomy_file = constax.file,
                         mothur_shared_file = shared.file)

# Convert stats data into phyloseq sample data (syntax of Isolate names changed to match mothur sample name syntax)
duplicate_stats <- split(stats_meanCoefs, duplicated(stats_meanCoefs$Isolate) | duplicated(stats_meanCoefs$Isolate, fromLast = TRUE), 
                         drop == "FALSE")
  
duplicates <-  rbind(duplicate_stats[["TRUE"]]) %>%
  filter(!Isolate == "AC")
  
sam.dataMiseqNames <- stats_meanCoefs |>
  mutate(Isolate = str_replace(Isolate, ",|\\.", "point"),
         Isolate = str_replace(Isolate, "DF", "D31"),
         seq_type = "Miseq") 

sam.dataJinny <- stats_meanCoefs |> 
  filter(str_detect(Isolate, "DF") & host_species == "chlorella") |>
  #Match up experimental isolates to the ones Jinny Sequenced that are missing
  mutate(Isolate = case_when(Isolate == "10DF" ~ "S10OMO",
                             Isolate == "11DF" ~ "S11",
                             Isolate == "12DF" ~ "S12",
                             Isolate == "13DF" ~ "S13",
                             Isolate == "14DF" ~ "S14",
                             Isolate == "16DF" ~ "S16",
                             Isolate == "17DF" ~ "S17A",
                             Isolate == "18DF" ~ "S18",
                             Isolate == "19DF" ~ "S19AB",
                             Isolate == "2DF" ~ "S2",
                             Isolate == "20DF" ~ "S20",
                             Isolate == "21DF" ~ "S21",
                             Isolate == "22DF" ~ "S22",
                             Isolate == "24DF" ~ "S24",
                             Isolate == "25DF" ~ "S25",
                             Isolate == "26DF" ~ "S26",
                             Isolate == "27DF" ~ "S27",
                             Isolate == "28DF" ~ "S28",
                             Isolate == "29DF" ~ "S29",
                             Isolate == "3DF" ~ "S3",
                             Isolate == "30,1DF" ~ "S30W",
                             Isolate == "31DF" ~ "S31",
                             Isolate == "32DF" ~ "S32AB",
                             Isolate == "33DF" ~ "S33",
                             Isolate == "34DF" ~ "S34",
                             Isolate == "31DF" ~ "S31",
                             Isolate == "5DF" ~ "S5",
                             Isolate == "6DF" ~ "S6",
                             Isolate == "7DF" ~ "S7",
                             Isolate == "8DF" ~ "S8",
                             Isolate == "9DF" ~ "S9W"
                             ),
         seq_type = "Miseq") |>
  filter(!is.na(Isolate))

sam.dataSangerNames <- stats_meanCoefs |>
  mutate(seq_type = "Sanger",
         Isolate = str_replace(Isolate, "D3", "_D3"),
         Isolate = str_replace(Isolate, ",|\\.", "point"),
         Isolate = str_replace(Isolate, "DF", "_DF"),
         Isolate = str_replace(Isolate, "D31", "_D31"))

sam.dataDF <- rbind(sam.dataMiseqNames, sam.dataSangerNames, sam.dataJinny)

sam.data <- sam.dataDF |>
  filter(! Isolate == "199point1_D3" & ! Isolate == "AC") |>
  group_by(Isolate, seq_type) |>
  distinct(Isolate) |>
  column_to_rownames(var = "Isolate")|>
  sample_data(.)

# verify that sample names match (there will be more samples in mo.data, because we included environmental samples in our analysis) - most important here is that syntax is the same (underscores between number and day, ex: 23D3 becomes 23_D3, 30,1DF becomes 30_1DF)
sample_names(mo.data)
sample_names(sam.data)
#diffobj::diffChr(sample_names(mo.data), sample_names(sam.data))

# change taxonomic "rank names" to recognizable classifications (rank 1-7 to KPCOFGS)
rank_names(mo.data)
colnames(tax_table(mo.data)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
rank_names(mo.data)

# full phyloseq object that combines the tax data from mothur with the coculture data
all.data <- merge_phyloseq(mo.data, sam.data)

#### Classify Isolates as Pure or Mixed Cultures ####
taxTable <- as.data.frame(all.data@tax_table)|>
  rownames_to_column(var = "asv")

sampleData <- sam.dataDF |>
  filter(! Isolate %in% c('199point1_D3','11_DF','S15','S16','S20','S25','34_DF','S3','6_DF','7_DF','8_DF','S9W') &
  ! Isolate == "AC") |>
  distinct(Isolate, .keep_all = TRUE)

asvTable <- as.data.frame(all.data@otu_table)|>
  rownames_to_column(var = "asv")|>
  pivot_longer(!asv, names_to = "Isolate", values_to = "count")|>
  filter(count != 0 & 
           ! Isolate %in% c('199point1_D3','11_DF','S15','S16','S20','S25','34_DF','S3','6_DF','7_DF','8_DF','S9W'))|>
  group_by(Isolate)|>
  mutate(totalReads = sum(count),
         asvMatches = n(),
         contamFlag = ifelse(count <= 0.1*totalReads, T, F))|> # 10% contam threshold
         #contamFlag = ifelse(count <= 0.05*totalReads, T, F))|> #5% contam threshold
  ungroup()

isolateTax <- asvTable |> 
## This filter pulls out only sample data (removes control, background, and pond asvs)
  filter(contamFlag == F, str_detect(Isolate, "^[:digit:]|^S|^_"))|> 
  # distinct()|> ## not needed
  group_by(Isolate)|>
  mutate(numasvs = n(),
         mixed = ifelse(numasvs > 1, T, F))|>
  left_join(sampleData)|>
  left_join(taxTable)

#### Write CSV Files ####
## Save Colors for Future Figures
write.csv(color_df, "./csv_files/colors.csv", row.names = F)

## Save Plate Reader Data joined to map Data
write.csv(data, "./csv_files/coculture_data.csv", row.names = F)

## Save All Isolate Info (Not Exactly Needed)
# write.csv(iKey, file = "./csv_files/isolateKey.csv", row.names = F)

## Just the Coefficients (includes controls)
# write.csv(coefs, file = "./csv_files/logistic_mod_coefs.csv", row.names = F)

## Just the Stats (does not include controls)
# write.csv(stats, file = "./csv_files/logistic_mod_stats.csv", row.names = F)

## Save all Stats and Coefficients (raw coefficiencts, not mean), Includes Controls
write.csv(stats_coefs, file = "./csv_files/logistic_mod_stats_coefs.csv", row.names = F)

## Save all Stats and Mean/Std. Error Coefficients, Does not include Controls
write.csv(stats_meanCoefs, file = "./csv_files/logistic_mod_stats_meanCoefs.csv", row.names = F)

## Save all Stats, Coefficients (raw) and normalized/log(normalized) coefficients. Does not Include controls
write.csv(stats_normCoefs, file = "./csv_files/logistic_mod_stats_normCoefs.csv")

## Save all data including taxonomy, impact on growth outcomes (p values, normalizations), host, plate, EVERYTHING
write.csv(isolateTax, file = "./csv_files/collection_tax_data.csv", row.names = F)

## Save the taxonomic information for only mixed cultures
mixed_cultures <- isolateTax |> filter(mixed == T)
write.csv(mixed_cultures, file = "./csv_files/collection_tax_data_mixedOnly.csv", row.names = F)

## Save the taxonomic information for only pure cultures
pure_cultures <- isolateTax |> filter(mixed == F)
write.csv(pure_cultures, file = "./csv_files/collection_tax_data_pureOnly.csv", row.names = F)

## Save the taxonomic information for the pond/natural community Data
#write.csv(pondTax, file = "./csv_files/natCom_tax_data.csv", row.names = F)

## Save the taxonomic information for the phycosphere data
#write.csv(phycosphereTax, file = "./csv_files/phycosphere_tax_data.csv", row.names = F)
