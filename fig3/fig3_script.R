## This script pulls plate reader data from the "Away From Home" series of Coculture experiments
## This script ultimately produces figure 3, and uses area under the growth curve to determine statistical differences in algal growth.

#### Libraries ####

library(tidyverse)
library(lubridate)
library(DescTools)
library(broom)
library(ggpubr)
library(growthrates)
# for Foreach/doParallel
library(foreach)
library(doParallel)
numCores <- detectCores()
registerDoParallel(numCores-1)

#### Read in Data From Master ####
colors <- read.csv("./csv_files/colors.csv")|>
  mutate(speciesCode = hostShort)|>
  select(speciesCode, color)

tax <- read.csv("./csv_files/collection_tax_data.csv")

extract_dfs <- function(x, full_results) { 
  tryCatch(
    expr = {
      num <- x
      name <- names(full_results[num]) 
      df <- as.data.frame(full_results[[num]]$coefficients) |>
        mutate(sample = name) |>
        dplyr::select(sample, `Std. Error`, Estimate) |>
        rename(growthrate = Estimate,
               std_error = `Std. Error`) |>
        slice(2) }, 
    finally = {
      return(df)
    }
  )
}

## Grubbs Test checks for outliers and gives statistical basis for exclusion
## Used just before normalization
grubbs <- function(df){
  data_grubbs <- split(df, df$Isolate, df$new_host)
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

#### Plate Data #####
afh_plate_list <- list.files(path = "./raw_data/AFH_data", pattern = "^Plate_AW(\\d\\d|\\d)_Day_(\\d\\d|\\d).csv",
                         full.names = TRUE, ignore.case = TRUE)

read_plate_afh <- function(afh_plate_list) {
  read_csv(afh_plate_list, show_col_types = F)|>
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
  afh_plate_list |> 
  map_df(~read_plate_afh(.))|>
  rename(ChlA_100 = "Mean RFU [ChlA_100:460,685]")|>
  select(Well, ChlA_100, read_time, read_day, plate_no)

#### Map Data ####
afh_map_list <- list.files(path = "./raw_data/AFH_data",
                       pattern = "^plate_AW(\\d\\d|\\d)_map.csv",
                       full.names = TRUE, ignore.case = TRUE)

read_map_afh <- function(afh_map_list) {
  read_csv(afh_map_list, col_types = cols(.default = col_character())) |>
    mutate(Sample = as.character(Sample),
           plate_no = str_extract(str_extract(afh_map_list, pattern = "plate_AW(\\d\\d|\\d)"),pattern = regex("(\\d\\d|\\d)")), 
           filename = afh_map_list) |>
    select(-filename)
}

afh_mapDataAll <- 
  afh_map_list |>
  map_df(~read_map_afh(.))


#### Key Files ####
isolate_key <- read_csv("./raw_data/AFH_data/plate_map_isolate_codes.csv")

species_key <- read_csv("./raw_data/AFH_data/plate_map_species_codes.csv")

#### Assemble Data ####
afhData <- inner_join(afh_plateDataAll, afh_mapDataAll,
                   by = c("Well", "plate_no")) |>
  mutate(Sample = ifelse(Sample == "MC", "MC_MC", Sample))|>
  separate(Sample, c("sampleCode", "speciesCode"), remove = F)|>
  left_join( isolate_key, by = "sampleCode") |>
  left_join( species_key, by = "speciesCode") |>
  mutate(plateCode = paste(plate_no, speciesCode, sep = "_"),
         collectionCode = ifelse(sampleCode == "AC", "AC", collectionCode),
         hostCode = ifelse(sampleCode == "AC", speciesCode, hostCode),
         sample_exact = paste("plate", plate_no, "well", Well, Sample, sep = "_"),
         afh = ifelse(sampleCode == "AC", "Axenic", 
                      ifelse(sampleCode == "MC", "Media",
                             ifelse(speciesCode == hostCode, "Native", "Non-Native"))))|>
  # calc time replacement
  group_by(sample_exact)|>
  mutate(begin = min(read_time),
         read_interval = begin %--% read_time,
         read_timeHours = as.numeric(as.duration((read_interval)/dhours(1))))|>
  select(-begin)|>
  group_by(plate_no, Sample, read_day)|>
  mutate(day_mean = mean(ChlA_100, na.rm=TRUE),
         day_sd = sd(ChlA_100, na.rm=TRUE))|>
  ungroup()|>
  mutate(read_day = as.numeric(read_day))|>
  filter(read_day <= 12) %>%
  arrange(read_day)
 
#### Calc AUC ####
gr_chlA <- afhData |>
  filter(collectionCode != "MC") |>
  all_easylinear(ChlA_100 ~ read_timeHours | sample_exact)
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
growthrates <- inner_join(results_chlA, std_err_chlA, by = c("sample_exact" = "sample")) |>
  dplyr::select(sample_exact, growthrate,
                r2, std_error) |>
  mutate(growthrate = ifelse(is.na(growthrate), 0, growthrate))

afhData_split <- split(afhData, afhData$sample_exact)

coefList <- foreach(i = 1:length(afhData_split), .packages = "tidyverse") %dopar%{
  tryCatch({
    df <- afhData_split[[i]]
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
    
    coefDF <- data.frame(sample_exact = unique(df$sample_exact),
                         asymptote = aMod,
                         intercept = iMod,
                         growthParam = xMod)
    return(coefDF)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

iKey <- afhData |>
  dplyr::select(sample_exact, collectionCode, speciesCode, hostCode, plate_no, afh)|>
  distinct()

coefs <- bind_rows(coefList)|>
  left_join(iKey)|>
  left_join(growthrates) |>
  mutate(growthrate = ifelse(is.na(growthrate), growthParam, growthrate))

ACcoefs <- coefs |> filter(collectionCode == "AC")

coefList_4t <- split(coefs, list(coefs$collectionCode, coefs$speciesCode, coefs$plate_no), drop = T)#remove empty lists from resulting split dataframe

statsList <- foreach(i = 1:length(coefList_4t), .packages = c("tidyverse", "broom")) %dopar%{
  tryCatch({
    tt_df <- coefList_4t[[i]]
    new_host <- unique(tt_df$speciesCode)
    orig_host <- unique(tt_df$hostCode)
    plate <- unique(tt_df$plate_no)
    Isolate <- unique(tt_df$collectionCode)
    afh <- unique(tt_df$afh)
    ac_df <- ACcoefs |> filter(hostCode == new_host & plate_no == plate)
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
             new_host = new_host,
             orig_host = orig_host,
             afh = afh)
    return(stats_df)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

stats <- bind_rows(statsList)|>
  dplyr::select(Isolate, new_host, orig_host, plate_no, pCC_greater, pCC_less, pGR_greater, pGR_less)|>
  mutate(ccEffect = ifelse(pCC_greater <= 0.05 & pCC_less > 0.05, "Increased CC",
                           ifelse(pCC_greater > 0.05 & pCC_less <= 0.05, "Decreased CC",
                                  ifelse(pCC_greater > 0.05 & pCC_less > 0.05,
                                         "No Significant CC Change", "Error"))),
         grEffect = ifelse(pGR_greater <= 0.05 & pGR_less > 0.05, "Increased GR",
                           ifelse(pGR_greater > 0.05 & pGR_less <= 0.05, "Decreased GR",
                                  ifelse(pGR_greater > 0.05 & pGR_less > 0.05,
                                         "No Significant GR Change", "Error"))),
         annotation = paste(ccEffect, grEffect, sep = "\n"))

library(outliers)
stats_coefs <- left_join(stats, coefs, by = c("Isolate" = "collectionCode", 
                                              "new_host" = "speciesCode",
                                              "orig_host" = "hostCode",
                                              "plate_no"))|>
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
  distinct(Isolate, plate_no, .keep_all = T)
#Calculate area under the curve for a timeseries for each well on each plate.
aucData_split <- as.list(1:length(afhData_split)) 

for(i in 1:length(afhData_split)){
  df <- afhData_split[[i]]
  auc <- AUC(x = df$read_time, y = df$ChlA_100, method = "spline")
  aucData_split[[i]] <- data.frame("sample_exact" = df$sample_exact,
                              "Sample" = df$Sample,
                              "sampleCode" = df$sampleCode,
                              "plateCode" = df$plateCode,
                              "auc" = auc)|>
    distinct()
}

#Combine area under the curve data into one list.
aucData <- bind_rows(aucData_split)

#### AUC T-Test ####
ctrl_aucData <- filter(aucData, sampleCode == "AC")

sample_aucData_split <- filter(aucData, sampleCode != "AC", sampleCode != "MC") %>%
  split( .$Sample)

tTestData_split <- as.list(1:length(sample_aucData_split))
for(i in 1:length(sample_aucData_split)){
  df <- sample_aucData_split[[i]]
  plate <- unique(df$plateCode)
  sample <- unique(df$Sample)
  ctrl <- filter(ctrl_aucData, plateCode == plate)
  stat_greater <- tidy(t.test(df$auc, ctrl$auc, alternative = "greater"))|>
    mutate(Sample = sample,
           p_greater = p.value) |>
    select(Sample, p_greater)
  stat_less <- tidy(t.test(df$auc, ctrl$auc, alternative = "less"))|>
    mutate(Sample = sample,
           p_less = p.value) |>
    select(Sample, p_less)
  tTestData_split[[i]] <- full_join(stat_greater, stat_less)|>
    mutate(Effect = ifelse(p_greater >= 0.05 & p_less >= 0.05, "Not Significant",
                    ifelse(p_greater <= 0.05, "Positive",
                    ifelse(p_less <= 0.05, "Negative",
                    ifelse(p_greater <= 0.05 & p_less <= 0.05, "Error", NA)))))
}

tTestData <- bind_rows(tTestData_split)

#### Normalization ####
acData_auc <- group_by(ctrl_aucData, plateCode)|>
  mutate(auc_mean = mean(auc),
         auc_sd = sd(auc))|>
  ungroup()|>
  select(Sample, sampleCode, plateCode, auc_mean, auc_sd)|>
  distinct()

acData <- filter(afhData, sampleCode == "AC")|>
  select(Sample, sampleCode, plateCode, read_day, read_timeHours, day_mean, day_sd)|>
  distinct()|>
  left_join(acData_auc)|>
  mutate(acDay_mean = day_mean,
         acDay_sd = day_sd,
         acAuc_mean = auc_mean,
         acAuc_sd = auc_sd)|>
  select(plateCode, read_day, read_timeHours, acDay_mean, acDay_sd, acAuc_mean, acAuc_sd)

afhData_all <- left_join(afhData, aucData)|>
  left_join( acData)|>
  mutate(ChlA_100norm = ChlA_100/acDay_mean,
         aucNorm = auc/acAuc_mean)
#Rename to current naming conventions to match up to taxonomy data
sam.dataMiseqNames <- afhData_all |>
  mutate(Isolate = str_replace(collectionCode, ",|\\.", "point"),
         Isolate = str_replace(collectionCode, "DF", "D31"),
         seq_type = "Miseq") 

sam.dataJinny <- afhData_all |> 
  filter(str_detect(collectionCode, "DF") & hostCode == "CS") |>
  #Match up experimental collectionCodes to the ones Jinny Sequenced that are missing
  mutate(Isolate = case_when(collectionCode == "10DF" ~ "S10OMO",
                             collectionCode == "11DF" ~ "S11",
                             collectionCode == "12DF" ~ "S12",
                             collectionCode == "13DF" ~ "S13",
                             collectionCode == "14DF" ~ "S14",
                             collectionCode == "16DF" ~ "S16",
                             collectionCode == "17DF" ~ "S17A",
                             collectionCode == "18DF" ~ "S18",
                             collectionCode == "19DF" ~ "S19AB",
                             collectionCode == "2DF" ~ "S2",
                             collectionCode == "20DF" ~ "S20",
                             collectionCode == "21DF" ~ "S21",
                             collectionCode == "22DF" ~ "S22",
                             collectionCode == "24DF" ~ "S24",
                             collectionCode == "25DF" ~ "S25",
                             collectionCode == "26DF" ~ "S26",
                             collectionCode == "27DF" ~ "S27",
                             collectionCode == "28DF" ~ "S28",
                             collectionCode == "29DF" ~ "S29",
                             collectionCode == "3DF" ~ "S3",
                             collectionCode == "30,1DF" ~ "S30W",
                             collectionCode == "31DF" ~ "S31",
                             collectionCode == "32DF" ~ "S32AB",
                             collectionCode == "33DF" ~ "S33",
                             collectionCode == "34DF" ~ "S34",
                             collectionCode == "31DF" ~ "S31",
                             collectionCode == "5DF" ~ "S5",
                             collectionCode == "6DF" ~ "S6",
                             collectionCode == "7DF" ~ "S7",
                             collectionCode == "8DF" ~ "S8",
                             collectionCode == "9DF" ~ "S9W"
  ),
  seq_type = "Miseq") |>
  filter(!is.na(collectionCode))

sam.dataSangerNames <- afhData_all |>
  mutate(seq_type = "Sanger",
         Isolate = str_replace(collectionCode, "D3", "_D3"),
         Isolate = str_replace(collectionCode, ",|\\.", "point"),
         Isolate = str_replace(collectionCode, "DF", "_DF"),
         Isolate = str_replace(collectionCode, "D31", "_D31"))

sam.dataDF <- rbind(sam.dataMiseqNames, sam.dataSangerNames, sam.dataJinny)

sam.data <- sam.dataDF |>
  filter(! collectionCode == "AC") |>
  group_by(sample_exact) |>
  distinct(Isolate, .keep_all = T)

renamed_data <- afhData_all |> 
  left_join(sam.data) |>
  left_join(tax)

#### Plots ####
## Faceted Bar Plot
### Idea for this plot: bars colored by algal host, stat/effect with annotation above/below bar
bar_chart <- afhData_all|>
  filter(sampleCode != "AC", sampleCode != "MC")|>
  group_by(Sample, Well)|>
  distinct(Sample, .keep_all =T) |>
  mutate(aucLogNorm_mean = log(mean(aucNorm)),
         aucLogNorm_sd = log(sd(aucNorm))|>
  left_join( tax)|>
  select(Sample, collectionCode, facet_title, algalSpecies, afh, plateCode,
         aucLogNorm_mean, aucLogNorm_sd)|>
  distinct()|>
  left_join( tTestData)|>
  # ggplot( aes(x = algalSpecies, y = aucLogNorm_mean, fill = Effect, color = afh))+
  ggplot( aes(x = algalSpecies, y = aucLogNorm_mean, fill = afh))+ ## use this for native/non-native fill
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = aucLogNorm_mean-(2*(aucLogNorm_sd/sqrt(3))), ymax = aucLogNorm_mean+((2*aucLogNorm_sd/sqrt(3)))))+ ## error bars are currently 2* std dev, to make std error divide by sqrt # of sample size (3)
  # ylim(-0.5, 1)+ ## Sets scale limits
  coord_flip()+ ## Flips coordinates to make axes more readable
  facet_wrap(vars(collectionCode))
  
#bar_chart
### NOTE: using facet_title as the facet var will group by assigned genus
### using collectionCode will group by isolate number


## Line Plot
line_plot <- afhData_all|>
  filter(sampleCode != "AC", sampleCode != "MC")|>
  mutate(logNormRFU = log(ChlA_100norm))|>
  left_join( tax)|>
  left_join( colors)|>
  # group_by(Sample, read_day)|>
  # mutate()
  select(Sample, collectionCode, facet_title, read_day, algalSpecies, color,
         plateCode, afh, read_timeHours, ChlA_100norm, logNormRFU)|>
  left_join( tTestData)|>
  ggplot( aes(x = read_timeHours, y = logNormRFU, color = algalSpecies, linetype = afh))+
  geom_hline(yintercept = 0)+
  # geom_point()+
  geom_smooth(alpha = 0)+
  # facet_wrap(vars(facet_title))+
  facet_wrap(vars(collectionCode))+
  scale_color_manual(values = colors$color)

#line_plot
### NOTE: using facet_title as the facet var will group by assigned genus
### using collectionCode will group by isolate number


