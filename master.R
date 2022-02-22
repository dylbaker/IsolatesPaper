## Master Script for CoCulture Manuscript

## This script pulls plate reader data and mothur outputs into R and combines them to make a variety of data files used to build figures

#### Directory Set ####

setwd("C:/Users/jlaue/OneDrive/Documents/Denef_Lab_Stuff/IsolatesPaper/")
getwd()

#### Libraries ####

## Definitely Needed
# Data Processing
library(tidyverse)
library(lubridate)

# Phyloseq Processing
library(phyloseq)

# Logistic Mod
library(broom)

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
  read_csv(plate_list, show_col_types = F)%>%
    rename(ChlA_100 = `Mean RFU [ChlA_100:460,685]`, 
           read_time = `Reading Date/Time`)%>%
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
                                  format = "%Y-%m-%d %H:%M:%S")) %>%
    select(Well, ChlA_100, read_time, read_day, plate_no, host_species)
}

read_map <- function(map_list) {
  read_csv(map_list, col_types = cols(.default = col_character())) %>%
    mutate(Isolate = as.character(Isolate),
           plate_no = str_extract(str_extract(map_list,
                                              pattern = "Plate_(\\d\\d|\\d)"),
                                  pattern = regex("(\\d\\d|\\d)")), 
           host_species = str_extract(map_list,
                                      pattern = "(chlorella|coelastrum|scenedesmus|monoraphidium|selenastrum)"),
           filename = map_list) %>%
    select(-filename)
}

## Grubbs Test checks for outliers and gives statistical basis for exclusion
## Used just before normalization
grubbs <- function(df){
  data_grubbs <- split(df, df$host_species)
  for(i in 1:length(data_grubbs)){
    x <- data_grubbs[[i]]
    data_grubbs[[i]] <- x %>%
      mutate(p_high = grubbs.test(x$asymptote, type = 10)[[3]],
             p_low = grubbs.test(x$asymptote, type = 10, opposite = T)[[3]],
             maxVal = max(x$asymptote),
             minVal = min(x$asymptote),
             outlier = ifelse(p_high <= 0.05 & maxVal == asymptote, T,
                              ifelse(p_low <= 0.05 & minVal == asymptote, T, F)))
    # %>%
    #   filter(outlier == F)%>%
    #   mutate(p_high = grubbs.test(x$asymptote, type = 10)[[3]],
    #          p_low = grubbs.test(x$asymptote, type = 10, opposite = T)[[3]],
    #          maxVal = max(x$asymptote),
    #          minVal = min(x$asymptote),
    #          outlier = ifelse(p_high <= 0.05 & maxVal == asymptote, T,
    #                           ifelse(p_low <= 0.05 & minVal == asymptote, T, F)))
  }
  df_final <- bind_rows(data_grubbs)%>%
    select(-p_high, -p_low, -maxVal, -minVal)
  return(df_final)
}


#### Set Up Colors ####
color_df <- data.frame(hostLong = c("Chlorella", "Coelastrum", "Monoraphidium",
                                    "Scenedesmus", "Selenastrum"),
                       hostShort = c("CS", "CM", "MM", "SA", "SC"),
                       color = c("#f8766d", "#a3a500", "#00bf7d", "#00b0f6", "#e76bf3"))
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
                   by = c("Well", "plate_no", "host_species"))%>%
  mutate(exact_isolate = paste("Plate", plate_no, Isolate, Well, sep = "_"),
         Day_isolated = ifelse(str_detect(Isolate, "D3")==TRUE, "D3",
                               ifelse(str_detect(Isolate, "DF")==TRUE, "D31", "Ctrl")))%>%
  group_by(exact_isolate)%>%
  mutate(begin = min(read_time),
         read_interval = begin %--% read_time,
         read_timeHours = as.numeric(as.duration((read_interval)/dhours(1))))%>%
  select(-begin)%>%
  group_by(plate_no, Isolate, read_day)%>%
  mutate(day_mean = mean(ChlA_100, na.rm=TRUE),
         day_sd = sd(ChlA_100, na.rm=TRUE))%>%
  ungroup()%>%
# This chunk flags outliers - defined as any point more than 2 std. dev. from day mean
  mutate(Outlier = ifelse((day_mean-(2*day_sd))>ChlA_100|(day_mean+(2*day_sd))<ChlA_100,
                          TRUE, FALSE))%>%
  select(Well, host_species, plate_no, Day_isolated, 
         Isolate, exact_isolate, read_day, read_time,
         read_timeHours, ChlA_100, day_mean, day_sd, Outlier)

#### Statistical Analysis of Plate Reader Data ####
algal_species <- unique(data$host_species)

## data_flagged appends a T/F flag to any wells that demonstrate unreasonably high Chlorophyll A fluorescence on day 0 (defined as more than 2* the axenic control wells) as these wells likely contain algal contaminants or cyanobacteria that will disrupt fluorescence measurement throughout the experiment
data_flagged <- as.list(1:length(algal_species))
for(i in 1:length(algal_species)){
  host <- algal_species[i]
  host_data <- data %>% filter(host_species == host)
  AC_day0_avg <- filter(host_data, Isolate == "AC" & read_day == 0)%>%
    summarise(ChlA_100_mean = mean(ChlA_100))%>%
    as.numeric()
  flags <- filter(host_data, read_day == 0)%>%
    mutate(algae_flag = ifelse(ChlA_100 > 2*AC_day0_avg, T, F))%>%
    select(Isolate, algae_flag)%>%
    distinct()
  data_flagged[[i]] <- left_join(host_data, flags)
}

data_split <- bind_rows(data_flagged) %>%
  mutate(Day_isolated = ifelse(Isolate == "AC", "AC", Day_isolated))%>%
  filter(algae_flag != T,
         Day_isolated != "Ctrl")%>%
  split(., .$exact_isolate)

coefList <- foreach(i = 1:length(data_split), .packages = "tidyverse") %dopar%{
  tryCatch({
    df <- data_split[[i]]
    aPara <- max(df$ChlA_100)+20
    parameters <- coef(lm(qlogis(ChlA_100/aPara) ~ read_timeHours, data = df))
    iPara <- parameters[[1]]
    xPara <- parameters[[2]]
    
    mod <- nls(ChlA_100 ~ A/(1 + exp(-(I + X * read_timeHours))),
               start = list(A = aPara, I = iPara, X = xPara),
               data = df, trace = T, nls.control(maxiter = 100, warnOnly = T))
    # modSummary <- summary(mod)
    
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

iKey <- data %>%
  select(exact_isolate, Isolate, host_species, plate_no)%>%
  distinct()

coefs <- bind_rows(coefList)%>%
  left_join(., iKey)%>%
  ungroup()

ACcoefs <- coefs %>% filter(Isolate == "AC")
coefList_4t <- split(coefs, coefs$Isolate)

statsList <- foreach(i = 1:length(coefList_4t), .packages = c("tidyverse", "broom")) %dopar%{
  tryCatch({
    tt_df <- coefList_4t[[i]]
    host <- unique(tt_df$host_species)
    plate <- unique(tt_df$plate_no)
    Isolate <- unique(tt_df$Isolate)
    ac_df <- ACcoefs %>% filter(host_species == host & plate_no == plate)
    cc_greater <- tidy(t.test(tt_df$asymptote, 
                              ac_df$asymptote, 
                              alternative = "greater"))%>%
      mutate(pCC_greater = p.value, Isolate = Isolate)%>%
      select(pCC_greater, Isolate)
    cc_less <- tidy(t.test(tt_df$asymptote, 
                           ac_df$asymptote, 
                           alternative = "less"))%>%
      mutate(pCC_less = p.value, Isolate = Isolate)%>%
      select(pCC_less, Isolate)
    gr_greater <- tidy(t.test(tt_df$growthParam, 
                              ac_df$growthParam, 
                              alternative = "greater"))%>%
      mutate(pGR_greater = p.value, Isolate = Isolate)%>%
      select(pGR_greater, Isolate)
    gr_less <- tidy(t.test(tt_df$growthParam, 
                           ac_df$growthParam, 
                           alternative = "less"))%>%
      mutate(pGR_less = p.value, Isolate = Isolate)%>%
      select(pGR_less, Isolate)
    cc_stat <- full_join(cc_greater, cc_less)
    gr_stat <- full_join(gr_greater, gr_less)
    stats_df <- full_join(cc_stat, gr_stat)%>%
      mutate(plate_no = plate,
             host_species = host)
    return(stats_df)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

stats <- bind_rows(statsList)%>%
  select(Isolate, host_species, plate_no, pCC_greater, pCC_less, pGR_greater, pGR_less)%>%
  mutate(ccEffect = ifelse(pCC_greater <= 0.05 & pCC_less > 0.05, "Increased CC",
                           ifelse(pCC_greater > 0.05 & pCC_less <= 0.05, "Decreased CC",
                                  ifelse(pCC_greater > 0.05 & pCC_less > 0.05,
                                         "No Significant CC Change", "Error"))),
         grEffect = ifelse(pGR_greater <= 0.05 & pGR_less > 0.05, "Increased GR",
                           ifelse(pGR_greater > 0.05 & pGR_less <= 0.05, "Decreased GR",
                                  ifelse(pGR_greater > 0.05 & pGR_less > 0.05,
                                         "No Significant GR Change", "Error"))),
         annotation = paste(ccEffect, grEffect, sep = "\n"))

stats_coefs <- full_join(stats, coefs)%>%
  grubbs() # This test added here when normalization was added 01/19/2022 -JL

#### Normalization of Plate Reader Data to Axenic Controls ####
ac_coefs <- stats_coefs %>% 
  filter(Isolate == "AC" & outlier != T)%>%
  group_by(host_species, plate_no)%>%
  mutate(nSamples = n(),
         meanCC = mean(asymptote),
         meanGR = mean(growthParam))%>%
  ungroup()%>%
  # These columns are empty, as there is no statistics data for axenic controls samples
  select(-pCC_greater, -pCC_less, -pGR_greater, -pGR_less,
         -ccEffect, -grEffect, -annotation)
# running grubbs() again doesnt make much sense to me
# %>%grubbs()

stats_coefs_split <- stats_coefs %>%
  filter(Isolate != "AC" & Isolate != "MC")%>%
  mutate(split_var = paste(Isolate, plate_no, sep = "_"))%>%
  split(., .$split_var)
for(i in 1:length(stats_coefs_split)){
  df <- stats_coefs_split[[i]]
  plate <- unique(df$plate_no)
  host <- unique(df$host_species)
  ac <- ac_coefs %>%
    filter(host_species == host & plate_no == plate & asymptote >= 0)%>%
    summarize(acCC = mean(asymptote),
              acGR = mean(growthParam))
  acCC <- as.numeric(ac$acCC)
  acGR <- as.numeric(ac$acGR)
  stats_coefs_split[[i]] <- df %>%
    mutate(normCC = asymptote/acCC,
           normGR = growthParam/acGR,
           logNormCC = log(normCC),
           logNormGR = log(normGR))
}
stats_normCoefs <- bind_rows(stats_coefs_split)

#### Calculate Triplicate Mean of Stats Coefficients ####
stats_meanCoefs <- coefs %>%
  group_by(Isolate, host_species, plate_no)%>%
  mutate(n = n())%>%
  summarize(meanCC = mean(asymptote),
            seCC = sd(asymptote)/sqrt(n),
            meanGR = mean(growthParam),
            seGR = sd(growthParam)/sqrt(n))%>%
  # NOTE: need to write in normalization code --> figure this out Tuesday
            # ,
            # mean_normCC = mean(normCC),
            # mean_normGR = mean(normGR),
            # mean_logNormCC = mean(logNormCC),
            # mean_logNormGR = mean(logNormGR))%>%
  distinct()%>%
  left_join(stats, .)
# full_join(stats, .) # keeps AC and other problem isolates (things that got rmvd from stats)
# %>% drop_na() # removes outliers and includes only the isolates included in the stats df
# "probs" contains everything that is removed by the use of left_join() -or- drop_na()
# probs <- filter(stats_meanCoefs, (Isolate =="105D3"|Isolate =="114D3"|Isolate =="129D3"|Isolate =="199,1D3"|Isolate =="224D3"|Isolate =="94D3"|Isolate =="96D3"|Isolate =="AC"))

#### Read in Mothur Outputs and Combine with Coculture Data and Stats ####
# Define mothur outputs
## NOTE: These are the outputs from the old pipeline
# tree.file <- "./raw_data/mothur_outputs/old_pipeline_outputs/fullTree.phylip.tre"
# constax.file <- "./raw_data/mothur_outputs/old_pipeline_outputs/fullTree.cons.taxonomy"
# shared.file <- "./raw_data/mothur_outputs/old_pipeline_outputs/fullTree.an.shared"
# list.file <- "./raw_data/mothur_outputs/old_pipeline_outputs/fullTree.an.list"

#### NEW MOTHUR OUTPUTS ####
# list.file <- "./raw_data/mothur_outputs/final.an.list"
# constax.file <- "./raw_data/mothur_outputs/final.an.0.01.cons.taxonomy"
# shared.file <- "./raw_data/mothur_outputs/final.an.shared"
# tree.file <- "./raw_data/mothur_outputs/final.an.unique.rep.otuRename.phylip.tre"
# 
# list.file <- "./raw_data/mothur_outputs/final.opti_mcc.list"
# constax.file <- "./raw_data/mothur_outputs/final.opti_mcc.0.03.cons.taxonomy"
# shared.file <- "./raw_data/mothur_outputs/final.opti_mcc.shared"
# tree.file <- "./raw_data/mothur_outputs/final.opti_mcc.0.03.rep.otuRename.phylip.tre"

list.file <- "./raw_data/mothur_outputs/final.an.list"
constax.file <- "./raw_data/mothur_outputs/final.an.unique.cons.taxonomy"
shared.file <- "./raw_data/mothur_outputs/final.an.shared"
tree.file <- "./raw_data/mothur_outputs/final.an.unique.rep.otuRename.phylip.tre"
# tree.file <- "./raw_data/mothur_outputs/final.an.jclass.unique.tre"

# Read mothur outputs into phyloseq/R
mo.data <- import_mothur(mothur_list_file = list.file,
                         mothur_constaxonomy_file = constax.file,
                         mothur_shared_file = shared.file,
                         mothur_tree_file = tree.file)
tree.data <- read_tree(treefile = tree.file)

# Convert stats data into phyloseq sample data (syntax of Isolate names changed to match mothur sample name syntax)
sam.dataMiseqNames <- stats_meanCoefs %>%
  mutate(Isolate = str_replace(Isolate, ",", "sep"),
         Isolate = str_replace(Isolate, "(?<!sep\\d)D", "sepD"))

sam.dataSangerNames <- sam.dataMiseqNames %>%
  mutate(Isolate = paste("_", Isolate, sep=""))

sam.dataPondNames <- data.frame(Isolate = c("P1sepD0sepR1sep022um", "P1sepD0sepR2sep022um",
                                            "P1sepD0sepR3sep022um", "P2sepD0sepR1sep022um",
                                            "P2sepD0sepR2sep022um", "P2sepD0sepR3sep022um",
                                            "P3sepD0sepR1sep022um", "P3sepD0sepR2sep022um",
                                            "P3sepD0sepR3sep022um"))%>%
  mutate(host_species = "naturalCommunity",
         plate_no = NA,
         pCC_greater = NA,
         pCC_less = NA,
         pGR_greater = NA,
         pGR_less = NA,
         ccEffect = NA,
         grEffect = NA,
         annotation = NA,
         meanCC = NA,
         seCC = NA,
         meanGR = NA,
         seGR = NA)

sam.data <- rbind(sam.dataMiseqNames, sam.dataSangerNames, sam.dataPondNames)%>%
  column_to_rownames(var = "Isolate")%>%
  sample_data(.)


# Here We Filter the Mothur outputs to only include samples of interest and pond (natural community) data
allSamples <- rbind(sam.dataMiseqNames, sam.dataSangerNames, sam.dataPondNames)%>%
  pull(Isolate)

# this pruning function should tell phyloseq to only return the samples we selected above
mo.dataPruned <- prune_samples(allSamples, mo.data)
sample_names(mo.dataPruned)

# verify that sample names match (there will be more samples in mo.data, because we included environmental samples in our analysis) - most important here is that syntax is the same (underscores between number and day, ex: 23D3 becomes 23_D3, 30,1DF becomes 30_1DF)
# sample_names(mo.data) <- str_replace(sample_names(mo.data), "_D", "sepD")
# sample_names(mo.data)
sample_names(sam.data)

# change taxonomic "rank names" to recognizable classifications (rank 1-7 to KPCOFGS)
rank_names(mo.dataPruned)
colnames(tax_table(mo.dataPruned)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
rank_names(mo.dataPruned)
# rank_names(mo.data)
# colnames(tax_table(mo.data)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
# rank_names(mo.data)

# plot a test tree using the mothur data (not needed right here)
# plot_tree(mo.data, "treeonly", ladderize = T)

# full phyloseq object that combines the tax data from mothur with the coculture data
all.data <- merge_phyloseq(mo.dataPruned, sam.data)
# all.data <- merge_phyloseq(mo.data, sam.data)

#### Classify Isolates as Pure or Mixed Cultures ####
taxTable <- as.data.frame(all.data@tax_table)%>%
  rownames_to_column(., var = "otu")

otuTable <- as.data.frame(all.data@otu_table)%>%
  rownames_to_column(., var = "otu")%>%
  pivot_longer(!otu, names_to = "Isolate", values_to = "count")%>%
  filter(count != 0)%>%
  group_by(Isolate)%>%
  mutate(totalReads = sum(count),
         otuMatches = n(),
         readType = ifelse(totalReads == 1, "sanger", "miseq"),
         contamFlag = ifelse(count <= 0.1*totalReads, T, F))%>% # 10% contam threshold
         # contamFlag = ifelse(count <= 0.05*totalReads, T, F))%>% #5% contam threshold
  ungroup()

isolateTax <- otuTable %>% 
## This filter pulls out only sample data (removes control, background, and pond OTUs)
  filter(contamFlag == F, str_detect(.$Isolate, "^[:digit:]|^S|^_"))%>% 
  # distinct()%>% ## not needed
  group_by(Isolate)%>%
  mutate(numOTUs = n(),
         mixed = ifelse(numOTUs > 1, T, F),
         Isolate = str_replace(Isolate, "(?<=\\d)sep(?=\\d)", ","), #these str_replace functions revert sample names to R naming convention (removes "sep", underscores, ect.)
         Isolate = str_replace(Isolate, "sep", ""),
         Isolate = str_replace(Isolate, "_", ""))%>%
  left_join(., taxTable)

#### Pull the Pond Data Only ####
pondOTUs <- as.data.frame(all.data@otu_table)%>%
  rownames_to_column(., var = "otu")%>%
  pivot_longer(!otu, names_to = "Isolate", values_to = "count")%>%
  filter(count != 0,
         Isolate == "P1sepD0sepR1sep022um" |
           Isolate == "P1sepD0sepR2sep022um"|
           Isolate == "P2sepD0sepR2sep022um"|
           Isolate == "P2sepD0sepR3sep022um"|
           Isolate == "P3sepD0sepR1sep022um"|
           Isolate == "P3sepD0sepR2sep022um"|
           Isolate == "P3sepD0sepR3sep022um")%>%
  group_by(Isolate)%>%
  # mutate(totalReads = sum(count),
  #        otuMatches = n())%>%
  ungroup()

pondTax <- pondOTUs %>%
  group_by(Isolate)%>%
  mutate(numOTUs = n(),
         mixed = ifelse(numOTUs > 1, T, F),
         Isolate = str_replace(Isolate, "(?<=\\d)sep(?=\\d)", ","), #these str_replace functions revert sample names to R naming convention (removes "sep", underscores, ect.)
         Isolate = str_replace(Isolate, "sep", ""),
         Isolate = str_replace(Isolate, "_", ""))%>%
  left_join(., taxTable)

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

## Save the taxonomic information for all cultures (mixed and pure)
write.csv(isolateTax, file = "./csv_files/collection_tax_data.csv", row.names = F)

## Save the taxonomic information for only mixed cultures
mixed_cultures <- isolateTax %>% filter(mixed == T)
write.csv(mixed_cultures, file = "./csv_files/collection_tax_data_mixedOnly.csv", row.names = F)

## Save the taxonomic information for only pure cultures
pure_cultures <- isolateTax %>% filter(mixed == F)
write.csv(pure_cultures, file = "./csv_files/collection_tax_data_pureOnly.csv", row.names = F)

## Save the taxonomic information for the pond/natural community Data
write.csv(pondTax, file = "./csv_files/natCom_tax_data.csv", row.names = F)

