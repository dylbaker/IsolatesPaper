### Eliminated/Old Fxns

calc_time <- function(df){
  data <- df%>%
    mutate(exact_isolate = paste("Plate", plate_no, Isolate, Well, sep = "_"))
  day_0 <- data %>%
    arrange(read_day)%>% 
    distinct(exact_isolate, .keep_all = TRUE)%>%
    mutate(begin = read_time) %>%
    select(exact_isolate, begin)
  data_with_time <- left_join(data, day_0)%>%
    mutate(read_interval = begin %--% read_time,
           read_timeHours = as.numeric(as.duration(read_interval)/dhours(1)))%>%
    select(-begin, -read_interval)
  return(data_with_time)
}

day_Isolated <- function(df){
  data <- df%>%
    mutate(Day_isolated = ifelse(str_detect(Isolate, "D3")==TRUE, "D3", ifelse(str_detect(Isolate, "DF")==TRUE, "D31", "Ctrl")))
  return(data)
}


## from phyloseq_processing
# stats_annotation <- select(stats, -host_species, -plate_no)
# # stats_annotation <- read.csv("logistic_mod_stats.csv")%>%
# #   select(-X, -host_species, -plate_no)
# 
# ## Pure cultures database info with these selections is just mean, se, meanNorm, meanLogNorm, coefficients
# sample_key <- read.csv("pureCultures_databaseInfo.csv")%>%
#   select(-kingdom, -phylum, -class, -order, -family, -genus, -species)%>%
#   distinct(.)%>%
#   left_join(., stats_annotation)%>%
#   mutate(Isolate = str_replace(Isolate, ",", "_"),
#          Isolate = str_replace(Isolate, "(?<!\\_\\d)D", "\\_D"))%>%
#   column_to_rownames(var = "Isolate")


## Function Flags Contaminants (less than 5% of total reads) Needed for Phyloseq Analysis
## not needed (can do this with pipes)
applyFlag <- function(df){
  df <- df %>%
    left_join(., numberMatches)%>%
    mutate(contamFlag = ifelse(count <= 0.05*totalReads, T, F))
  return(df)
}

### Cut from getTaxa_sort_pure_mixed.R and phloyseq_processing_treeBuilding.R

## Perhaps unneeded (no filtration of contaminants)
# numberMatches <- otuTable %>% group_by(Isolate) %>%
#   summarize(otuMatches = n(), totalReads = sum(count))

# potentially unneeded
# hitEval <- otuTable %>%
#   group_by(Isolate)%>%
#   mutate(otuMatches = n(),
#          totalReads = sum(count),
#          contamFlag = ifelse(count <= 0.05*totalReads, T, F))%>%
#   ungroup()

## not sure why this is here, gonna proceed without it
# difNames <- hitEval %>% filter(str_detect(.$Isolate, "[^[:digit:]]"))

## This chunk determines filters out contaminent OTUs, leaving only OTUs present in each isolate
# otuKey <- otuTable %>%
#   group_by(Isolate)%>%
#   mutate(otuMatches = n(),
#          totalReads = sum(count),
#          contamFlag = ifelse(count <= 0.05*totalReads, T, F))%>%
#   ungroup() %>% 
#   filter(contamFlag == F, str_detect(.$Isolate, "^[:digit:]|^S|^_"))%>% 
#   distinct() %>%
#   left_join(., taxTable)

