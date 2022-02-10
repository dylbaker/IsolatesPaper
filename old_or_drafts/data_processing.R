# Writes coculture_data.csv

library(tidyverse)
library(lubridate)

setwd("./raw_data")

#### Functions ####
read_plate <- function(plate_list) {
  read_csv(plate_list,
           col_types = cols(.default = col_double(),
                            Well = col_character(),
                            `Well ID` = col_character(),
                            # Format no longer an issue
                            `Reading Date/Time` = col_datetime(format = "%Y-%m-%d %H:%M:%S"))) %>%
    rename(ChlA_100 = `Mean RFU [ChlA_100:460,685]`, 
           read_time = `Reading Date/Time`)%>%
    #Extract plate information and species from each file to aid in isolate mapping.
    mutate(plate_no = str_extract(str_extract(plate_list, pattern = "Plate_(\\d\\d|\\d)"),pattern = regex("(\\d\\d|\\d)")), 
           host_species = str_extract(plate_list, pattern = "(chlorella|coelastrum|scenedesmus|monoraphidium|selenastrum)"),
           filename = plate_list,
           read_day = str_extract(str_extract(plate_list, pattern = "(\\d\\d|\\d).csv"), pattern = regex("(\\d\\d|\\d)"))) %>%
    select(Well, ChlA_100, read_time, read_day, plate_no, host_species)
}

read_map <- function(map_list) {
  read_csv(map_list, col_types = cols(.default = col_character())) %>%
    mutate(Isolate = as.character(Isolate),
           plate_no = str_extract(str_extract(map_list, pattern = "Plate_(\\d\\d|\\d)"),pattern = regex("(\\d\\d|\\d)")), 
           host_species = str_extract(map_list, pattern = "(chlorella|coelastrum|scenedesmus|monoraphidium|selenastrum)"),
           filename = map_list) %>%
    select(-filename)
}

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

flag_outlier <- function(df){
  data <- df %>%
    group_by(plate_no, Isolate, read_day)%>%
    mutate(day_mean = mean(ChlA_100, na.rm=TRUE),
           day_sd = sd(ChlA_100, na.rm=TRUE))%>%
    ungroup()%>%
    mutate(Outlier = ifelse((day_mean-(2*day_sd))>ChlA_100|(day_mean+(2*day_sd))<ChlA_100, TRUE, FALSE))
  return(data)
}

day_Isolated <- function(df){
  data <- df%>%
    mutate(Day_isolated = ifelse(str_detect(Isolate, "D3")==TRUE, "D3", ifelse(str_detect(Isolate, "DF")==TRUE, "D31", "Ctrl")))
  return(data)
}

# flag_algae <- function(df){
#   data_raw <- df
#   AC_day0_avg <- df %>%
#     filter(Isolate == "AC" & read_day == 0)%>%
#     summarize(ChlA_100_mean = mean(ChlA_100))%>%
#     as.numeric()
#   data_flagged <- df %>%
#     filter(read_day == 0)%>%
#     mutate(algae_flag = ifelse(ChlA_100 > 2*AC_day0_avg, TRUE, FALSE))%>%
#     select(Isolate, algae_flag)%>%
#     distinct()%>%
#     left_join(data_raw, .)
#   return(data_flagged)
# }

#### Plate Data ####
folder_list <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE) %>%
  ifelse(str_detect(. , pattern = regex("plate_(\\d|\\d\\d)_(chlorella|coelastrum|scenedesmus|monoraphidium|selenastrum)$", ignore_case = TRUE)), .,NA) %>%
  na.omit()

plate_list <- list.files(path = folder_list, 
                         pattern = "^Plate_(\\d\\d|\\d)_Day_(\\d\\d|\\d).csv",
                         full.names = TRUE, 
                         ignore.case = TRUE)

plate_data_all <- map_df(plate_list, ~read_plate(.))

#### Map Data ####
map_list <- list.files(path = folder_list,
                       pattern = "^plate_(\\d\\d|\\d)_map.csv",
                       full.names = TRUE,
                       ignore.case = TRUE)

map_data_all <- map_df(map_list, ~read_map(.))

#### Combine Plate and Map ####
data <- inner_join(plate_data_all, map_data_all,
                   by = c("Well", "plate_no", "host_species"))%>%
  calc_time()%>%
  day_Isolated()%>%
  flag_outlier()%>%
  # mutate(Isolate_Plate = paste(Isolate, plate_no, sep = "_"))%>%
  select(Well, host_species, plate_no, Day_isolated, Isolate, exact_isolate, read_day, read_time, read_timeHours, ChlA_100, day_mean, day_sd, Outlier)

#### Write CSV with all data ####
write.csv(data, "./coculture_data.csv")