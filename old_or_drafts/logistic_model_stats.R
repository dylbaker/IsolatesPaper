### Writes logistic_mod_coefs.csv
### Writes logistic_mod_stats.csv
### Writes isolateKey.csv
### requires coculture_data.csv (full raw data spreadsheet)

library(tidyverse)
library(broom)
# for Foreach/doParallel
library(foreach)
library(doParallel)
numCores <- detectCores()
registerDoParallel(numCores-1)

flag_algae <- function(df){
  data_raw <- df
  AC_day0_avg <- df %>%
    filter(Isolate == "AC" & read_day == 0)%>%
    summarize(ChlA_100_mean = mean(ChlA_100))%>%
    as.numeric()
  data_flagged <- df %>%
    filter(read_day == 0)%>%
    mutate(algae_flag = ifelse(ChlA_100 > 2*AC_day0_avg, T, F))%>%
    select(Isolate, algae_flag)%>%
    distinct()%>%
    left_join(data_raw, .)
  return(data_flagged)
}

data <- read.csv("./raw_data/coculture_data.csv") %>% 
  select(-X)

algal_species <- unique(data$host_species)

data_flagged <- as.list(1:length(algal_species))
system.time(for(i in 1:length(algal_species)){
  host <- algal_species[i]
  host_data <- data %>% filter(host_species == host)
  data_flagged[[i]] <- flag_algae(host_data)
}) ## System Time elapsed = 0.06

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


# # Temporarily Commented to make the script run faster
# plotList <- foreach(i = 1:length(data_split), .packages = "tidyverse") %dopar%{
#   tryCatch({
#     df <- data_split[[i]]
#     coef_df <- coefList[[i]]
# 
#     aMod <- coef_df[["asymptote"]]
#     iMod <- coef_df[["intercept"]]
#     xMod <- coef_df[["growthParam"]]
# 
#     x <- c(min(df$read_timeHours):max(df$read_timeHours))
#     y <- aMod/(1 + exp(-(iMod + xMod * x)))
#     predict <- data.frame(x, y)
#     plotTitle <- paste(unique(df$exact_isolate), "Logistic Model", sep = " ")
# 
#     plot <- ggplot(data = df, aes(x = read_timeHours, y = ChlA_100))+
#       geom_point(color = "blue")+
#       theme_bw()+
#       labs(x = "Read Time (Hours)",
#            y = "Chlorophyll A Fluorescence (RFU)",
#            title = plotTitle)+
#       geom_line(data=predict, aes(x = x, y = y), size = 1)
# 
#     return(plot)
#     }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
#   }

iKey <- data %>%
  select(exact_isolate, Isolate, host_species, plate_no)%>%
  distinct()

write.csv(iKey, file = "./raw_data/isolateKey.csv")

coefs <- bind_rows(coefList)%>%
  left_join(., iKey)

write.csv(coefs, file = "./raw_data/logistic_mod_coefs.csv")

ACcoefs <- coefs %>% filter(Isolate == "AC")
coefList_4t <- split(coefs, coefs$Isolate)

stats <- foreach(i = 1:length(coefList_4t), .packages = c("tidyverse", "broom")) %dopar%{
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

stats_all <- bind_rows(stats)%>%
  select(Isolate, host_species, plate_no, pCC_greater, pCC_less, pGR_greater, pGR_less)%>%
  mutate(ccEffect = ifelse(pCC_greater <= 0.05 & pCC_less > 0.05, "Increased CC",
                           ifelse(pCC_greater > 0.05 & pCC_less <= 0.05, "Decreased CC",
                                  ifelse(pCC_greater > 0.05 & pCC_less > 0.05, "No Significant CC Change", "Error"))),
         grEffect = ifelse(pGR_greater <= 0.05 & pGR_less > 0.05, "Increased GR",
                           ifelse(pGR_greater > 0.05 & pGR_less <= 0.05, "Decreased GR",
                                  ifelse(pGR_greater > 0.05 & pGR_less > 0.05, "No Significant GR Change", "Error"))),
         annotation = paste(ccEffect, grEffect, sep = "\n"))

write.csv(stats_all, file = "./raw_data/logistic_mod_stats.csv")

#### Triplicate Plots ####
acData <- bind_rows(data_flagged) %>%
  mutate(Day_isolated = ifelse(Isolate == "AC", "AC", Day_isolated))%>%
  filter(algae_flag != T,
         Day_isolated == "AC")%>%
  left_join(coefs)

data_split_trip <- bind_rows(data_flagged) %>%
  mutate(Day_isolated = ifelse(Isolate == "AC", "AC", Day_isolated))%>%
  filter(algae_flag != T,
         Day_isolated != "Ctrl")%>%
  split(., .$Isolate)

plotList_trip <- foreach(i = 1:length(data_split_trip),
                         .packages = c("tidyverse", "foreach")) %dopar%{
                           df <- data_split_trip[[i]]
                           isolate <- unique(as.character(df$Isolate))
                           coef_df <- coefs[coefs$Isolate==isolate,]
                           
                           aMod <- coef_df[["asymptote"]]
                           iMod <- coef_df[["intercept"]]
                           xMod <- coef_df[["growthParam"]]
                           
                           x <- c(min(df$read_timeHours):max(df$read_timeHours))
                           predict <- foreach(j = 1:nrow(coef_df),
                                              .combine = "rbind") %do% {
                             predict <- data.frame(exact_isolate = coef_df$exact_isolate[j],
                                                   x = c(min(df$read_timeHours):max(df$read_timeHours)),
                                                   y = aMod[j]/(1 + exp(-(iMod[j] + xMod[j] * x))))
                             return(predict)
                           }
                           
                           meanCoef_df <- coef_df %>%
                             summarize(aMod = mean(asymptote),
                                       iMod = mean(intercept),
                                       xMod = mean(growthParam))
                           
                           meanPredict <- data.frame(x = c(min(df$read_timeHours):max(df$read_timeHours)),
                                                     y = meanCoef_df$aMod/(1 + exp(-(meanCoef_df$iMod + meanCoef_df$xMod * x))))
                           
                           plotTitle <- paste(unique(df$Isolate), "Logistic Model", sep = " ")
                           
                           fileName <- paste("./CoCulture_Figures/isolate_plots/mod_plots/",
                                             unique(df$Isolate), "_logModel_plot.png",
                                             sep = "")
                           
                           plot <- ggplot(data = df, aes(x = read_timeHours, y = ChlA_100))+
                             geom_point(color = "forestgreen")+
                             theme_bw()+
                             labs(x = "Read Time (Hours)",
                                  y = "Chlorophyll A Fluorescence (RFU)",
                                  title = plotTitle)+
                             geom_line(data=predict, 
                                       aes(x = x, y = y, group = exact_isolate), 
                                       size = 0.5, color = "green", linetype = 2)+
                             geom_line(data = meanPredict, aes(x = x, y = y),
                                       size = 1, color = "forestgreen")
                           
                           ggsave(filename = fileName, plot = plot, device = "png", height = 4, width = 3)
                           
                           return(plot)
                         }

#### Something Here Breaks Everything, work this out later ####
# stats_bar <- drop_na(stats_all)%>%
#   ggplot(., aes(x = grEffect, fill = grEffect))+
#   geom_bar(stat = "count")+
#   facet_grid(vars(ccEffect))
# stats_bar
# 
# stats_barHosts <- data %>%
#   select(Isolate, host_species)%>%
#   left_join(stats_all, .)%>%
#   drop_na()%>%
#   ggplot(., aes(x = grEffect, fill = grEffect))+
#   geom_bar(stat = "count")+
#   facet_grid(host_species ~ ccEffect)
# stats_barHosts


#### Interesting, Revisit this ####
# line_fig_data <- left_join(data, stats_all)%>%
#   select(Isolate, host_species, read_timeHours, ChlA_100, ccEffect, grEffect, annotation)%>%
#   split(., .$host_species)
# 
# plots <- as.list(1:length(line_fig_data))
# for(i in 1:length(line_fig_data)){
#   df <- line_fig_data[[i]]%>%
#     drop_na()
#   ac_df <- filter(line_fig_data[[i]], Isolate == "AC")
# 
#   plots[[i]] <- ggplot(df, aes(x = read_timeHours, y = ChlA_100, color = annotation))+
#     geom_point()+
#     geom_smooth()+
#     geom_point(data = ac_df, color = "black")+
#     geom_smooth(data = ac_df, color = "black")
# }
