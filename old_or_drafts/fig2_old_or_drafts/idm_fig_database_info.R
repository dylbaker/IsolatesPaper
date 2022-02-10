### Requires collection_tax_data.csv
### Requires logistic_mod_stats.csv
### Requires logistic_mod_coefs.csv

library(tidyverse)
library(outliers)

#### Functions (Grubbs Test) ####
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

#### READ CSV HERE ####

tax <- read.csv("./CoCulture_data/collection_tax_data.csv")%>%
  select(-X)

stats <- read.csv("./CoCulture_Data/logistic_mod_stats.csv")%>%
  select(-X)%>%
  mutate(Isolate = ifelse(str_detect(Isolate, "\\.") == T,
                          str_replace(Isolate, "\\.", ","), Isolate))

coefs <- read.csv("./CoCulture_Data/logistic_mod_coefs.csv")%>%
  select(-X)%>%
  mutate(Isolate = ifelse(str_detect(Isolate, "\\.") == T,
         str_replace(Isolate, "\\.", ","), Isolate))%>% 
  grubbs()

#### Script Continues ####

ac_coefs <- coefs %>% filter(Isolate == "AC")%>%
  group_by(host_species, plate_no)%>%
  mutate(nSamples = n(),
         meanCC = mean(asymptote),
         meanGR = mean(growthParam))%>%
  grubbs()

tax_stats <- left_join(stats, tax, by = "Isolate")

tax_coefs <- left_join(coefs, tax, by = "Isolate")

data_list <- coefs %>%
  filter(Isolate != "AC", Isolate != "MC")%>%
  split(., .$Isolate)
for(i in 1:length(data_list)){
  df <- data_list[[i]]
  plate <- unique(df$plate_no)
  host <- unique(df$host_species)
  ac <- ac_coefs %>%
    filter(host_species == host & plate_no == plate & outlier == F & asymptote >= 0)%>%
    summarize(acCC = mean(asymptote),
              acGR = mean(growthParam))
  acCC <- as.numeric(ac$acCC)
  acGR <- as.numeric(ac$acGR)
  data_list[[i]] <- df %>%
    mutate(normCC = asymptote/acCC,
           normGR = growthParam/acGR,
           logNormCC = log(asymptote/acCC),
           logNormGR = log(growthParam/acGR))
}
normCoefs <- bind_rows(data_list)

#### Data for the Database ####

data_4jacob <- normCoefs %>%
  left_join(., tax, by = "Isolate")%>%
  mutate(algalHost = host_species)%>%
  filter(mixed != T, outlier == F)%>%
  group_by(Isolate, algalHost, kingdom, phylum, class, order, family, genus, species)%>%
  mutate(n = n())%>%
  summarize(meanCC_RFU = mean(asymptote),
            seCC_RFU = sd(asymptote)/sqrt(n),
            mean_GR = mean(growthParam),
            se_GR = sd(growthParam)/sqrt(n),
            mean_normCC = mean(normCC),
            mean_normGR = mean(normGR),
            mean_logNormCC = mean(logNormCC),
            mean_logNormGR = mean(logNormGR))%>%
  distinct()

write.csv(data_4jacob, "pureCultures_databaseInfo.csv")

#### Plot Options ####

CCboxplot <- normCoefs %>%
  left_join(., tax)%>%
  ggplot(., aes(x = family, y = logNormCC))+
  geom_boxplot()+ 
  geom_hline(yintercept = 0, linetype = 2)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 315, hjust = 0))+
  labs(x = "Bacterial Family", y = "Log Fold Change in Algal Carrying Capacity", title = "Carrying Capacity Change by Bacterial Family")
CCboxplot
 
GRboxplot <- normCoefs %>%
  left_join(., tax)%>%
  ggplot(., aes(x = family, y = logNormGR))+
  geom_boxplot()+ 
  geom_hline(yintercept = 0, linetype = 2)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 315, hjust = 0))+
  labs(x = "Bacterial Family", y = "Log Fold Change in Algal Growth Rate", title = "Growth Rate Change by Bacterial Family")
GRboxplot

#### Plot Based on Vincent Paper ####
grPlot <- stats %>%
  select(Isolate, plate_no, pGR_greater, pGR_less, grEffect)%>%
  mutate(significant = ifelse(pGR_greater <= 0.05| pGR_less <= 0.05, T, F))%>%
  left_join(normCoefs, .)%>%
  left_join(., tax)%>%
  ggplot(., aes(y = family, x = logNormGR, color = phylum))+
  geom_vline(aes(xintercept = 0), color = "grey", linetype = "longdash")+
  geom_point(aes(shape = significant), size = 2)+
  xlim(-4, 4)
grPlot

ccPlot <- stats %>%
  select(Isolate, plate_no, pCC_greater, pCC_less, ccEffect)%>%
  mutate(significant = ifelse(pCC_greater <= 0.05 | pCC_less <= 0.05, T, F))%>%
  left_join(normCoefs, .)%>%
  left_join(., tax)%>%
ggplot(., aes(y = family, x = logNormCC, color = phylum))+
  geom_vline(aes(xintercept = 0), color = "grey", linetype = "longdash")+
  geom_point(aes(shape = significant), size = 2)+
  xlim(-4, 4)
ccPlot

mean_normCoefs <- normCoefs %>%
  group_by(Isolate, host_species, plate_no)%>%
  mutate()

grMeanPlot <- stats %>%
  select(Isolate, plate_no, pGR_greater, pGR_less, grEffect)%>%
  mutate(significant = ifelse(pGR_greater <= 0.05| pGR_less <= 0.05, T, F))%>%
  left_join(normCoefs, .)%>%
  group_by(Isolate)%>%
  mutate(n = n(),
         mean_normGR = mean(logNormGR, na.rm = T),
         se_normGR = sd(logNormGR, na.rm = T)/sqrt(n))%>%
  select(-exact_isolate, -asymptote, -intercept, -growthParam,
         -logNormCC, -logNormGR, -n)%>%
  distinct()%>%
  left_join(., tax)%>%
  ggplot(., aes(y = family, x = mean_normGR, color = host_species))+
  geom_vline(aes(xintercept = 0), color = "grey", linetype = "longdash")+
  geom_point(aes(shape = significant), size = 2)+
  geom_errorbar(aes(xmin = mean_normGR - se_normGR,
                    xmax = mean_normGR + se_normGR))
grMeanPlot

ccMeanPlot <- stats %>%
  select(Isolate, plate_no, pCC_greater, pCC_less, ccEffect)%>%
  mutate(significant = ifelse(pCC_greater <= 0.05| pCC_less <= 0.05, T, F))%>%
  left_join(normCoefs, .)%>%
  group_by(Isolate)%>%
  mutate(mean_normCC = mean(logNormCC, na.rm = T))%>%
  select(-exact_isolate, -asymptote, -intercept, -growthParam, -logNormCC, -logNormGR)%>%
  distinct()%>%
  left_join(., tax)%>%
  ggplot(., aes(y = family, x = mean_normCC, color = host_species))+
  geom_vline(aes(xintercept = 0), color = "grey", linetype = "longdash")+
  geom_point(aes(shape = significant), size = 2)
ccMeanPlot

## Significant Only ##
grPlot <- stats %>%
  select(Isolate, plate_no, pGR_greater, pGR_less, grEffect)%>%
  mutate(significant = ifelse(pGR_greater <= 0.05| pGR_less <= 0.05, T, F))%>%
  left_join(normCoefs, .)%>%
  drop_na()%>%
  # filter(significant == T)%>%
  left_join(., tax)%>%
  # filter(significant == T)%>%
  ggplot(., aes(y = family, x = logNormGR, color = host_species, shape = significant))+
  geom_vline(aes(xintercept = 0), color = "grey", linetype = "longdash")+
  geom_point(size = 2, alpha = 0.5, position = "jitter")+
  scale_shape_manual(values = c(1, 16))+
  xlim(-4, 4)+
  labs(x = "log(Coculture Growth Rate / Axenic Growth Rate)",
       y = "Bacterial Family",
       title = "Family Level Impact on Algal Growth Rate")
grPlot

ggsave("./CoCulture_Figures/cc_gr_plots/growthRate_plot_sig.png", grPlot, device = "png", height = 6, width = 8)

ccPlot <- stats %>%
  select(Isolate, plate_no, pCC_greater, pCC_less, ccEffect)%>%
  mutate(significant = ifelse(pCC_greater <= 0.05 | pCC_less <= 0.05, T, F))%>%
  left_join(normCoefs, .)%>%
  drop_na()%>%
  # filter(significant == T)%>%
  left_join(., tax)%>%
  # filter(significant == T)%>%
  ggplot(., aes(y = family, x = logNormCC, color = host_species, shape = significant))+
  geom_vline(aes(xintercept = 0), color = "grey", linetype = "longdash")+
  geom_point(size = 2, alpha = 0.5, position = "jitter")+
  scale_shape_manual(values = c(1, 16))+
  xlim(-4, 4)+
  labs(x = "log(Coculture Carrying Capacity / Axenic Carrying Capacity)",
       y = "Bacterial Family",
       title = "Family Level Impact on Algal Carrying Capacity")
ccPlot

ggsave("./CoCulture_Figures/cc_gr_plots/carryingCapacity_plot_sig.png", ccPlot, device = "png", height = 6, width = 8)

#### Faceted Plots ####
## Significant Only ##
grPlot_facet <- stats %>%
  select(Isolate, plate_no, pGR_greater, pGR_less, grEffect)%>%
  mutate(significant = ifelse(pGR_greater <= 0.05| pGR_less <= 0.05, T, F))%>%
  left_join(normCoefs, .)%>%
  drop_na()%>%
  # filter(significant == T)%>%
  left_join(., tax)%>%
  # filter(significant == T)%>%
  ggplot(., aes(y = family, x = logNormGR, color = host_species, shape = significant))+
  geom_vline(aes(xintercept = 0), color = "grey", linetype = "longdash")+
  geom_point(size = 2, alpha = 0.5, position = "jitter")+
  scale_shape_manual(values = c(1, 16))+
  xlim(-4, 4)+
  labs(x = "log(Coculture Growth Rate / Axenic Growth Rate)",
       y = "Bacterial Family",
       title = "Family Level Impact on Algal Growth Rate")+
  facet_grid(vars(host_species), scales = "free_y")
grPlot_facet

ggsave("./CoCulture_Figures/cc_gr_plots/growthRate_plot_sig_facet.png", grPlot_facet, device = "png", height = 10, width = 8)

ccPlot_facet <- stats %>%
  select(Isolate, plate_no, pCC_greater, pCC_less, ccEffect)%>%
  mutate(significant = ifelse(pCC_greater <= 0.05 | pCC_less <= 0.05, T, F))%>%
  left_join(normCoefs, .)%>%
  drop_na()%>%
  # filter(significant == T)%>%
  left_join(., tax)%>%
  # filter(significant == T)%>%
  ggplot(., aes(y = family, x = logNormCC, color = host_species, shape = significant))+
  geom_vline(aes(xintercept = 0), color = "grey", linetype = "longdash")+
  geom_point(size = 2, alpha = 0.5, position = "jitter")+
  scale_shape_manual(values = c(1, 16))+
  xlim(-4, 4)+
  labs(x = "log(Coculture Carrying Capacity / Axenic Carrying Capacity)",
       y = "Bacterial Family",
       title = "Family Level Impact on Algal Carrying Capacity")+
  facet_grid(vars(host_species), scales = "free_y")
ccPlot_facet

ggsave("./CoCulture_Figures/cc_gr_plots/carryingCapacity_plot_sig_facet.png", ccPlot_facet, device = "png", height = 10, width = 8)

#### Chlorella Test Plot ####
chlorella <- tax_stats %>%
  # mutate(ccPN = ifelse(ccEffect == "Increased CC", 1,
  #                      ifelse(ccEffect == "Decreased CC", -1, 0)),
  #        grPN = ifelse(grEffect == "Increased GR", 1,
  #                      ifelse(grEffect == "Decreased GR", -1, 0)),
  #        pn = ccPN + grPN,
  #        color = ifelse(pn == 2, "Positive",
  #                       ifelse(pn == 1, "Slightly Positive",
  #                              ifelse(pn == -1, "Slightly Negative",
  #                                     ifelse(pn == -2, "Negative",
  #                                            ifelse(pn == 0 & ccPN == 0 & grPN == 0, "No Effect",
  #                                                   ifelse(pn == 0 & (ccPN != 0 | grPN != 0), "Conflicting", "ERROR")))))))%>%
  # select(-ccPN, -grPN, -pn) %>%
  filter(host_species == "chlorella") %>%
  ggplot(data = ., aes(x = family, fill = grEffect))+
  geom_bar(stat = "count", position = "dodge")+
  facet_grid(ccEffect ~ .)+
  theme(axis.text.x = element_text(angle = 315, hjust = 0),
        strip.text.y = element_text(angle = 0))
chlorella

speciesCC_plot <- function(chr){
  plot <- tax_stats %>%
    filter(host_species == chr) %>%
    ggplot(data = ., aes(x = family, fill = grEffect))+
    geom_bar(stat = "count", position = "dodge")+
    facet_grid(ccEffect ~ .)+
    labs(title = chr)+
    theme(axis.text.x = element_text(angle = 315, hjust = 0),
          strip.text.y = element_text(angle = 0))
  plotName <- paste("./CoCulture_Figures/cc_gr_plots/", chr, "_ccPlot.png", sep = "")
  ggsave(plotName, plot, device = "png", height = 6, width = 8)
  return(plot)
} 

speciesGR_plot <- function(chr){
  plot <- tax_stats %>%
    filter(host_species == chr) %>%
    ggplot(data = ., aes(x = family, fill = ccEffect))+
    geom_bar(stat = "count", position = "dodge")+
    facet_grid(grEffect ~ .)+
    labs(title = chr)+
    theme(axis.text.x = element_text(angle = 315, hjust = 0),
          strip.text.y = element_text(angle = 0))
  plotName <- paste("./CoCulture_Figures/cc_gr_plots/", chr, "_grPlot.png", sep = "")
  ggsave(plotName, plot, device = "png", height = 6, width = 8)
  return(plot)
} 

names <- as.list(unique(tax_stats$host_species))
lapply(names, speciesCC_plot)
lapply(names, speciesGR_plot)


