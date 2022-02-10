#### Writes collection_pure_data.csv (needed for fig 1)
#### Writes collection_tax_data.csv (needed for fig 2)
#### rename both of these to make things a little clearer

library(tidyverse)
library(phyloseq)

setwd("C:/Users/jlaue/OneDrive/Documents/Denef_Lab_Stuff/R_scripts_git")

load("./phyloseq_data/Phyloseq.RData")

# taxTable_kosher <- as.data.frame(tax_table(phy.otu))

taxTable <- as.data.frame(phy.otu@tax_table)%>%
  rownames_to_column(., var = "otu")%>%
  rename(kingdom = Rank1,
         phylum = Rank2,
         class = Rank3,
         order = Rank4,
         family = Rank5,
         genus = Rank6,
         species = Rank7)

otuTable <- as.data.frame(phy.otu@otu_table) %>%
  rownames_to_column(., var = "otu")%>%
  pivot_longer(!otu, names_to = "Isolate", values_to = "count") %>%
  filter(count != 0)

numberMatches <- otuTable %>% group_by(Isolate) %>%
  summarize(otuMatches = n(), totalReads = sum(count))

## Change Criteria Here pending more info on proper way to filter
applyFlag <- function(df){
  df <- df %>%
    left_join(., numberMatches)%>%
    mutate(contamFlag = ifelse(count <= 0.05*totalReads, T, F))
  return(df)
}

hitEval <- otuTable %>%
  split(., .$Isolate)%>%
  lapply(., applyFlag)%>%
  bind_rows(.)

difNames <- hitEval %>% filter(str_detect(.$Isolate, "[^[:digit:]]"))

otuKey <- hitEval %>% 
  filter(contamFlag == F, str_detect(.$Isolate, "^[:digit:]|^S|^_"))%>% 
  distinct() %>%
  left_join(., taxTable)

howMany <- length(unique(otuKey$otu))

isolatesPerOTU <- otuKey %>% select(Isolate, otu, genus) %>% distinct() %>%
  ggplot(., aes(y = fct_infreq(genus)))+
  geom_bar(stat = "count")
# + theme(axis.text.x = element_text(angle = 315, hjust = 0))
isolatesPerOTU

OTUsPerIsolate <- otuKey %>%
  select(Isolate, otu)%>%
  distinct()%>%
  group_by(Isolate)%>%
  summarize(numOTUs  = n())%>%
  mutate(mixed = ifelse(numOTUs > 1, T, F))%>%
  left_join(., otuKey)

otus_in_multiCultures <- OTUsPerIsolate %>%
  filter(mixed == T)%>%
  ggplot(., aes(y = fct_infreq(genus)))+
  geom_bar(stat = "count")+
  labs(title = "Bacterial Families found in Cultures containing Multiple OTUs",
       x = "Number of Cultures OTU was Found in",
       y = "genus")
otus_in_multiCultures

howManyPure <- OTUsPerIsolate%>%
  filter(mixed == F)%>%
  select(otu, Isolate)%>%
  distinct()%>%
  group_by(otu)%>%
  summarize(numIsolates = n())

#### Test Plots ####
otus_in_pureCultures <- OTUsPerIsolate %>%
  filter(mixed == F)%>%
  ggplot(., aes(y = fct_infreq(family)))+
  geom_bar(stat = "count")+
  labs(title = "Bacterial Families found in Cultures containing Single OTUs",
       x = "Number of Cultures OTU was Found in",
       y = "family")
otus_in_pureCultures

contaminants <- hitEval %>%
  filter(contamFlag == T, str_detect(.$Isolate, "^[:digit:]|^S|^_"))%>%
  mutate(percentage = count/totalReads *100)%>%
  distinct()

commonContams <- left_join(contaminants, taxTable)%>%
  select(otu, Isolate, family)%>%
  group_by(otu, family)%>%
  summarize(numIsolates = n())%>%
  filter(numIsolates >= 6)%>%
  ggplot(., aes(y = reorder(family, numIsolates), x = numIsolates))+
  geom_bar(stat = "Identity")+
  labs(title = "Families Commonly found as Contaminants",
       x = "Number of Cultures OTU was Found in",
       y = "family")
commonContams

#### End Test Plots ####
isolateTax <- OTUsPerIsolate %>%
  mutate(Isolate = ifelse(str_detect(Isolate, "^_") == T, 
                          str_extract(Isolate, "(\\d+)(_|)D."),
                          ifelse(str_detect(Isolate, "^S") == T, paste(str_extract(Isolate, "(\\d+)"), "DF", sep = ""), Isolate)),
         Isolate = ifelse(str_detect(Isolate, "\\d_D") == T,
                          str_replace(Isolate, "_", ""),
                          ifelse(str_detect(Isolate, "\\d_\\d") == T,
                                 str_replace(Isolate, "_", ","), Isolate)))

write.csv(isolateTax, file = "./CoCulture_data/collection_tax_data_raw.csv")

tax4stats <- isolateTax %>%
  mutate(kingdom = ifelse(mixed == T, "mixed", kingdom),
         phylum = ifelse(mixed == T, "mixed", phylum),
         class = ifelse(mixed == T, "mixed", class),
         order = ifelse(mixed == T, "mixed", order),
         family = ifelse(mixed == T, "mixed", family),
         genus = ifelse(mixed == T, "mixed", genus),
         species = ifelse(mixed == T, "mixed", species))%>%
  select(Isolate, numOTUs, mixed, contamFlag, kingdom, phylum, class, order, family, genus, species)%>%
  distinct()

write.csv(tax4stats, file = "./CoCulture_data/collection_tax_data.csv")

mixed_cultures <- isolateTax %>% filter(mixed == T)
write.csv(mixed_cultures, file = "./CoCulture_data/collection_mixed_cultures.csv")

pure_cultures <- isolateTax %>% filter(mixed == F)
write.csv(pure_cultures, file = "./CoCulture_data/collection_pure_cultures.csv")

write.csv(isolateTax, file = "./CoCulture_data/collection_all_cultures.csv")

