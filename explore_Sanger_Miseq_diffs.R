### source("./master.R")
## checking out the differences between Sanger and Miseq reads of the same Isolate(s)


## Quick Check

SM_duplicates <- isolateTax %>%
  select(Isolate, readType)%>%
  distinct()%>%
  group_by(Isolate)%>%
  mutate(duplicate = ifelse(n() > 1, T, F))%>%
  filter(duplicate == T)%>%
  left_join(., isolateTax)%>%
  ungroup()%>%
  group_by(Isolate, readType)%>%
  mutate(maxCount = max(count))%>%
  ungroup()%>%
  filter(count == maxCount)%>%
  ungroup()%>%
  select(Isolate, readType, otu, Family, Genus)%>%
  pivot_wider(names_from = readType, values_from = c(otu, Family, Genus))%>%
  mutate(match_fam = ifelse(Family_miseq == Family_sanger, T, F),
         match_gen = ifelse(Genus_miseq == Genus_sanger, T, F))
