library(tidyverse)
library(phyloseq)

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

list.file <- "./raw_data/mothur_outputs/final.opti_mcc.list"
constax.file <- "./raw_data/mothur_outputs/final.opti_mcc.0.03.cons.taxonomy"
shared.file <- "./raw_data/mothur_outputs/final.opti_mcc.shared"
# tree.file <- "./raw_data/mothur_outputs/final.opti_mcc.0.03.rep.otuRename.phylip.tre"

tree.file <- "./raw_data/mothur_outputs/final.opti_mcc.0.03.rep.otuRename.tree.nwk"


# Read mothur outputs into phyloseq/R
mo.data <- import_mothur(mothur_list_file = list.file,
                         mothur_constaxonomy_file = constax.file,
                         mothur_shared_file = shared.file)
                         # ,
                         # mothur_tree_file = tree.file)



# Convert stats data into phyloseq sample data (syntax of Isolate names changed to match mothur sample name syntax)
sam.data <- stats_meanCoefs %>%
  mutate(Isolate = str_replace(Isolate, ",", "sep"),
         Isolate = str_replace(Isolate, "(?<!sep\\d)D", "sepD"))%>%
  column_to_rownames(var = "Isolate")%>%
  sample_data(.)

# verify that sample names match (there will be more samples in mo.data, because we included environmental samples in our analysis) - most important here is that syntax is the same (underscores between number and day, ex: 23D3 becomes 23_D3, 30,1DF becomes 30_1DF)
sample_names(mo.data)
sample_names(sam.data)

# change taxonomic "rank names" to recognizable classifications (rank 1-7 to KPCOFGS)
rank_names(mo.data)
colnames(tax_table(mo.data)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
rank_names(mo.data)

# plot a test tree using the mothur data (not needed right here)
# plot_tree(mo.data, "treeonly", ladderize = T)

# full phyloseq object that combines the tax data from mothur with the coculture data
all.data <- merge_phyloseq(mo.data, sam.data)

#### Classify Isolates as Pure or Mixed Cultures ####
taxTable <- as.data.frame(all.data@tax_table)%>%
  rownames_to_column(., var = "otu")

otuTable <- as.data.frame(all.data@otu_table)%>%
  rownames_to_column(., var = "otu")%>%
  pivot_longer(!otu, names_to = "Isolate", values_to = "count")%>%
  filter(count != 0)%>%
  group_by(Isolate)%>%
  mutate(totalReads = sum(count),
         otuMatches = n())%>%
  ungroup()

isolateTax <- otuTable %>%
  mutate(contamFlag = ifelse(count <= 0.05*totalReads, T, F))%>%
  ## This filter pulls out only sample data (removes control, background, and pond OTUs)
  filter(contamFlag == F, str_detect(.$Isolate, "^[:digit:]|^S|^_"))%>% 
  # distinct()%>% ## not needed
  group_by(Isolate)%>%
  mutate(numOTUs = n(),
         mixed = ifelse(numOTUs > 1, T, F),
         Isolate = str_replace(Isolate, "(?<=\\d)sep(?=\\d)", ","),
         Isolate = str_replace(Isolate, "sep", ""))%>%
  # ,
  # ## This ifelse() removes "_" from the beginning of some sample names
  #          Isolate = ifelse(str_detect(Isolate, "^_") == T,
  #                           str_extract(Isolate, "(\\d+)(_|)D."),
  # ## some samples used a different naming convention ("S" means "DF")
  # ## This ifelse() changes these sample names to be uniform
  #                           ifelse(str_detect(Isolate, "^S") == T,
  #                                  paste(str_extract(Isolate, "(\\d+)"),
  #                                        "DF", sep = ""), Isolate)),
  # ## This ifelse removes mothur naming syntax and reverts Isolate names back to coculture naming conventions (removes underscores mainly)
  #          Isolate = ifelse(str_detect(Isolate, "\\d_D") == T,
#                           str_replace(Isolate, "_", ""),
#                           ifelse(str_detect(Isolate, "\\d_\\d") == T,
#                                  str_replace(Isolate, "_", ","), Isolate)))%>%

left_join(., taxTable)

#### Tree Files ####
treeGR <- plot_tree(all.data, ladderize=T, color = "grEffect", justify = "left", label.tips = "Genus")
treeGR

treeCC <- plot_tree(all.data, ladderize=T, color = "ccEffect", justify = "left")
treeCC

treeHost <- plot_tree(all.data, ladderize = T, color = "host_species", justify = "left")
treeHost


#### Additional Troubleshooting ####

sample_names <- data.frame(names = sample_names(all.data))%>%
  mutate(names = str_replace(names, "sep(?=\\d)", "."),
         names = str_replace(names, "sep(?=D)", ""),
         number = str_extract(names, ".+(?=D)"),
         day = str_extract(names, "D.+"))%>%
  distinct()

sample_names_MO <- data.frame(names = sample_names(mo.data))%>%
  mutate(names = str_replace(names, "sep(?=\\d)", "."),
         names = str_replace(names, "sep(?=D)", ""),
         number = str_extract(names, ".+(?=D)"),
         day = str_extract(names, "D.+"),
         mo.sample=T)%>%
  distinct()

sample_names_SAM <- data.frame(names = sample_names(sam.data))%>%
  mutate(names = str_replace(names, "sep(?=\\d)", "."),
         names = str_replace(names, "sep(?=D)", ""),
         number = str_extract(names, ".+(?=D)"),
         day = str_extract(names, "D.+"),
         co.sample=T)%>%
  distinct()

sample_names_all <- full_join(sample_names_MO, sample_names_SAM)%>%
  mutate(both = ifelse(mo.sample == T & co.sample == T, T, F))%>%
  filter(both == T)%>%
  distinct()

how_many_mixed <- select(mixed_cultures, Isolate)%>% distinct()

how_many_pure <- select(isolateTax, Isolate, mixed)%>% distinct()

numbers_onlyMO <- data.frame(number = as.numeric(unique(sample_names_MO$number))) # 274
numbers_onlySAM <- data.frame(number = as.numeric(unique(sample_names_SAM$number))) # 285

numbersAntiMO <- anti_join(numbers_onlyMO, numbers_onlySAM)
numbersAntiSAM <- anti_join(numbers_onlySAM, numbers_onlyMO)

#### Pruning ####
otu_df <- as.data.frame(all.data@otu_table) %>%
  rownames_to_column(., var = "otu") %>%
  pivot_longer(., cols = 2:250, names_to = "Isolate", values_to = "reads")

not_in_any_isolate <- otu_df %>%
  group_by(otu) %>%
  summarize(total_reads = sum(reads)) %>%
  filter(total_reads == 0)

otu_in_isolate <- select(otu_df, otu) %>%
  distinct()%>%
  anti_join(., not_in_any_isolate)%>%
  pull(otu)

pruned.data <- prune_taxa(otu_in_isolate, all.data)

pruned.tree <- plot_tree(pruned.data, label.tips="taxa_names", ladderize="left", color="Class")
pruned.tree

pruned.treeGR <- plot_tree(pruned.data, ladderize = T, color = "grEffect", justify = "left")
pruned.treeGR

pruned.treeCC <- plot_tree(pruned.data, ladderize = T, color = "ccEffect", justify = "left")
pruned.treeCC
#### Is anything quality?####
quality <- otuTable%>%
  mutate(percent = count/totalReads*100)%>%
  group_by(Isolate)%>%
  mutate(major = max(percent))
