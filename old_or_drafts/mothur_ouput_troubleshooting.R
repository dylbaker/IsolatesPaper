# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.16")

BiocManager::install("phyloseq")


library(tidyverse)
library(phyloseq)
library(ape)

####Newest outputs from mothur 10/19/22
fasttree <- "./mothur_outputs/GRBC_seqs/final.asvRename.rep.tree.nwk"
clearcut <- "./mothur_outputs/GRBC_seqs/final.an.rep.asvRename.phylip.tre"
constax.file <- "./mothur_outputs/GRBC_seqs/final.pick.asv.ASV.cons.taxonomy"
shared.file <- "./mothur_outputs/GRBC_seqs/final.pick.asv.shared"
list.file <- "./mothur_outputs/GRBC_seqs/final.pick.asv.list"

# Read mothur outputs into phyloseq/R
mo.data <- import_mothur(mothur_constaxonomy_file = constax.file,
                         mothur_shared_file = shared.file,
                         mothur_list_file = list.file,
                         mothur_tree_file = fasttree)
# change taxonomic "rank names" to recognizable classifications (rank 1-7 to KPCOFGS)
rank_names(mo.data)
colnames(tax_table(mo.data)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
rank_names(mo.data)
#Change sample names to desired ones.
otuTable <- as.data.frame(mo.data@otu_table)|>
  rownames_to_column(var = "otu")|>
  pivot_longer(!otu, names_to = "Isolate", values_to = "count")|>
  filter(count != 0)|>
  group_by(Isolate)|>
  mutate(totalReads = sum(count),
         otuMatches = n())|>
  ungroup()
sample_names <- as_tibble(unique(otuTable$Isolate)) |>
  mutate(Isolate = case_when(str_detect(value, "^S") ~ paste(str_extract(value, "([A-Z0-9]+)"),"D31", sep = "_"),
                             TRUE ~ value),
         Isolate = str_replace(Isolate, "^S",""),
         Isolate = str_replace(Isolate, "F","31"),
         Isolate = str_replace(Isolate, "point", "_"))

# Convert stats data into phyloseq sample data (syntax of Isolate names changed to match mothur sample name syntax)
sam.data <- sample_names |>
  column_to_rownames(var = "value")|>
  sample_data()
# full phyloseq object that combines the tax data from mothur with the coculture data
all.data <- merge_phyloseq(mo.data, sam.data)
rel_abund <- transform_sample_counts(all.data, function(x) x /sum(x)) |>
  filter_taxa(function(x) sum(x) > .05, TRUE)
sample_names(all.data) <- all.data@sam_data[["value"]]
sample_names(all.data) <- sample_data(all.data)$Isolate

# Agglomerate similar taxa
tip_glom <- tip_glom(all.data, h = 0.1) #Cluster based on cophenetic distances of <0.2
genus_glom <- tax_glom(all.data, taxrank = rank_names(mo.data)[6])



# verify that sample names match (there will be more samples in mo.data, because we included environmental samples in our analysis) - most important here is that syntax is the same (underscores between number and day, ex: 23D3 becomes 23_D3, 30,1DF becomes 30_1DF)
sample_names(mo.data) <- sample_names[2]
sample_names(sam.data)


# plot a test tree using the mothur data (not needed right here)
#plot_tree(mo.data, "treeonly", ladderize = T)
phylum.tree <- plot_tree(rel_abund, 
                         label.tips="Genus", 
                         ladderize = "left",
                         color = "Phylum")
order.tree <- plot_tree(rel_abund, 
                         label.tips="Genus", 
                         ladderize = "left",
                         color = "Order")
plot_tree(tip_glom(mo.data, h=0.2),label.tips="taxa_names", size="abundance", title="After tip_glom()")


#### Classify Isolates as Pure or Mixed Cultures ####
taxTable <- as.data.frame(all.data@tax_table)|>
  rownames_to_column(., var = "otu")
genustaxTable <- as.data.frame(genus_glom@tax_table) |>
  rownames_to_column(., var = "otu")
genusotuTable <- as.data.frame(genus_glom@otu_table)|>
  rownames_to_column(., var = "otu")|>
  pivot_longer(!otu, names_to = "Isolate", values_to = "count")|>
  filter(count != 0)|>
  group_by(Isolate)|>
  mutate(totalReads = sum(count),
         otuMatches = n())|>
  ungroup()
isolateTax <- genusotuTable |>
  mutate(contamFlag = ifelse(count <= 0.1*totalReads, T, F))|>
  filter(contamFlag == F) |>
  group_by(Isolate)|>
  mutate(numOTUs = n(),
         mixed = ifelse(numOTUs > 1, T, F),
         Isolate = ifelse(str_detect(Isolate, "^S"),
                                   paste(str_extract(Isolate, "([A-Z0-9]+)"),"D31", sep = "_"), Isolate),
         Isolate = str_replace(Isolate, "^S",""),
         Isolate = str_replace(Isolate, "F","31"))|>
  # ## This ifelse removes mothur naming syntax and reverts Isolate names back to coculture naming conventions (removes underscores mainly)
  #          Isolate = ifelse(str_detect(Isolate, "\\d_D") == T,
#                           str_replace(Isolate, "_", ""),
#                           ifelse(str_detect(Isolate, "\\d_\\d") == T,
#                                  str_replace(Isolate, "_", ","), Isolate)))|>

left_join(., taxTable)

#### Tree Files ####
treeGR <- plot_tree(all.data, ladderize=T, color = "grEffect", justify = "left", label.tips = "Genus")
treeGR

treeCC <- plot_tree(all.data, ladderize=T, color = "ccEffect", justify = "left")
treeCC

treeHost <- plot_tree(all.data, ladderize = T, color = "host_species", justify = "left")
treeHost


#### Additional Troubleshooting ####

sample_names <- data.frame(names = sample_names(all.data))|>
  mutate(names = str_replace(names, "sep(?=\\d)", "."),
         names = str_replace(names, "sep(?=D)", ""),
         number = str_extract(names, ".+(?=D)"),
         day = str_extract(names, "D.+"))|>
  distinct()

sample_names_MO <- data.frame(names = sample_names(mo.data))|>
  mutate(names = str_replace(names, "sep(?=\\d)", "."),
         names = str_replace(names, "sep(?=D)", ""),
         number = str_extract(names, ".+(?=D)"),
         day = str_extract(names, "D.+"),
         mo.sample=T)|>
  distinct()

sample_names_SAM <- data.frame(names = sample_names(sam.data))|>
  mutate(names = str_replace(names, "sep(?=\\d)", "."),
         names = str_replace(names, "sep(?=D)", ""),
         number = str_extract(names, ".+(?=D)"),
         day = str_extract(names, "D.+"),
         co.sample=T)|>
  distinct()

sample_names_all <- full_join(sample_names_MO, sample_names_SAM)|>
  mutate(both = ifelse(mo.sample == T & co.sample == T, T, F))|>
  filter(both == T)|>
  distinct()

how_many_mixed <- select(mixed_cultures, Isolate)|> distinct()

how_many_pure <- select(isolateTax, Isolate, mixed)|> distinct()

numbers_onlyMO <- data.frame(number = as.numeric(unique(sample_names_MO$number))) # 274
numbers_onlySAM <- data.frame(number = as.numeric(unique(sample_names_SAM$number))) # 285

numbersAntiMO <- anti_join(numbers_onlyMO, numbers_onlySAM)
numbersAntiSAM <- anti_join(numbers_onlySAM, numbers_onlyMO)

#### Pruning ####
otu_df <- as.data.frame(all.data@otu_table) |>
  rownames_to_column(., var = "otu") |>
  pivot_longer(., cols = 2:250, names_to = "Isolate", values_to = "reads")

not_in_any_isolate <- otu_df |>
  group_by(otu) |>
  summarize(total_reads = sum(reads)) |>
  filter(total_reads == 0)

otu_in_isolate <- select(otu_df, otu) |>
  distinct()|>
  anti_join(., not_in_any_isolate)|>
  pull(otu)

pruned.data <- prune_taxa(otu_in_isolate, all.data)

pruned.tree <- plot_tree(pruned.data, label.tips="taxa_names", ladderize="left", color="Phylum")
pruned.tree

pruned.treeGR <- plot_tree(pruned.data, ladderize = T, color = "grEffect", justify = "left")
pruned.treeGR

pruned.treeCC <- plot_tree(pruned.data, ladderize = T, color = "ccEffect", justify = "left")
pruned.treeCC
#### Is anything quality?####
quality <- otuTable|>
  mutate(percent = count/totalReads*100)|>
  group_by(Isolate)|>
  mutate(major = max(percent))

