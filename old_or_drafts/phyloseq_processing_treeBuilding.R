library(phyloseq)
library(ape)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(ggtree)

setwd(".")

#### Tree Setup ####
tree.file <- "fullTree.phylip.tre"
constax.file <- "fullTree.cons.taxonomy"
shared.file <- "fullTree.an.shared"
list.file <- "fullTree.an.list"

stats_annotation <- read.csv("logistic_mod_stats.csv")%>%
  select(-X, -host_species, -plate_no)

sample_key <- read.csv("pureCultures_databaseInfo.csv")%>%
  select(-kingdom, -phylum, -class, -order, -family, -genus, -species)%>%
  distinct(.)%>%
  left_join(., stats_annotation)%>%
  mutate(Isolate = str_replace(Isolate, ",", "_"),
         Isolate = str_replace(Isolate, "(?<!\\_\\d)D", "\\_D"))%>%
  column_to_rownames(var = "Isolate")

sam.data <- sample_data(sample_key)

mo.data <- import_mothur(mothur_list_file = list.file,
                         mothur_constaxonomy_file = constax.file,
                         mothur_shared_file = shared.file,
                         mothur_tree_file = tree.file)
sample_names(mo.data)
sample_names(sam.data)

rank_names(mo.data)
colnames(tax_table(mo.data)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rank_names(mo.data)

plot_tree(mo.data, "treeonly", ladderize = T)

all.data <- merge_phyloseq(mo.data, sam.data)

treeGR <- plot_tree(all.data, ladderize=T, color = "grEffect", justify = "left", label.tips = "Genus")
treeGR

treeCC <- plot_tree(all.data, ladderize=T, color = "ccEffect", justify = "left")
treeCC

treeHost <- plot_tree(all.data, ladderize = T, color = "algalHost", justify = "left")
treeHost

#### Testing Fasttree ####

fasttree.file <- "fullTree.fasttree.nwk"

ft.data <- import_mothur(mothur_list_file = list.file,
                         mothur_constaxonomy_file = constax.file,
                         mothur_shared_file = shared.file,
                         mothur_tree_file = fasttree.file)%>%
  merge_phyloseq(., sam.data)

rank_names(ft.data)
colnames(tax_table(ft.data)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rank_names(ft.data)

fasttree <- plot_tree(ft.data, "treeonly", ladderize = T) + coord_polar(theta = "y")
fasttree

fasttreeGR <- plot_tree(ft.data, ladderize=T, color = "grEffect")
fasttreeGR

fasttreeCC <- plot_tree(ft.data, ladderize=T, color = "ccEffect", justify = "left")
fasttreeCC

fasttreeHost <- plot_tree(ft.data, ladderize = T, color = "algalHost", justify = "left")
fasttreeHost


#### Filter OTUs ####

otu_df <- as.data.frame(all.data@otu_table) %>%
  rownames_to_column(., var = "otu") %>%
  pivot_longer(., cols = 2:177, names_to = "Isolate", values_to = "reads")

not_in_any_isolate <- otu_df %>%
  group_by(otu) %>%
  summarize(total_reads = sum(reads)) %>%
  filter(total_reads == 0)

otu_in_isolate <- select(otu_df, otu) %>%
  distinct()%>%
  anti_join(., not_in_any_isolate)%>%
  pull(otu)

# in_multiple_isolates <- otu_df %>%
#   filter(reads != 0)%>%
#   group_by(otu) %>%
#   summarize(num_Isolates = n())

#### Pruned Trees ####

pruned.data <- prune_taxa(otu_in_isolate, all.data)

pruned.ft.data <- prune_taxa(otu_in_isolate, ft.data)

pruned.tree <- plot_tree(pruned.data, label.tips="taxa_names", ladderize="left", color="Class")
pruned.tree

pruned.ft.tree <- plot_tree(pruned.ft.data, nodelabf=nodeplotboot(), label.tips="taxa_names", ladderize="left", color="Class")
pruned.ft.tree

pruned.treeGR <- plot_tree(pruned.ft.data, ladderize = T, color = "grEffect", justify = "left")
pruned.treeGR

pruned.treeCC <- plot_tree(pruned.data, ladderize = T, color = "ccEffect", justify = "left")
pruned.treeCC


#### Useful Functions for Algal Host Level Analysis ####

filter_host <- function(host){
  species_filtered <- filter(sample_key, algalHost == host)%>%
    rownames_to_column(., var = "Isolate")%>%
    select(Isolate)%>%
    left_join(otu_df)%>%
    group_by(Isolate)%>%
    mutate(total_reads = sum(reads),
           percent_total = reads/sum(reads)*100)%>%
    filter(percent_total >= 80)
  
  otus <- pull(species_filtered, otu)
  
  isolates <- select(species_filtered, Isolate)
  
  sample.data <- sample_key %>%
    rownames_to_column(., var="Isolate")%>%
    left_join(isolates, .)%>%
    column_to_rownames(., var="Isolate")%>%
    sample_data(.)
  
  final.data <- prune_taxa(otus, mo.data)%>%
    merge_phyloseq(., sample.data)
  
  return(final.data)
}

grTree <- function(data, host){
  title <- paste(host, "Growth Rate", sep = " ")
  tree <- plot_tree(data, ladderize = T, color = "mean_logNormGR", shape = "Class", justify = "left", title = title)+
    lims(color = c(-1.5, 1.5))
  
  return <- tree
}

ccTree <- function(data, host){
  title <- paste(host, "Carrying Capacity", sep = " ")
  tree <- plot_tree(data, ladderize = T, color = "mean_logNormCC", shape = "Class", justify = "left", title = title)+
    lims(color = c(-1.5, 1.5))
  
  return <- tree
}

#### Chlorella Only ####
cs.data <- filter_host(host = "chlorella")

cs.treeGR <- grTree(cs.data, "Chlorella")
cs.treeGR

cs.treeCC <- ccTree(cs.data, "Chlorella")
cs.treeCC

#### Coelastrum ####
cm.data <- filter_host(host = "coelastrum")

cm.treeGR <- grTree(cm.data, "Coelastrum")
cm.treeGR

cm.treeCC <- ccTree(cm.data, "Coelastrum")
cm.treeCC

#### Scenedesmus ####
sa.data <- filter_host(host = "scenedesmus")

sa.treeGR <- grTree(sa.data, "Scenedesmus")
sa.treeGR

sa.treeCC <- ccTree(sa.data, "Scenedesmus")
sa.treeCC

#### Monoraphidium ####
mm.data <- filter_host(host = "monoraphidium")

mm.treeGR <- grTree(mm.data, "Monoraphidium")
mm.treeGR

mm.treeCC <- ccTree(mm.data, "Monoraphidium")
mm.treeCC

#### Selenastrum ####
sc.data <- filter_host(host = "selenastrum")

sc.treeGR <- grTree(sc.data, "Selenastrum")
sc.treeGR

sc.treeCC <- ccTree(sc.data, "Selenastrum")
sc.treeCC

#### Patchwork ####

patchwork <- (cs.treeGR|cs.treeCC)/(cm.treeGR|cm.treeCC)/(sa.treeGR|sa.treeCC)/(mm.treeGR|mm.treeCC)/(sc.treeGR|sc.treeCC) + plot_layout(guides = "collect")

patchwork

