# loading and pre-process of the phyloseq object produced by Jay Jiang

library(tidyverse) # v2.0.0
library(phyloseq) # v1.41.1
library(ape) # v5.8-1
library(copiome) # v0.1.0

# load the physeq object
phy_10y <- readRDS("phy_10y_share2.rds")

otu <- otu_df(phy_10y) # taxa are columns, samples are rows
tax <- tax_df(phy_10y)
samp <- sample_df(phy_10y)
colnames(samp)
# [1] "SampleID"          "Timepoint"         "abcno"             "Sequencer_ID"      "Flow_Cell_ID"      "Flow_Cell_Barcode"
# [7] "Sample_Barcode"    "Library_Size"      "Avg_Read_Length"   "Total_Bases" 

-----------------------
  
# ROOTING THE PHYLOSEQ OBJECT - code from Jay Jiang
# Extract the phylogenetic tree from your phyloseq object
tree <- phy_tree(phy_10y)

# Extract the taxonomic table from the phyloseq object, check the number of phyla present, 
# and determine which phylum is the most distant from the others.
# In this case, Patescibacteria is the most distant phylum.
tax <- tax_df(phy_10y)

# Identify the tip labels (OTUs/ASVs) that belong to the phylum "Patescibacteria."
# We use rownames(tax_info) to get the taxa names and filter based on the "Phylum" column.
Patescibacteria_tips <- rownames(tax)[tax[, "Phylum"] == "Patescibacteria"]

# Find the most recent common ancestor (MRCA) of the selected Patescibacteria taxa.
# The MRCA is the lowest node in the tree that includes all the selected tips.
Patescibacteria_node <- getMRCA(tree, tip = Patescibacteria_tips)

# Print the node number corresponding to the MRCA of Patescibacteria.
print(Patescibacteria_node) # Example output: 8882 # same output here

# Root the tree at the MRCA of Patescibacteria.
# The 'resolve.root = TRUE' argument ensures that the root is fully bifurcating (i.e., no polytomies remain).
rooted_tree <- root(tree, node = 8882, resolve.root = TRUE)

# Replace the original tree in the phyloseq object with the newly rooted tree.
phy_tree(phy_10y) <- rooted_tree

# check if it worked:
# re-extract the tree
tree <- phy_tree(phy_10y)
is.rooted(tree) 

------------------------------------
  
# trimming physeq object to remove abcno's without food/nutr data
# REQUIRES RUNING "FOOD NUTR MANIPS" AND "VAR10Y AND DIAGNOSES" FIRST
# Identify samples to keep
samples_to_keep <- rownames(samp)[samp$abcno %in% food10y2.0$abcno]
# should show n = 309
# Subset the phyloseq object
phy_trimmed <- prune_samples(samples_to_keep, phy_10y)
# Verify the new phyloseq object
phy_trimmed
rm(samples_to_keep)

----------------------

# prune the data so you only keep the taxa present in at least 10 samples
phy_trimmed_pruned <- filter_taxa(phy_trimmed, function(x) sum(x > 0)>10, TRUE)

# Agglomorate to Genus level for later differential abundance analysis
phy_genus <- phyloseq::tax_glom(phy_trimmed_pruned, "Genus")
# Agglomorate to species level for later differential abundance analysis
phy_species <- phyloseq::tax_glom(phy_trimmed_pruned, "Species") 

--------------------
  
# transforming counts to relative abundance (aka TSS normalization): used for beta diversity
phy_rel <- transform_sample_counts(phy_trimmed_pruned, function(x) x/sum(x))
otu_table(phy_rel)[1:5, 1:10] # check how it looks
sample_sums(phy_rel) # each sample should have a sum of 1

otu_rel <- otu_df(phy_rel)

# mean and sd of species rel abundance across the dataset
taxa_stats <- otu_rel %>%
  as_tibble() %>%
  summarise(across(everything(), list(Mean = mean, SD = sd))) %>%
  pivot_longer(cols = everything(),
               names_to = c("Taxon", ".value"),
               names_sep = "_")

--------------------------
  
# SAMPLE DATA
  
samp_trimmed_pruned <- sample_df(phy_trimmed_pruned)

# first make the dfs compatible for merging with samp data
samp_trimmed_pruned$abcno <- as.character(samp_trimmed_pruned$abcno)
variablesmerged$abcno <- as.character(variablesmerged$abcno)
ADHDRS_10y$abcno <- as.character(ADHDRS_10y$abcno)
food_pc$abcno <- as.character(food_pc$abcno)
picky_13y$abcno <- as.character(picky_13y$abcno)

# create new df for sample data + variables
samp_trimmed_pruned_var <- samp_trimmed_pruned %>%
  left_join(variablesmerged, by = "abcno") %>%
  select(-Sequencer_ID:-Total_Bases) %>%
  left_join(ADHDRS_10y, by = "abcno") %>%
  left_join(food_pc, by = "abcno") %>%
  left_join(picky_13y, by = "abcno") %>%
  select(-s_331:-s_336)
row.names(samp_trimmed_pruned_var) <- samp_trimmed_pruned_var[,1]

# replace old sample data in physeq with the new one
sample_data(phy_rel) <- sample_data(samp_trimmed_pruned_var)
sample_data(phy_genus) <- sample_data(samp_trimmed_pruned_var)
sample_data(phy_species) <- sample_data(samp_trimmed_pruned_var)

-----------------------------
  
# ADHDRS PHYLOSEQ: trimmed to only include samples with ADHDRS data
# Identify samples to keep, remove df afterwards
samples_to_keep <- rownames(samp_trimmed_pruned_var)[]
samples_to_keep <- rownames(samp_trimmed_pruned_var)[!is.na(samp_trimmed_pruned_var$total_score)]
phy_genus_adhdrs <- prune_samples(samples_to_keep, phy_genus)
# Verify the new phyloseq object
phy_genus_adhdrs
rm(samples_to_keep)

------------------------------
  
# creating dfs with general alpha + beta diversity measures for future analysis

# ALPHA DIVERSITY- DO NOT USE PRUNED DATA HERE
# create otu df
otu_trimmed <- otu_df(phy_trimmed)
# relabeling trimmed otu table with abcno
row.names(otu_trimmed) <- samp_trimmed[,3]

# observed richness
obs_rich <- apply(otu_trimmed, MARGIN = 1, function(x) sum(x > 0)) 

# shannon diversity
shannon <- apply(otu_trimmed, MARGIN = 1, function(x) -sum(x/sum(x)*log(x/sum(x)), na.rm = TRUE))
# na.rm = T ensures taxa with 0 abundance are ignored

# create data frame with shannon diversity with sample names
alpha <- merge(shannon, obs_rich, by = "row.names") 
alpha <- setNames(alpha, c("abcno", "shannon", "obs.richness"))

# insert all variables into the alpha df
alpha <- merge(alpha, samp_trimmed_pruned_var, by = "abcno")

------------------
  
## BETA DIVERSITY- USE phy_rel
# jaccard: presence/absence only
jac <- phyloseq::distance(phy_rel, method="binary")
as.matrix(jac)[1:5, 1:5]

# bray curtis: rel abundance
braycurtis <- phyloseq::distance(phy_rel, method="bray")
as.matrix(braycurtis)[1:5, 1:5]

# unifrac distances: based on how phylogenetically similar the present taxa are in one sample
# compared to those in another sample (aka another abcno)
# unweighted: presence/absence only
unifrac_un <- UniFrac(phy_rel, weighted=FALSE) 
as.matrix(unifrac_un)[1:5, 1:5]

# weighted: considers rel abundance
# smallest values with weighted unifrac compared to other metrics ->
# most abundant taxa in the samples are closely related, and the samples are therefore phylogenetically similar
unifrac_we <- UniFrac(phy_rel, weighted=TRUE) 
as.matrix(unifrac_we)[1:5, 1:5]








