## Using PC algorithm for causal discovery 
  # learns causal structure from the data > infers structure from statistical dependencies
  # attempts to infer direction

# first model: western diet + picky eating score + sex + socialcircumstances + ADHD-RS total score
  # + 6 species that represent the GM component

# second model: without the GM component

# n = 269

library(pcalg) # v2.7-12
library(igraph) # v2.1.4
library(tidyr) # v1.3.1
library(tibble) # v3.2.1
library(phyloseq) # v1.41.1
library(dplyr) # v1.1.4
library(purrr) # v1.0.4

----------------------------

## SETTING UP A GM VARIABLE: identify specific bacterial species that are 
  # increased/decreased with both western diet and ADHD symptoms
# use dfs created from the species-level DA analysis
pears_species_Diet_2FG_trend <- pears_species_Diet_2FG %>% filter(pval <= 0.05)
pears_species_adhdrs_trend <- pears_species_adhdrs %>% filter(pval <= 0.05)

# pull names of the ones in common
common_clusters <- intersect(pears_species_Diet_2FG_trend$Feature, pears_species_adhdrs_trend$Feature)
# 6/7 of these common clusters have associations going in the same directions -> these will be used

# create an otu table with the abundances of only these taxa
# transform species-level physeq object to relative abundances
phy_species_rel <- transform_sample_counts(phy_species, function(x) x/sum(x))
otu_6species <- otu_df(phy_species_rel) %>% select(Cluster11693, Cluster14747, Cluster1493, Cluster2217, Cluster23561, Cluster31850)
otu_6species$SampleID <- rownames(otu_6species) %>% as.character(otu_6species$SampleID)

otu_6species <- otu_6species %>%
  merge(samp_trimmed_pruned_var %>% mutate(SampleID = as.character(SampleID)) %>% select(SampleID, abcno), by = "SampleID") %>%
  select(abcno, everything(), -SampleID)

--------------------
  
# CREATE INPUT DFs

# input df: GM species individually + covariates
causdis_input_cov1 <- samp_trimmed_pruned_var %>%
  select(abcno, sex, socialcircumstances,
         total_score, Diet_2FG, picky_score) %>%
  mutate(sex = ifelse(sex=="male",1,0)) %>%
  left_join(otu_6species, by = "abcno") %>%
  rename(Sex = sex,
         Social_circumstances = socialcircumstances,
         ADHDRS_score = total_score,
         Western_diet = Diet_2FG,
         Picky_score = picky_score,
         Odoribacter_splanchnicus = Cluster2217,
         Pauljensenia_sp001838165 = Cluster11693,
         Allisonella_histaminiformans = Cluster14747,
         CAG_273_sp000437855 = Cluster1493,
         Enterenecus_sp900549885 = Cluster23561,
         UBA4285_sp900545225 = Cluster31850) %>% 
  na.omit() %>%
  mutate(across(-c(abcno, Sex), ~ scale(.)[,1])) %>%
  {rownames(.) <- .$abcno; select(., -abcno)}

# input df: no gm included
causdis_input_cov2 <- causdis_input_cov %>% select(-GM_6Species)


---------------------------------
  
# CREATE BOOTSTRAP FUNCTION: 
bootstrap_func <- function(data, R = 1000, alpha = 0.05, verbose = FALSE) {
  n <- nrow(data)
  p <- ncol(data)
  labels <- colnames(data)
  run_bootstrap <- function(i) {
    if (verbose && i %% 50 == 0) message("Bootstrap sample: ", i)
    # Bootstrap sample
    boot_data <- data[sample(n, n, replace = TRUE), ]
    suffStat <- list(C = cor(boot_data), n = n)
    # PC algorithm
    pc_fit <- pcalg::pc(suffStat, indepTest = gaussCItest, alpha = alpha, labels = labels, verbose = FALSE)
    adj <- as(pc_fit@graph, "matrix")
    # Partial correlation weights for existing edges
    inds <- which(upper.tri(adj) & adj == 1, arr.ind = TRUE)
    weight_mat <- matrix(0, p, p)
    if (nrow(inds) > 0) {
      map2(inds[, 1], inds[, 2], ~{
        val <- pcorOrder(.x, .y, k = NULL, C = suffStat$C)
        weight_mat[.x, .y] <<- val
        weight_mat[.y, .x] <<- val})}
    list(adj = adj, weights = weight_mat)}
  # Run bootstraps
  results <- map(1:R, run_bootstrap)
  # Stack results
  adj_array <- simplify2array(map(results, "adj"))       # p x p x R
  weight_array <- simplify2array(map(results, "weights"))# p x p x R
  # Compute edge frequencies
  edge_frequency <- apply(adj_array, c(1, 2), mean)
  # Compute means and CIs
  get_stats <- function(i, j) {
    vals <- weight_array[i, j, ]
    vals <- vals[vals != 0]
    if (length(vals) > 1) {
      tibble(mean = mean(vals), ci_low = quantile(vals, 0.025), ci_high = quantile(vals, 0.975))} 
    else {tibble(mean = 0, ci_low = NA, ci_high = NA)}}
  edge_stats <- expand_grid(from = 1:p, to = 1:p) %>%
    mutate(stats = map2(from, to, get_stats)) %>%
    unnest_wider(stats)
  list(edge_frequency = edge_frequency,
       edge_stats_tibble = edge_stats)}

-----------------------------

# MODEL 8: GM SPECIES INDIVIDUALLY    
set.seed(00)
model8 <- bootstrap_func(causdis_input_cov1, R = 1000, alpha = 0.05, verbose = TRUE)
model8_edgeweights <- as.data.frame(model8$edge_stats_tibble) %>% mutate(from = colnames(causdis_input_cov1)[from],
                                                                         to = colnames(causdis_input_cov1)[to])
model8_edgefreqs <- as.data.frame(model8$edge_frequency) 
model8_edgefreqs_filtered <- as.data.frame(as.table(model8$edge_frequency)) %>%
  rename(from = Var1, to = Var2, frequency = Freq) %>%
  filter(from != to, frequency >= 0.5) # edge is present in at least 50% of runs

# Join weights and CIs from edge_stats_tibble
edge_df8 <- model8_edgefreqs_filtered %>%
  left_join(model8_edgeweights, by = c("from", "to")) %>%
  mutate(from = as.character(from), to = as.character(to), mean = as.numeric(mean)) %>%
  mutate(across(where(is.character), ~ recode(.x, "ADHDRS_score" = "ADHD-RS\nTotal\nScore",
                                              "Sex" = "Male\nSex",
                                              "Picky_score" = "Picky\nEating",
                                              "Social_circumstances" = "Social\nCircumstances",
                                              "Western_diet" = "Western\nDiet\nPattern",
                                              "CAG_273_sp000437855" = "CAG-273\nsp000437855",
                                              "UBA4285_sp900545225" = "UBA4285\nsp900545225")))

# Create igraph object
model8_plotmetrics <- graph_from_data_frame(edge_df8, directed = TRUE)

# Edge labels
E(model8_plotmetrics)$label <- paste0(
  round(E(model8_plotmetrics)$mean, 2),
  " [", round(E(model8_plotmetrics)$ci_low, 2),
  ", ", round(E(model8_plotmetrics)$ci_high, 2), "]")
# Vertex colors
V(model8_plotmetrics)$color[V(model8_plotmetrics)$name %in% c("Male\nSex", "Social\nCircumstances")] <- "darkolivegreen3"
V(model8_plotmetrics)$color[V(model8_plotmetrics)$name %in% c("Western\nDiet\nPattern", "Picky\nEating", "ADHD-RS\nTotal\nScore")] <- "darkorchid1"
V(model8_plotmetrics)$color[V(model8_plotmetrics)$name %in% c("CAG-273\nsp000437855", "UBA4285\nsp900545225")] <- "cadetblue"


plot(
  model8_plotmetrics,
  edge.label = E(model8_plotmetrics)$label, 
  edge.label.cex = 0.7,   
  edge.label.color = "black", 
  edge.width = 1.8,  
  vertex.label.cex = 0.8,  
  vertex.label.color = "black",
  vertex.size = 50) 

-----------------------------------------

# model 10: no gm included at all
set.seed(00)
model10 <- bootstrap_func(causdis_input_cov2, R = 1000, alpha = 0.05, verbose = TRUE)
model10_edgeweights <- as.data.frame(model10$edge_stats_tibble) %>% mutate(from = colnames(causdis_input_cov2)[from],
                                                                         to = colnames(causdis_input_cov2)[to])
model10_edgefreqs <- as.data.frame(model10$edge_frequency) 
model10_edgefreqs_filtered <- as.data.frame(as.table(model10$edge_frequency)) %>%
  rename(from = Var1, to = Var2, frequency = Freq) %>%
  filter(from != to, frequency >= 0.5) # edge is present in at least 50% of runs

# Join weights and CIs from edge_stats_tibble
edge_df10 <- model10_edgefreqs_filtered %>%
  left_join(model10_edgeweights, by = c("from", "to")) %>%
  mutate(from = as.character(from), to = as.character(to), mean = as.numeric(mean)) %>%
  mutate(across(where(is.character), ~ recode(.x, "ADHDRS_score" = "ADHD-RS\nTotal\nScore",
                                              "Sex" = "Male\nSex",
                                              "Picky_score" = "Picky\nEating",
                                              "Social_circumstances" = "Social\nCircumstances",
                                              "Western_diet" = "Western\nDiet\nPattern")))

# Create igraph object
model10_plotmetrics <- graph_from_data_frame(edge_df10, directed = TRUE)

# Edge labels
E(model10_plotmetrics)$label <- paste0(
  round(E(model10_plotmetrics)$mean, 2),
  " [", round(E(model10_plotmetrics)$ci_low, 2),
  ", ", round(E(model10_plotmetrics)$ci_high, 2), "]")
# Vertex colors
V(model10_plotmetrics)$color[V(model10_plotmetrics)$name %in% c("Male\nSex", "Social\nCircumstances")] <- "darkolivegreen3"
V(model10_plotmetrics)$color[V(model10_plotmetrics)$name %in% c("Western\nDiet\nPattern", "Picky\nEating", "ADHD-RS\nTotal\nScore")] <- "darkorchid1"

plot(
  model10_plotmetrics,
  layout = layout_as_tree(model10_plotmetrics, root = "ADHD-RS\nTotal\nScore", mode = "out"),
  edge.label = E(model10_plotmetrics)$label,
  edge.label.cex = 0.7, 
  edge.label.color = "black",
  edge.width = 2,    
  vertex.label.cex = 0.8, 
  vertex.label.color = "black", 
  vertex.size = 65)    

------------------------------

# dfs for tables for report
edge_df8_expanded <- model8_edgefreqs %>%
  as.data.frame() %>%
  rownames_to_column("from") %>%
  pivot_longer(-from, names_to = "to", values_to = "edge_frequency") %>% 
  left_join(model8_edgeweights, by = c("from", "to")) %>% filter(!is.na(ci_low))
edge_df8_expanded$label <- paste0(
  round(edge_df8_expanded$mean, 2),
  " [", round(edge_df8_expanded$ci_low, 2),
  ", ", round(edge_df8_expanded$ci_high, 2), "]")
edge_df8_expanded <- edge_df8_expanded %>% select(-mean:-ci_high)

edge_df10_expanded <- model10_edgefreqs %>%
  as.data.frame() %>%
  rownames_to_column("from") %>%
  pivot_longer(-from, names_to = "to", values_to = "edge_frequency") %>% 
  left_join(model10_edgeweights, by = c("from", "to")) %>% filter(!is.na(ci_low))
edge_df10_expanded$label <- paste0(
  round(edge_df10_expanded$mean, 2),
  " [", round(edge_df10_expanded$ci_low, 2),
  ", ", round(edge_df10_expanded$ci_high, 2), "]")
edge_df10_expanded <- edge_df10_expanded %>% select(-mean:-ci_high)






