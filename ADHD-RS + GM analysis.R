## analysis of the ADHD symptom measures with the gut microbiome

# includes alpha + beta diversity and differential abundance

# n = 298

## REQUIRES DF FROM "Phyloseq pre-process + diversity"

library(phyloseq) # v1.41.1
library(dplyr) # v1.1.4
library(ggplot2) # v3.5.1
library(ggpubr) # v0.6.0
library(DAtest) # v2.8.0
library(ggrepel) # v0.9.6
library(vegan) # v2.6-10
library(purrr) # v1.0.4

## ALPHA DIVERSITY 
# input df: "alpha" from Phyloseq pre-process + diversity

# anova for shannon alpha diversity:
model <- aov(shannon ~ attention, data = alpha)
summary(model) 
model <- aov(shannon ~ implus_hyperact, data = alpha)
summary(model)
model <- aov(shannon ~ behavior_subscale, data = alpha)
summary(model) 
model <- aov(shannon ~ total_score, data = alpha)
summary(model)

# anova for observed richness alpha diversity:
model <- aov(obs.richness ~ attention, data = alpha)
summary(model) 
model <- aov(obs.richness ~ implus_hyperact, data = alpha)
summary(model)
model <- aov(obs.richness ~ behavior_subscale, data = alpha)
summary(model) 
model <- aov(obs.richness ~ total_score, data = alpha)
summary(model)

-------------------------
  
## BETA DIVERSITY 
# REQUIRES dfs FROM "Phyloseq pre-process + diversity"
  # run PERMANOVA to estimate the effect ADHD-RS total_score has on the entire community
  # use adonis func to run PERMANOVA
  
# bray-curtis
braycurtis <- distance(phy_rel, method="bray")

adonis2(braycurtis ~ total_score, na.action = na.omit, data = samp_trimmed_pruned_var)
adonis2(braycurtis ~ attention, na.action = na.omit, data = samp_trimmed_pruned_var)
adonis2(braycurtis ~ implus_hyperact, na.action = na.omit, data = samp_trimmed_pruned_var)
adonis2(braycurtis ~ behavior_subscale, na.action = na.omit, data = samp_trimmed_pruned_var)

# weighted unifrac
unifrac_we <- UniFrac(phy_rel, weighted=TRUE) 
adonis2(unifrac_we ~ total_score, na.action = na.omit, data = samp_trimmed_pruned_var)
adonis2(unifrac_we ~ attention, na.action = na.omit, data = samp_trimmed_pruned_var)
adonis2(unifrac_we ~ implus_hyperact, na.action = na.omit, data = samp_trimmed_pruned_var)
adonis2(unifrac_we ~ behavior_subscale, na.action = na.omit, data = samp_trimmed_pruned_var)

-------------------------------
  
# DIFFERENTIAL ABUNDANCE: 

# trimming phyloseq to include only those with adhdrs info
# Identify samples to keep, remove df afterwards
samples_to_keep <- rownames(samp_trimmed_pruned_var)[!is.na(samp_trimmed_pruned_var$total_score)]
# should show n = 298
# Subset the phyloseq object, remove df afterwards
phy_genus_adhdrs <- prune_samples(samples_to_keep, phy_genus)
phy_species_adhdrs <- prune_samples(samples_to_keep, phy_species)
# Verify the new phyloseq object
phy_genus_adhdrs
phy_species_adhdrs
rm(samples_to_keep)

--------------------------------
  
## GENUS LEVEL
# test which differential abundance tests will work best for this dataset + predictor
res <- testDA(phy_genus_adhdrs, predictor = "total_score", relative=F)
summary(res) # llm, lrm, and pea are the top 
res <- testDA(phy_genus_adhdrs, predictor = "attention", relative=F)
summary(res) # lli, llm, pea, lrm, and lim are the top
res <- testDA(phy_genus_adhdrs, predictor = "implus_hyperact", relative=F)
summary(res) # llm, lrm, and pea are the top
res <- testDA(phy_genus_adhdrs, predictor = "behavior_subscale", relative=F)
summary(res) # llm, lli, lrm, and pea are the top

# genus-level DA analysis using pearson correlation
pears_genus_adhdrs <- c("attention", "implus_hyperact", "behavior_subscale", "total_score") %>% set_names() %>%
  map_dfr(~ DA.pea(phy_genus, predictor = .x, relative = FALSE, p.adj = "none") %>%
            select(-Method, -Domain, -pval.adj, -Species, -Cluster, -name), .id = "adhdrs") %>%
  group_by(adhdrs) %>% mutate(pval_adj = p.adjust(pval, method = "fdr")) %>% ungroup() %>% 
  relocate(pval_adj, .after = pval)
  # Apply FDR correction within each ADHDRS score individually

---------------------------
  
## SPECIES LEVEL
res <- testDA(phy_species_adhdrs, predictor = "total_score", relative=F)
summary(res) # llm, pea, and lrm are the top
res <- testDA(phy_species_adhdrs, predictor = "attention", relative=F)
summary(res) # llm, pea, and lrm are the top
res <- testDA(phy_species_adhdrs, predictor = "implus_hyperact", relative=F)
summary(res) # llm, lli, lrm, and pea are the top
res <- testDA(phy_species_adhdrs, predictor = "behavior_subscale", relative=F)
summary(res) # pea and lrm are the top

# species-level DA analysis using pearson correlation
pears_species_adhdrs <- c("attention", "implus_hyperact", "behavior_subscale", "total_score") %>% set_names() %>%
  map_dfr(~ DA.pea(phy_species, predictor = .x, relative = FALSE, p.adj = "none") %>%
            select(-Method, -Domain, -pval.adj, -Cluster, -name), .id = "adhdrs") %>%
  group_by(adhdrs) %>% mutate(pval_adj = p.adjust(pval, method = "fdr")) %>% ungroup() %>% 
  relocate(pval_adj, .after = pval)
# Apply FDR correction within each ADHDRS score individually

# plotting species-level DA of pearson results
p <- ggplot(pears_species_adhdrs, aes(x = cor, y = -log10(pval), color = adhdrs)) +
  theme_bw() +
  geom_point() +
  geom_hline(yintercept = 4) +
  geom_text_repel(data = pears_species_adhdrs_sig,
                  aes(x = cor, y = -log10(pval), label = Species),
                  size = 3) +
  xlab("Pearson correlation estimate") +
  ylab("-log10(P-value)") +
  scale_color_manual(
    values = c(
      "attention" = "brown1", 
      "behavior_subscale" = "chartreuse3",
      "implus_hyperact" = "darkmagenta", 
      "total_score" = "darkgoldenrod1"),  # Customize colors
    labels = c(
      "attention" = "Inattention", 
      "behavior_subscale" = "Oppositional Behavior",
      "implus_hyperact" = "Impulsivity-Hyperactivity", 
      "total_score" = "ADHD-RS Total Score"),
    name = "ADHD-RS\nSymptom\nMeasure")
p

# combine with the volcano plot from diet
pp <- ggarrange(p1,p, labels = c("A", "B"),
                ncol = 2, nrow = 1)
pp
