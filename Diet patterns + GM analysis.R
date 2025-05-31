## analysis of the diet patterns created in "FG level PCs" with the gut microbiome
## REQUIRES DF FROM "Phyloseq pre-process + diversity"

# includes alpha + beta diversity and differential abundance

# n = 309

library(DAtest) # v2.8.0
library(ggrepel) # v0.9.6
library(vegan) # v2.6-10
library(purrr) # v1.0.4

-------------------------

## ALPHA DIVERSITY + DIET PATTERNS
# create input df: "alpha" from physeq manips

# anova for shannon alpha diversity:
model <- aov(shannon ~ Diet_1FG, data = alpha)
summary(model) # Pr(>F) 0.0422 *
model <- aov(shannon ~ Diet_2FG, data = alpha)
summary(model)

# anova for observed richness alpha diversity:
model <- aov(obs.richness ~ Diet_1FG, data = alpha)
summary(model)
model <- aov(obs.richness ~ Diet_2FG, data = alpha)
summary(model)

----------------------------------
  
## BETA DIVERSITY + DIET PATTERNS
# bray curtis: none signif
adonis2(braycurtis ~ Diet_1FG, na.action = na.exclude, data = samp_trimmed_pruned_var) # varied
adonis2(braycurtis ~ Diet_2FG, na.action = na.exclude, data = samp_trimmed_pruned_var) # unhealthy diet

# weighted unifrac: none signif
adonis2(unifrac_we ~ Diet_1FG, na.action = na.exclude, data = samp_trimmed_pruned_var) # varied
adonis2(unifrac_we ~ Diet_2FG, na.action = na.exclude, data = samp_trimmed_pruned_var) # unhealthy diet
              
-----------------------------

# DIFFERENTIAL ABUNDANCE: GENUS LEVEL
# test which differential abundance tests will work best for this dataset + predictor
res <- testDA(phy_genus, predictor = "Diet_1FG", relative=F)
summary(res)
res <- testDA(phy_genus, predictor = "Diet_2FG", relative=F)
summary(res)

# using pearson because it is top 3 test for both diets (at genus level)
pears_genus_diet <- c("Diet_1FG", "Diet_2FG") %>% set_names() %>%
  map_dfr(~ DA.pea(phy_genus, predictor = .x, relative = FALSE, p.adj = "none") %>%
            select(-Method, -Domain, -pval.adj, -Species, -Cluster, -name), .id = "diet") %>%
  group_by(diet) %>% mutate(pval_adj = p.adjust(pval, method = "fdr")) %>% ungroup() %>% 
  relocate(pval_adj, .after = pval)
# Apply FDR correction within each diet individually

pears_genus_diet_subset <- pears_genus_diet %>% select(-Phylum:-Family) %>% 
  relocate(pval_adj, .after = pval) %>% filter(pval_adj < 0.05)

# DIFFERENTIAL ABUNDANCE: SPECIES LEVEL
# test which differential abundance tests will work best for this dataset + predictor
res <- testDA(phy_species, predictor = "Diet_1FG", relative=F)
summary(res) # llm, lrm, and pea are the top
res <- testDA(phy_species, predictor = "Diet_2FG", relative=F)
summary(res) # spe, lli, llm, lrm, and pea are the top

# using pearson because it is top 3 test for both diets (at species level)
pears_species_diet <- c("Diet_1FG", "Diet_2FG") %>% set_names() %>%
  map_dfr(~ DA.pea(phy_species, predictor = .x, relative = FALSE, p.adj = "none") %>%
            select(-Method, -Domain, -pval.adj, -Cluster, -name), .id = "diet") %>%
  group_by(diet) %>% mutate(pval_adj = p.adjust(pval, method = "fdr")) %>% ungroup() %>% 
  relocate(pval_adj, .after = pval)
    # Apply FDR correction within each diet individually

# plotting DA 
pears_species_Diet_2FG <- pears_species_diet %>% filter(diet != "Diet_1FG")
pears_species_Diet_2FG_sig <- pears_species_Diet_2FG %>% filter(pval_adj <= 0.05)

p1 <- ggplot(pears_species_diet, aes(x = cor, y = -log10(pval), color = diet)) +
  theme_bw() +
  geom_point() +
  geom_hline(yintercept = 3.552842) +
  geom_text_repel(data = pears_species_Diet_2FG_sig,
                  aes(x = cor, y = -log10(pval), label = Species),
                  size = 3) +
  xlab("Pearson correlation estimate") +
  ylab("-log10(P-value)") +
  scale_color_manual(
    values = c("Diet_2FG" = "darkturquoise",
               "Diet_1FG" = "deeppink3"),  
    labels = c("Diet_2FG" = "Western Diet",
               "Diet_1FG" = "Varied Diet"),
    name = "Diet\nPattern") +
  theme(axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 10), # labels of the tick marks
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10))

p1

