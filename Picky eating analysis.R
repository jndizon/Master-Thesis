# analyzing how the picky eating total score correlates with the ADHD-RS scores and the diet patterns

library(dplyr) # v1.1.4
library(tidyr) # v1.3.1
library(ggplot2) # v3.5.1
library(broom) # v1.0.7
library(vegan) # v2.6-10
library(phyloseq) # v1.41.1

## ANALYZING PICKY EATING SCORE WITH ADHD-RS
# create input df
picky_adhdrs <- merge(ADHDRS_10y, picky_13y, by = "abcno") 
picky_adhdrs <- picky_adhdrs %>% select(-STARTDATE, -s_331:-s_336) %>% filter(!is.na(total_score))

# spearman correlation 
spear_picky_adhdrs <- picky_adhdrs %>%
  gather(adhdrs, yy, `attention`:`total_score`) %>%  # endpoints > col names for all ADHD-RS scores
  gather(picky, y, picky_score) %>%  # predictors > nutrient intakes
  group_by(adhdrs, picky) %>%
  summarise(
    cor_result = list(tidy(cor.test(yy, y, method = "spearman", exact = FALSE))),
    .groups = "drop"
  ) %>%
  unnest(cor_result) %>%
  select(adhdrs, picky, estimate, p.value)
spear_picky_adhdrs <- spear_picky_adhdrs %>% # Adjust p-values for multiple comparisons
  mutate(p_adjusted = p.adjust(p.value, method = "fdr"))

------------------------------------
  
# plot the results combined with the picky eating spear cors with the diet patterns
  # requires df from "diet pattern analysis"
spear_picky_diet <- FGdietpatterns_cor %>% filter(vars == "picky_score") %>% 
  select(-exposure_1, -exposure_1.1) %>% rename(picky = vars) %>% rename(adhdrs = diet)

spear_picky <- rbind(spear_picky_adhdrs, spear_picky_diet)

# reorder to the order the x-axis should appear
spear_picky$exposure_1 <- factor(spear_picky$adhdrs, 
                                 levels=rev(c("Diet_1FG", "Diet_2FG", "behavior_subscale", 
                                              "implus_hyperact", "attention", "total_score")))
spear_picky$exposure_1.1   <- order(spear_picky$exposure_1, spear_picky$estimate, decreasing=TRUE)
spear_picky <- spear_picky[order(spear_picky$estimate),]

ggplot(spear_picky, aes(x = exposure_1, y = estimate)) +
  geom_segment(aes(xend = adhdrs, y = 0, yend = estimate), color = "darkgoldenrod2", linewidth = 0.75) +
  geom_point(size = 4, color = "darkgoldenrod2") +
  geom_text(aes(label = case_when(
    p_adjusted < 0.001 ~ "***",
    p_adjusted < 0.01  ~ "**",
    p_adjusted < 0.05  ~ "*",
    TRUE         ~ "")), color = "black", nudge_y = 0.01) +
  theme_minimal() +
  scale_y_continuous(limits = c(-0.25, 0.25)) +
  scale_x_discrete(labels = c("total_score" = "Total ADHD-RS Score", 
                              "implus_hyperact" = "Impulsivity/Hyperactivity", 
                              "behavior_subscale" = "Oppositional Behavior", 
                              "socialcircumstances" = "Social Circumstances", 
                              "attention" = "Inattention",
                              "Diet_1FG" = "Varied Diet Pattern",
                              "Diet_2FG" = "Western Diet Pattern")) + # Renaming y-axis labels
  labs(y = "Spearman Correlation Estimate", x = "", title = "Spearman Correlation with Y") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

-------------------------------
  
## ANALYZING PICKY EATING + GM 

sum(!is.na(samp_trimmed_pruned_var$picky_score))
# analysis will include n = 279

## ALPHA DIVERSITY 
# input df: "alpha" from physeq manips
# anova for shannon alpha diversity:
model <- aov(shannon ~ picky_score, data = alpha)
summary(model) 
# anova for observed richness alpha diversity:
model <- aov(obs.richness ~ picky_score, data = alpha)
summary(model) 

-------------------------
  
## BETA DIVERSITY 
# REQUIRES dfs FROM "physeq manips"

# bray-curtis (already done in physeq manips)
BC_ord <- ordinate(phy_rel, "PCoA", "bray") # perform the PCoA ordination
braycurtis <- distance(phy_rel, method="bray")

# run PERMANOVA to estimate the effect picky score has on the entire community
# use adonis func to run PERMANOVA
adonis2(braycurtis ~ picky_score, na.action = na.omit, data = samp_trimmed_pruned_var)
# weighted unifrac
adonis2(unifrac_we ~ picky_score, na.action = na.omit, data = samp_trimmed_pruned_var)

