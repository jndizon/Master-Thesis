# covariate analysis on beta diversity 
# requires samp_trimmed_pruned_var from "physeq manips"

library(phyloseq) # v1.41.1
library(vegan) # v2.6-10
library(ggplot2) # v3.5.1
library(ggpubr) # v0.6.0

# BASELINE CHARACTERISTICS
# sex 
sum(samp_trimmed_pruned_var$sex == "male")
(165/309)*100 
# race
sum(samp_trimmed_pruned_var$race == "caucasian")
(300/309)*100
# social circumstances
mean(samp_trimmed_pruned_var$socialcircumstances, na.rm = TRUE)
sd(samp_trimmed_pruned_var$socialcircumstances, na.rm = TRUE)
# bmi
mean(samp_trimmed_pruned_var$zbmi, na.rm = TRUE)
sd(samp_trimmed_pruned_var$zbmi, na.rm = TRUE)
# total daily energy consumption
mean(samp_trimmed_pruned_var$energi, na.rm = TRUE)
sd(samp_trimmed_pruned_var$energi, na.rm = TRUE)

------------------------
  
# covariate analysis with PERMANOVA
braycurtis <- distance(phy_rel, method="bray")
# sex:
adonis2(braycurtis ~ sex, data = samp_trimmed_pruned_var)
# race
adonis2(braycurtis ~ race, data = samp_trimmed_pruned_var)
#socialcircumstances
adonis2(braycurtis ~ socialcircumstances, data = samp_trimmed_pruned_var)
# zbmi 
adonis2(braycurtis ~ zbmi, na.action = na.exclude, data = samp_trimmed_pruned_var)
# total energy intake 
## none of the five tested covariates affect the beta diversity -> don't need to correct/adjust for

--------------------------------------

##QC ANALYSIS OF THE MICROBIOME DATA

# testing alpha diversity + antibiotics at birth to mother. this should show that 
  # alpha div is lower in those whose mother received antibiotics
alphaOR_abmotheratbirth <- ggplot(alpha, aes(x = abmotheratbirth, y = obs.richness)) +
  geom_point()
alphaOR_abmotheratbirth
alphaOR_abmotheratbirth <- ggplot(alpha, aes(x = abmotheratbirth, y = obs.richness)) +
  theme_bw() +
  geom_boxplot() +
  xlab("abmotheratbirth") +
  ylab("Observed Richness") +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", 
                     size = 5)
alphaOR_abmotheratbirth

ORanova <- aov(obs.richness ~ abmotheratbirth, data = alpha)
summary(ORanova) # signficant

Sanova <- aov(shannon ~ abmotheratbirth, data = alpha)
summary(Sanova) # significant



