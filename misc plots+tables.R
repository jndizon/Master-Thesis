# additional script with the code for misc. tables and plots for the thesis report


# table: DA results for the species included in the PC-DAG
# input: the two dfs of trending DA results from the boot func pcalg script
pears_species_Diet_2FG_trend <- pears_species_Diet_2FG_trend %>% rename(adhdrs = diet)

# set up table:
pears_6species <- rbind(
  pears_species_Diet_2FG_trend[pears_species_Diet_2FG_trend$Feature %in% 
                                 intersect(unique(pears_species_Diet_2FG_trend$Feature), 
                                           unique(pears_species_adhdrs_trend$Feature)), ],
  pears_species_adhdrs_trend[pears_species_adhdrs_trend$Feature %in% 
                               intersect(unique(pears_species_Diet_2FG_trend$Feature), 
                                         unique(pears_species_adhdrs_trend$Feature)), ])
pears_6species <- pears_6species %>% filter(Feature != "Cluster96658") 

# create new table for making the labels that combine multiple columns
pears_6specieslabels <- pears_6species
pears_6specieslabels <- pears_6specieslabels %>% rename(cor2.5 = `cor_2.5%`,
                                                        cor97.5 = `cor_97.5%`)
pears_6specieslabels$label <- paste0(
  round(pears_6specieslabels$cor, 2),
  " [", round(pears_6specieslabels$cor2.5, 2),
  ", ", round(pears_6specieslabels$cor97.5, 2), "]")
pears_6specieslabels <- pears_6specieslabels %>% select(-cor:-cor97.5)

# table: signficant DA results for ADHD-RS
pears_species_adhdrs_sig_labels <- pears_species_adhdrs_sig
pears_species_adhdrs_sig_labels <- pears_species_adhdrs_sig_labels %>% rename(cor2.5 = `cor_2.5%`,
                                                        cor97.5 = `cor_97.5%`)
pears_species_adhdrs_sig_labels$label <- paste0(
  round(pears_species_adhdrs_sig_labels$cor, 2),
  " [", round(pears_species_adhdrs_sig_labels$cor2.5, 2),
  ", ", round(pears_species_adhdrs_sig_labels$cor97.5, 2), "]")
pears_species_adhdrs_sig_labels <- pears_species_adhdrs_sig_labels %>% select(-cor:-cor97.5)

# table: signficant DA results for diet patterns
pears_species_Diet_2FG_sig_labels <- pears_species_Diet_2FG_sig
pears_species_Diet_2FG_sig_labels <- pears_species_Diet_2FG_sig_labels %>% rename(cor2.5 = `cor_2.5%`,
                                                                              cor97.5 = `cor_97.5%`)
pears_species_Diet_2FG_sig_labels$label <- paste0(
  round(pears_species_Diet_2FG_sig_labels$cor, 2),
  " [", round(pears_species_Diet_2FG_sig_labels$cor2.5, 2),
  ", ", round(pears_species_Diet_2FG_sig_labels$cor97.5, 2), "]")
pears_species_Diet_2FG_sig_labels <- pears_species_Diet_2FG_sig_labels %>% select(-cor:-cor97.5)

---------------------------

# alpha/beta div plots using gradients with the continuous variable 
BC_ord <- ordinate(phy_rel, "PCoA", "bray") # perform the PCoA ordination
wunifrac_ord <- ordinate(phy_rel, "PCoA", "UniFrac", weighted = TRUE) # perform the PCoA ordination
  
# total score
alphaS_scatter_tot <- ggplot(alpha, aes(x = total_score, y = shannon)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "darkgoldenrod1") +
  theme_minimal() +
  labs(x = "Total ADHD-RS score", 
       y = "Shannon Diversity")
alphaS_scatter_tot
alphaOR_scatter_tot <- ggplot(alpha, aes(x = total_score, y = obs.richness)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "darkgoldenrod1") +
  theme_minimal() +
  labs(x = "Total ADHD-RS score", 
       y = "Observed Richness")
alphaOR_scatter_tot
BC_plot_tot <- plot_ordination(phy_rel, BC_ord, type = "samples", color = "total_score") %>%
  { .data <- .; .data$data <- na.omit(.data$data); .data } %>%
  `+`(scale_color_gradient2(low = "darkgoldenrod2", mid = "antiquewhite", high = "darkmagenta", midpoint = 25,
                            name = "Total ADHD-RS Score")) %>%
  `+`(ggtitle("Bray-Curtis distances")) %>%
  `+`(theme_minimal())
print(BC_plot_tot)
wunifrac_plot_tot <- plot_ordination(phy_rel, wunifrac_ord, type = "samples", color = "total_score") %>%
  { .data <- .; .data$data <- na.omit(.data$data); .data } %>%
  `+`(scale_color_gradient2(low = "darkgoldenrod2", mid = "antiquewhite", high = "darkmagenta", midpoint = 25,
                            name = "Total ADHD-RS Score")) %>%
  `+`(ggtitle("Weighted UniFrac distances")) %>%
  `+`(theme_minimal())
print(wunifrac_plot_tot)

# inattention
alphaS_scatter_att <- ggplot(alpha, aes(x = attention, y = shannon)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "darkgoldenrod1") +
  theme_minimal() +
  labs(x = "Inattention", 
       y = "Shannon Diversity")
alphaS_scatter_att
alphaOR_scatter_att <- ggplot(alpha, aes(x = attention, y = obs.richness)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "darkgoldenrod1") +
  theme_minimal() +
  labs(x = "Inattention", 
       y = "Observed Richness")
alphaOR_scatter_att
BC_plot_att <- plot_ordination(phy_rel, BC_ord, type = "samples", color = "attention") %>%
  { .data <- .; .data$data <- na.omit(.data$data); .data } %>%
  `+`(scale_color_gradient2(low = "darkgoldenrod2", mid = "antiquewhite", high = "darkmagenta", midpoint = 12,
                            name = "Inattention")) %>%
  `+`(ggtitle("Bray-Curtis distances")) %>%
  `+`(theme_minimal())
print(BC_plot_att)
wunifrac_plot_att <- plot_ordination(phy_rel, wunifrac_ord, type = "samples", color = "attention") %>%
  { .data <- .; .data$data <- na.omit(.data$data); .data } %>%
  `+`(scale_color_gradient2(low = "darkgoldenrod2", mid = "antiquewhite", high = "darkmagenta", midpoint = 12,
                            name = "Inattention")) %>%
  `+`(ggtitle("Weighted UniFrac distances")) %>%
  `+`(theme_minimal())
print(wunifrac_plot_att)

# impuls hyperact
alphaS_scatter_imphyp <- ggplot(alpha, aes(x = implus_hyperact, y = shannon)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "darkgoldenrod1") +
  theme_minimal() +
  labs(x = "Impulsivity/Hyperactivity", 
       y = "Shannon Diversity")
alphaS_scatter_imphyp
alphaOR_scatter_imphyp <- ggplot(alpha, aes(x = implus_hyperact, y = obs.richness)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "darkgoldenrod1") +
  theme_minimal() +
  labs(x = "Impulsivity/Hyperactivity", 
       y = "Observed Richness")
alphaOR_scatter_imphyp
BC_plot_imphyp <- plot_ordination(phy_rel, BC_ord, type = "samples", color = "implus_hyperact") %>%
  { .data <- .; .data$data <- na.omit(.data$data); .data } %>%
  `+`(scale_color_gradient2(low = "darkgoldenrod2", mid = "antiquewhite", high = "darkmagenta", midpoint = 12,
                            name = "Impulsivity/Hyperactivity")) %>%
  `+`(ggtitle("Bray-Curtis distances")) %>%
  `+`(theme_minimal())
print(BC_plot_imphyp)
wunifrac_plot_imphyp <- plot_ordination(phy_rel, wunifrac_ord, type = "samples", color = "implus_hyperact") %>%
  { .data <- .; .data$data <- na.omit(.data$data); .data } %>%
  `+`(scale_color_gradient2(low = "darkgoldenrod2", mid = "antiquewhite", high = "darkmagenta", midpoint = 12,
                            name = "Impulsivity/Hyperactivity")) %>%
  `+`(ggtitle("Weighted UniFrac distances")) %>%
  `+`(theme_minimal())
print(wunifrac_plot_imphyp)

# opp beh
alphaS_scatter_beh <- ggplot(alpha, aes(x = behavior_subscale, y = shannon)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "darkgoldenrod1") +
  theme_minimal() +
  labs(x = "Oppositional Behavior", 
       y = "Shannon Diversity")
alphaS_scatter_beh
alphaOR_scatter_beh <- ggplot(alpha, aes(x = behavior_subscale, y = obs.richness)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "darkgoldenrod1") +
  theme_minimal() +
  labs(x = "Oppositional Behavior", 
       y = "Observed Richness")
alphaOR_scatter_beh
BC_plot_beh <- plot_ordination(phy_rel, BC_ord, type = "samples", color = "behavior_subscale") %>%
  { .data <- .; .data$data <- na.omit(.data$data); .data } %>%
  `+`(scale_color_gradient2(low = "darkgoldenrod2", mid = "antiquewhite", high = "darkmagenta", midpoint = 10, 
                            name = "Oppositional Behavior")) %>%
  `+`(ggtitle("Bray-Curtis distances")) %>%
  `+`(theme_minimal())
print(BC_plot_beh)
wunifrac_plot_beh <- plot_ordination(phy_rel, wunifrac_ord, type = "samples", color = "behavior_subscale") %>%
  { .data <- .; .data$data <- na.omit(.data$data); .data } %>%
  `+`(scale_color_gradient2(low = "darkgoldenrod2", mid = "antiquewhite", high = "darkmagenta", midpoint = 10,
                            name = "Oppositional Behavior")) %>%
  `+`(ggtitle("Weighted UniFrac distances")) %>%
  `+`(theme_minimal())
print(wunifrac_plot_beh)

-------------------------------------

# alpha/beta div plots for diet patterns
# Diet_1FG
alphaS_scatter_1FG <- ggplot(alpha, aes(x = Diet_1FG, y = shannon)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "darkgoldenrod1") +
  theme_minimal() +
  labs(x = "Varied Diet Pattern", 
       y = "Shannon Diversity")
alphaS_scatter_1FG
alphaOR_scatter_1FG <- ggplot(alpha, aes(x = Diet_1FG, y = obs.richness)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "darkgoldenrod1") +
  theme_minimal() +
  labs(x = "Varied Diet Pattern", 
       y = "Observed Richness")
alphaOR_scatter_1FG
BC_plot_1FG <- plot_ordination(phy_rel, BC_ord, type="samples", color="Diet_1FG") +
  scale_color_gradient2(low = "darkgoldenrod2", mid = "antiquewhite", high = "deeppink2", midpoint = 0,
                        name = "Varied Diet") +
  ggtitle("Bray-Curtis distances") +
  theme_minimal() 
print(BC_plot_1FG)
wunifrac_plot_1FG <- plot_ordination(phy_rel, wunifrac_ord, type="samples", color="Diet_1FG") +
  scale_color_gradient2(low = "darkgoldenrod2", mid = "antiquewhite", high = "deeppink2", midpoint = 0,
                        name = "Varied Diet") +
  ggtitle("Weighted UniFrac distances") +
  theme_minimal() 
print(wunifrac_plot_1FG)

# Diet_2FG
alphaS_scatter_2FG <- ggplot(alpha, aes(x = Diet_2FG, y = shannon)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "darkgoldenrod1") +
  theme_minimal() +
  labs(x = "Western Diet Pattern", 
       y = "Shannon Diversity")
alphaS_scatter_2FG
alphaOR_scatter_2FG <- ggplot(alpha, aes(x = Diet_2FG, y = obs.richness)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "darkgoldenrod1") +
  theme_minimal() +
  labs(x = "Western Diet Pattern", 
       y = "Observed Richness")
alphaOR_scatter_2FG
BC_plot_2FG <- plot_ordination(phy_rel, BC_ord, type="samples", color="Diet_2FG") +
  scale_color_gradient2(low = "darkgoldenrod2", mid = "antiquewhite", high = "darkturquoise", midpoint = -0.5,
                        name = "Western Diet") +
  ggtitle("Bray-Curtis distances") +
  theme_minimal() 
print(BC_plot_2FG)
wunifrac_plot_2FG <- plot_ordination(phy_rel, wunifrac_ord, type="samples", color="Diet_2FG") +
  scale_color_gradient2(low = "darkgoldenrod2", mid = "antiquewhite", high = "darkturquoise", midpoint = -0.5,
                        name = "Western Diet") +
  ggtitle("Weighted UniFrac distances") +
  theme_minimal() 
print(wunifrac_plot_2FG)

------------------------------

# combined alpha div plots
combined_alpha <- ggarrange(alphaOR_scatter_1FG, alphaS_scatter_1FG, 
                            alphaOR_scatter_2FG, alphaS_scatter_2FG,
                            alphaOR_scatter_tot, alphaS_scatter_tot, 
                                  labels = c("A", "B", "C", "D", "E", "F"),
                                  ncol = 2, nrow = 3)
print(combined_alpha)

combined_alpha_supp <- ggarrange(alphaOR_scatter_att, alphaS_scatter_att, 
                            alphaOR_scatter_imphyp, alphaS_scatter_imphyp,
                            alphaOR_scatter_beh, alphaS_scatter_beh, 
                            labels = c("A", "B", "C", "D", "E", "F"),
                            ncol = 2, nrow = 3)
print(combined_alpha_supp)

-------------------------------

# combined beta div plots
combined_beta_diet_1FG <- ggarrange(BC_plot_1FG, wunifrac_plot_1FG, 
                                      labels = c("A", "B"),
                                      ncol = 2, nrow = 1,
                                      common.legend = TRUE, legend = "bottom")
print(combined_beta_diet_1FG)

combined_beta_diet_2FG <- ggarrange(BC_plot_2FG, wunifrac_plot_2FG, 
                          labels = c("C", "D"),
                          ncol = 2, nrow = 1,
                          common.legend = TRUE, legend = "bottom")
print(combined_beta_diet_2FG)


combined_beta_tot <- ggarrange(BC_plot_tot, wunifrac_plot_tot,  
                                    labels = c("E", "F"),
                                    ncol = 2, nrow = 1,
                                    common.legend = TRUE, legend = "bottom")
print(combined_beta_tot)

combined_beta_att <- ggarrange(BC_plot_att, wunifrac_plot_att,
                                    labels = c("A", "B"),
                                    ncol = 2, nrow = 1,
                                    common.legend = TRUE, legend = "bottom")
print(combined_beta_att)

combined_beta_imphyp <- ggarrange(BC_plot_imphyp, wunifrac_plot_imphyp,
                                    labels = c("C", "D"),
                                    ncol = 2, nrow = 1,
                                    common.legend = TRUE, legend = "bottom")
print(combined_beta_imphyp)

combined_beta_beh <- ggarrange(BC_plot_beh, wunifrac_plot_beh, 
                               labels = c("E", "F"),
                               ncol = 2, nrow = 1,
                               common.legend = TRUE, legend = "bottom")
print(combined_beta_beh)

