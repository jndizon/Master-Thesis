## analysis of associations between ADHD symptoms and diet patterns
# input data: dfs from "FG level PCs", "10y vars and diagnoses", and "ADHDRS analysis"


# bootstrapped Quasi-Poisson models
# outcome = ADHDRS // predictor = diet PCs
# model A: univariate, n = 575
# model B: univariate, removing ADHD+HCD diagnoses, n = 533
# model C: multivariate, adjusted for sex, age at ADHD-RS response, 
  # total energy intake, and social circumstances, n = 575

# robustness analysis
# Quasi-Poisson models with train-test split
  # n = 403 training // n = 172 testing
# spearman correlation

library(dplyr) # v1.1.4
library(ggplot2) # v3.5.1
library(ggpubr) # v0.6.0
library(boot) # v1.3-31
library(broom) # v1.0.7
library(tidyr) # v1.3.1
library(caret) # v7.0-1

--------------------------------

# create input data df: also made in "diet pattern analysis"
FGdietpatterns_adhdrs <- merge(food_pc, ADHDRS_10y, by = "abcno") %>% select(-total_wo_behavior) %>% filter(!is.na(total_score))
FGdietpatterns_adhdrs_covar <- merge(FGdietpatterns_adhdrs, variablesmerged, by = "abcno") %>% select(-bmi, -Overweight, -abmotheratbirth)
FGdietpatterns_adhdrs_covar <- FGdietpatterns_adhdrs_covar %>% mutate(across(starts_with("Diet_"), as.numeric))

------------------------------------

# raw data plots:
# varied diet + adhd-rs:
# Diet_1FG + total_score
vtot <- ggplot(FGdietpatterns_adhdrs, aes(x = Diet_1FG, y = total_score)) +
  geom_point() +
  theme_minimal() +
  geom_smooth(method = "lm", se = FALSE, color = "deeppink2") +
  labs(x = "Varied Diet Pattern", 
       y = "ADHD-RS Total Score") +
  theme(axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 10), # labels of the tick marks
        axis.text.y = element_text(size = 10))
# Diet_1FG + attention  
vatt <- ggplot(FGdietpatterns_adhdrs, aes(x = Diet_1FG, y = attention)) +
  geom_point() +
  theme_minimal() +
  geom_smooth(method = "lm", se = FALSE, color = "deeppink2") +
  labs(x = "Varied Diet Pattern", 
       y = "Inattention") +
  ylim(0,27) +
  theme(axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 10), # labels of the tick marks
        axis.text.y = element_text(size = 10))
# Diet_1FG + implus_hyperact  
vimphyp <- ggplot(FGdietpatterns_adhdrs, aes(x = Diet_1FG, y = implus_hyperact)) +
  geom_point() +
  theme_minimal() +
  geom_smooth(method = "lm", se = FALSE, color = "deeppink2") +
  labs(x = "Varied Diet Pattern", 
       y = "Impulsivity/Hyperactivity") +
  ylim(0,27) +
  theme(axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 10), # labels of the tick marks
        axis.text.y = element_text(size = 10))
# Diet_1FG + behavior_subscale
vbeh <- ggplot(FGdietpatterns_adhdrs, aes(x = Diet_1FG, y = behavior_subscale)) +
  geom_point() +
  theme_minimal() +
  geom_smooth(method = "lm", se = FALSE, color = "deeppink2") +
  labs(x = "Varied Diet Pattern", 
       y = "Oppositional Behavior") +
  ylim(0,27) +
  theme(axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 10), # labels of the tick marks
        axis.text.y = element_text(size = 10))

# western diet + adhd-rs:
# Diet_2FG + total_score
wtot <- ggplot(FGdietpatterns_adhdrs, aes(x = Diet_2FG, y = total_score)) +
  geom_point() +
  theme_minimal() +
  geom_smooth(method = "lm", se = FALSE, color = "darkturquoise") +
  labs(x = "Western Diet Pattern", 
       y = "ADHD-RS Total Score") +
  theme(axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 10), # labels of the tick marks
        axis.text.y = element_text(size = 10))
# Diet_2FG + attention  
watt <- ggplot(FGdietpatterns_adhdrs, aes(x = Diet_2FG, y = attention)) +
  geom_point() +
  theme_minimal() +
  geom_smooth(method = "lm", se = FALSE, color = "darkturquoise") +
  labs(x = "Western Diet Pattern", 
       y = "Inattention") +
  ylim(0,27) +
  theme(axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 10), # labels of the tick marks
        axis.text.y = element_text(size = 10))
# Diet_2FG + implus_hyperact  
wimphyp <- ggplot(FGdietpatterns_adhdrs, aes(x = Diet_2FG, y = implus_hyperact)) +
  geom_point() +
  theme_minimal() +
  geom_smooth(method = "lm", se = FALSE, color = "darkturquoise") +
  labs(x = "Western Diet Pattern", 
       y = "Impulsivity/Hyperactivity") +
  ylim(0,27) +
  theme(axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 10), # labels of the tick marks
        axis.text.y = element_text(size = 10))
# Diet_2FG + behavior_subscale
wbeh <- ggplot(FGdietpatterns_adhdrs, aes(x = Diet_2FG, y = behavior_subscale)) +
  geom_point() +
  theme_minimal() +
  geom_smooth(method = "lm", se = FALSE, color = "darkturquoise") +
  labs(x = " Western Diet Pattern", 
       y = "Oppositional Behavior") +
  ylim(0,27) +
  theme(axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 10), # labels of the tick marks
        axis.text.y = element_text(size = 10))

# combined plots
combinedplots_rawtot <- ggarrange(vtot, wtot, 
                                   labels = c("A", "B"),
                                   ncol = 2, nrow = 1)
print(combinedplots_rawtot)

combinedplots_rawsupp1.1 <- ggarrange(vatt, watt, vimphyp, wimphyp,
                                    labels = c("A", "B", "C", "D"),
                                    ncol = 2, nrow = 2)
print(combinedplots_rawsupp1.1)

combinedplots_rawsupp1.2 <- ggarrange(vimphyp, wimphyp,vbeh, wbeh,
                                      labels = c("C", "D", "E", "F"),
                                      ncol = 2, nrow = 2)
print(combinedplots_rawsupp1.2)



combinedplots_rawsupp1 <- ggarrange(vatt, watt,
                                     labels = c("A", "B"),
                                     ncol = 2, nrow = 1)
print(combinedplots_rawsupp1)
combinedplots_rawsupp2 <- ggarrange(vimphyp, wimphyp,
                                   labels = c("C", "D"),
                                   ncol = 2, nrow = 1)
print(combinedplots_rawsupp2)
combinedplots_rawsupp3 <- ggarrange(vbeh, wbeh, 
                                   labels = c("E", "F"),
                                   ncol = 2, nrow = 1)
print(combinedplots_rawsupp3)


----------------------------

## QUASI POISSON MODELS FOR DIET PATTERNS + ADHDRS
  # bootstrapping allows us to have CI's around the rate ratio and gives more certainty to findings
  
# Define outcomes and diet predictors
outcomes <- c("attention", "implus_hyperact", "total_score", "behavior_subscale")
diet_patterns <- c("Diet_1FG", "Diet_2FG")

# create bootstrapping function 
boot_fn <- function(data, indices, formula) {
  train_data <- data[indices, ]  # Resample data from full dataset
  model <- glm(formula, family = quasipoisson, data = train_data)
  return(coef(model))}  # Return coefficients


# MODEL A: UNIVARIATE MODEL
# Initialize results dataframe
bootstrap_results <- data.frame(Outcome = character(),
                                       DietPattern = character(),
                                       Term = character(),
                                       Estimate = numeric(),
                                       RateRatio = numeric(),
                                       Conf.Low = numeric(),
                                       Conf.High = numeric(),
                                       pseudor2 = numeric(),
                                       pval = numeric(),
                                       stringsAsFactors = FALSE)

# Loop through each outcome and diet predictor
set.seed(11)
for (outcome in outcomes) {
  for (diet in diet_patterns) {
    # Define formula dynamically
    formula <- as.formula(paste(outcome, "~", diet))
    # Run bootstrap
    model <- boot(FGdietpatterns_adhdrs_covar, 
                        function(data, indices) boot_fn(data, indices, formula), 
                        R = 1000)
    # Fit model on full dataset to get p-value
    full_model <- glm(formula, family = quasipoisson, data = FGdietpatterns_adhdrs_covar)
    p_value <- summary(full_model)$coefficients[diet, "Pr(>|t|)"]
    glance <- glance(full_model)
    pseudo_r2 <- 1 - glance$deviance / glance$null.deviance
    # Get the coefficient names (Intercept + Diet predictor)
    term_names <- names(coef(full_model))
    # Extract bootstrap confidence intervals
    boot_ci <- t(sapply(1:length(coef(full_model)), function(i) {
      ci <- boot.ci(model, index = i, type = "perc")
      if (!is.null(ci)) {return(ci$percent[4:5])}  # Extract lower and upper CI
      else {return(c(NA, NA))}}))  # Handle cases where CI is not computed
    # Store results
    for (i in 1:length(model$t0)) {
      bootstrap_results <- rbind(bootstrap_results, data.frame(
        Outcome = outcome,
        DietPattern = diet,
        Term = term_names[i],  
        Estimate = model$t0[i],  
        RateRatio = exp(model$t0[i]),  
        Conf.Low = ifelse(!is.na(boot_ci[i, 1]), exp(boot_ci[i, 1]), NA),
        Conf.High = ifelse(!is.na(boot_ci[i, 2]), exp(boot_ci[i, 2]), NA),
        pseudor2 = pseudo_r2,
        pval = ifelse(i == 2, p_value, NA)))}}} # only save pval for the diet variable

bootstrap_results <- bootstrap_results %>% filter(Term != "(Intercept)") %>% group_by(DietPattern) %>% 
  mutate(adj.pval = p.adjust(pval, method = "fdr")) %>% ungroup()

# Forest plots
modelAplot <- ggplot(bootstrap_results, 
                     aes(x = RateRatio, y = Outcome, color = DietPattern)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +  # Plot estimates
  geom_errorbarh(aes(xmin = Conf.Low, xmax = Conf.High), height = 0.2, position = position_dodge(width = 0.3)) +  # CI bars
  geom_text(aes(label = case_when(
    adj.pval < 0.001 ~ "***",
    adj.pval < 0.01  ~ "**",
    adj.pval < 0.05  ~ "*",
    TRUE         ~ "")),
    nudge_y = 0.2) +  # Position asterisks slightly above data points
  scale_color_manual(values = c("deeppink2", "darkturquoise"),  # Custom colors
                     labels = c("Diet_1FG" = "Varied Diet", "Diet_2FG" = "Western Diet")) +  # Rename predictors in legend
  scale_y_discrete(labels = c("attention" = "Inattention", 
                              "implus_hyperact" = "Impulsivity/\nHyperactivity",
                              "total_score" = "Total ADHD-RS\nScore",
                              "behavior_subscale" = "Oppositional\nBehavior")) +  # Rename y-axis labels
  scale_x_continuous(limits = c(0.8, 1.5)) + 
  theme_minimal() 
  labs(color = "Diet Pattern") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 17),
        axis.text.y = element_text(size = 16),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 16)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black")  # Reference line at RR = 1

---------------------------

# MODEL B: UNIVARIATE MODEL WITHOUT DIAGNOSED ADHD PARTICIPANTS
FGdietpatterns_adhdrs_covar_noF90 <- FGdietpatterns_adhdrs_covar %>% filter(F900 != "F900") %>% filter(F901 != "F901")

# Initialize results dataframe
bootstrap_results_noF90<- data.frame(Outcome = character(),
                                DietPattern = character(),
                                Term = character(),
                                Estimate = numeric(),
                                RateRatio = numeric(),
                                Conf.Low = numeric(),
                                Conf.High = numeric(),
                                pseudor2 = numeric(),
                                pval = numeric(),
                                stringsAsFactors = FALSE)

# Loop through each outcome and diet predictor
set.seed(11)
for (outcome in outcomes) {
  for (diet in diet_patterns) {
    # Define formula dynamically
    formula <- as.formula(paste(outcome, "~", diet))
    # Run bootstrap
    model <- boot(FGdietpatterns_adhdrs_covar_noF90, 
                  function(data, indices) boot_fn(data, indices, formula), 
                  R = 1000)
    # Fit model on full dataset to get p-value
    full_model <- glm(formula, family = quasipoisson, data = FGdietpatterns_adhdrs_covar_noF90)
    p_value <- summary(full_model)$coefficients[diet, "Pr(>|t|)"]
    glance <- glance(full_model)
    pseudo_r2 <- 1 - glance$deviance / glance$null.deviance
    # Get the coefficient names (Intercept + Diet predictor)
    term_names <- names(coef(full_model))
    # Extract bootstrap confidence intervals
    boot_ci <- t(sapply(1:length(coef(full_model)), function(i) {
      ci <- boot.ci(model, index = i, type = "perc")
      if (!is.null(ci)) {return(ci$percent[4:5])}  # Extract lower and upper CI
      else {return(c(NA, NA))}}))  # Handle cases where CI is not computed
    # Store results
    for (i in 1:length(model$t0)) {
      bootstrap_results_noF90 <- rbind(bootstrap_results_noF90, data.frame(
        Outcome = outcome,
        DietPattern = diet,
        Term = term_names[i],  
        Estimate = model$t0[i],  
        RateRatio = exp(model$t0[i]),  
        Conf.Low = ifelse(!is.na(boot_ci[i, 1]), exp(boot_ci[i, 1]), NA),
        Conf.High = ifelse(!is.na(boot_ci[i, 2]), exp(boot_ci[i, 2]), NA),
        pseudor2 = pseudo_r2,
        pval = ifelse(i == 2, p_value, NA)))}}} # only save pval for the diet variable

bootstrap_results_noF90 <- bootstrap_results_noF90 %>% filter(Term != "(Intercept)") %>% group_by(DietPattern) %>% 
  mutate(adj.pval = p.adjust(pval, method = "fdr")) %>% ungroup()

# Forest plots
modelBplot <- ggplot(bootstrap_results_noF90, 
                     aes(x = RateRatio, y = Outcome, color = DietPattern)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +  # Plot estimates
  geom_errorbarh(aes(xmin = Conf.Low, xmax = Conf.High), height = 0.2, position = position_dodge(width = 0.3)) +  # CI bars
  geom_text(aes(label = case_when(
    adj.pval < 0.001 ~ "***",
    adj.pval < 0.01  ~ "**",
    adj.pval < 0.05  ~ "*",
    TRUE         ~ "")),
    nudge_y = 0.2) +  # Position asterisks slightly above data points
  scale_color_manual(values = c("deeppink2", "darkturquoise"),  # Custom colors
                     labels = c("Diet_1FG" = "Varied Diet", "Diet_2FG" = "Western Diet")) +  # Rename predictors in legend
  scale_y_discrete(labels = c("attention" = "Inattention", 
                              "implus_hyperact" = "Impulsivity/\nHyperactivity",
                              "total_score" = "Total ADHD-RS\nScore",
                              "behavior_subscale" = "Oppositional\nBehavior")) +  # Rename y-axis labels
  scale_x_continuous(limits = c(0.8, 1.5)) + 
  theme_minimal() +
  labs(color = "Diet Pattern") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 17),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 16)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black")  # Reference line at RR = 1


--------------------------------------

## MODEL C: MULTIVARIATE MODEL
# Initialize results dataframe
bootstrap_results_multi3 <- data.frame(Outcome = character(),
                                       DietPattern = character(),
                                       Term = character(),
                                       Estimate = numeric(),
                                       RateRatio = numeric(),
                                       Conf.Low = numeric(),
                                       Conf.High = numeric(),
                                       pseudor2 = numeric(),
                                       pval = numeric(),
                                       stringsAsFactors = FALSE)

# create bootstrapping function 
boot_fn_multi <- function(data, indices, formula) {
  train_data <- data[indices, ]  # Resample data from full dataset
  model_multi <- glm(formula_multi, family = quasipoisson, data = train_data)
  return(coef(model_multi))}  # Return coefficients

# Loop through each outcome and diet predictor
set.seed(11)
for (outcome in outcomes) {
  for (diet in diet_patterns) {
    # Define formula dynamically
    formula_multi <- as.formula(paste(outcome, "~", diet, "+ sex + socialcircumstances + energi + ageatADHDRS"))
    # Run bootstrap
    model_multi <- boot(FGdietpatterns_adhdrs_covar, 
                        function(data, indices) boot_fn_multi(data, indices, formula_multi), 
                        R = 1000)
    # Fit model on full dataset to get p-value
    full_model <- glm(formula_multi, family = quasipoisson, data = FGdietpatterns_adhdrs_covar)
    p_value <- summary(full_model)$coefficients[diet, "Pr(>|t|)"]
    glance <- glance(full_model)
    pseudo_r2 <- 1 - glance$deviance / glance$null.deviance
    # Get the coefficient names (Intercept + Diet predictor)
    term_names <- names(coef(full_model))
    # Extract bootstrap confidence intervals
    boot_ci <- t(sapply(1:length(coef(full_model)), function(i) {
      ci <- boot.ci(model_multi, index = i, type = "perc")
      if (!is.null(ci)) {return(ci$percent[4:5])}  # Extract lower and upper CI
      else {return(c(NA, NA))}}))  # Handle cases where CI is not computed
    # Store results
    for (i in 1:length(model_multi$t0)) {
      bootstrap_results_multi3 <- rbind(bootstrap_results_multi3, data.frame(
        Outcome = outcome,
        DietPattern = diet,
        Term = term_names[i],  
        Estimate = model_multi$t0[i],  
        RateRatio = exp(model_multi$t0[i]),  
        Conf.Low = ifelse(!is.na(boot_ci[i, 1]), exp(boot_ci[i, 1]), NA),
        Conf.High = ifelse(!is.na(boot_ci[i, 2]), exp(boot_ci[i, 2]), NA),
        pseudor2 = pseudo_r2,
        pval = ifelse(i == 2, p_value, NA)))}}} # only save pval for the diet variable

bootstrap_results_multi3 <- bootstrap_results_multi3 %>% filter(Term %in% c("Diet_1FG", "Diet_2FG")) %>% group_by(DietPattern) %>% 
  mutate(adj.pval = p.adjust(pval, method = "fdr")) %>% ungroup()

# forest plot multivariate model
modelCplot1 <- ggplot(bootstrap_results_multi3, 
                     aes(x = RateRatio, y = Outcome, color = DietPattern)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +  # Plot estimates
  geom_errorbarh(aes(xmin = Conf.Low, xmax = Conf.High), height = 0.2, position = position_dodge(width = 0.3)) +  # CI bars
  geom_text(aes(label = case_when(
    adj.pval < 0.001 ~ "***",
    adj.pval < 0.01  ~ "**",
    adj.pval < 0.05  ~ "*",
    TRUE         ~ "")),
    nudge_y = 0.2) +  # Position asterisks slightly above data points
  scale_color_manual(values = c("deeppink2", "darkturquoise"),  # Custom colors
                     labels = c("Diet_1FG" = "Varied Diet", "Diet_2FG" = "Western Diet")) +  # Rename predictors in legend
  scale_color_manual(values = c("deeppink2", "darkturquoise"),  # Custom colors
                     labels = c("Diet_1FG" = "Varied Diet", "Diet_2FG" = "Western Diet")) +  # Rename predictors in legend
  scale_y_discrete(labels = c("attention" = "Inattention", 
                              "implus_hyperact" = "Impulsivity/\nHyperactivity",
                              "total_score" = "Total ADHD-RS\nScore",
                              "behavior_subscale" = "Oppositional\nBehavior")) +  # Rename y-axis labels
  scale_x_continuous(limits = c(0.8, 1.45)) + 
  theme_minimal() + xlim(0.8,1.5) +
  labs(color = "Diet Pattern") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 17),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 16)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black")  # Reference line at RR = 1

---------------------
  
# create side-by-side figure
combinedplots <- ggarrange(modelAplot,modelBplot, modelCplot1, labels = c("A", "B", "C"),
                           ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom", widths = c(1.3, 1, 1))
combinedplots <- annotate_figure(combinedplots, bottom = text_grob("Rate Ratio (95% CIs)"))
print(combinedplots)

------------------------

## ROBUSTNESS ANALYSIS 1: SPEARMAN CORRELATIONS
FGspeardietpatterns_adhdrs <- FGdietpatterns_adhdrs %>%
  gather(diet, yy, `Diet_1FG`:`Diet_2FG`) %>%  # endpoints > the three diets (three main PCs fromPCA of nutrients)
  gather(adhdrs, y, `attention`:`total_score`) %>%  # predictors > various variables
  group_by(diet, adhdrs) %>%
  summarise(
    cor_result = list(tidy(cor.test(yy, y, method = "spearman", exact = FALSE))),
    .groups = "drop"
  ) %>%
  unnest(cor_result) %>%
  select(diet, adhdrs, estimate, p.value)
FGspeardietpatterns_adhdrs <- FGspeardietpatterns_adhdrs %>% # Adjust p-values for multiple comparisons
  mutate(p_adjusted = p.adjust(p.value, method = "fdr"))

FGspeardietpatterns_adhdrs$estimate <- as.numeric(FGspeardietpatterns_adhdrs$estimate)

# lollipop graph
# Create numeric x and dodged version of data for plotting 
FGspeardietpatterns_adhdrs$x_numeric <- as.numeric(factor(FGspeardietpatterns_adhdrs$adhdrs,
                                                          levels = c("total_score", "implus_hyperact", 
                                                                     "behavior_subscale", "attention")))
FGspeardietpatterns_adhdrs$x_dodged <- with(FGspeardietpatterns_adhdrs,
                                            ifelse(diet == "Diet_1FG",
                                                   x_numeric - 0.05,
                                                   x_numeric + 0.05))
# Plot lollipop graph
ggplot(FGspeardietpatterns_adhdrs, aes(x = x_dodged, y = estimate, color = diet)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_segment(aes(x = x_dodged, xend = x_dodged, y = 0, yend = estimate),
               color = "grey60", show.legend = FALSE) +
  geom_point(size = 5) +
  geom_text(aes(label = case_when(
    p_adjusted < 0.001 ~ "***",
    p_adjusted < 0.01  ~ "**",
    p_adjusted < 0.05  ~ "*",
    TRUE         ~ "")),
    vjust = -1, color = "black") +
  scale_color_manual(values = c("Diet_1FG" = "deeppink2", "Diet_2FG" = "darkturquoise"),
                     labels = c("Diet_1FG" = "Varied Diet", "Diet_2FG" = "Western Diet"),
                     name = "Diet Pattern") +
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("Total ADHD-RS\nScore", 
                                                        "Impulsivity/\nHyperactivity", 
                                                        "Oppositional\nBehavior", 
                                                        "Inattention")) +
  scale_y_continuous(limits = c(-0.25, 0.25)) +
  labs(x = "ADHD Measure",
       y = "Spearman Correlation Estimate") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))

-------------------------------------
  
# ROBUSTNESS ANALYSIS 2: TRAIN-TEST SPLIT OF Q-P MODELS

# create training and testing dfs
set.seed(00)
inTraining <- createDataPartition(FGdietpatterns_adhdrs_covar$F900, p = .7, list = FALSE)
training_data <- FGdietpatterns_adhdrs_covar[ inTraining,]
test_data  <- FGdietpatterns_adhdrs_covar[-inTraining,]
rm(inTraining)

## UNIVARIATE MODEL
# Initialize results dataframe
results <- data.frame(Outcome = character(),
                      DietPattern = character(),
                      Estimate = numeric(),
                      RateRatio = numeric(),
                      DevianceReduction = numeric(),
                      PseudoR2 = numeric(),
                      RMSE = numeric(),
                      pval = numeric(),
                      stringsAsFactors = FALSE)

# Loop through each outcome and diet predictor
for (outcome in outcomes) {
  for (diet in diet_patterns) {
    # Fit the model
    formula <- as.formula(paste(outcome, "~", diet))
    model <- glm(formula, family = quasipoisson, data = training_data)
    # Get model performance metrics
    perf <- broom::glance(model)
    deviance_reduction <- perf$null.deviance - perf$deviance
    pseudo_r2 <- 1 - perf$deviance / perf$null.deviance
    # Extract p-value for the diet predictor
    model_summary <- summary(model)
    p_value <- model_summary$coefficients[diet, "Pr(>|t|)"]
    # Extract estimate for the diet predictor (excluding intercept)
    estimate <- model_summary$coefficients[diet, "Estimate"]
    # Transform estimate to RateRatio (exp(estimate))
    rate_ratio <- exp(estimate)
    # Predict on test data
    predicted <- predict(model, type = "response", newdata = test_data)
    # Compute RMSE
    actual <- test_data[[outcome]]
    rmse <- sqrt(mean((actual - predicted)^2, na.rm = TRUE))
    # Store results
    results <- rbind(results, data.frame(Outcome = outcome,
                                         DietPattern = diet,
                                         Estimate = estimate,
                                         RateRatio = rate_ratio,
                                         DevianceReduction = deviance_reduction,
                                         PseudoR2 = pseudo_r2,
                                         RMSE = rmse,
                                         pval = p_value,
                                         stringsAsFactors = FALSE))}}

results$RRgraph <- results$RateRatio -1

## MULTIVARIATE MODEL
# adjust for sex, socialcircumstances, and energy intake
# Initialize results dataframe
results_multivariate <- data.frame(Outcome = character(),
                                   DietPattern = character(),
                                   Estimate = numeric(),
                                   RateRatio = numeric(),
                                   DevianceReduction = numeric(),
                                   PseudoR2 = numeric(),
                                   RMSE = numeric(),
                                   pval = numeric(),
                                   stringsAsFactors = FALSE)

# Loop through each outcome and diet predictor
set.seed(11)
for (outcome in outcomes) {
  for (diet in diet_patterns) {
    # Fit the model
    formula <- as.formula(paste(outcome, "~", diet, "+ sex + socialcircumstances + energi + ageatADHDRS"))
    model <- glm(formula, family = quasipoisson, data = training_data)
    # Get model performance metrics
    perf <- broom::glance(model)
    deviance_reduction <- perf$null.deviance - perf$deviance
    pseudo_r2 <- 1 - perf$deviance / perf$null.deviance
    # Extract p-value for the diet predictor
    model_summary <- summary(model)
    p_value <- model_summary$coefficients[diet, "Pr(>|t|)"]
    # Extract estimate for the diet predictor (excluding intercept)
    estimate <- model_summary$coefficients[diet, "Estimate"]
    # Transform estimate to RateRatio (exp(estimate))
    rate_ratio <- exp(estimate)
    # Predict on test data
    predicted <- predict(model, type = "response", newdata = test_data)
    # Compute RMSE
    actual <- test_data[[outcome]]
    rmse <- sqrt(mean((actual - predicted)^2, na.rm = TRUE))
    # Store results
    results_multivariate <- rbind(results_multivariate, data.frame(Outcome = outcome,
                                                                   DietPattern = diet,
                                                                   Estimate = estimate,
                                                                   RateRatio = rate_ratio,
                                                                   DevianceReduction = deviance_reduction,
                                                                   PseudoR2 = pseudo_r2,
                                                                   RMSE = rmse,
                                                                   pval = p_value,
                                                                   stringsAsFactors = FALSE))}}

results_multivariate$RRgraph <- results_multivariate$RateRatio -1

