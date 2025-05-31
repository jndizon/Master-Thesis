## ANALYSIS ON THE (FG-BASED) DIET PATTERNS

# covariates, potential confounders, ADHDRS, and picky eating analysis

# full dataset: n = 597
# microbiome data: n = 309 of the 597
# picky eating data: n = 494 of the 598
# ADHD-RS data: n = 575 of the 598

library(ggplot2) # v3.5.1
library(dplyr) # v1.1.4
library(ggpubr) # v0.6.0
library(broom) # v1.0.7

--------------------

# histograms to visualize data distribution
hist(food_pc$Diet_1FG, 15, xlim = c(-4, 4))
hist(food_pc$Diet_2FG, 15, xlim = c(-10, 5))

----------------------------------
  
## spearman correlations between FG-based dietary patterns and other factors
# requires dfs from "FG level PCs" and many others

# create input df
FGdietpatterns_analysis <- food_pc # add the diet patterns
FGdietpatterns_analysis <- merge(FGdietpatterns_analysis, variablesmerged, by = "abcno") %>% 
  select(-bmi, -abmotheratbirth, -Overweight, -F900, -F901) # add variables
FGdietpatterns_analysis <- merge(FGdietpatterns_analysis, picky_13y, by = "abcno", all.x = TRUE) %>% 
  select(-s_331:-s_336) # add picky eating score
FGdietpatterns_analysis <- FGdietpatterns_analysis %>%
  mutate(sex = if_else(sex == "male", 1, 0)) %>% # create numeric variables for sex and race 
  mutate(race = if_else(race == "caucasian", 1, 0)) # male + causaian = 1 // female + non-causaian = 0
FGdietpatterns_analysis$sex <- as.numeric(as.character(FGdietpatterns_analysis$sex))
FGdietpatterns_analysis$race <- as.numeric(as.character(FGdietpatterns_analysis$race))

# run spearman cor
FGdietpatterns_cor <- FGdietpatterns_analysis %>%
  gather(diet, yy, `Diet_1FG`:`Diet_2FG`) %>%  # endpoints > the diets
  gather(vars, y, `sex`:`picky_score`) %>%  # predictors > various variables
  group_by(diet, vars) %>%
  summarise(
    cor_result = list(tidy(cor.test(yy, y, method = "spearman", exact = FALSE))),
    .groups = "drop") %>%
  unnest(cor_result) %>%
  select(diet, vars, estimate, p.value)
FGdietpatterns_cor <- FGdietpatterns_cor %>% # Adjust p-values for multiple comparisons
  mutate(p_adjusted = p.adjust(p.value, method = "fdr"))

# manipulating the df for plotting
# reorder 
FGdietpatterns_cor$exposure_1 <- factor(FGdietpatterns_cor$vars, 
                                        levels=rev(c("picky_score", "energi", "zbmi", "socialcircumstances", "race", "sex")))
FGdietpatterns_cor$exposure_1.1   <- order(FGdietpatterns_cor$exposure_1, FGdietpatterns_cor$estimate, decreasing=TRUE)
FGdietpatterns_cor <- FGdietpatterns_cor[order(FGdietpatterns_cor$estimate),]
# relevel so that Western is plotted on the bottom
FGdietpatterns_cor$diet <- factor(FGdietpatterns_cor$diet, levels = c("Diet_2FG", "Diet_1FG"))
FGdietpatterns_cor <- FGdietpatterns_cor %>% filter(vars != "picky_score") 

# create heatmap
ggplot(FGdietpatterns_cor, aes(x = exposure_1, y = diet, fill = estimate)) +
  geom_tile(color = "white") +  # Creates heatmap tiles
  scale_fill_gradient2(low = "darkgoldenrod2", mid = "white", high = "darkorchid1", midpoint = 0, 
                       name = "Spearman\nCorrelation") +
  geom_text(aes(label = case_when(
    p_adjusted < 0.001 ~ "***",
    p_adjusted < 0.01  ~ "**",
    p_adjusted < 0.05  ~ "*",
    TRUE         ~ "")), color = "black", size = 5) +
  labs(title = "Diet Patterns via Spearman Correlation",
       x = "Variables",
       y = "Diet patterns") +
  scale_x_discrete(labels = c("energi" = "Total Energy Intake", 
                              "zbmi" = "BMI", 
                              "socialcircumstances" = "Social Circumstances", 
                              "race" = "Race", 
                              "sex" = "Sex")) +  # Renaming x-axis labels
  scale_y_discrete(labels = c("Diet_1FG" = "Varied Diet Pattern", 
                              "Diet_2FG" = "Western Diet Pattern")) +  # Renaming y-axis labels (adjust accordingly)
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

-------------------------------------
  
# COVARIATE ANALYSIS with diet patterns
# Linear model for the covariates
# use same input df (FGdietpatterns_analyss)

# diet_1FG:
food_pc_covars_fullmodel <- lm(Diet_1FG ~ sex + race + socialcircumstances + zbmi + energi, data = FGdietpatterns_analysis)
anova(food_pc_covars_fullmodel) 
nobs(food_pc_covars_fullmodel)
# [1] 589 -> n = 589 of 597 were included in the analysis
food_pc_covars_fullmodel1 <- lm(Diet_1FG ~ race + socialcircumstances + zbmi + energi, data = FGdietpatterns_analysis)
anova(food_pc_covars_fullmodel1) 

food_pc_covars_fullmodel2 <- lm(Diet_1FG ~ socialcircumstances + zbmi + energi, data = FGdietpatterns_analysis)
anova(food_pc_covars_fullmodel2) 
# social circumstances, bmi z-score, and total energy intake are significant

# diet_2FG
food_pc_covars_fullmodel <- lm(Diet_2FG ~ sex + race + socialcircumstances + zbmi + energi, data = FGdietpatterns_analysis)
anova(food_pc_covars_fullmodel) 
nobs(food_pc_covars_fullmodel)
# [1] 589 -> n = 589 of 597 were included in the analysis
food_pc_covars_fullmodel1 <- lm(Diet_2FG ~ race + socialcircumstances + zbmi + energi, data = FGdietpatterns_analysis)
anova(food_pc_covars_fullmodel1) 

food_pc_covars_fullmodel2 <- lm(Diet_2FG ~ socialcircumstances + race + energi, data = FGdietpatterns_analysis)
anova(food_pc_covars_fullmodel2) 
# social circumstances, race, and total energy intake are significant

rm(list = ls(, pattern = "^food_pc_covars"))

---------------------------------

## nutrient intakes represented by the diet patterns
# create input df with participants' diet pattern component and z-scores of nutrient intakes
FGdietpatterns_nutr <- merge(food_pc, nutr_zscores, by = "abcno")

spearFGdiet_nutr <- FGdietpatterns_nutr %>%
  gather(diet, yy, `Diet_1FG`:`Diet_2FG`) %>%  # endpoints > the diets
  gather(nutr, y, `suktil`:`serin`) %>%  # predictors > nutrient intakes
  group_by(diet, nutr) %>%
  summarise(
    cor_result = list(tidy(cor.test(yy, y, method = "spearman", exact = FALSE))),
    .groups = "drop") %>%
  unnest(cor_result) %>%
  select(diet, nutr, estimate, p.value)
spearFGdiet_nutr <- spearFGdiet_nutr %>% # Adjust p-values for multiple comparisons
  mutate(p_adjusted = p.adjust(p.value, method = "fdr"))

# df manips for plotting purposes
# new labels 
nutrient_labels <- c(
  "nikkel"   = "Nickel",
  "selen"    = "Selenium",
  "chrom"    = "Chromium",
  "mangan"   = "Manganese",
  "jod"      = "Iodine",
  "zink"     = "Zinc",
  "kobber"   = "Copper",
  "jern"     = "Iron",
  "fosfor"   = "Phosphor",
  "magnes"   = "Magnesium",
  "calciu"   = "Calcium",
  "kaliu"    = "Potassium",
  "natriu"   = "Sodium",
  "vand"     = "Water",
  "aske"     = "Ashes",
  "ldehydr"  = "L-dehydroascorbic acid",
  "lascorb"  = "L-ascorbic acid",
  "cvit"     = "Vitamin C",
  "b12vit"   = "Vitamin B12",
  "folaci"   = "Folate",
  "biotin"   = "Biotin",
  "pantot"   = "Pantothenic acid",
  "b6vit"    = "Vitamin B6",
  "niacin"   = "Niacin",
  "neniac"   = "Niacin equivalents",
  "ribofl"   = "Vitamin B2",
  "thiami"   = "Vitamin B1",
  "kvit"     = "Vitamin K",
  "alfato"   = "Alphatocoferol",
  "evit"     = "Vitamin E",
  "v5hydr"   = "25-hydroxycholecalciferol",
  "d3cho"    = "Vitamin D3",
  "dvit"     = "Vitamin D",
  "betaca"   = "Betacarotene",
  "retino"   = "Retinol",
  "avit"     = "Vitamin A",
  "serin"    = "Serine",
  "prolin"   = "Proline",
  "glycin"   = "Glycine",
  "glutam"   = "Glutamic acid",
  "aspara"   = "Asparagine acid",
  "alanin"   = "Alanine",
  "histid"   = "Histidine",
  "argini"   = "Arginine",
  "valin"    = "Valine",
  "trypto"   = "Tryptophan",
  "threon"   = "Threonine",
  "thyros"   = "Tyrosine",
  "phenyl"   = "Phenylalanine",
  "cystin"   = "Cystine",
  "methio"   = "Methionine",
  "lysin"    = "Lysine",
  "leucin"   = "Leucine",
  "isoleu"   = "Isoleucine",
  "ntrypto"  = "N-Tryptophan",
  "choles"   = "Cholesterol",
  "transf"   = "Trans fatty acids, total",
  "andrfs"   = "Other fatty acids",
  "c22x6n3"  = "Docosahexaenoic acid DHA",
  "c22x5n3"  = "Docosapentaenoic acid DPA",
  "c20x5n3"  = "Eicosapentaenoic acid EPA",
  "c20x4n6"  = "Arachidonic acid AA",
  "c18x4n3"  = "Stearidonic acid SDA",
  "c18x3n3"  = "Alpha-linolenic acid ALA",
  "c18x2n6"  = "Linoleic acid",
  "c24x1n9"  = "Nervonic acid",
  "c22x1n11" = "Cetoleic acid",
  "c22x1n9"  = "Erucic acid",
  "c20x1n11" = "Gadoleic acid",
  "c18x1n7"  = "Vaccenic acid",
  "c18x1n9"  = "Oleic acid",
  "c16x1"    = "Palmitoleic acid",
  "c14x1"    = "Myristoleic acid",
  "c24x0"    = "Lignoceric acid",
  "c22x0"    = "Behenic acid",
  "c20x0"    = "Arachidic acid",
  "c18x0"    = "Stearic acid",
  "c17x0"    = "Margaric acid",
  "c16x0"    = "Palmitic acid",
  "c15x0"    = "Pentadecylic acid",
  "c14x0"    = "Myristic acid",
  "c12x0"    = "Lauric acid",
  "c10x0"    = "Capric acid",
  "c8x0"     = "Caprylic acid",
  "c6x0"     = "Caproic acid",
  "c4x0"     = "Butyric acid",
  "stivel"   = "Starch",
  "saccha"   = "Sucrose",
  "maltos"   = "Maltose",
  "laktos"   = "Lactose",
  "glukos"   = "Glucose",
  "frukto"   = "Fructose",
  "alko"     = "Alcohol",
  "kfibre"   = "Dietary fiber",
  "suktil"   = "Added sugar")
spearFGdiet_nutr <- spearFGdiet_nutr %>%
  mutate(label = recode(nutr, !!!nutrient_labels))

# reorder the nutr variable so each group of nutrients will be plotted tg
spearFGdiet_nutr$exposure_1 <- factor(spearFGdiet_nutr$label, 
                                      levels = rev(c(
                                        "Nickel", "Selenium", "Chromium", "Manganese", "Iodine", "Zinc", "Copper", "Iron", "Phosphor", "Magnesium", 
                                        "Calcium", "Potassium", "Sodium", "Water", "Ashes", "L-dehydroascorbic acid", "L-ascorbic acid", 
                                        "Vitamin C", "Vitamin B12", "Folate", "Biotin", "Pantothenic acid", "Vitamin B6", "Niacin", 
                                        "Niacin equivalents", "Vitamin B2", "Vitamin B1", "Vitamin K", "Alphatocoferol", 
                                        "Vitamin E", "25-hydroxycholecalciferol", "Vitamin D3", "Vitamin D", "Betacarotene", 
                                        "Retinol", "Vitamin A", "Serine", "Proline", "Glycine", "Glutamic acid", "Asparagine acid", "Alanine", 
                                        "Histidine", "Arginine", "Valine", "Tryptophan", "Threonine", "Tyrosine", "Phenylalanine", "Cystine", 
                                        "Methionine", "Lysine", "Leucine", "Isoleucine", "N-Tryptophan", "Cholesterol", "Trans fatty acids, total", 
                                        "Other fatty acids", "Docosahexaenoic acid DHA", "Docosapentaenoic acid DPA", "Eicosapentaenoic acid EPA", 
                                        "Arachidonic acid AA", "Stearidonic acid SDA", "Alpha-linolenic acid ALA", "Linoleic acid", 
                                        "Nervonic acid", "Cetoleic acid", "Erucic acid", "Gadoleic acid", "Vaccenic acid", "Oleic acid", 
                                        "Palmitoleic acid", "Myristoleic acid", "Lignoceric acid", "Behenic acid", "Arachidic acid", 
                                        "Stearic acid", "Margaric acid", "Palmitic acid", "Pentadecylic acid", "Myristic acid", "Lauric acid", 
                                        "Capric acid", "Caprylic acid", "Caproic acid", "Butyric acid", "Starch", "Sucrose", "Maltose", 
                                        "Lactose", "Glucose", "Fructose", "Alcohol", "Dietary fiber", "Added sugar")))
spearFGdiet_nutr$exposure_1.1   <- order(spearFGdiet_nutr$exposure_1, spearFGdiet_nutr$estimate, decreasing=TRUE)
spearFGdiet_nutr <- spearFGdiet_nutr[order(spearFGdiet_nutr$estimate),]

# relevel so that Western is plotted on the bottom
spearFGdiet_nutr$diet <- factor(spearFGdiet_nutr$diet, levels = c("Diet_2FG", "Diet_1FG"))

# heatmap of correlations of nutr z-scores and diet patterns
ggplot(spearFGdiet_nutr, aes(x = exposure_1, y = diet, fill = estimate)) +
  geom_tile(color = "white") +  # Creates heatmap tiles
  scale_fill_gradient2(low = "darkgoldenrod2", mid = "white", high = "darkorchid1", midpoint = 0, 
                       name = "Spearman\nCorrelation") +
  labs(title = "Spearman Correlation Heatmap (Nutrient intakes represented in the diet patterns)",
       x = "Nutrient",
       y = "Dietary pattern") +
  scale_y_discrete(labels = c("Diet_1FG" = "Varied Diet Pattern", 
                              "Diet_2FG" = "Western Diet Pattern")) +  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

---------------------------

# characterizing group-level nutrient intake with the diet patterns
# create input df
FGdietpatterns_nutr_macros<- merge(food_pc, nutr_macros_zscores, by = "abcno")

spearFGdiet_nutr_macros <- FGdietpatterns_nutr_macros %>%
  gather(diet, yy, `Diet_1FG`:`Diet_2FG`) %>%  # endpoints > the diets
  gather(nutr, y, `prot`:`kfibre`) %>%  # predictors > nutrient intakes
  group_by(diet, nutr) %>%
  summarise(
    cor_result = list(tidy(cor.test(yy, y, method = "spearman", exact = FALSE))),
    .groups = "drop") %>%
  unnest(cor_result) %>%
  select(diet, nutr, estimate, p.value)
spearFGdiet_nutr_macros <- spearFGdiet_nutr_macros %>% # Adjust p-values for multiple comparisons
  mutate(p_adjusted = p.adjust(p.value, method = "fdr"))

# df manips for plotting purposes
nutr_macros_labels <- c(
  "prot"    = "Protein",
  "fedt"    = "Fat",
  "sfa"   = "SFAs",
  "mufa"      = "MUFAs",
  "pufa"     = "PUFAs",
  "kulhy"   = "Carbohydrates",
  "suktil"     = "Added sugar",
  "kfibre" = "Fiber",
  "sumn3" = "Omega-3 FAs")
spearFGdiet_nutr_macros <- spearFGdiet_nutr_macros %>%
  mutate(label = recode(nutr, !!!nutr_macros_labels))

# reorder the nutr variable so each group of nutrients will be plotted tg
spearFGdiet_nutr_macros$exposure_1 <- factor(spearFGdiet_nutr_macros$label, 
                                      levels = rev(c("Omega-3 FAs", "PUFAs","SFAs","MUFAs", "Fat",
                                                     "Protein", "Fiber", "Added sugar", "Carbohydrates")))
spearFGdiet_nutr_macros$exposure_1.1   <- order(spearFGdiet_nutr_macros$exposure_1, spearFGdiet_nutr_macros$estimate, decreasing=TRUE)
spearFGdiet_nutr_macros <- spearFGdiet_nutr_macros[order(spearFGdiet_nutr_macros$estimate), ]

# relevel so that Western is plotted on the bottom
spearFGdiet_nutr_macros$diet <- factor(spearFGdiet_nutr_macros$diet, levels = c("Diet_2FG", "Diet_1FG"))

# heatmap of correlations macros and diet patterns
macrosheatmap <- ggplot(spearFGdiet_nutr_macros, aes(x = exposure_1, y = diet, fill = estimate)) +
  geom_tile(color = "white") +  # Creates heatmap tiles
  scale_fill_gradient2(low = "darkgoldenrod2", mid = "white", high = "darkorchid1", midpoint = 0, 
                       name = "Spearman\nCorrelation") +
  labs(x = "Nutrient Group",
       y = "Dietary pattern") +
  scale_y_discrete(labels = c("Diet_1FG" = "Varied Diet Pattern", 
                              "Diet_2FG" = "Western Diet Pattern")) +  # Renaming y-axis labels (adjust accordingly)
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.title.x = element_text(size = 19), 
        axis.title.y = element_text(size = 19),
        axis.text.y = element_text(size = 15),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 14))


# label the heatmap to add it to the combined figure with the factor loadings
ggarrange( macrosheatmap, labels = c("C"),
          ncol = 1, nrow = 1)

