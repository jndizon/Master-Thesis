# creating diet patterns via PCA of the food group intake data
  # code inspiration from David Horner

library(ggplot2) # v3.5.1
library(dplyr) # v1.1.4
library(ggpubr) # v0.6.0
library(ggbiplot) # v0.6.2

## RUNNING PCA
mdl2 <- prcomp(food10y2.0 %>% select(-1:-3), scale. = T) # Scales the variables 
summary(mdl2)
# have 95 PCs but the first 4 cover ~30% of variance in samples
#                           PC1     PC2     PC3     PC4     
# Standard deviation     1.9614 1.49215 1.34647 1.25411
# Proportion of Variance 0.1202 0.06958 0.05666 0.04915  
# Cumulative Proportion  0.1202 0.18979 0.24645 0.29560

# visualizations
ggbiplot::ggbiplot(mdl2)
screeplot(mdl2, npcs = 24, type = "lines") 

------------
# CREATING DF OF THE LOADINGS FOR FUTURE ANALYSIS
food_pc <- food10y2.0 %>% left_join(data.frame(abcno = food10y2.0$abcno, Newervalue = -mdl2$x[,1])) # inverts the sign of the loading
food_pc <- food_pc %>% dplyr::rename(Diet_1FG = Newervalue)

food_pc <- food_pc %>% left_join(data.frame(abcno = food10y2.0$abcno, Newervalue = -mdl2$x[,2]))
food_pc <- food_pc %>% dplyr::rename(Diet_2FG = Newervalue)

food_pc <- food_pc %>% select(abcno, Diet_1FG, Diet_2FG)
food_pc <- food_pc %>% mutate(Diet_1FG = scale(Diet_1FG), Diet_2FG = scale(Diet_2FG))
  # scale the components 

------------------------------

## MAKING FIGURES OF THE PC LOADINGS

# diet_1 fig
# Extract PC1 loadings and create a data frame
pc1_loadings <- data.frame(
  Variable = rownames(mdl2$rotation),  
  PC1_Loading = mdl2$rotation[, 1])  %>% 
  mutate(Variable = case_when(
    rownames(.) == "breakfast_cereals" ~ "Breakfast Cereals",
    rownames(.) == "water" ~ "Water",
    rownames(.) == "hig_fat_milk" ~ "High Fat Milk",
    rownames(.) == "low_fat_milk" ~ "Low Fat Milk",
    rownames(.) == "beans" ~ "Beans",
    rownames(.) == "other_vegetables" ~ "Other Vegetables",
    rownames(.) == "tomatoes" ~ "Tomatoes",
    rownames(.) == "soy_products" ~ "Soy Products",
    rownames(.) == "green_leafy_vegetables" ~ "Leafy Greens",
    rownames(.) == "fish" ~ "Fish",
    rownames(.) == "animal_fat" ~ "Animal Fat",
    rownames(.) == "cheese" ~ "Cheese",
    rownames(.) == "fruits" ~ "Fruits",
    rownames(.) == "whole_grains" ~ "Whole Grains",
    rownames(.) == "vegetable_fat" ~ "Vegetable Fat",
    rownames(.) == "tea" ~ "Tea",
    rownames(.) == "eggs" ~ "Eggs",
    rownames(.) == "dried_fruits" ~ "Dried Fruits",
    rownames(.) == "nuts" ~ "Nuts",
    rownames(.) == "poultry" ~ "Poultry",
    rownames(.) == "processed_meat" ~ "Processed Meat",
    rownames(.) == "refined_grains" ~ "Refined Grains",
    rownames(.) == "fruit_juice" ~ "Fruit Juice",
    rownames(.) == "sweets_desserts" ~ "Sweets & Desserts",
    rownames(.) == "fruit_syrup_marm" ~ "Fruit Syrup & Marmalade",
    rownames(.) == "red_meat" ~ "Red Meat",
    rownames(.) == "ice_cream" ~ "Ice Cream",
    rownames(.) == "margarine" ~ "Margarine",
    rownames(.) == "high_energi_drinks" ~ "High Energy Drinks",
    rownames(.) == "potatoes_and_products" ~ "Potatoes & Products",
    rownames(.) == "low_energi_drinks" ~ "Low Energy Drinks",
    rownames(.) == "snacks" ~ "Snacks"))

pc1_loadings <- pc1_loadings %>% mutate(PC1_Loading=(PC1_Loading*-1)) # sign inversion to depict varied diet

# Assign fill color based on loading value
pc1_loadings <- pc1_loadings %>%
  mutate(
    FILL = case_when(
      PC1_Loading >= 0.09 ~ "High",          # High positive loadings
      PC1_Loading < 0.09 & PC1_Loading > -0.09 ~ "Mid",  # Near zero
      PC1_Loading <= -0.09 ~ "Low"),           # High negative loadings
    FILL = as.factor(FILL))

# Order variables by loading values
pc1_loadings <- pc1_loadings[order(pc1_loadings$PC1_Loading), ]  
# Convert Variable to a factor for ordered plotting
pc1_loadings$Variable <- factor(pc1_loadings$Variable, levels = pc1_loadings$Variable)

# Create the bar chart
food1 <- ggplot(pc1_loadings %>% filter(FILL != "Mid"), aes(y = Variable, x = PC1_Loading, fill = FILL)) +
  geom_bar(stat = 'identity', position = 'dodge', color = 'black') +  # Bars
  xlim(c(-0.5, 0.5)) +
  scale_fill_manual(values = c("deeppink2", "darkgoldenrod2", "blueviolet")) + 
  guides(fill = "none") +
  xlab("PC1 Loadings") +
  ylab(NULL) +
  geom_vline(xintercept = 0, color = "black", size = 0.1) +
  ggtitle("Varied Diet Pattern (12%)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.text.x = element_text(size = 15),
    axis.title = element_text(size = 18),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) +
  geom_text(data = subset(pc1_loadings, FILL == "Low"),
            aes(label = Variable, x = 0.01),
            hjust = 0, size = 6) +
  geom_text(data = subset(pc1_loadings, FILL == "High"),
            aes(label = Variable, x = -0.01),
            hjust = 1, size = 6)
food1
rm(pc1_loadings)

# diet_2 fig
# Extract PC2 loadings and create a data frame
pc2_loadings <- data.frame(
  Variable = rownames(mdl2$rotation), 
  PC2_Loading = mdl2$rotation[, 2]) %>% 
  mutate(Variable = case_when(
    rownames(.) == "breakfast_cereals" ~ "Breakfast Cereals",
    rownames(.) == "water" ~ "Water",
    rownames(.) == "hig_fat_milk" ~ "High Fat Milk",
    rownames(.) == "low_fat_milk" ~ "Low Fat Milk",
    rownames(.) == "beans" ~ "Beans",
    rownames(.) == "other_vegetables" ~ "Other Vegetables",
    rownames(.) == "tomatoes" ~ "Tomatoes",
    rownames(.) == "soy_products" ~ "Soy Products",
    rownames(.) == "green_leafy_vegetables" ~ "Leafy Greens",
    rownames(.) == "fish" ~ "Fish",
    rownames(.) == "animal_fat" ~ "Animal Fat",
    rownames(.) == "cheese" ~ "Cheese",
    rownames(.) == "fruits" ~ "Fruits",
    rownames(.) == "whole_grains" ~ "Whole Grains",
    rownames(.) == "vegetable_fat" ~ "Vegetable Fat",
    rownames(.) == "tea" ~ "Tea",
    rownames(.) == "eggs" ~ "Eggs",
    rownames(.) == "dried_fruits" ~ "Dried Fruits",
    rownames(.) == "nuts" ~ "Nuts",
    rownames(.) == "poultry" ~ "Poultry",
    rownames(.) == "processed_meat" ~ "Processed Meat",
    rownames(.) == "refined_grains" ~ "Refined Grains",
    rownames(.) == "fruit_juice" ~ "Fruit Juice",
    rownames(.) == "sweets_desserts" ~ "Sweets & Desserts",
    rownames(.) == "fruit_syrup_marm" ~ "Fruit Syrup & Marmalade",
    rownames(.) == "red_meat" ~ "Red Meat",
    rownames(.) == "ice_cream" ~ "Ice Cream",
    rownames(.) == "margarine" ~ "Margarine",
    rownames(.) == "high_energi_drinks" ~ "High Energy Drinks",
    rownames(.) == "potatoes_and_products" ~ "Potatoes & Products",
    rownames(.) == "low_energi_drinks" ~ "Low Energy Drinks",
    rownames(.) == "snacks" ~ "Snacks"))

pc2_loadings <- pc2_loadings %>% mutate(PC2_Loading=(PC2_Loading*-1)) # sign inversion to depict Western diet

# Assign fill color based on loading value
pc2_loadings <- pc2_loadings %>%
  mutate(
    FILL = case_when(
      PC2_Loading >= 0.09 ~ "High",          # High positive loadings
      PC2_Loading < 0.09 & PC2_Loading > -0.09 ~ "Mid",  # Near zero
      PC2_Loading <= -0.09 ~ "Low"),           # High negative loadings
    FILL = as.factor(FILL))

# Order variables by loading values
pc2_loadings <- pc2_loadings[order(pc2_loadings$PC2_Loading), ]  
# Convert Variable to a factor for ordered plotting
pc2_loadings$Variable <- factor(pc2_loadings$Variable, levels = pc2_loadings$Variable)

# Create the ggplot bar chart
food2 <- ggplot(pc2_loadings %>% filter(FILL != "Mid"), aes(y = Variable, x = PC2_Loading, fill = FILL)) +
  geom_bar(stat = 'identity', position = 'dodge', color = 'black') +  # Bars
  xlim(c(-0.5, 0.5)) +
  scale_fill_manual(values = c("darkturquoise", "darkgoldenrod2", "darkolivegreen3")) +
  guides(fill = "none") +
  xlab("PC2 Loadings") +
  ylab(NULL) +
  geom_vline(xintercept = 0, color = "black", size = 0.1) +
  ggtitle("Western Diet Pattern (7%)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.text.x = element_text(size = 15),
    axis.title = element_text(size = 18),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) +
  geom_text(data = subset(pc2_loadings, FILL == "Low"),
            aes(label = Variable, x = 0.01),
            hjust = 0, size = 6) +
  geom_text(data = subset(pc2_loadings, FILL == "High"),
            aes(label = Variable, x = -0.01),
            hjust = 1, size = 6)
food2
rm(pc2_loadings)


# create side-by-side figure
ggarrange(food1,food2, labels = c("A", "B"),
          ncol = 2, nrow = 1)

