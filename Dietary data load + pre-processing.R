# load and pre-process files with food group and nutrient intake data,
  # plus load and pre-process picky eating data

library(dplyr) # v1.1.4
library(readr) # v2.1.5

## FOOD GROUP DATA PRE-PROCESSING

# foodgroups10y is pre-loaded into this project
colnames(foodgroups10y)
# [1] "ABCNO"                  "energi"                 "total_amount_consumed"  "low_fat_milk"          
# [5] "hig_fat_milk"           "ice_cream"              "breakfast_cereals"      "whole_grains"          
# [9] "refined_grains"         "fruits"                 "offal"                  "processed_meat"        
# [13] "red_meat"               "fish"                   "shelfish"               "poultry"               
# [17] "eggs"                   "animal_fat"             "vegetable_fat"          "margarine"             
# [21] "sweets_desserts"        "wine"                   "bear"                   "tea"                   
# [25] "coffe"                  "high_energi_drinks"     "low_energi_drinks"      "water"                 
# [29] "spirits"                "snacks"                 "vegetable_juice"        "fruit_juice"           
# [33] "fruit_syrup_marm"       "nuts"                   "other_vegetables"       "potatoes_and_products" 
# [37] "green_leafy_vegetables" "tomatoes"               "soy_products"           "dried_fruits"          
# [41] "dressings"              "beans"                  "cheese"       

## DEALING WITH DUPLICATES:
# df of only duplicates 
food10ydups <- foodgroups10y[foodgroups10y$ABCNO %in% foodgroups10y$ABCNO[duplicated(foodgroups10y$ABCNO)], ]

# averaging the data from each abcno with duplicate entries
Fdups1 <- foodgroups10y[c("275","449"), ]
Faveraged1 <- colMeans(Fdups1)
Fdups2 <- foodgroups10y[c("67","319"), ]
Faveraged2 <- colMeans(Fdups2)
Fdups3 <- foodgroups10y[c("361","460"), ]
Faveraged3 <- colMeans(Fdups3)
Fdups4 <- foodgroups10y[c("263","613"), ]
Faveraged4 <- colMeans(Fdups4)
Fdups5 <- foodgroups10y[c("302","528"), ]
Faveraged5 <- colMeans(Fdups5)

FoodAveragedRows <- data.frame(Faveraged1, Faveraged2, Faveraged3, 
                               Faveraged4, Faveraged5)
food10y1.0 <- foodgroups10y
food10y1.0 <- rbind(food10y1.0, t(FoodAveragedRows))
# adding the rows of averaged/merged data but still need to remove the original rows

# creating the df without the original entries since they're now accounted 
# for in the merged entries
food10y2.0 <- food10y1.0[!(rownames(food10y1.0) %in% 
                             c("275", "449", "67", "319","361", "460", 
                               "263", "613", "302", "528")), ]

# labeling the rows of the new df by abcno
row.names(food10y2.0) <- food10y2.0[,1]

# remove all of the created objects
rm(list = ls(pattern = "^Fdups"))
rm(food10ydups)
rm(list = ls(pattern = "^Fav"))
rm(FoodAveragedRows)
rm(food10y1.0)

---------------
  
# tidying df
food10y2.0 <- food10y2.0 %>% rename(abcno = ABCNO)
food10y2.0 <- food10y2.0[sapply(food10y2.0, function(x) length(unique(x))  > 1)]

# EXCLUDING FROM ANALYSIS
# exclusion based on unrealistic total energy intakes:
summary(foodgroups10y$energi)
sd(foodgroups10y$energi)
food10y2.0 <- food10y2.0 %>% filter(energi < 19700 & energi > 2000)

## remove those with celiac's, intestinal malabsorption, or eating disorders
# requires df from "var10y and diagnoses"
food10y2.0 <- food10y2.0 [!(food10y2.0$abcno %in% K90F50$abcno), ]

# final number included in diet-related analysis: n = 597

--------------------------------

## NUTRIENT INTAKE DATA PRE-PROCESSING

# nutrients10y is already pre-loaded in the project
colnames(nutrients10y)
# [1] "ABCNO"    "energi"   "prot"     "totaln"   "fedt"     "sfa"      "mufa"     "pufa"     "kulhtil"  "kulhy"   
# [11] "suktil"   "kfibre"   "alko"     "aske"     "vand"     "avit"     "retino"   "betaca"   "dvit"     "d3cho"   
# [21] "d2erca"   "v5hydr"   "evit"     "alfato"   "kvit"     "thiami"   "ribofl"   "neniac"   "niacin"   "ntrypto" 
# [31] "b6vit"    "pantot"   "biotin"   "folaci"   "b12vit"   "cvit"     "lascorb"  "ldehydr"  "natriu"   "kaliu"   
# [41] "calciu"   "magnes"   "fosfor"   "jern"     "kobber"   "zink"     "jod"      "mangan"   "chrom"    "selen"   
# [51] "nikkel"   "frukto"   "glukos"   "laktos"   "maltos"   "saccha"   "sukker"   "stivel"   "totkfib"  "c4x0"    
# [61] "c6x0"     "c8x0"     "c10x0"    "c12x0"    "c14x0"    "c15x0"    "c16x0"    "c17x0"    "c18x0"    "c20x0"   
# [71] "c22x0"    "c24x0"    "fedtms"   "c14x1"    "c16x1"    "c18x1n9"  "c18x1n7"  "c20x1n11" "c22x1n9"  "c22x1n11"
# [81] "c24x1n9"  "fedtus"   "c18x2n6"  "c18x3n3"  "c18x4n3"  "c20x4n6"  "c20x5n3"  "c22x5n3"  "c22x6n3"  "andrfs"  
# [91] "fedtps"   "sumn3"    "sumn6"    "transf"   "choles"   "isoleu"   "leucin"   "lysin"    "methio"   "cystin"  
# [101] "phenyl"   "thyros"   "threon"   "trypto"   "valin"    "argini"   "histid"   "alanin"   "aspara"   "glutam"  
# [111] "glycin"   "prolin"   "serin"    "bsmbgrp1" "prote"    "fedte"    "kulhye"   "alkoe"    "energybr" "energxbr"
# [121] "vaegt"   
  
# Removing duplicate entries from nturients10y df
# df of only duplicates 
nutrients10ydups <- nutrients10y[nutrients10y$ABCNO %in% nutrients10y$ABCNO[duplicated(nutrients10y$ABCNO)], ]

# averaging the data from each abcno with duplicate entries
dups1 <- nutrients10y[c("275","449"), ]
averaged1 <- colMeans(dups1)
dups2 <- nutrients10y[c("67","319"), ]
averaged2 <- colMeans(dups2)
dups3 <- nutrients10y[c("361","460"), ]
averaged3 <- colMeans(dups3)
dups4 <- nutrients10y[c("263","613"), ]
averaged4 <- colMeans(dups4)
dups5 <- nutrients10y[c("302","528"), ]
averaged5 <- colMeans(dups5)

NutrientsAveragedRows <- data.frame(averaged1, averaged2, averaged3, 
                                    averaged4, averaged5)

nutrients10y1.0 <- nutrients10y
nutrients10y1.0 <- rbind(nutrients10y1.0, t(NutrientsAveragedRows))

# creating the df without the original entries since they're now accounted 
# for in the merged entries
nutrients10y2.0 <- nutrients10y1.0[!(rownames(nutrients10y1.0) %in% 
                             c("275", "449", "67", "319","361", "460", 
                               "263", "613", "302", "528")), ]

# labeling the rows of the new df by abcno
row.names(nutrients10y2.0) <- nutrients10y2.0[,1]

# remove all of the created objects
rm(list = ls(pattern = "^dups"))
rm(nutrients10ydups)
rm(list = ls(pattern = "^av"))
rm(NutrientsAveragedRows)
rm(nutrients10y1.0)

-------------------------

# tidying df
nutrients10y2.0 <- nutrients10y2.0 %>% rename(abcno = ABCNO)
nutrients10y2.0 <- nutrients10y2.0[sapply(nutrients10y2.0, function(x) length(unique(x))  > 1)]  

# EXCLUDING FROM ANALYSIS
# exclusion based on unrealistic total energy intakes:
nutrients10y2.0 <- nutrients10y2.0 %>% filter(energi < 19700 & energi > 2000)

## remove those with celiac's, intestinal malabsorption, or eating disorders
# requires df from "var10y and diagnoses"
nutrients10y2.0 <- nutrients10y2.0 [!(nutrients10y2.0$abcno %in% K90F50$abcno), ]

# final number included in diet-related analysis: n = 597

--------------------------------------

## REMOVING COLUMNS WHICH ARE TOTALS AND ENERGY CALCULATIONS BASED ON NUTRIENT CONSTITUENTS
nutrients10y3.0 <- nutrients10y2.0 %>% select(-prot,-fedt,-sfa,-mufa,-pufa, -kulhy, -kulhtil, 
                                              -sukker, -fedtms , -fedtus ,-fedtps ,- sumn3 ,- sumn6, -totaln, -bsmbgrp1:-vaegt)
# creating df with only the relevant grouped data 
nutrients10y4.0 <- nutrients10y2.0 %>% select(abcno, prot, fedt, sfa, mufa, pufa, sumn3, kulhy, suktil, kfibre)


# raw nutrient intakes (same as with scale(df), center = T)
nutr_zscores <- cbind(
  nutrients10y3.0[, 1, drop = FALSE],  # keep first column as-is
  as.data.frame(apply(nutrients10y3.0[, -1], 2, function(x) (x - mean(x)) / sd(x))))  # z-score the rest
nutr_macros_zscores <- cbind(
  nutrients10y4.0[, 1, drop = FALSE],  # keep first column as-is
  as.data.frame(apply(nutrients10y4.0[, -1], 2, function(x) (x - mean(x)) / sd(x))))  # z-score the rest

------------------------------------------

# PICKY EATING DATA PRE-PROCESSING AND ANALYSIS
answers_13y <- read_delim("13y food behavior/dataset.csv", 
                          delim = ";", escape_double = FALSE, trim_ws = TRUE)
# remove all columns outside the questions on picky eating
picky_13y <- answers_13y%>% select(s_1, s_331:s_336)
picky_13y <- picky_13y %>% dplyr::rename(abcno = `s_1`)
picky_13y <- picky_13y %>% filter(!if_any(everything(), is.na))

# remove duplicates
picky_13y <- picky_13y[-which(duplicated(picky_13y$abcno)), ]

## processing the answers for analysis 
# invert answers for s_334 to s_336 so all high scores = picky
# max_value = 5 // min_value =1  -> max val + min val - actual entry -> inverted val 
picky_13y <- picky_13y %>%
  mutate(across(c(s_334:s_336), ~ 6 - .))

# subtract 1 from all entries so the lowest score will be zero instead of 1
picky_13y <- picky_13y %>%
  mutate(across(c(s_331:s_336), ~ . - 1))

# create a new column with a total picky eating score
picky_13y$picky_score <- rowSums(picky_13y[, -1])
# 0 = not picky at all, 24 = completely picky

# characterizing the data
# visualize distribution of picky eating in the cohort
hist(picky_13y$picky_score, 24, xlim = c(0,25))

# characterize based solely on the respondents with diet data available
picky_13y_trimmed <- picky_13y %>% merge(food_pc, by = "abcno")

mean(picky_13y_trimmed$picky_score) # 9.25
median(picky_13y_trimmed$picky_score) # 8
IQR(picky_13y_trimmed$picky_score) # 10
quantile(picky_13y_trimmed$picky_score) 
# 0%  25%  50%  75% 100% 
# 0    4    8   14   24 
