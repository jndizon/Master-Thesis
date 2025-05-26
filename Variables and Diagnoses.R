# read-in files with covariates and relevant diagnoses for the cohort

library(dplyr) # v1.1.4
library(tidyr) # v1.3.1
library(zscorer) # v0.3.1
library(readxl) # v1.4.5

# read in file
var_10y <- readRDS("var_10y.rds")
#sex
  # male/female
#race
  # cacasian, non-caucasian, or NA
#age_when_growth_measurements
#height
#weight
#bmi
#oldchild: number of older siblings
#furanimal
# Exposure to any kind of furred animal within proband's first year of life.
# Exposure include being in the same room as animal, or touching animal if outdoors.
# Exposure can be anywhere, not limited to own home.
#furredanimalcat
# Exposure to cat at home within proband's first year of life.
#numberofcats
# How many cats were living in proband's home during first year of life?
#furredanimaldog
# Exposure to dog at home within proband's first year of life.
#numberofdogs
# How many dogs were living in proband's home during first year of life?
#furredanimalother
# Exposure to furred animal other than cat or dog at home within proband's first year of life.
#income:
#household income at child's primary home
# 1: <400.000 dKr.
# 2: 400.000-600.000 dkr.
# 3: 600.000-800.000 dKr.
# 4: 800.000-1.000.000 dKr.
# 5: >1.000.000 dKr.
#m_edu
#mother's highest completed education level
# 1: College or less. 
# 2: Medium length academical education OR tradesman's certification.
# 3: University or higher.
#mother_age_2yrs
# mother's age at child's 2nd birthday
#socialcircumstances
# PCA based on income, mother's age and mother's education at child's 2nd birthday
#passivesmoke
# Number of days exposed to passive smoking during 10th year of life.
# All indoor smoking, and smoking outside if smoker comes back inside.
#event_j45_10yr
# Indicator of event (J45 start) or censoring
#eventage_j45_10yr:
# Age in days at the event of J45
# OR not full follow
# OR age 10yrs
# whichever comes first.
# j45_10yr_ever
# Any J45-diagnosis before age 10
# j45_10yr_cross
# Any J45-diagnosis between 9th and 10th birthday.
# j45_10yr_cross_type: asthma/no asthma/intermittent
# j45_10yr_transient
# Had a J45-diagnosis earlier than 9th birthday, but no open diagnosis between 9th and 10th birthday. I.e. had ever=Yes, but cross=No.
#delivery
#vaginal delivery/emergency CS/ planned CS
#age_days
#gestational days
#antibioticsbirth
#antibiotics to mother/child at birth
#ab_child_05yr_ever
#antibiotics to children at 5 years of age
#ab_child_5yr_ever
#antibiotics to children in the first 5 years of life
#ab_child_number_05y
#number of antibiotics to children at 5 years of age
#ab_child_number_5y
#number of antibiotics to children in the first 5 years of age

-------------------
  

# create a new df containing only the variables from var_10y which I want for analyses
variables <- var_10y %>% select(abcno, sex, race, socialcircumstances, bmi, antibioticsbirth) 
variables <- variables %>%
  mutate(abmotheratbirth = case_when(
    antibioticsbirth %in% c("mother", "both") ~ "yes",  # Convert specific values to "yes"
    antibioticsbirth %in% c("proband", "no") ~ "no"))   # Convert others to "no"
variables <- variables %>% select(-antibioticsbirth)

# re-create the bmi column with z-scores instead of raw bmi

bmi <- var_10y[, -c(3, 8:33)]
bmi <- bmi %>% mutate(AGE=(age_when_growth_measurements*365.25))
bmi <- bmi %>% mutate(Sex = ifelse(sex=="male",1,2))
bmi <- bmi %>% dplyr::rename(WEIGHT=weight, HEIGHT=height)
bmi <- addWGSR(data = bmi, sex = "Sex", firstPart = "WEIGHT",
               secondPart = "HEIGHT", thirdPart = "AGE", index = "bfa",
               output = "zbmi", digits = 2)
bmi <- bmi %>% dplyr::rename(height=HEIGHT, weight=WEIGHT) %>% mutate(Sex=ifelse(Sex==1,1,0))
bmi <- bmi %>% mutate(Overweight = ifelse(zbmi >= 1.04, 1, ifelse(zbmi < 1.04, 0, NA))) %>%
  mutate(Overweight=as.factor(Overweight))
  
bmi1 <- bmi %>% select(abcno, zbmi, Overweight) 
rm(bmi)

# add zbmi + overweight to variables
variables <- merge(variables, bmi1, by = "abcno")


-------------------------
  
## DIAGNOSES: all diagnoses within the K and E series, F84, F90, and F50
# read in file
icd10diagnoses <- read_excel("20241220_C2010_diagnoses_Jacqueline.xlsx")
View(icd10diagnoses)

# rename column of names to match other dfs
icd10diagnoses <- icd10diagnoses %>% rename(abcno = copsac_id)

# subset the ADHD diagnoses in the dataset
f900diagnoses <- subset(icd10diagnoses, grepl("^F900", code)) # ADHD
f901diagnoses <- subset(icd10diagnoses, grepl("^F901", code)) # hyperkinetic conduct disorder 
 
# create a df with diagnoses to exclude: K900 (celiac), K909 (intestinal 
  # malabsorption), and ^F50 (eating disorders)
K90F50 <- icd10diagnoses[grepl("^(K90|F50)", icd10diagnoses$code), ]

----------------------------

# create a df with the var10y and diagnoses data together
variablesmerged <- variables

# insert column for those with F900/F901 diagnosis vs not
variablesmerged <- variablesmerged %>%
  mutate(F900 = if_else(abcno %in% f900diagnoses$abcno, "F900", "noF900"))
variablesmerged <- variablesmerged %>%
  mutate(F901 = if_else(abcno %in% f901diagnoses$abcno, "F901", "noF901"))
variablesmerged$allF90 <- ifelse(
  variablesmerged$F900 == "F900" | variablesmerged$F901 == "F901", "F90", "noF90")

# insert column for total energy intake and trim to only those with diet data
energy <- food10y2.0 %>% select(abcno, energi)
variablesmerged <- merge(variablesmerged, energy, by = "abcno")
rm(energy)

# exclude those with celiac, eating disorders, and intestinal malabsorption in analysis
variablesmerged <- variablesmerged [!(variablesmerged$abcno %in% K90F50$abcno), ]

# remove duplicates
variablesmerged <- variablesmerged[-which(duplicated(variablesmerged$abcno)), ]
## DF CAN NOW BE USED FOR ANALYSIS

----------------------

# BASELINE CHARACTERISTICS OF TOTAL  N = 597 COHORT

# histograms to visualize
hist(variablesmerged$zbmi, 15, xlim = c(-6, 5))
hist(variablesmerged$socialcircumstances, 15, xlim = c(-6, 5))
hist(variablesmerged$energi, 30, xlim = c(2000, 25000))

# sex 
sum(variablesmerged$sex == "male")
(309/597)*100
# race
sum(variablesmerged$race == "caucasian")
(570/597)*100

# socialcircumstances
mean(variablesmerged$socialcircumstances, na.rm = TRUE)
sd(variablesmerged$socialcircumstances, na.rm = TRUE)

# bmi
mean(variablesmerged$zbmi, na.rm = TRUE)
sd(variablesmerged$zbmi, na.rm = TRUE)

# total daily energy consumption
mean(variablesmerged$energi, na.rm = TRUE)
sd(variablesmerged$energi, na.rm = TRUE)




