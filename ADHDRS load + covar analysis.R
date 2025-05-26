# loading ADHD rating scales data, pre-processing, and covariate analysis

library(dplyr) # v1.1.4
library(readr) # v2.1.5
library(readxl) # v1.4.5

# read in file
ADHDRS_10y <- read_csv("ADHDRS_10y.csv")
View(ADHDRS_10y)

# rename columns
ADHDRS_10y <- ADHDRS_10y %>% rename(abcno = ABCNO,
                                    attention = ADHDRS_opmaerksomhed_sc_F,
                                    implus_hyperact = ADHDRS_impuls_hyperakt_sc_F,
                                    behavior_subscale = ADHDRS_adfaerdsforstyrrelse_sc_F,
                                    total_score = ADHDRS_ADHD_CD_symptoms_sc_F,
                                    total_wo_behavior = ADHDRS_ADHD_symptoms_sc_F)

----------------------------
  
# creating age variable: age at which they filled in the ADHDRS
birthdates <- read_excel("birth.xlsx") %>% rename(abcno = ABCNO)
ADHDRS_10y <- merge(ADHDRS_10y, birthdates, by = "abcno") 
ADHDRS_10y$STARTDATE <- as.Date(ADHDRS_10y$STARTDATE, format = "%m/%d/%y")
ADHDRS_10y$BIRTHDATE <- as.Date(ADHDRS_10y$BIRTHDATE)
ADHDRS_10y$ageatADHDRS <- as.integer((ADHDRS_10y$STARTDATE - ADHDRS_10y$BIRTHDATE) / 365.25)
ADHDRS_10y <- ADHDRS_10y %>% select(-STARTDATE, - BIRTHDATE)

--------------------------

## COVARIATE ANALYSIS WITH ADHD-RS
# REQUIRES DF FROM "var10y and diagnoses"
# dataset: n = 575

# create df with covariates and the ADHDRS
adhdrs_covar <- merge(variablesmerged, ADHDRS_10y, by = "abcno") %>% filter(!is.na(total_score))

# BASELINE CHARACTERISTICS

# histograms to visualize distribution
hist(adhdrs_covar$attention, 10, xlim = c(0,30))
hist(adhdrs_covar$implus_hyperact, 10, xlim = c(0,30))
hist(adhdrs_covar$behavior_subscale, 15, xlim = c(0,30))
hist(adhdrs_covar$total_score, 15, xlim = c(0,70))

# number of participants with ADHD or HCD diagnoses
sum(adhdrs_covar$F900 == "F900", na.rm = TRUE)
sum(adhdrs_covar$F901 == "F901", na.rm = TRUE)
# age at time of ADHD-RS collection
mean(adhdrs_covar$ageatADHDRS, na.rm = TRUE)
sd(adhdrs_covar$ageatADHDRS, na.rm = TRUE)
# sex 
sum(adhdrs_covar$sex == "male")
(301/575)*100 
# race
sum(adhdrs_covar$race == "caucasian")
(549/575)*100
# social circumstances
mean(adhdrs_covar$socialcircumstances, na.rm = TRUE)
sd(adhdrs_covar$socialcircumstances, na.rm = TRUE)
# bmi
mean(adhdrs_covar$zbmi, na.rm = TRUE)
sd(adhdrs_covar$zbmi, na.rm = TRUE)
# total daily energy consumption
mean(adhdrs_covar$energi, na.rm = TRUE)
sd(adhdrs_covar$energi, na.rm = TRUE)

# analyze based ONLY on total_score
adhdrs_covars_fullmodel <- lm(total_score ~ sex + race + socialcircumstances + zbmi + energi + ageatADHDRS, data = adhdrs_covar)
anova(adhdrs_covars_fullmodel)
nobs(adhdrs_covars_fullmodel)
# [1] 572 -> n = 572 were included in the analysis
adhdrs_covars_fullmodel0 <- lm(total_score ~ sex + race + socialcircumstances + zbmi + energi, data = adhdrs_covar) # remove age
anova(adhdrs_covars_fullmodel0)
# no change in remaining signif values
adhdrs_covars_fullmodel1 <- lm(total_score ~ sex + race + socialcircumstances + zbmi, data = adhdrs_covar) # remove energi
anova(adhdrs_covars_fullmodel1)
  # no change in remaining signif values
adhdrs_covars_fullmodel2 <- lm(total_score ~ sex + socialcircumstances + zbmi, data = adhdrs_covar) 
anova(adhdrs_covars_fullmodel2)
  # sex, social circumstances, and zbmi are significantly changing with ADHD-RS total score




