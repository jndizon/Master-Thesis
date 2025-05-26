# Master-Thesis

This repository presents the code used for the masters thesis project, "Investigating correlations between the Western diet pattern, ADHD symptoms, and the gut microbiome composition in 10-year-old Danish children (COPSAC2010)." This thesis was carried out from Novemeber 2024 to June 2025 for the completion of the MSc Human Nutrition from the University of Copenhagen. This project is a data exploration study investigating correlations between ADHD, diet patterns, and the gut microbiome composition in 10-year-old Danish children. The data is collected from the Copenhagen Prospective Studies on Asthma in Childhood (COPSAC2010) project. 

The code requires multiple input dataframes not publically avaialable. For all dataframes, samples are rows.
These are:
foodgroups10y and nutrients10y: food group and nutrient intakes, respectively, in g/day.
picky_13y: raw scores from questionnaire responses
var_10y: misc. variables relating to the cohort, including sex, race, etc.
icd10diagnoses: includes relevant diagnoses for participants (ICD-10 diagnoses within the K and E series, F84, F90, and F50)
ADHDRS_10y: raw sub-scores from ADHD-RS questionnaire
birthdates: the birthdates of the participants
phy_10y: phyloseq object generated prior to the present project

The scripts build on one another and are therefore intended to be run in the following order:
1. Variables and Diagnoses
2. Dietary data load + pre-processing
3. ADHDRS load + covar analysis
4. Phyloseq pre-process + diversity
5. Phyloseq QC + covar analysis
6. Creating diet patterns
7. Diet pattern characterization + covar analysis
8. Diet patterns + ADHD-RS analysis
9. Diet patterns + GM analysis
10. Picky eating analysis
11. ADHD-RS + GM analysis
12. PC algorithm + CPDAGs
13. misc plots+tables
