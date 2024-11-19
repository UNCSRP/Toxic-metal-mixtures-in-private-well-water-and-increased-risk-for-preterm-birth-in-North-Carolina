# Toxic-metal-mixtures-in-private-well-water-and-increased-risk-for-preterm-birth-in-North-Carolina

Code and data for manuscript: Toxic metal mixtures in private well water and increased risk for preterm birth in North Carolina, which was published in Environmental Health in 2024 (PMID: 37845729). doi: 10.1186/s12940-023-01021-7

# Abstract 

Prenatal exposure to metals in private well water may increase the risk of preterm birth (PTB) (delivery <37 weeks’ gestation). In this study, we estimated associations between arsenic, manganese, lead, cadmium, chromium, copper, and zinc concentrations in private well water and PTB incidence in North Carolina (NC). Birth certificates from 2003-2015 (n = 1,329,071) were obtained and pregnancies were assigned exposure using the mean concentration and the percentage of tests above the maximum contaminant level (MCL) for the census tract of each individuals’ residence at the time of delivery using the NC-WELL database (117,960 well water tests from 1998-2019). We evaluated associations between single metals and PTB using adjusted logistic regression models. Metals mixtures were assessed using quantile-based g-computation. Compared with participants in other tracts, individuals residing in census tracts where >25% of tests exceeded the maximum contaminant level (MCL) for lead (aOR 1.10, 95%CI 1.02,1.18) or cadmium (aOR 1.11, 95% CI 1.00,1.23) had an increased odds of PTB. Conversely, those residing in areas with >25% MCL for zinc and copper had a reduced odds of PTB. A quartile increase in the concentrations of a mixture of lead, cadmium, and chromium was associated with a small increased odds for PTB (aOR 1.02, 95% CI 1.01, 1.03).  This metal mixture effect was most pronounced among American Indian individuals (aOR per quartile increase in all metals: 1.19 (95% CI 1.06,1.34)), highlighting racial health disparities and the need for more awareness and action to protect vulnerable populations from private well metal contamination.

# Files 

Data: 
Some data used in this analysis is not publicaly available. Specifically, geocoded birth certificate data contain sensitive information and thus are not made public in order to protect human subject data. Additionally, to protect the privacy of private well owners, the dataset that includes individual well level data, including their point locations, is not publicly available online. Please contact the authors directly for information about these datasets. 

Data uploaded in this repo:
1. Tract-level average concentrations of metals in well water following imputation: Metals_wellwater_tractlevel_distributions
2. Codebook for (1): CODEBOOK_Metals_wellwater_tractlevel_distributions

Code: 
1. Script used for imputation of metals data: 0_metals_imputation 
2. Script used to evaluate single metal exposures: 1_singlemetalmodels
3. Script used to evaluate metal mixtures exposures: 2_mixturesmodels
4. script used to repeat partial effect quantile-based g-computation: 3_partialeffects_repeats 
5. Script used for sensitivity assessment based on tract-levle percentage of well users: 4_singlemetalmodels_wwsensitivity and 5_mixturesmodels_wwsensitivity  
