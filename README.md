# Achieving-the-EHE-incidence-reduction-goals-in-the-US-South
This repository contains the EpiModelHIV version used for the analysis reported in “Achieving the “Ending the HIV Epidemic in the U.S.” incidence reduction goals among at-risk populations in the South.” The R package can be downloaded from the SouthComb branch of this repository.

Hamilton, D.T., Hoover, K.W., Smith, D.K. et al. Achieving the “Ending the HIV Epidemic in the U.S.” incidence reduction goals among at-risk populations in the South. BMC Public Health 23, 716 (2023). https://doi.org/10.1186/s12889-023-15563-5

and “Potential contribution of PrEP uptake by adolescents 15-17 years old to achieving the EHE incidence reduction goals in the US South” [Cite PAPER2]. 

Before proceeding it is highly recommended that researchers thoroughly read both papers and their associated technical appendix. 

# Installation

You can install the version of `EpiModelHIV` used for this analysis in R using `remotes`:
```
remotes::install_github("dth2/Achieving-the-EHE-incidence-reduction-goals-in-the-US-South",ref="SouthComb")
```

The versions for all other R packages used for this analysis can be found in the renv.lock file


# Achieving EHE goals in the South: Routine opt-out HIV screening can facilitate entry to the HIV treatment and prevention continuums of care
This repository contains the EpiModelHIV version used for the analysis reported in "Achieving EHE goals in the South: Routine opt-out HIV screening can facilitate entry to the HIV treatment and prevention continuums of care" The R package can be downloaded from the SouthComb branch of this repository.

# Installation

You can install the version of `EpiModelHIV` used for this analysis in R using `remotes`:
```
remotes::install_github("dth2/Achieving-the-EHE-incidence-reduction-goals-in-the-US-South",ref="Routine-Opt-out-testing")
```

The versions for all other R packages used for this analysis can be found in the renv.lock file

# Data requirements
The primary data sources for this analysis are ARTnet and the NSFG. The ARTnet data can be downloaded from the ARTnet repository after completing a data use agreement https://github.com/EpiModel/ArtNet. The NSFG is publicly available from the Centers for Disease Control and Prevention website https://www.cdc.gov/nchs/nsfg/nsfg_questionnaires.htm. The National Center for Health Statistics provides detailed guides for working with these data and options to download the data using SPSS, SAS, or STATA. For this study we opted to use the SPSS versions so the scripts for cleaning variable construction are written in SPSS code.     
