# CIMBAL

### Description
CIMBAL is a new statistical approach for meta-analyzing cohorts with confounder imbalance such as those found in cohort collaborations. This requires two sets of cohorts: one with unadjusted or crude estimate of the exposure-outcome association, and another with fully adjusted estimate of the same association. This software is based on the following manuscript in progress: 

Ray et al. "Meta-analysis under Imbalance in Measurement of Confounders in Cohort Collaborations". *In progress*.

**Key Words:** Collective analysis; Confounders; Confounder imbalance; Meta-analysis

### Requirements
R (>= 3.0.1)


### How to Install within R
```{r}
require(devtools)
source_url("https://github.com/RayDebashree/metaUSAT/blob/master/CIMBAL_v0.6.R?raw=TRUE")
```
It is recommended to download/copy the stand-alone R program in this repository, save it in your local directory of choice and `source()` it from your local directory. When a new version of the software is available, older versions may be removed from this repository, and the above `devtools::source_url()` line may not work.


### Changes
Version 0.6 - March 31, 2021
> Software release intended for internal ECHO use only.
