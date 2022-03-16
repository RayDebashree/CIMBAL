# CIMBAL

### Description
CIMBAL is a new statistical approach for meta-analyzing cohorts with confounder imbalance such as those found in cohort collaborations. This requires two sets of cohorts: one with unadjusted or crude estimate of the exposure-outcome association, and another with fully adjusted estimate of the same association. This software is based on the following manuscript in progress: 

Ray et al. (2022+) "Meta-analysis under imbalance in measurement of confounders in cohort studies using only summary-level data". *In revision*.

**Key Words:** Collective analysis; Confounders; Confounder imbalance; Data integration; Meta-analysis; Omitted variable bias

### Requirements
R (>= 3.0.1)


### How to Install within R
```{r}
require(devtools)
source_url("https://github.com/RayDebashree/CIMBAL/blob/master/CIMBAL_v0.7.R?raw=TRUE")
```
It is recommended to download/copy the stand-alone R program in this repository, save it in your local directory of choice and `source()` it from your local directory. When a new version of the software is available, older versions may be removed from this repository, and the above `devtools::source_url()` line may not work.


### Changes
Version 0.7 - March 17, 2022
> First public release of the software.
Version 0.6 - March 31, 2021
> Software release intended for internal ECHO use only.

### Notes
1. Check the manuscript to understand when CIMBAL may or may not be appropriate to use, and how to interpret its results.
2. Requires two sets of estimates (i.e., effect estimate and its standard error) of exposure-outcome association:
    a. crude or unadjusted estimates
    b. fully adjusted estimates
2. Current implementation requires a few cohorts with fully adjusted estimates. We recommend having at least 20 such cohorts.

Contact **dray@jhu.edu** for any question on CIMBAL or to report any issue/feedback.
