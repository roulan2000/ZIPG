# ZIPG
## Overview

We provide R code for Zero-inflated Poisson-Gamma Model (ZIPG) witn an application to longitudinal microbiome count data. You can download R code from this github, and see more details in `ZIPG_demo.Rmd`.

## Usage

```R
ZIPG_main(W,M,X,X_star,bWald_list = NULL,pbWald_list = NULL)
```

**Input :**

W : Observed count data.

M : Sequencing depth.

X, X\_star : Covariats of interesting with intercept.

pbWald_list, bWald_list: Setting of  non-parameteric  bootstrap and parameteric bootstrap.

**Output list:**

init : pscl results, used as initialization.

res : ZIPG output evaluate at last EM iteration.

res\$par : ZIPG estimation $(\beta,\beta^*,\gamma)$.

wald_test : ZIPG Wald test with Wald statistics, SE, pvalue for each parameter.

logli : ZIPG log-likelihood

