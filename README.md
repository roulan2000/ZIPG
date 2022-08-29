
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ZIPG

<!-- badges: start -->
<!-- badges: end -->

We provide R code for Zero-inflated Poisson-Gamma Model (ZIPG) with an
application to longitudinal microbiome count data.

## Installation

You can install the development version of ZIPG like so:

``` r
devtools::install_github("roulan2000/ZIPG")
```

## Example

### Load Data

Load dietary data. Complete Dietary data can be found in “Daily sampling
reveals personalized diet-microbiome associations in humans.” (Johnson
et al. 2019)

``` r
library(ZIPG)
library(ggplot2)
data("Dietary")
dat = Dietary
taxa_num = 100
dat$taxa_name[taxa_num] # taxa name
#>                             OTU100 
#> "Burkholderiales bacterium 1_1_47"
W = dat$OTU[,taxa_num] # taxa count
M = dat$M # sequencing depth
ggplot(NULL)+
  geom_boxplot(aes(
  x = as.factor(dat$COV$ALC01),y=log((W+0.5)/M)))+
  labs(title = dat$taxa_name[taxa_num],
       x = 'ALC',y='Relative Abundance')+
  theme_bw()
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

## ZIPG Wald

Use function `ZIPG_main()` to run our ZIPG model.

**Input :**

`W` : Observed taxa count data.

`M` : Sequencing depth, ZIPG use log(M) as offset by default.

`X`, `X_star` : Covariates of interesting of differential abundance and
differential varibility, input as formula.

**Output list:**

`ZIPG_res$init` : pscl results, used as initialization.

`ZIPG_res$res` : ZIPG output evaluated at last EM iteration.

`ZIPG_res$res$par` : ZIPG estimation for
$\Omega = (\beta,\beta^*,\gamma)$.

`ZIPG_res$wald_test` : ZIPG Wald test

`ZIPG_res$logli` : ZIPG log-likelihood

``` r
ZIPG_res <- ZIPG_main(data = dat$COV,
                      X = ~ALC01+nutrPC1+nutrPC2, X_star = ~ ALC01,
                      W = W, M = M)
res = ZIPG_summary(ZIPG_res)
#>           ZIPG Wald  
#>        Estimation     SE     pval    
#> beta0      -7.371 0.1512 0.00e+00 ***
#> beta1       0.121 0.1985 5.41e-01    
#> beta2       0.106 0.0188 1.41e-08 ***
#> beta3      -0.118 0.0287 4.17e-05 ***
#> beta0*      0.525 0.1199 1.20e-05 ***
#> beta1*      0.606 0.1406 1.63e-05 ***
#> gamma      -2.080 0.1460 4.93e-46 ***
```

## ZIPG bWald

Set the bootstrap replicates `B` in `bWald_list` to conduct ZIPG-bWald,
results and covariance matrix can be find in `ZIPG_res$bWald`.

``` r
# Set bootstrap replicates B
bWald_list = list(B = 100)
# Wait for a wile
ZIPG_res1 = ZIPG_main(
  data = dat$COV,
  X = ~ALC01+nutrPC1+nutrPC2, X_star = ~ ALC01,
  W = W, M = M,
  bWald_list = bWald_list)
#> Running non-parametric bootstrap wald test 
#> Finish
res = ZIPG_summary(ZIPG_res1,type = 'bWald')
#>           ZIPG bWald   
#>        Estimation     SE     pval    
#> beta0      -7.371 0.1757 0.00e+00 ***
#> beta1       0.121 0.2170 5.76e-01    
#> beta2       0.106 0.0186 9.77e-09 ***
#> beta3      -0.118 0.0308 1.36e-04 ***
#> beta0*      0.525 0.1440 2.68e-04 ***
#> beta1*      0.606 0.1847 1.04e-03  **
#> gamma      -2.080 0.3168 5.23e-11 ***
res = ZIPG_CI(ZIPG_res1,type='bWald',alpha = 0.05)
#>         ZIPG Wald Confidence interval 
#>        Estimation      lb      ub
#> beta0      -7.371 -7.7150 -7.0262
#> beta1       0.121 -0.3038  0.5468
#> beta2       0.106  0.0701  0.1428
#> beta3      -0.118 -0.1780 -0.0572
#> beta0*      0.525  0.2427  0.8073
#> beta1*      0.606  0.2439  0.9680
#> gamma      -2.080 -2.7009 -1.4588
```

To test more complicated hypothesis, you may use the covariance matirx
driven from bootstrap.

``` r
round(ZIPG_res1$bWald$vcov,3)
#>        [,1]   [,2]   [,3]   [,4]   [,5]   [,6]   [,7]
#> [1,]  0.031 -0.033  0.001  0.003  0.009 -0.008  0.006
#> [2,] -0.033  0.047 -0.001 -0.003 -0.010  0.006  0.003
#> [3,]  0.001 -0.001  0.000  0.000  0.000  0.000  0.000
#> [4,]  0.003 -0.003  0.000  0.001  0.002 -0.002  0.001
#> [5,]  0.009 -0.010  0.000  0.002  0.021 -0.020 -0.009
#> [6,] -0.008  0.006  0.000 -0.002 -0.020  0.034 -0.016
#> [7,]  0.006  0.003  0.000  0.001 -0.009 -0.016  0.100
```

## ZIPG pbWald

Set bootstrap replicates `B` and the null hypothesis by formula `X0` and
`X_star0` in `pbWald_list` to conduct ZIPG-pbWald, results can be find
in ZIPG_res\$pbWald

``` r
# test beta1star, the 6th parameter
# 
pbWald_list = list(
  X0 = ~ALC01 + nutrPC1+nutrPC2,
  X_star0 = ~ 1,
  B = 100
)

ZIPG_res2 = ZIPG_main(
  data = dat$COV,
  X = ~ALC01+nutrPC1+nutrPC2, X_star = ~ ALC01,
  W = W, M = M,
  pbWald_list= pbWald_list)
#> Running parametric bootstrap wald test 
#> Finish

res = ZIPG_summary(ZIPG_res2,type ='pbWald')
#>    ZIPG pbWald 
#>  H0: beta1* = 0 
#>  pvalue =  0.0099
```
