# sBIC Package

## Purpose

This package allows you to compute the sinuglar bayesian information criterion
as described in

Mathias Drton, and Martyn Plummer. "A Bayesian information criterion for
singular models." arXiv preprint arXiv:1309.0911 (2015).

for collections of the following model types

1. Binomial mixture
2. Gaussian mixture
3. Latent class analysis
4. Gaussian latent forest
5. Reduced rank regression
6. Factor analysis

All of these models, excluding gaussian latent forests, are described in the
above paper. For details regardings the use of the sBIC iwht gaussian latent
forests see 

Mathias Drton, Shaowei Lin, Luca Weihs, and Piotr Zwiernik. "Marginal likelihood
and model selection for Gaussian latent tree and forest models." arXiv preprint
arXiv:1412.8285 (2014).

## Example

This package makes extensive use of the R.oo package which allows for the use
of some object oriented principles in R. While not strictly necessary to use
this package it may be helpful to read sections 1 and 2 of
http://www1.maths.lth.se/help/R/R.oo/ which serve as an introduction to R.oo.

Each collection of models is defined as its own class, to see all available
classes check the help for the sBIC function. As an example for how to use the
package we will compute the sBIC for a collection of gaussian mixture models
with at most 6 components:

```
> library(sBIC)
>
> # Set a seed for reproducibility
> set.seed(123)
> 
> # Create an object representing a collection of Gaussian mixture models
> # with at most 6 components.
> gms = GaussianMixtures(maxNumComponents = 6)
>
> # Generate some simulation data, a mixture of 2 normals
> class = sample(1:2, 100, replace = T)
> X = (class == 1) * rnorm(100) + (class == 2) * rnorm(100, 3)
> 
> # Compute the sBIC on the mixture models with the randomly generated data,
> # notice that the BIC too strongly penalizes the (true) model with 2
> # components.
> sBIC(X, gms)
$logLike
[1] -196.6984 -191.3533 -190.1739 -189.0820 -188.0583 -188.7303

$sBIC
[1] -201.3036 -200.5109 -202.8495 -205.2173 -207.6479 -211.7738

$BIC
[1] -201.3036 -202.8662 -208.5946 -214.4105 -220.2945 -227.8742

$modelPoset
[1] "GaussianMixtures: 0x10c2954b0"

$modelPoset
[1] "GaussianMixtures: 0x10c2954b0"
```