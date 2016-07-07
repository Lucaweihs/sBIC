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
above paper. For details regardings the use of the sBIC with gaussian latent
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
> # with at most 6 components in two dimensions.
> gms = GaussianMixtures(maxNumComponents = 6, dim = 2)
>
> # Generate some simulation data, a mixture of 2 normals
> library(MASS)
> class = sample(0:1, 100, replace = T)
> X = (class == 0) * mvrnorm(100, mu = c(0,0), Sigma = diag(2)) +
>     (class == 1) * mvrnorm(100, mu = c(2.5,2.5), Sigma = diag(2) + t(.3 * diag(2)))
> 
> # Compute the sBIC on the mixture models with the randomly generated data,
> # notice that the BIC too strongly penalizes the (true) model with 2
> # components.
> sBIC(X, gms)
$logLike
[1] -353.7890 -341.3404 -336.1208 -328.5612 -319.4646 -309.0157

$sBIC
[1] -367.0348 -364.8352 -367.6173 -368.0392 -366.9820 -364.9812

$BIC
[1] -367.0348 -370.4811 -381.1565 -389.4919 -396.2902 -401.7362

$modelPoset
[1] "GaussianMixtures: 0x10b852200"
```