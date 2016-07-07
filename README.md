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
with at most 8 components:

```
> library(sBIC)
>
> # Set a seed for reproducibility
> set.seed(123)
> 
> # Create an object representing a collection of Gaussian mixture models
> # with at most 8 components in 2 dimensions.
> gms = GaussianMixtures(maxNumComponents = 8, dim = 2, restarts = 100)
>
> # Generate some simulation data, a mixture of 3 normals.
> library(MASS)
> n = 175
> class = sample(0:2, n, replace = T)
> X = (class == 0) * mvrnorm(n, mu = c(0, 0), Sigma = diag(2)) +
>     (class == 1) * mvrnorm(n, mu = c(2.5, 2.5), Sigma = diag(2) + t(.3 * diag(2))) +
>     (class == 2) * mvrnorm(n, mu = c(-3, 2.5), Sigma = diag(2) + t(.2 * diag(2)))
> 
> # Compute the sBIC on the mixture models with the randomly generated data,
> # notice that the BIC too strongly penalizes the (true) model with 3
> # components.
> sBIC(X, gms)
$logLike
[1] -732.9610 -697.5850 -683.9564 -676.2722 -668.5807 -661.2474 -653.9223 -646.2799

$sBIC
[1] -747.6058 -729.8036 -727.4262 -728.7474 -729.9027 -731.3697 -732.8346 -733.9800

$BIC
[1] -747.6058 -729.8037 -733.7488 -743.6385 -753.5207 -763.7613 -774.0100 -783.9413

$modelPoset
[1] "GaussianMixtures: 0x10e24eca8"
```