# Integrated Completed Likelihood (ICL) for any mixture model

ICL = BIC + 2 \* Entropy, on the "lower is better" (-2LL) scale, where
Entropy = -sum_ik z_ik \* log(z_ik) (Biernacki, Celeux & Govaert, 2000).

## Usage

``` r
icl(z, loglik, npar, n = nrow(z))
```

## Arguments

- z:

  n x K matrix of posterior class probabilities (rows sum to 1).

- loglik:

  Log-likelihood of the fitted model.

- npar:

  Number of free parameters.

- n:

  Number of observations (default: nrow(z)).

## Value

Named numeric vector: loglik, npar, n, BIC, entropy, ICL.
