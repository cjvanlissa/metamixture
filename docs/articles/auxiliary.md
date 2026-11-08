# Auxiliary Analyses for Random-Effects Mixture Meta-Analysis

Random-effects mixture meta-analysis can reveal unobserved subgroups of
studies with different underlying effects or heterogeneity. Once a
[`remix()`](https://cjvanlissa.github.io/metamixture/reference/remix.md)
model has been estimated, researchers may wish to investigate whether
those latent classes differ on an *auxiliary variable*: a study
characteristic not used to fit the original mixture model. An auxiliary
variable can be a predictor of class membership, or an outcome of class
membership, or even an entire auxiliary model whose parameters differ
between classes.

This vignette illustrates three approaches available through the
`tidySEM` ecosystem: Bolck–Croon–Hagenaars (BCH) weighting, pseudo-class
draws, and two-step maximum likelihood. These methods account for
class-membership uncertainty in different ways. There are important
differences between their underlying assumptions.

## Simulating Data

We will simulate data from two latent subpopulations with different true
effect sizes `yi` and their sampling variances `vi`. An auxiliary
variable `z` is a study characteristic associated with class membership
but was not included in the original mixture fit.

``` r
library(metamixture)
#> Loading required package: OpenMx
#> Warning: package 'OpenMx' was built under R version 4.6.1
#> Loading required package: tidySEM
library(tidySEM)
library(OpenMx)

set.seed(956)
n <- 200
class_true <- sample(1:2, n, replace = TRUE, prob = c(0.55, 0.45))
vi <- runif(n, 0.02, 0.12)
yi <- rnorm(n, mean = c(0, 2)[class_true], sd = sqrt(vi + 0.3))
z <- rnorm(n, mean = c(0, 1)[class_true], sd = 1)

dat <- data.frame(yi = yi, vi = vi, z = z)
```

Next, we fit the `remix` model:

``` r
fit <- remix(
  yi = dat$yi,
  vi = dat$vi,
  classes = 2,
  tau2 = "equal"
  )
```

### Inspect the class solution

Before conducting auxiliary analyses, let’s inspect model fit, the
estimated class proportions, and posterior class-membership
probabilities. Poor separation makes it particularly important not to
treat modal assignments as observed classes.

``` r
table_fit(fit)
#>   Minus2LogLikelihood   n Parameters observedStatistics  df RMSEASquared
#> 1            579.0184 200          4                  0 196            0
#>   RMSEANull   modelName      AIC      BIC    saBIC Classes   Entropy  prob_min
#> 1      0.05 RE_Mixture2 587.0184 600.2116 587.5392       2 0.8832082 0.9673427
#>    prob_max n_min n_max        LL
#> 1 0.9695407 0.455 0.545 -289.5092
```

For a two-class solution, we can glean a lot of information from the fit
table. We see that the classes are approximately equally sized: `n_min`
is the proportion of the sample assigned to the smallest class, and
`n_max` the proportion assigned to the largest class.

We also see that cases belonging to each true class have a high
posterior probability of being assigned to the matching latent class
(`prob_min` and `prob_max`). These two quantities are derived from the
classification probabilities table:

``` r
class_prob(fit, "mostlikely.class")
#> $mostlikely.class
#>           assigned.1 assigned.2
#> avgprob.1 0.96734272 0.03265728
#> avgprob.2 0.03045934 0.96954066
```

If C is the true class of an observation (in rows), and N is the
assigned class based on highest posterior class probability (in
columns), then this table shows the probability that cases that belong
to true class C are assigned to latent class N. The diagonal thus
represents the probability that observations in each class will be
correctly classified (classification accuracy); values near 1 are
excellent.

The fit table also gives an overall metric of class separation: the
entropy criterion. Values over $`0.80`$ indicate good separation. If
separation was excellent - say $`0.95`$ or greater - then we could
essentially tread predicted class as an observed variable, with little
measurement error:

``` r
metamixture:::predict_class_remix(fit)
#>   [1] 2 2 1 1 2 1 2 2 2 1 1 2 1 1 1 1 2 2 1 1 1 2 2 2 2 1 2 2 2 2 1 1 2 2 1 1 1
#>  [38] 1 1 2 1 1 1 1 2 1 2 1 1 1 1 2 2 1 2 1 1 2 2 1 1 1 1 1 2 2 1 2 1 1 1 1 2 1
#>  [75] 2 2 2 2 2 2 1 1 2 1 1 1 1 2 2 1 2 1 1 1 1 1 2 2 2 1 2 2 1 2 1 1 1 1 1 1 2
#> [112] 2 1 2 2 2 1 1 1 1 1 2 2 1 2 1 1 2 1 2 1 2 1 1 2 2 2 2 1 1 2 2 1 2 1 1 2 2
#> [149] 2 1 1 1 2 2 2 1 1 2 2 1 1 1 2 1 2 2 1 2 1 1 2 1 1 1 1 2 1 1 2 1 1 1 1 2 2
#> [186] 2 2 2 2 1 2 1 1 2 1 2 1 1 2 2
```

Yet even with such excellent correct classification probabilities,
errors will occur:

``` r
table(true = class_true, pred = metamixture:::predict_class_remix(fit))
#>     pred
#> true   1   2
#>    1 106   9
#>    2   3  82
```

Therefore, we must account for classification error somehow.

## BCH: comparing a distal outcome across classes

The BCH method uses classification-error-adjusted weights to estimate
class-specific auxiliary models without re-estimating the measurement
model. Note that the auxiliary data must have exactly the same number of
cases in the same order as the observations used to estimate
[`remix()`](https://cjvanlissa.github.io/metamixture/reference/remix.md).
Be careful: prepare one `data.frame` for your
[`remix()`](https://cjvanlissa.github.io/metamixture/reference/remix.md)
analysis and the subsequent auxiliary analyses, then select your
variables from the same `data.frame` to prevent mix-ups. Do not
independently sort or filter the auxiliary data after fitting the
measurement model.

For a single continuous auxiliary variable, the simplest call is:

``` r
bch_fit <- BCH(fit, data = dat$z)
table_results(bch_fit)
#> Warning: Unavailable columns requested: class
#>                label est_sig   se pval       confint  group
#> 1     Means.y.class1   -0.02 0.09 0.82 [-0.20, 0.16] class1
#> 2 Variances.y.class1 0.91*** 0.12 0.00  [0.67, 1.15] class1
#> 3     Means.y.class2 0.96*** 0.11 0.00  [0.74, 1.17] class2
#> 4 Variances.y.class2 1.10*** 0.16 0.00  [0.78, 1.43] class2
```

We see that the mean of class 2 is higher than that of class 1. Note
that the standard errors of BCH analyses are not valid. Do not interpret
the resulting p-values as if they properly account for classification
error.

We can test the significance of the mean difference with a likelihood
ratio test:

``` r
bch_test <- lr_test(bch_fit, compare = "M")
bch_test
#> BCH test for equality across classes
#> 
#> Overall likelihood ratio test:
#>  LL_baseline LL_restricted   LL_dif df            p
#>     566.1196      607.5235 41.40394  1 1.238077e-10
#> 
#> Pairwise comparisons using likelihood ratio tests:
#>  Model1 Model2 LL_baseline LL_restricted   LL_dif df            p
#>  class1 class2    566.1196      607.5235 41.40394  1 1.238077e-10
```

Note that, for a single auxiliary variable, a significant between-class
difference can be interpreted either as evidence that class causes
differences in this auxiliary variable - or as evidence that the
auxiliary variable causes class membership. The interpretation depends
on theoretical assumptions, not statistical results.

## Two-step maximum likelihood

Two-step maximum likelihood starts from a fitted measurement model,
fixes its parameters and mixing proportions, and then estimates
class-specific auxiliary parameters using the joint likelihood. Unlike
modal assignment, this likelihood integrates over uncertain class
membership. It can therefore produce valid standard errors that account
for classification uncertainty.

``` r
fit_two <- two_step(
  x = fit,
  model = "z ~ 1",
  data = dat["z"]
)
```

We can get the results as follows:

``` r
table_results(fit_two)
#> Warning: Unavailable columns requested: class
#>              label est_sig   se pval       confint group
#> 1     Means.z.aux1   -0.03 0.09 0.72 [-0.22, 0.15]  aux1
#> 2 Variances.z.aux1 0.96*** 0.13 0.00  [0.71, 1.22]  aux1
#> 3     Means.z.aux2 0.99*** 0.11 0.00  [0.77, 1.21]  aux2
#> 4 Variances.z.aux2 0.99*** 0.15 0.00  [0.68, 1.29]  aux2
#> 5  weights[1,1].NA    0.55 0.00 <NA>          <NA>  <NA>
#> 6  weights[1,2].NA    0.45 0.00 <NA>          <NA>  <NA>
#> 
#> Note: adjusted standard errors were used for: Variances.z.aux1, Means.z.aux1, Variances.z.aux2, Means.z.aux2.
```

Notice the message informing you that adjusted standard errors were used
for the auxiliary parameters. These adjusted standard errors are valid
for inference on the auxiliary parameters, as they take classification
uncertainty into account.

As with
[`BCH()`](https://cjvanlissa.github.io/tidySEM/reference/BCH.html), it
is possible to perform likelihood ratio tests on parameter matrices, for
example, the `M` matrix (means):

``` r
lr_test(fit_two, compare = "M")
#> BCH test for equality across classes
#> 
#> Overall likelihood ratio test:
#>  LL_baseline LL_restricted   LL_dif df            p
#>     1145.927      1187.533 41.60663  1 1.116141e-10
#> 
#> Pairwise comparisons using likelihood ratio tests:
#>  Model1 Model2 LL_baseline LL_restricted   LL_dif df            p
#>    aux1   aux2    1145.927      1187.533 41.60663  1 1.116141e-10
```

## One-step estimation

Another approach is the one-step solution, which jointly estimates the
measurement and auxiliary models. In a one-step model, the inclusion of
the auxiliary model can influence the estimated class solution. This is
appropriate when the goal is to estimate a single joint model, but may
be less desirable when researchers want to perform follow-up analyses on
a previously established class solution.

The one-step method has the same interface as
[`BCH()`](https://cjvanlissa.github.io/tidySEM/reference/BCH.html) and
[`two_step()`](https://cjvanlissa.github.io/tidySEM/reference/two_step.html):

``` r
fit_one <- one_step(
  x = fit,
  data = dat$z
)
table_results(fit_one,  columns = c("label", "est", "se", "confint", "group"))
#>               label   est   se       confint group
#> 1      Means.y.aux1 -0.03 0.09 [-0.22, 0.15]  aux1
#> 2  Variances.y.aux1  0.96 0.13  [0.71, 1.22]  aux1
#> 3      Means.y.aux2  1.00 0.11  [0.78, 1.22]  aux2
#> 4  Variances.y.aux2  0.98 0.15  [0.68, 1.29]  aux2
#> 5           eta1.NA  0.24 0.15 [-0.05, 0.54]  <NA>
#> 6        Means.1.NA  0.01 0.06 [-0.11, 0.12]  <NA>
#> 7        Means.2.NA  2.05 0.07  [1.92, 2.19]  <NA>
#> 8    Variances.1.NA  0.25 0.04  [0.18, 0.32]  <NA>
#> 9   weights[1,1].NA  0.56 0.04  [0.49, 0.63]  <NA>
#> 10  weights[1,2].NA  0.44 0.04  [0.37, 0.51]  <NA>
```

Note that the results include *all* model parameters, including the
latent class parameters, as they are all free in this model. Compare
these results to those of the measurement model, and the two-step
solution. Note that the class parameters are not reported in the results
of `fit_two`, as they are fixed and excluded from model optimization:

``` r
table_results(fit, columns = c("label", "est", "se", "confint"))
#>          label   est   se       confint
#> 1         eta1  0.20 0.15 [-0.10, 0.50]
#> 2      Means.1 -0.02 0.06 [-0.13, 0.10]
#> 3      Means.2  2.03 0.07  [1.89, 2.16]
#> 4  Variances.1  0.25 0.04  [0.19, 0.33]
#> 5 weights[1,1]  0.55 0.04  [0.47, 0.62]
#> 6 weights[1,2]  0.45 0.04  [0.38, 0.53]
table_results(fit_two, columns = c("label", "est", "se", "confint", "group"))
#>              label   est   se       confint group
#> 1     Means.z.aux1 -0.03 0.09 [-0.22, 0.15]  aux1
#> 2 Variances.z.aux1  0.96 0.13  [0.71, 1.22]  aux1
#> 3     Means.z.aux2  0.99 0.11  [0.77, 1.21]  aux2
#> 4 Variances.z.aux2  0.99 0.15  [0.68, 1.29]  aux2
#> 5  weights[1,1].NA  0.55 0.00          <NA>  <NA>
#> 6  weights[1,2].NA  0.45 0.00          <NA>  <NA>
#> 
#> Note: adjusted standard errors were used for: Variances.z.aux1, Means.z.aux1, Variances.z.aux2, Means.z.aux2.
```

Notice that the estimates for the measurement model differ in the
one-step model, as the auxiliary variable influences their estimation.
Two-step estimation instead preserves the original measurement model.

## Pseudo-class analysis

A different solution for accounting for classification uncertainty is
pseudo-class analysis. It generates multiple datasets; in each datasets,
individual studies’ latent class membership is randomly drawn with
probabilities corresponding to its posterior class probabilities. Thus,
a study with class probabilities `c(0.97, 0.03)` will be assigned to
class 1 in 97% of datasets, whereas a study with class probabilities
`c(0.45, 0.55)` will be assigned to class 1 in 45% of datasets.

The auxiliary analysis is fitted in each dataset, and estimates are
pooled across datasets, thus incorporating the classification
uncertainty.

### Predicting an Outcome from Class Membership

The following example estimates class-specific means using an ordinary
linear model. The `class` variable is generated by the pseudo-class
procedure. Using the model specification `~ factor(class)` means that we
treat class as a two-group factor (otherwise it will incorrectly be
treated as numeric!), and `-1` means the intercept should be removed
from the model, so we get two class-specific means in the output:

``` r
pseudo_fit <- tidySEM::pseudo_class(
  x = fit,
  model = lm(z ~ factor(class) - 1, data = data),
  data = dat["z"],
  m = 20
)
#> Registered S3 method overwritten by 'lme4':
#>   method           from
#>   na.action.merMod car
pseudo_fit
#>             term    estimate  std.error   statistic       df      p.value
#> 1 factor(class)1 0.000166715 0.09841807 0.001693947 186.0779 9.986502e-01
#> 2 factor(class)2 0.936606202 0.11272493 8.308775861 161.8635 3.691711e-14
```

Unlike BCH and two-step estimation, pseudo-class analysis can
accommodate auxiliary analyses implemented using any ordinary R modeling
functions, as long as they are compatible with the pooling methods
supported by `tidySEM`.

### Predicting Class Membership from an Auxiliary Variable

So far, we have examined whether latent classes differ on auxiliary
variables. We can also reverse the question: Can auxiliary variables
predict latent class membership?

For example, suppose we expect that the true population effect is zero,
but that when research on this effect was first introduced, early
studies all reported large positive effects due to publication bias
(and/or questionable research practices). In such a case, study
characteristic `z` might represent “year of publication”. Our theory
implies the hypothesis that `z` causes class membership.

To test this hypothesis, we can use
[`pseudo_class()`](https://cjvanlissa.github.io/tidySEM/reference/pseudo_class.html)
to repeatedly draw plausible class assignments from the posterior
probabilities, fit a logistic regression in each dataset, and pool the
regression estimates.

``` r
# Predict latent class membership from an auxiliary variable
set.seed(783)
pred <- pseudo_class(
  x = fit,
  model = glm(
    I(as.integer(class) == 2L) ~ z,
    data = data,
    family = binomial()
  ),
  data = data.frame(z = dat$z),
  m = 10
)
pred
#>          term   estimate std.error statistic       df      p.value
#> 1 (Intercept) -0.6089003 0.1921151 -3.169455 123.7412 1.924902e-03
#> 2           z  0.8992796 0.1891752  4.753687  86.3312 7.913785e-06
```

The regression estimates the association between `z` and the log odds of
belonging to class 2 rather than class 1. A negative coefficient
indicates that higher values of `z` are associated with greater odds of
class 1 membership. We see that this effect is significant; thus, it is
plausible that publication year explains the different classes.

The present model can use logistic regression because there are only two
classes. If there are more than two classes, multinomial logistic
regression would be more appropriate, e.g.:

    pred <- pseudo_class(
      x = fit,
      model = function(data) {
        nnet::multinom(
          factor(class) ~ z,
          data = data,
          trace = FALSE
        )
      },
      data = data.frame(z = dat$z),
      m = 10
    )

## Reporting auxiliary analyses

When reporting auxiliary analyses, describe the measurement model and
class enumeration, the auxiliary variable and its measurement scale, the
estimation method, and the way classification uncertainty was handled.
Report class-specific estimates with uncertainty intervals and specify
whether standard errors are valid (as in
[`one_step()`](https://cjvanlissa.github.io/tidySEM/reference/one_step.html)
and
[`pseudo_class()`](https://cjvanlissa.github.io/tidySEM/reference/pseudo_class.html),
but not
[`BCH()`](https://cjvanlissa.github.io/tidySEM/reference/BCH.html)), and
whether they are adjusted for first-stage estimation (as in
[`two_step()`](https://cjvanlissa.github.io/tidySEM/reference/two_step.html)).
