# CondIntCal

This repository provides `R` functions to assess the conditional calibration of interval forecasts (i.e. prediction intervals) via decompositions of the interval score.

### Features

- Two theoretically-appealing methods to calculate the decomposition terms
- Additional functions to assess the average interval length and unconditional (i.e. marginal) coverage 
- Miscalibration-Discrimination plots to visualise the decomposition terms

---

## Interval Forecasts and the Interval Score

Interval forecasts (or prediction intervals) are comprised of a lower bound $\ell$ and an upper bound $u$, with $\ell < u$. The forecast is made such that the observation $y$ is predicted to fall within the interval with a given coverage level $1 - \alpha$. Typically, interval forecasts are _central_ prediction intervals, where it is assumed that $\mathrm{P}(y < l) = \mathrm{P}(y > u) = \alpha/2$.

Competing interval forecasts can be assessed and compared using the _interval score_ (also sometimes called the _Winkler score_),

$$\mathrm{IS}_{\alpha}([\ell, u], y) = |u - \ell| + \frac{2}{\alpha} 1 [ y < \ell ] (\ell - y) + \frac{2}{\alpha} 1 [ y > u ] (y - u), $$

where $1[\cdotp]$ is the indicator function.

In practice, we observe several interval forecasts $[\ell_i, u_i]$ and corresponding observations $y_i$, for $i = 1, \dots, n$, and we wish to compare forecasters or forecast methods based on the average interval score

$$ \mathrm{\bar{S}} = \frac{1}{n} \sum_{i=1}^{n} \mathrm{IS}_{\alpha}([\ell_i, u_i], y_i). $$

While the average score $\mathrm{\bar{S}}$ provides a single value that can be used to rank different forecasters, it can also be additively decomposed into terms quantifying different aspects of forecast performance. Most commonly, these terms measure how much _uncertainty_ there is in the forecast problem, the forecast's _discrimination_ ability or information content, and the degree of forecast _miscalibration_. These terms can be calculated using the following decomposition:

$$ \mathrm{\bar{S}} = \mathrm{\bar{S}^R} - (\mathrm{\bar{S}^R} - \mathrm{\bar{S}^C}) + (\mathrm{\bar{S}} - \mathrm{\bar{S}^C}), $$

where

$$ \mathrm{\bar{S}^R} = \frac{1}{n} \sum_{i=1}^{n} \mathrm{IS}_{\alpha}(R, y_i) $$

is the average interval score for a reference prediction interval $R$, and

$$ \mathrm{\bar{S}^C} = \frac{1}{n} \sum_{i=1}^{n} \mathrm{IS}_{\alpha}(C_i, y_i) $$

is the average interval score for predictions $C_i$, which correspond to recalibrated versions of the original prediction intervals $[\ell_i, u_i].$

The first term of the decomposition, $\mathrm{\bar{S}^R}$, provides a baseline measure of forecast performance, thereby quantifying the inherent uncertainty or unpredictability of the observations; the second term, $\mathrm{\bar{S}^R} - \mathrm{\bar{S}^C}$, represents the improvement in accuracy that is obtained from the recalibrated predictions compared to the (uninformative) reference predictions, thus providing a measure of how much information is contained within the original interval forecasts; the third term, $\mathrm{\bar{S}} - \mathrm{\bar{S}^C}$, is the improvement in accuracy that is obtained by recalibrating the original interval forecasts, which quantifies the degree of miscalibration in the original predictions. Note that the miscalibration corresponds to conditional calibration, where the conditioning is on the original forecasts themselves.

This decomposition holds for any choice of the reference prediction $R$ and recalibration method to produce $C_i$. However, the interpretation of the terms requires specific choices. In particular, it is desirable that the decomposition terms are non-negative; in this case, we can say that a miscalibration of zero corresponds to a calibrated prediction interval, for example. One method that has been proposed in the literature to calculate the decomposition terms with desirable theoretical guarantees corresponds to recalibrating the original prediction intervals using _isotonic distributional regression_ (see references below).

This repository provides the functionality to perform these theoretically-appealing interval score decompositions, facilitating an assessment of the conditional calibration of interval forecasts in practice.

---

## Example

```r
set.seed(744)

n <- 1000 # sample size
alpha <- 0.1 # 90% prediction intervals
mu <- rnorm(n)
y <- rnorm(n, mean = mu, sd = 1) # simulate observations

# Ideal forecaster: F = N(mu, 1)
int_id <- data.frame(Lower = qnorm(alpha/2, mu), Upper = qnorm(1 - alpha/2, mu))

out_iso <- is_decomp(y, int_id, level = 1 - alpha)                     # get isotonic decomposition
out_lin <- is_decomp(y, int_id, level = 1 - alpha, method = "linear")  # get linear decomposition

print(round(out_iso, 3)) # print the isotonic decomposition terms
#   IS   UNC   DSC   MCB 
# 4.130 5.882 2.037 0.286 
 
print(round(out_lin, 3)) # print the linear decomposition terms
#   IS   UNC   DSC   MCB 
# 4.130 5.882 1.757 0.005

plot_mcbdsc(list(Linear = out_lin, Isotonic = out_iso), MCB_lim = 2.5) # plot the decomposition terms
```

---

## Installation and development

This package has not been submitted to `CRAN`, but can be installed in `R` using `devtools`
```r
# install.packages("devtools")
devtools::install_github("sallen12/CondIntCal")
```
The package is still in active development. Comments, suggestions, and input are more than welcome.

---

## Citation

The interval score decomposition based on isotonic regression is detailed in the following arXiv preprint:

```
@article{AllenEtAl2025,
  title={Assessing the conditional calibration of interval forecasts using decompositions of the interval score},
  author={Allen, Sam and Burnello, Julia and Ziegel, Johanna},
  journal={arXiv preprint arXiv:2508.18034},
  year={2025}
}
```


