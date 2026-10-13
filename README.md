# mgcv `tw()` family reuse

Why the p-value for a parametric term in `mgcv::gam()` changes depending on
whether `family = tw()` is written inline or passed as a pre-built object
(`family_val <- tw()`), and how to get a more stable Tweedie fit.

Tested with mgcv 1.9-3, R 4.5.1.

## Cause

`tw()` stores the Tweedie power parameter (`.Theta`) in an environment, and
`gam()` overwrites it in place via `family$putTheta()`. Environments are
reference objects, so a reused family object carries the previous fit's
estimate and the next fit starts its search from there. Its first use matches
an inline `tw()` exactly; later uses do not.

## Scripts

- `R/tw_family_reuse.R` reproduces the effect and runs a small simulation:
  - n = 200, well-behaved: 25/100 fits differ, max difference in the `method` p-value about 2e-6.
  - n = 50, 60–96% zeros: 58/150 differ, 4 by more than 0.01, 1 changes significance at 0.05.
- `R/tw_stable_p.R` compares multi-start (`tw(theta = -start)`) with profiling
  REML over a grid of fixed Tweedie powers (`tw(theta = p)`) on the hard datasets.
  In 2 of 4 datasets, none of the 10 starts reached the profile minimum. In one,
  the `method` p-value ranged from 0.05 to 1.0 across powers within 2 REML units
  of the best.

## Recommendations

1. Write `family = tw()` inline in every call. Don't reuse a family object.
2. If the data are mostly zeros, consider a hurdle model: a binomial GAM for
   presence plus a Gamma GAM for the positive values.
3. Otherwise, profile the Tweedie power on a grid (`profile_p()` in
   `R/tw_stable_p.R`), fit at the minimum, and report how sensitive the
   `method` test is across the plausible range.

The Tweedie power "p" is unrelated to the significance p-value.

## Status

- Not a documented issue as far as we can find. mgcv has no GitHub issue
  tracker, so reports go by email to the maintainer (Simon Wood).
- Waiting on the student's data (proportion of zeros, `.RData` + model call)
  before choosing the hurdle or profiling approach. Data files are gitignored.
