## More stable estimation of the Tweedie power p in mgcv::gam
## (a) multi-start: fit with several initial p values, keep lowest REML
## (b) profile: fit over a grid of *fixed* p, keep lowest REML
## Applied to the "hard" datasets from tw_family_reuse.R where fresh vs reused
## starts disagreed.

suppressPackageStartupMessages(library(mgcv))

hard_data <- function(i) {
  set.seed(5000 + i)
  n <- 50
  contvar <- runif(n)
  method  <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  eta <- -1 + sin(2 * pi * contvar) + c(A = 0, B = 0.5, C = -0.3)[as.character(method)]
  data.frame(response = rTweedie(exp(eta), p = runif(1, 1.1, 1.9), phi = 3),
             contvar, method)
}

fit <- function(d, fam) {
  m <- gam(response ~ s(contvar) + method, family = fam, data = d, method = "REML")
  c(p = m$family$getTheta(TRUE), reml = unname(m$gcv.ubre),
    pv = summary(m)$pTerms.pv[["method"]])
}

## (a) negative theta = starting value for p (estimated)
multistart <- function(d, starts = seq(1.05, 1.95, by = 0.1)) {
  r <- t(sapply(starts, function(s) fit(d, tw(theta = -s))))
  cbind(start = starts, r)
}

## (b) positive theta = p held fixed
profile_p <- function(d, grid = seq(1.02, 1.98, by = 0.02)) {
  r <- t(sapply(grid, function(p) fit(d, tw(theta = p))))
  cbind(p_fixed = grid, r)
}

for (i in c(83, 48, 118, 20, 122)) {
  d <- hard_data(i)
  ms <- multistart(d)
  pr <- profile_p(d)
  best <- pr[which.min(pr[, "reml"]), ]
  ## p values within 2 REML units of the best: a rough "plausible range" for p
  ok <- pr[pr[, "reml"] - min(pr[, "reml"]) < 2, , drop = FALSE]
  cat(sprintf("\n=== dataset %d (%.0f%% zeros) ===\n", i, 100 * mean(d$response == 0)))
  cat("Multi-start results (estimated p, REML, method p-value):\n")
  print(round(ms, 4))
  cat(sprintf("Profile minimum: p = %.2f, REML = %.3f, method p-value = %.4f\n",
              best["p_fixed"], best["reml"], best["pv"]))
  cat(sprintf("p within 2 REML units of min: %.2f to %.2f; method p-value over that range: %.3f to %.3f\n",
              min(ok[, "p_fixed"]), max(ok[, "p_fixed"]), min(ok[, "pv"]), max(ok[, "pv"])))
}
