## Reproduce: p-value for a parametric term in mgcv::gam differs depending on
## whether family = tw() is created inline or passed as a pre-built object.
##
## Hypothesis: tw() stores the Tweedie power parameter in an environment that
## gam() updates in place (family$putTheta). R environments are reference
## objects, so a family object reused across fits remembers the previous fit's
## estimate and the next fit starts its optimisation from there.

suppressPackageStartupMessages(library(mgcv))
cat("mgcv", as.character(packageVersion("mgcv")), "\n\n")

sim_data <- function(n = 200, seed = 1, p = 1.6, phi = 2) {
  set.seed(seed)
  contvar <- runif(n)
  method  <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  eta <- 1 + sin(2 * pi * contvar) + c(A = 0, B = 0.2, C = -0.1)[as.character(method)]
  response <- mgcv::rTweedie(exp(eta), p = p, phi = phi)
  data.frame(response, contvar, method)
}

fit_one <- function(dat, fam) {
  m <- gam(response ~ s(contvar) + method, family = fam, data = dat, method = "REML")
  s <- summary(m)
  c(p_tweedie = m$family$getTheta(TRUE),
    pval_methodB = s$p.table["methodB", "Pr(>|t|)"],
    pval_methodC = s$p.table["methodC", "Pr(>|t|)"],
    method_anova = s$pTerms.pv[["method"]],
    reml = unname(m$gcv.ubre))
}

## ---- 1. Minimal demonstration ---------------------------------------------
dat1 <- sim_data(seed = 1)
dat2 <- sim_data(seed = 2, p = 1.2)

family_val <- tw()
cat("Stored theta in family_val before any fit:", family_val$getTheta(), "\n")

r_inline_1 <- fit_one(dat1, tw())
r_object_1 <- fit_one(dat1, family_val)       # first use: identical to inline
cat("Stored theta in family_val after fit 1:   ", family_val$getTheta(), "\n")

fit_one(dat2, family_val)                       # fit some other model/data
cat("Stored theta in family_val after fit on other data:", family_val$getTheta(), "\n\n")

r_object_2 <- fit_one(dat1, family_val)       # same data as fit 1 again
r_inline_2 <- fit_one(dat1, tw())

print(rbind(inline_first  = r_inline_1,
            object_first  = r_object_1,
            inline_again  = r_inline_2,
            object_reused = r_object_2), digits = 10)

## ---- 2. Small simulation study --------------------------------------------
## For each simulated dataset, fit (a) with a fresh tw() and (b) with one
## family object created once outside the loop (as in a typical script).
nsim <- 100
family_val <- tw()
res <- do.call(rbind, lapply(seq_len(nsim), function(i) {
  d <- sim_data(seed = 1000 + i, p = runif(1, 1.1, 1.9))
  fresh  <- fit_one(d, tw())
  reused <- fit_one(d, family_val)
  data.frame(sim = i,
             p_fresh = fresh["p_tweedie"], p_reused = reused["p_tweedie"],
             pv_fresh = fresh["method_anova"], pv_reused = reused["method_anova"],
             reml_fresh = fresh["reml"], reml_reused = reused["reml"])
}))
rownames(res) <- NULL

res$abs_diff_pv <- abs(res$pv_fresh - res$pv_reused)
res$rel_diff_pv <- res$abs_diff_pv / res$pv_fresh
res$abs_diff_p  <- abs(res$p_fresh - res$p_reused)

cat("\nSimulation (", nsim, " datasets)\n", sep = "")
cat("Fits whose method p-value differs (rel diff > 1e-6):",
    sum(res$rel_diff_pv > 1e-6), "of", nsim, "\n")
cat("Fits that flip significance at 0.05:",
    sum((res$pv_fresh < 0.05) != (res$pv_reused < 0.05)), "\n")
cat("Max |diff| in method p-value:", signif(max(res$abs_diff_pv), 3), "\n")
cat("Max |diff| in estimated Tweedie p:", signif(max(res$abs_diff_p), 3), "\n")
cat("Max |diff| in REML score:",
    signif(max(abs(res$reml_fresh - res$reml_reused)), 3), "\n\n")
print(head(res[order(-res$rel_diff_pv), ], 8), digits = 6)

## ---- 2b. Harder data: small n, many zeros (flat / multimodal REML in p) ----
## Lower REML score (gcv.ubre) is better; note neither start wins consistently.
family_val <- tw()
hard <- do.call(rbind, lapply(seq_len(150), function(i) {
  set.seed(5000 + i)
  n <- 50
  contvar <- runif(n)
  method  <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  eta <- -1 + sin(2 * pi * contvar) + c(A = 0, B = 0.5, C = -0.3)[as.character(method)]
  d <- data.frame(response = rTweedie(exp(eta), p = runif(1, 1.1, 1.9), phi = 3),
                  contvar, method)
  fresh  <- try(fit_one(d, tw()), silent = TRUE)
  reused <- try(fit_one(d, family_val), silent = TRUE)
  if (inherits(fresh, "try-error") || inherits(reused, "try-error")) return(NULL)
  data.frame(sim = i, prop_zero = mean(d$response == 0),
             p_fresh = fresh["p_tweedie"], p_reused = reused["p_tweedie"],
             pv_fresh = fresh["method_anova"], pv_reused = reused["method_anova"],
             reml_fresh = fresh["reml"], reml_reused = reused["reml"])
}))
rownames(hard) <- NULL
hard$abs_diff_pv <- abs(hard$pv_fresh - hard$pv_reused)
cat("\nHard scenario (n = 50, mostly zeros):", nrow(hard), "fits\n")
cat("Method p-value differs (rel diff > 1e-6):",
    sum(hard$abs_diff_pv / hard$pv_fresh > 1e-6), "\n")
cat("Differs by > 0.01:", sum(hard$abs_diff_pv > 0.01), "\n")
cat("Flips significance at 0.05:",
    sum((hard$pv_fresh < 0.05) != (hard$pv_reused < 0.05)), "\n")
print(head(hard[order(-hard$abs_diff_pv), ], 6), digits = 5)

## ---- 3. Fix: a fresh family per fit gives order-independent results -------
make_fam <- function() tw()
a <- fit_one(dat1, make_fam()); fit_one(dat2, make_fam()); b <- fit_one(dat1, make_fam())
cat("\nFresh family each time, identical results:", isTRUE(all.equal(a, b, tolerance = 0)), "\n")
