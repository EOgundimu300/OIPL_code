## Deterministic checks of Table S13, the intercept-only remainder in Theorem 2 and the
## observed-information Wald interval in Section 5.4. Run from the directory containing
## obsinfo_core.R and Endometrial.dat: Rscript final_checks.R
source("obsinfo_core.R")
a <- 0.5

## Table S13: exact two-group means and first-order approximations
fo1 <- function(link, eta, m) {
  L <- make_link(link); b <- L$base(eta)
  L$tau(b) / (2 * m * exp(L$lomega(b)$l))
}
dm <- function(link, eta, m)
  exact_bias_intercept(link, eta, m, "O")$bias - exact_bias_intercept(link, eta, m, "E")$bias
for (lk in c("probit", "cloglog")) for (st in list(c(-0.5, 1), c(-1.5, 2))) for (n in c(25, 50)) {
  n1 <- floor(0.3 * n + 0.5); n0 <- n - n1
  mn <- dm(lk, sum(st), n1) - dm(lk, st[1], n0)
  f  <- fo1(lk, sum(st), n1) - fo1(lk, st[1], n0)
  cat(sprintf("%-7s (%g, %g) n=%d: mean %.4f, first order %.4f, ratio %.2f\n",
              lk, st[1], st[2], n, mn, f, mn / f))
}

## Theorem 2 in intercept-only models: n^(3/2) E|R_n|, R_n = estimate - eta0 - Delta
remainder <- function(link, eta0, n, type) {
  L <- make_link(link); b <- L$base(eta0)
  p1 <- L$psi1(b, 1L); p0 <- L$psi0(b, 1L); G <- L$prob(eta0)
  i <- n * exp(L$lomega(b)$l)
  ex <- exact_bias_intercept(link, eta0, n, type)
  t <- 0:n
  S <- t * b$r1 - (n - t) * b$r0
  j <- t * p1$v + (n - t) * p0$v
  Tt <- -(t * p1$d1 + (n - t) * p0$d1)          # second derivative of the score
  Tbar <- -n * (G * p1$d1 + (1 - G) * p0$d1)
  B <- if (type == "E") a * L$lomega(b, p1$v, p0$v, 1L)$d1 else a * (-Tt) / j
  s <- S / i
  Delta <- s + B / i - (j - i) * s / i + 0.5 * Tbar * s^2 / i
  n^1.5 * sum(ex$weights * abs(ex$estimates - eta0 - Delta))
}
for (lk in list(c("probit", 0.5), c("cloglog", -0.5), c("cloglog", 0.5)))
  for (ty in c("E", "O"))
    cat(lk[1], lk[2], ty, round(sapply(c(25, 100, 400, 1600), function(n)
      remainder(lk[1], as.numeric(lk[2]), n, ty)), 3), "\n")

## Table 6 and the observed-information Wald interval for beta_v
dat <- read.table("Endometrial.dat", header = TRUE)
fit <- obsinfo_fit(HG ~ NV + PI + EH, data = dat, link = "probit", method = "O")
gam <- fit$coefficients * fit$scale
inf <- bp_info(fit$model, gam)
se_j <- sqrt(diag(solve(inf$observed))) / fit$scale
print(round(fit$coefficients[["NV"]] + c(-1, 1) * qnorm(0.975) * se_j[["NV"]], 3))
