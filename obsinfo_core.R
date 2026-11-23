## Binomial-regression likelihoods, determinant penalties and adjusted scores.
## The expected-information penalty uses i(beta). The observed-information
## penalty uses j(beta) on its positive-definite domain. The search returns
## the highest accepted local maximum among its starting values.
## Load with source("obsinfo_core.R"). Only standard R packages are required.
##
## Quick start
##   fit <- obsinfo_fit(y ~ x1 + x2, data = d, link = "probit", method = "O")
##   fit; summary(fit)               # estimates, Wald tests, diagnostics
##   confint(fit)                    # penalized-profile confidence sets
##   confint(fit, type = "wald")     # Wald intervals
##   predict(fit, newdata, type = "response")
## The response is binary (0/1, logical or a two-level factor) or
## cbind(successes, failures). method: "O" (observed-information penalty,
## default), "E" (expected-information penalty), "A" (Firth's adjusted score),
## "ML" (maximum likelihood; only when the data overlap). See Section 10.
##
## Sections
##   1. Numerical utilities
##   2. Links: logit, probit, cloglog, loglog, cauchit
##   3. Model, penalized criteria and Firth's adjusted score
##   4. Optimization: trust region, quasi-Fisher scoring, Newton refinement
##   5. Protocol design, starts, multistart fitting, selection, inference
##   6. First-order biases, control variate, exact intercept-only means
##   7. Separation by linear programming
##   8. Penalized-profile confidence sets
##   9. Independent re-evaluation by singular value decompositions
##  10. Data analysis interface: obsinfo_fit() and its methods
##
## Main functions
##   obsinfo_fit(formula, data, link, method)  fit from a formula; methods
##                                        print, summary, coef, vcov, predict,
##                                        confint
##   bp_model(X, y, m, link)              data and inverse link
##   bp_eval(model, beta, type, a, level) value, gradient (level 1) and Hessian
##                                        (level 2) of a criterion:
##       "ML" l;  "E" l + a log det i;  "O" l + a log det j (on j > 0);
##       "Ri", "Rj" l - tr(j^{-1} e)/2 + a log det i or j, with a = 1/N
##   bp_adjscore(model, beta, jacobian)   Firth's S + A and its Jacobian
##   fit_criterion_protocol()             protocol fit of a criterion; $history
##                                        holds every candidate of every stage
##   fit_criterion_search()               the same plus, under separation, starts
##                                        along the separating direction
##   fit_adjscore_protocol()              protocol solution of S + A = 0
##   fit_criterion(), fit_adjscore()      multistart fits without rescue
##   profile_stat(), profile_set()        profile statistic; confidence set
##   wald_ci(), stable_reevaluation()     Wald intervals; SVD re-evaluation
##   separation_check()                   complete, quasi-complete or overlap
##   first_order_bias()                   b_ML, b^E_a, b^O_a, displacement
##   exact_bias_intercept()               exact means, intercept-only models
##   protocol_design(), protocol_truth, protocol_starts()
##   logdet_cauchy_binet()                independent log-determinant check


## =============================================================================
## 1. Numerical utilities
## =============================================================================
## Numerical utilities used throughout. Base R only.

## log(1 - exp(-x)) for x >= 0. The expm1 branch is used for x <= log 2 and
## the log1p branch otherwise, so that both small and large x are accurate.
log1mexp <- function(x) {
  out <- numeric(length(x))
  s <- x <= log(2)
  out[s] <- log(-expm1(-x[s]))
  out[!s] <- log1p(-exp(-x[!s]))
  out
}

## Product c * v with the convention 0 * v = 0. Used to combine contributions
## of response categories: a category absent from the data contributes
## nothing, whatever the value of its (possibly infinite) term.
cmul <- function(c, v) {
  r <- c * v
  r[c == 0] <- 0
  r
}

## Row-wise Kronecker product: row i is A[i, ] %x% B[i, ].
rowkron <- function(A, B) {
  pa <- ncol(A)
  pb <- ncol(B)
  A[, rep(seq_len(pa), each = pb), drop = FALSE] *
    B[, rep(seq_len(pb), times = pa), drop = FALSE]
}

## Cholesky factorization of a symmetric matrix M after symmetric diagonal
## equilibration: S M S = R'R with S = diag(diag(M)^(-1/2)). Returns NULL when
## M has non-finite entries, a non-positive diagonal entry, or the Cholesky
## factorization fails (M not numerically positive definite). Equilibration
## changes neither positive definiteness nor the determinant, which is
## recovered exactly: log det M = log det(R'R) - 2 sum(log diag S).
chol_equil <- function(M) {
  if (!all(is.finite(M))) return(NULL)
  d <- diag(M)
  if (any(d <= 0)) return(NULL)
  s <- 1 / sqrt(d)
  R <- tryCatch(chol(M * tcrossprod(s)), error = function(e) NULL)
  if (is.null(R)) return(NULL)
  list(R = R, scale = s, logdet = 2 * sum(log(diag(R))) - 2 * sum(log(s)))
}

## Weighted cross-product K = X' diag(wt) X with its equilibrated Cholesky
## factor. With need_Z = TRUE also returns Zt = (X K^{-1/2})' in the sense
## Z Z' = X K^{-1} X', and the leverages q_i = x_i' K^{-1} x_i.
wcp_factor <- function(X, wt, need_Z = TRUE) {
  K <- crossprod(X, wt * X)
  f <- chol_equil(K)
  if (is.null(f)) return(NULL)
  f$K <- K
  if (need_Z) {
    f$Zt <- backsolve(f$R, t(X) * f$scale, transpose = TRUE)
    f$q <- colSums(f$Zt^2)
  }
  f
}

## Gradient and Hessian of log det K(beta), K = sum_i k_i(eta_i) x_i x_i',
## given the factor f of K (from wcp_factor), c1 = dk_i/deta_i and
## c2 = d^2 k_i/deta_i^2:
##   grad = X' (c1 * q),
##   hess = X' diag(c2 * q) X - X' diag(c1) (Q o Q) diag(c1) X,
## with Q = X K^{-1} X' and (Q o Q) = W W', W = rowkron(Z, Z).
logdet_derivs <- function(X, f, c1, c2 = NULL) {
  out <- list(grad = drop(crossprod(X, c1 * f$q)), hess = NULL)
  if (!is.null(c2)) {
    Z <- t(f$Zt)
    G <- crossprod(rowkron(Z, Z), c1 * X)
    out$hess <- crossprod(X, (c2 * f$q) * X) - crossprod(G)
  }
  out
}

## Independent check of log det(X' diag(w) X) by the Cauchy-Binet formula
## (Supplementary Section D.3): the determinant is the sum over p-subsets S of rows
## of det(X_S)^2 prod_{i in S} w_i, accumulated on the log scale with signs, so
## that signed (cauchit) weights are allowed. Returns NA if the sum is not
## positive. A positive sum does not establish positive definiteness when some
## weights are negative; check that separately (chol_equil). choose(n, p)
## subsets: for small designs only.
logdet_cauchy_binet <- function(X, w) {
  X <- as.matrix(X)
  S <- utils::combn(nrow(X), ncol(X))
  ld <- numeric(ncol(S)); sg <- numeric(ncol(S))
  for (k in seq_len(ncol(S))) {
    idx <- S[, k]
    ld[k] <- 2 * as.numeric(determinant(X[idx, , drop = FALSE], logarithm = TRUE)$modulus) +
      sum(log(abs(w[idx])))
    sg[k] <- prod(sign(w[idx]))
  }
  keep <- is.finite(ld)
  if (!any(keep)) return(NA_real_)
  mx <- max(ld[keep])
  tot <- sum(sg[keep] * exp(ld[keep] - mx))
  if (!isTRUE(tot > 0)) return(NA_real_)
  mx + log(tot)
}

## Symmetric inverse square root of a symmetric positive definite matrix;
## NULL if the matrix is not finite or not positive definite.
isqrt_sym <- function(M) {
  if (!all(is.finite(M))) return(NULL)
  e <- eigen(M, symmetric = TRUE)
  if (any(e$values <= 0)) return(NULL)
  e$vectors %*% (t(e$vectors) / sqrt(e$values))
}

symmetrize <- function(M) (M + t(M)) / 2

## Solve M x = v and evaluate v' M^{-1} v from the equilibrated Cholesky
## factor f of M (chol_equil or wcp_factor): M^{-1} = S (R'R)^{-1} S.
chol_solve <- function(f, v) {
  f$scale * backsolve(f$R, backsolve(f$R, f$scale * v, transpose = TRUE))
}
chol_quadform <- function(f, v) {
  sum(backsolve(f$R, f$scale * v, transpose = TRUE)^2)
}


## =============================================================================
## 2. Links: logit, probit, cloglog, loglog, cauchit
## =============================================================================
## Inverse links G for binomial-response GLMs and the per-observation
## quantities used by the likelihood, the two information matrices and
## their derivatives. Notation follows the paper:
##   r1 = g/G, r0 = g/Gbar,
##   psi1 = -(log G)'', psi0 = -(log Gbar)''   (observed weights),
##   omega = g^2/(G Gbar)                       (expected weight),
##   tau = d/deta log{g/(G Gbar)}.
## Primes denote derivatives with respect to eta.
##
## Each link is a list with
##   base(eta)        list: eta, lG = log G, lGb = log Gbar, lg = log g,
##                    d1, d2, d3 = first three derivatives of log g, r1, r0
##                    (plus link-specific auxiliaries);
##   psi1(b, deriv)   list(v, d1, d2): psi1 and its first two derivatives;
##   psi0(b, deriv)   the same for psi0;
##   lomega(b, p1v, p0v, deriv)
##                    list(l, d1, d2): log omega and its first two
##                    derivatives (d2 requires psi1 and psi0 values);
##   tau(b)           tau;
##   prob(eta)        G(eta), for data generation;
##   glm_family       quasibinomial family for the adjusted-response start
##                    (NULL for loglog).
## Derivatives are computed for deriv = 0, 1, 2; unrequested entries are NULL.
##
## Generic identities (all links), with l_k the k-th derivative of log g:
##   psi1   = r1 (r1 - l1)
##   psi1'  = -psi1 (2 r1 - l1) - r1 l2
##   psi1'' = -psi1' (2 r1 - l1) + 2 psi1^2 + 2 psi1 l2 - r1 l3
##   psi0   = r0 (r0 + l1)
##   psi0'  = psi0 (2 r0 + l1) + r0 l2
##   psi0'' = psi0' (2 r0 + l1) + 2 psi0^2 + 2 psi0 l2 + r0 l3
##   (log omega)'  = 2 l1 - r1 + r0,   (log omega)'' = 2 l2 + psi1 + psi0,
##   tau = l1 - r1 + r0.
## Links override these where cancellation would otherwise lose accuracy.

psi1_generic <- function(b, deriv = 0L) {
  k <- 2 * b$r1 - b$d1
  v <- b$r1 * (b$r1 - b$d1)
  d1 <- d2 <- NULL
  if (deriv >= 1L) d1 <- -v * k - b$r1 * b$d2
  if (deriv >= 2L) d2 <- -d1 * k + 2 * v * v + 2 * v * b$d2 - b$r1 * b$d3
  list(v = v, d1 = d1, d2 = d2)
}

psi0_generic <- function(b, deriv = 0L) {
  k <- 2 * b$r0 + b$d1
  v <- b$r0 * (b$r0 + b$d1)
  d1 <- d2 <- NULL
  if (deriv >= 1L) d1 <- v * k + b$r0 * b$d2
  if (deriv >= 2L) d2 <- d1 * k + 2 * v * v + 2 * v * b$d2 + b$r0 * b$d3
  list(v = v, d1 = d1, d2 = d2)
}

lomega_generic <- function(b, p1v = NULL, p0v = NULL, deriv = 0L) {
  list(l = 2 * b$lg - b$lG - b$lGb,
       d1 = if (deriv >= 1L) 2 * b$d1 - b$r1 + b$r0 else NULL,
       d2 = if (deriv >= 2L) 2 * b$d2 + p1v + p0v else NULL)
}

tau_generic <- function(b) b$d1 - b$r1 + b$r0

## ---------------------------------------------------------------------------
## Logit: psi1 = psi0 = omega = g = G Gbar, tau = 0.
link_logit <- function() {
  base <- function(eta) {
    lG <- plogis(eta, log.p = TRUE)
    lGb <- plogis(eta, lower.tail = FALSE, log.p = TRUE)
    lg <- lG + lGb
    g <- exp(lg)
    d1 <- -tanh(eta / 2)                      # 1 - 2G
    list(eta = eta, lG = lG, lGb = lGb, lg = lg, g = g,
         d1 = d1, d2 = -2 * g, d3 = -2 * g * d1,
         r1 = exp(lGb), r0 = exp(lG))
  }
  psi <- function(b, deriv = 0L) {
    list(v = b$g,
         d1 = if (deriv >= 1L) b$g * b$d1 else NULL,
         d2 = if (deriv >= 2L) b$g * (b$d1 * b$d1 + b$d2) else NULL)
  }
  list(name = "logit", base = base, psi1 = psi, psi0 = psi,
       lomega = lomega_generic,
       tau = function(b) numeric(length(b$eta)),
       prob = function(eta) plogis(eta),
       glm_family = function() quasibinomial(link = "logit"))
}

## ---------------------------------------------------------------------------
## Probit. With u the signed predictor and v = -u > 0 in the lower tail,
## nu(v) = phi(v)/{1 - Phi(v)} - v = 1/(v + 2/(v + 3/(v + ...))) (Laplace's
## continued fraction). By Lemma 1, psi(u) = 1 - var(Z | Z > v) and, since the
## excess Z - v given Z > v forms an exponential family in -v with mean nu(v),
## psi(u) = 1 + nu'(v),  psi'(u) = -nu''(v),  psi''(u) = nu'''(v),
## and r(u) = phi(u)/Phi(u) = v + nu(v). For u < -3 these are evaluated by a
## backward recursion T_k = k/(v + T_{k+1}) differentiated three times in v,
## started at the fixed point of T = (K + 1)/(v + T); otherwise the direct
## formulas psi = r(u + r), psi' = r - psi(u + 2r), psi'' = -psi'(u + 2r) +
## 2 psi(psi - 1) are used. Against 1200-bit references the continued-fraction
## branch is accurate to about 1e-14 for u <= -3; the direct branch loses
## accuracy in psi' and psi'' as u decreases (relative error 7e-11 at u = -4,
## 1.3e-10 at u = -4.5), which fixes the switch at u = -3. The two branches
## agree to within 7.2e-11 on 3 < |eta| <= 4.

norm_cf_group <- function(v, K) {
  k1 <- K + 1
  s <- sqrt(v * v + 4 * k1)
  T0 <- (s - v) / 2
  T1 <- (v / s - 1) / 2
  T2 <- 2 * k1 / s^3
  T3 <- -6 * k1 * v / s^5
  for (k in K:1) {
    a <- v + T0
    ra <- (1 + T1) / a
    f <- k / a
    f1 <- -f * ra
    f2 <- f * (2 * ra * ra - T2 / a)
    f3 <- f * (-T3 / a + 6 * ra * T2 / a - 6 * ra^3)
    T0 <- f; T1 <- f1; T2 <- f2; T3 <- f3
  }
  list(T0 = T0, T1 = T1, T2 = T2, T3 = T3)
}

## nu(v) and its first three derivatives for v >= 3; the number of terms
## depends on v (60, 30 and 15 terms on [3, 5), [5, 10) and [10, Inf)).
norm_cf <- function(v) {
  n <- length(v)
  out <- list(T0 = numeric(n), T1 = numeric(n), T2 = numeric(n), T3 = numeric(n))
  grp <- list(v >= 10, v >= 5 & v < 10, v < 5)
  Ks <- c(15L, 30L, 60L)
  for (k in 1:3) {
    s <- grp[[k]]
    if (any(s)) {
      r <- norm_cf_group(v[s], Ks[k])
      out$T0[s] <- r$T0; out$T1[s] <- r$T1; out$T2[s] <- r$T2; out$T3[s] <- r$T3
    }
  }
  out
}

link_probit <- function(switch_at = 3) {
  base <- function(eta) {
    n <- length(eta)
    lG <- pnorm(eta, log.p = TRUE)
    lGb <- pnorm(eta, lower.tail = FALSE, log.p = TRUE)
    lg <- dnorm(eta, log = TRUE)
    lo <- eta < -switch_at                    # lower tail of psi1, r1
    hi <- eta > switch_at                     # lower tail of psi0, r0
    r1 <- r0 <- numeric(n)
    r1[!lo] <- exp(lg[!lo] - lG[!lo])
    r0[!hi] <- exp(lg[!hi] - lGb[!hi])
    cf1 <- cf0 <- NULL
    if (any(lo)) {
      cf1 <- norm_cf(-eta[lo])
      r1[lo] <- -eta[lo] + cf1$T0
    }
    if (any(hi)) {
      cf0 <- norm_cf(eta[hi])
      r0[hi] <- eta[hi] + cf0$T0
    }
    list(eta = eta, lG = lG, lGb = lGb, lg = lg,
         d1 = -eta, d2 = rep(-1, n), d3 = numeric(n),
         r1 = r1, r0 = r0, lo = lo, hi = hi, cf1 = cf1, cf0 = cf0)
  }
  ## psi1(eta) = psi(eta)
  psi1 <- function(b, deriv = 0L) {
    n <- length(b$eta)
    v <- numeric(n); d1 <- if (deriv >= 1L) numeric(n); d2 <- if (deriv >= 2L) numeric(n)
    ok <- !b$lo
    if (any(ok)) {
      u <- b$eta[ok]; r <- b$r1[ok]; k <- u + 2 * r
      pv <- r * (u + r)
      v[ok] <- pv
      if (deriv >= 1L) {
        pd1 <- r - pv * k
        d1[ok] <- pd1
        if (deriv >= 2L) d2[ok] <- -pd1 * k + 2 * pv * (pv - 1)
      }
    }
    if (any(b$lo)) {
      v[b$lo] <- 1 + b$cf1$T1
      if (deriv >= 1L) d1[b$lo] <- -b$cf1$T2
      if (deriv >= 2L) d2[b$lo] <- b$cf1$T3
    }
    list(v = v, d1 = d1, d2 = d2)
  }
  ## psi0(eta) = psi(-eta); psi0'(eta) = -psi'(-eta); psi0''(eta) = psi''(-eta)
  psi0 <- function(b, deriv = 0L) {
    n <- length(b$eta)
    v <- numeric(n); d1 <- if (deriv >= 1L) numeric(n); d2 <- if (deriv >= 2L) numeric(n)
    ok <- !b$hi
    if (any(ok)) {
      u <- -b$eta[ok]; r <- b$r0[ok]; k <- u + 2 * r
      pv <- r * (u + r)
      v[ok] <- pv
      if (deriv >= 1L) {
        pd1 <- r - pv * k
        d1[ok] <- -pd1
        if (deriv >= 2L) d2[ok] <- -pd1 * k + 2 * pv * (pv - 1)
      }
    }
    if (any(b$hi)) {
      v[b$hi] <- 1 + b$cf0$T1
      if (deriv >= 1L) d1[b$hi] <- b$cf0$T2
      if (deriv >= 2L) d2[b$hi] <- b$cf0$T3
    }
    list(v = v, d1 = d1, d2 = d2)
  }
  ## tau = -eta - r1 + r0; in the tails r0 - eta = nu(eta) or r1 + eta = nu(-eta)
  tau <- function(b) {
    t <- -b$eta - b$r1 + b$r0
    if (any(b$hi)) t[b$hi] <- b$cf0$T0 - b$r1[b$hi]
    if (any(b$lo)) t[b$lo] <- b$r0[b$lo] - b$cf1$T0
    t
  }
  list(name = "probit", base = base, psi1 = psi1, psi0 = psi0,
       lomega = lomega_generic, tau = tau,
       prob = function(eta) pnorm(eta),
       glm_family = function() quasibinomial(link = "probit"))
}

## ---------------------------------------------------------------------------
## Complementary log-log, G = 1 - exp(-z), z = e^eta.
##   log Gbar = -z, l1 = 1 - z, l2 = l3 = -z, r1 = z/(e^z - 1), r0 = z,
##   psi0 = psi0' = psi0'' = z,
##   psi1 = r1 A with A = r1 - 1 + z = N(z)/(e^z - 1), N(z) = (z - 1)e^z + 1,
##   psi1' = z r1 - psi1 (r1 + A),
##   psi1'' = -psi1' (r1 + A) + 2 psi1^2 - 2 z psi1 + z r1,
##   tau = 1 - r1 = (e^z - 1 - z)/(e^z - 1).
## For z < 1/2, A = z S_N(z) r1 and tau = z S_E(z) r1 with the positive-term
## series S_N(z) = sum_{k>=2} (k-1) z^(k-2)/k! and S_E(z) = sum_{k>=2} z^(k-2)/k!,
## which avoid the cancellation in r1 - 1 + z and 1 - r1.

.cll_kk <- 2:22
.cll_cN <- (.cll_kk - 1) / factorial(.cll_kk)
.cll_cE <- 1 / factorial(.cll_kk)
horner <- function(coef, z) {
  s <- coef[length(coef)]
  for (k in (length(coef) - 1L):1L) s <- coef[k] + z * s
  s
}

link_cloglog <- function() {
  ## Tails: for z = e^eta < 1e-5, log G = eta + log{(1 - e^{-z})/z} =
  ## eta + log1p(z^2/6 - z/2) + O(z^3), which stays exact when z is subnormal
  ## or underflows; r1 = z/(e^z - 1) is 0 when z overflows. Quantities equal to
  ## -e^eta beyond eta = log(.Machine$double.xmax) are -Inf.
  base <- function(eta) {
    z <- exp(eta)
    r1 <- z / expm1(z)
    r1[z == 0] <- 1
    r1[z == Inf] <- 0
    lG <- log1mexp(z)
    s <- z < 1e-5
    if (any(s)) lG[s] <- eta[s] + log1p(z[s] * (z[s] / 6 - 0.5))
    list(eta = eta, z = z, lG = lG, lGb = -z, lg = eta - z,
         d1 = 1 - z, d2 = -z, d3 = -z, r1 = r1, r0 = z)
  }
  psi0 <- function(b, deriv = 0L) {
    list(v = b$z, d1 = if (deriv >= 1L) b$z, d2 = if (deriv >= 2L) b$z)
  }
  psi1 <- function(b, deriv = 0L) {
    z <- b$z; r1 <- b$r1
    A <- r1 - 1 + z
    s <- z < 0.5
    if (any(s)) A[s] <- z[s] * horner(.cll_cN, z[s]) * r1[s]
    big <- !is.finite(z)                       # eta > log(.Machine$double.xmax)
    v <- r1 * A
    v[big] <- 0
    d1 <- d2 <- NULL
    if (deriv >= 1L) {
      d1 <- z * r1 - v * (r1 + A)
      d1[big] <- 0
      if (deriv >= 2L) {
        d2 <- -d1 * (r1 + A) + 2 * v * v - 2 * z * v + z * r1
        d2[big] <- 0
      }
    }
    list(v = v, d1 = d1, d2 = d2)
  }
  tau <- function(b) {
    z <- b$z
    t <- 1 - b$r1
    s <- z < 0.5
    if (any(s)) t[s] <- z[s] * horner(.cll_cE, z[s]) * b$r1[s]
    t
  }
  ## log omega = 2 log g - log G - log Gbar = 2 eta - z - log G, written so that
  ## no infinite terms cancel; derivatives 2 - z - r1 and psi1 - z.
  lomega <- function(b, p1v = NULL, p0v = NULL, deriv = 0L) {
    list(l = 2 * b$eta - b$z - b$lG,
         d1 = if (deriv >= 1L) 2 - b$z - b$r1 else NULL,
         d2 = if (deriv >= 2L) p1v - b$z else NULL)
  }
  list(name = "cloglog", base = base, psi1 = psi1, psi0 = psi0,
       lomega = lomega, tau = tau,
       prob = function(eta) -expm1(-exp(eta)),
       glm_family = function() quasibinomial(link = "cloglog"))
}

## ---------------------------------------------------------------------------
## Log-log, G(eta) = exp(-e^{-eta}) = 1 - Gc(-eta), Gc the cloglog inverse link.
## Every quantity is obtained from the cloglog quantities at -eta:
##   lG <-> lGb, r1 <-> r0, psi1 <-> psi0, l1, l3 and odd derivatives change
##   sign, tau changes sign.
link_loglog <- function() {
  L <- link_cloglog()
  base <- function(eta) {
    bi <- L$base(-eta)
    list(eta = eta, lG = bi$lGb, lGb = bi$lG, lg = bi$lg,
         d1 = -bi$d1, d2 = bi$d2, d3 = -bi$d3, r1 = bi$r0, r0 = bi$r1,
         inner = bi)
  }
  flip <- function(p) list(v = p$v, d1 = if (!is.null(p$d1)) -p$d1, d2 = p$d2)
  list(name = "loglog", base = base,
       psi1 = function(b, deriv = 0L) flip(L$psi0(b$inner, deriv)),
       psi0 = function(b, deriv = 0L) flip(L$psi1(b$inner, deriv)),
       lomega = function(b, p1v = NULL, p0v = NULL, deriv = 0L) {   # omega(eta) = omega_c(-eta)
         lc <- L$lomega(b$inner, p0v, p1v, deriv)
         list(l = lc$l, d1 = if (!is.null(lc$d1)) -lc$d1, d2 = lc$d2)
       },
       tau = function(b) -L$tau(b$inner),
       prob = function(eta) exp(-exp(-eta)),
       glm_family = NULL)
}

## ---------------------------------------------------------------------------
## Cauchit, G = 1/2 + arctan(eta)/pi. pcauchy evaluates both tails through
## arctan(1/eta), so lG and lGb are accurate. Generic psi formulas; psi1 is
## negative on (-Inf, eta_c), eta_c ~ -0.4290 (Proposition 1(e)).
link_cauchit <- function() {
  base <- function(eta) {
    lG <- pcauchy(eta, log.p = TRUE)
    lGb <- pcauchy(eta, lower.tail = FALSE, log.p = TRUE)
    lg <- dcauchy(eta, log = TRUE)
    e2 <- 1 + eta * eta
    list(eta = eta, lG = lG, lGb = lGb, lg = lg,
         d1 = -2 * eta / e2, d2 = -2 * (1 - eta * eta) / e2^2,
         d3 = 4 * eta * (3 - eta * eta) / e2^3,
         r1 = exp(lg - lG), r0 = exp(lg - lGb))
  }
  list(name = "cauchit", base = base, psi1 = psi1_generic, psi0 = psi0_generic,
       lomega = lomega_generic, tau = tau_generic,
       prob = function(eta) pcauchy(eta),
       glm_family = function() quasibinomial(link = "cauchit"))
}

make_link <- function(name = c("probit", "logit", "cloglog", "loglog", "cauchit")) {
  name <- match.arg(name)
  switch(name,
         probit = link_probit(), logit = link_logit(), cloglog = link_cloglog(),
         loglog = link_loglog(), cauchit = link_cauchit())
}


## =============================================================================
## 3. Model, penalized criteria and Firth's adjusted score
## =============================================================================
## Binomial-response GLM y_i ~ Bin(m_i, G(x_i' beta)): log-likelihood, score,
## information matrices, the penalized criteria of the paper and Firth's
## adjusted score, with analytic first and second derivatives.
##
## Criteria (type), with multiplier a:
##   "ML" : l(beta)
##   "E"  : l(beta) + a log det i(beta)                      (l^E(.; a))
##   "O"  : l(beta) + a log det j(beta),   on D = {j > 0}    (l^O(.; a))
##   "Ri" : l(beta) - tr{j^{-1} e}/2 + a log det i(beta), on D   (l^R_i, a = 1/N)
##   "Rj" : l(beta) - tr{j^{-1} e}/2 + a log det j(beta), on D   (l^R_j, a = 1/N)
## "Ri" with a = 0 is the unpenalized empirical bias-reducing objective.
##
## Per-observation derivatives with respect to eta_i:
##   score       u_i  = y_i r1 - (m_i - y_i) r0,         u_i' = -w_i
##   observed    w_i  = y_i psi1 + (m_i - y_i) psi0
##   expected    k_i  = m_i omega
##   empirical   u_i^2, with (u^2)' = -2 u w, (u^2)'' = 2 w^2 - 2 u w'.

bp_model <- function(X, y, m = NULL, link = "probit") {
  X <- as.matrix(X)
  if (!is.numeric(X) || !all(is.finite(X))) stop("X must be a finite numeric matrix")
  storage.mode(X) <- "double"
  n <- nrow(X); p <- ncol(X)
  if (is.null(m)) m <- rep(1, n)
  y <- as.numeric(y); m <- as.numeric(m)
  if (length(y) != n || length(m) != n) stop("length(y) and length(m) must equal nrow(X)")
  if (any(!is.finite(y)) || any(!is.finite(m))) stop("y and m must be finite")
  if (any(m < 1) || any(m != round(m))) stop("m must contain positive integers")
  if (any(y < 0) || any(y > m) || any(y != round(y))) stop("y must be integers in [0, m]")
  if (qr(X)$rank < p) stop("X must have full column rank")
  L <- if (is.character(link)) make_link(link) else link
  list(X = X, y = y, m = m, n0 = m - y, link = L, n = n, p = p)
}

## Gradient and Hessian of T(beta) = tr{j^{-1} e}/2, given the factor fj of j
## (with Zt, q), u, w = observed weights, w1 = w', w2 = w''. With
## P = X j^{-1} X' = Z Z', M = X j^{-1} e j^{-1} X' = Z We Z', We = Z' diag(u^2) Z:
##   grad T = X'(-w' diag(M) + e' q)/2,                         e' = -2 u w
##   hess T = {2 A'(P o M)A - X' diag(w'' diag M)X - A'(P o P)B - B'(P o P)A
##             + X' diag(e'' q) X}/2,  A = diag(w')X, B = diag(e')X,
##   e'' = 2 w^2 - 2 u w', P o P = W1 W1', P o M = W1 W2',
##   W1 = rowkron(Z, Z), W2 = rowkron(Z, Z We).
empirical_T_derivs <- function(X, fj, u, w, w1, w2 = NULL) {
  Z <- t(fj$Zt)
  We <- crossprod(Z, (u * u) * Z)
  ZWe <- Z %*% We
  dM <- rowSums(ZWe * Z)
  e1 <- -2 * u * w
  out <- list(grad = 0.5 * drop(crossprod(X, -w1 * dM + e1 * fj$q)), hess = NULL)
  if (!is.null(w2)) {
    W1 <- rowkron(Z, Z)
    W2 <- rowkron(Z, ZWe)
    Am <- w1 * X
    Bm <- e1 * X
    e2 <- 2 * w * w - 2 * u * w1
    G1A <- crossprod(W1, Am)
    G2A <- crossprod(W2, Am)
    G1B <- crossprod(W1, Bm)
    H <- 2 * crossprod(G1A, G2A) - crossprod(X, (w2 * dM) * X) -
      crossprod(G1A, G1B) - crossprod(G1B, G1A) + crossprod(X, (e2 * fj$q) * X)
    out$hess <- symmetrize(0.5 * H)
  }
  out
}

## Evaluate a criterion. level = 0 (value), 1 (+ gradient), 2 (+ Hessian).
## Returns a list with ok (FALSE if beta lies outside the domain or any
## requested quantity is not finite), value, loglik, penalty = value - loglik,
## grad, hess, for "Ri" and "Rj" trace_term = -tr{j^{-1} e}/2, and, for
## level >= 1, iexp = i(beta) (used to standardize gradients).
bp_eval <- function(model, beta, type = c("ML", "E", "O", "Ri", "Rj"),
                    a = 0.5, level = 2L) {
  type <- match.arg(type)
  fail <- list(ok = FALSE, value = -Inf)
  X <- model$X; y <- model$y; n0 <- model$n0; L <- model$link
  eta <- drop(X %*% beta)
  if (!all(is.finite(eta))) return(fail)
  b <- L$base(eta)
  ll <- sum(cmul(y, b$lG)) + sum(cmul(n0, b$lGb))
  if (!is.finite(ll)) return(fail)

  obs_pen <- type %in% c("O", "Rj")
  emp <- type %in% c("Ri", "Rj")
  exp_pen <- type %in% c("E", "Ri")
  need_psi <- obs_pen || emp || level >= 2L
  pd <- if (obs_pen || emp) level else 0L
  if (need_psi) {
    p1 <- L$psi1(b, pd)
    p0 <- L$psi0(b, pd)
    w <- cmul(y, p1$v) + cmul(n0, p0$v)
    w1 <- if (pd >= 1L) cmul(y, p1$d1) + cmul(n0, p0$d1)
    w2 <- if (pd >= 2L) cmul(y, p1$d2) + cmul(n0, p0$d2)
  }
  out <- list(ok = TRUE, value = ll, loglik = ll, penalty = 0,
              grad = NULL, hess = NULL, iexp = NULL)

  if (level >= 1L || exp_pen) {
    lo <- L$lomega(b, if (need_psi) p1$v, if (need_psi) p0$v,
                   deriv = if (exp_pen) level else 0L)
    kexp <- model$m * exp(lo$l)
  }
  if (level >= 1L) {
    u <- cmul(y, b$r1) - cmul(n0, b$r0)
    out$grad <- drop(crossprod(X, u))
    out$iexp <- crossprod(X, kexp * X)
  }
  if (level >= 2L) out$hess <- -crossprod(X, w * X)

  if (exp_pen) {
    fi <- wcp_factor(X, kexp, need_Z = level >= 1L)
    if (is.null(fi)) return(fail)
    pen <- a * fi$logdet
    if (level >= 1L) {
      k1 <- cmul(kexp, lo$d1)                  # 0 where omega underflows
      k2 <- if (level >= 2L) cmul(kexp, lo$d2 + lo$d1 * lo$d1)
      dd <- logdet_derivs(X, fi, k1, k2)
      out$grad <- out$grad + a * dd$grad
      if (level >= 2L) out$hess <- out$hess + a * dd$hess
    }
    out$penalty <- out$penalty + pen
  }
  if (obs_pen || emp) {
    fj <- wcp_factor(X, w, need_Z = emp || level >= 1L)
    if (is.null(fj)) return(fail)                       # outside D
    if (obs_pen) {
      out$penalty <- out$penalty + a * fj$logdet
      if (level >= 1L) {
        dd <- logdet_derivs(X, fj, w1, if (level >= 2L) w2)
        out$grad <- out$grad + a * dd$grad
        if (level >= 2L) out$hess <- out$hess + a * dd$hess
      }
    }
    if (emp) {
      if (level < 1L) u <- cmul(y, b$r1) - cmul(n0, b$r0)
      out$trace_term <- -0.5 * sum(u * u * fj$q)
      out$penalty <- out$penalty + out$trace_term
      if (level >= 1L) {
        dT <- empirical_T_derivs(X, fj, u, w, w1, if (level >= 2L) w2)
        out$grad <- out$grad - dT$grad
        if (level >= 2L) out$hess <- out$hess - dT$hess
      }
    }
  }
  out$value <- ll + out$penalty
  if (level >= 2L) out$hess <- symmetrize(out$hess)
  if (!is.finite(out$value) ||
      (level >= 1L && (!all(is.finite(out$grad)) || !all(is.finite(out$iexp)))) ||
      (level >= 2L && !all(is.finite(out$hess)))) return(fail)
  out
}

## Firth's mean bias-reducing adjusted score U = S + A, with
## A_r = tr{i^{-1}(P_r + Q_r)}/2 = sum_i h_i (g_i'/g_i) x_ir / 2,
## h_i = m_i omega_i x_i' i^{-1} x_i. With jacobian = TRUE also returns the
## analytic Jacobian dU/dbeta = -j + dA/dbeta, where
##   dA/dbeta = {X' diag(h l2) X + X' diag(l1) dh/dbeta}/2,
##   dh/dbeta = diag(k' q) X - diag(k) (Q o Q) diag(k') X,  k = m omega.
bp_adjscore <- function(model, beta, jacobian = FALSE) {
  fail <- list(ok = FALSE)
  X <- model$X; y <- model$y; n0 <- model$n0; L <- model$link
  eta <- drop(X %*% beta)
  if (!all(is.finite(eta))) return(fail)
  b <- L$base(eta)
  if (jacobian) {
    p1 <- L$psi1(b, 0L)
    p0 <- L$psi0(b, 0L)
  }
  lo <- L$lomega(b, if (jacobian) p1$v, if (jacobian) p0$v, deriv = if (jacobian) 2L else 0L)
  kexp <- model$m * exp(lo$l)
  fi <- wcp_factor(X, kexp)
  if (is.null(fi)) return(fail)
  u <- cmul(y, b$r1) - cmul(n0, b$r0)
  h <- kexp * fi$q
  A <- 0.5 * drop(crossprod(X, cmul(h, b$d1)))
  U <- drop(crossprod(X, u)) + A
  out <- list(ok = TRUE, U = U, A = A, iexp = fi$K, ifac = fi, J = NULL)
  if (jacobian) {
    w <- cmul(y, p1$v) + cmul(n0, p0$v)
    k1 <- cmul(kexp, lo$d1)
    Z <- t(fi$Zt)
    W <- rowkron(Z, Z)
    Hb <- (k1 * fi$q) * X - kexp * (W %*% crossprod(W, k1 * X))
    d1h <- b$d1; d1h[kexp == 0] <- 0            # rows of Hb are 0 where omega underflows
    JA <- 0.5 * (crossprod(X, cmul(h, b$d2) * X) + crossprod(X, d1h * Hb))
    out$J <- -crossprod(X, w * X) + JA
  }
  if (!all(is.finite(U)) || !all(is.finite(fi$K)) ||
      (jacobian && !all(is.finite(out$J)))) return(fail)
  out
}

## Score S, observed information j, expected information i and empirical
## information e = sum_i S_i S_i' at beta.
bp_info <- function(model, beta) {
  X <- model$X; y <- model$y; n0 <- model$n0; L <- model$link
  b <- L$base(drop(X %*% beta))
  p1 <- L$psi1(b, 0L); p0 <- L$psi0(b, 0L)
  u <- cmul(y, b$r1) - cmul(n0, b$r0)
  w <- cmul(y, p1$v) + cmul(n0, p0$v)
  kexp <- model$m * exp(L$lomega(b)$l)
  list(score = drop(crossprod(X, u)),
       observed = crossprod(X, w * X),
       expected = crossprod(X, kexp * X),
       empirical = crossprod(X, (u * u) * X),
       weights = cbind(observed = w, expected = kexp, score = u))
}


## =============================================================================
## 4. Optimization: trust region, quasi-Fisher scoring, Newton refinement
## =============================================================================
## Trust-region maximization of the penalized criteria, quasi-Fisher scoring
## and Newton refinement for Firth's adjusted score, acceptance checks.
## Nothing here modifies a criterion: shifts in the trust-region subproblem
## only choose steps.

## Criterion restricted to free coordinates. fixed: NULL or a length-p vector
## with NA for free coordinates and the fixed values elsewhere.
## info_scale: optional fixed matrix (free block) used in place of i(theta)
## to standardize gradients (ev$gscale); ev$iexp stays i at theta.
bp_criterion <- function(model, type = "O", a = 0.5, fixed = NULL, info_scale = NULL) {
  p <- model$p
  if (!is.null(fixed) && length(fixed) != p) stop("fixed must have length p")
  free <- if (is.null(fixed)) seq_len(p) else which(is.na(fixed))
  if (!length(free)) stop("no free coordinates")
  full <- function(theta) {
    if (is.null(fixed)) return(theta)
    b <- fixed
    b[free] <- theta
    b
  }
  fn <- function(theta, level = 2L) {
    ev <- bp_eval(model, full(theta), type, a, level)
    if (!ev$ok) return(ev)
    if (!is.null(fixed)) {
      if (level >= 1L) {
        ev$grad <- ev$grad[free]
        ev$iexp <- ev$iexp[free, free, drop = FALSE]
      }
      if (level >= 2L) ev$hess <- ev$hess[free, free, drop = FALSE]
    }
    if (level >= 1L && !is.null(info_scale)) ev$gscale <- info_scale
    ev
  }
  list(fn = fn, full = full, free = free, type = type, a = a, fixed = fixed)
}

## Standardized gradient i^{-1/2} U (symmetric inverse square root).
std_grad <- function(ev) {
  S <- isqrt_sym(if (is.null(ev$gscale)) ev$iexp else ev$gscale)
  if (is.null(S)) return(rep(NA_real_, length(ev$grad)))
  drop(S %*% ev$grad)
}

## Eigenvalues of the expected information (free block for constrained fits)
## and its condition number; recorded for every candidate.
info_spectrum <- function(M) {
  v <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  list(values = v, cond = if (min(v) > 0) max(v) / min(v) else Inf)
}

## Curvature check for a maximum: lambda_max(H) < -ctol * max(1, ||H||_2).
curvature_check <- function(H, ctol = 1e-10) {
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  nrm <- max(abs(ev))
  list(ok = max(ev) < -ctol * max(1, nrm), lambda_max = max(ev), norm2 = nrm,
       eigenvalues = ev)
}

## Trust-region subproblem: maximize g's + s'Hs/2 subject to ||s|| <= Delta.
## Solved through the eigendecomposition of B = -H and the secular equation
## 1/||s(sigma)|| = 1/Delta, s(sigma) = -(B + sigma I)^{-1}(-g), including the
## hard case (Moré and Sorensen, 1983).
tr_subproblem <- function(g, H, Delta) {
  p <- length(g)
  e <- eigen(-H, symmetric = TRUE)
  ord <- rev(seq_len(p))
  lam <- e$values[ord]                       # ascending
  V <- e$vectors[, ord, drop = FALSE]
  gh <- drop(crossprod(V, -g))               # gradient of the minimization form
  sc <- max(1, abs(lam))
  ltol <- 1e-12 * sc
  snorm <- function(sig) sqrt(sum((gh / (lam + sig))^2))
  step <- function(sig) -drop(V %*% (gh / (lam + sig)))
  if (lam[1] > ltol) {
    sN <- step(0)
    if (sqrt(sum(sN^2)) <= Delta) return(list(s = sN, boundary = FALSE, sigma = 0, hard = FALSE))
  }
  I1 <- which(lam - lam[1] <= ltol)
  g1 <- sqrt(sum(gh[I1]^2))
  gn <- sqrt(sum(gh^2))
  if (lam[1] <= ltol && g1 <= 1e-12 * max(gn, .Machine$double.xmin)) {
    sig <- -lam[1]
    idx <- setdiff(seq_len(p), I1)
    sh <- if (length(idx)) -drop(V[, idx, drop = FALSE] %*% (gh[idx] / (lam[idx] + sig))) else numeric(p)
    ns <- sqrt(sum(sh^2))
    if (ns <= Delta) {
      v1 <- V[, I1[1]]
      v1 <- v1 * sign(v1[which.max(abs(v1))])
      return(list(s = sh + sqrt(Delta^2 - ns^2) * v1, boundary = TRUE, sigma = sig, hard = TRUE))
    }
  }
  lo <- max(0, -lam[1])
  hi <- lo + gn / Delta                      # ||s(hi)|| <= Delta
  sig <- hi
  for (it in 1:200) {
    ns <- snorm(sig)
    phi <- 1 / ns - 1 / Delta
    if (abs(ns - Delta) <= 1e-12 * Delta) break
    if (phi < 0) lo <- sig else hi <- sig
    dphi <- sum(gh^2 / (lam + sig)^3) / ns^3
    nxt <- sig - phi / dphi
    sig <- if (is.finite(nxt) && nxt > lo && nxt < hi) nxt else (lo + hi) / 2
    if (hi - lo <= 1e-15 * max(1, hi)) break
  }
  list(s = step(sig), boundary = TRUE, sigma = sig, hard = FALSE)
}

## Trust-region maximization from one start. Radius starts at Delta0 = 1,
## is capped at Delta_max = 100, a trial step is accepted when the ratio of
## actual to predicted increase exceeds 0.1, the radius is quartered below
## ratio 0.25 and doubled above ratio 0.75 when the step reaches the boundary.
## A non-finite value or a point outside the domain counts as ratio -Inf.
## When the predicted increase is below 100 eps max(1, |Q|), the computed
## change in the criterion is dominated by its rounding error (which can exceed
## that bound when the summed terms span many orders of magnitude), so neither
## the ratio nor the computed change is informative; the step is then accepted
## when it reduces the standardized gradient, as in Newton refinement.
## Convergence: max|i^{-1/2} U| < gtol and the curvature check.
tr_maximize <- function(fn, theta0, maxit = 500L, gtol = 1e-8, ctol = 1e-10,
                        Delta0 = 1, Delta_max = 100) {
  theta <- as.numeric(theta0)
  ev <- fn(theta, 2L)
  res <- function(status, it, ev, theta, gs = NULL, cc = NULL) {
    if (!ev$ok) return(list(status = status, accepted = FALSE, theta = theta,
                            value = -Inf, iterations = it))
    if (is.null(gs)) gs <- std_grad(ev)
    if (is.null(cc)) cc <- curvature_check(ev$hess, ctol)
    gmax <- max(abs(gs))
    isp <- info_spectrum(ev$iexp)
    list(status = status,
         accepted = status == "converged" || (is.finite(gmax) && gmax < gtol && cc$ok),
         theta = theta, value = ev$value, loglik = ev$loglik, penalty = ev$penalty,
         grad = ev$grad, hess = ev$hess, std_grad_max = gmax,
         lambda_max = cc$lambda_max, hess_norm = cc$norm2,
         hess_eigenvalues = cc$eigenvalues, info_eigenvalues = isp$values,
         info_cond = isp$cond, iterations = it)
  }
  if (!ev$ok) return(res("start_infeasible", 0L, ev, theta))
  Delta <- Delta0
  eps_f <- 100 * .Machine$double.eps
  for (it in 0:maxit) {
    gs <- std_grad(ev)
    cc <- curvature_check(ev$hess, ctol)
    if (all(is.finite(gs)) && max(abs(gs)) < gtol && cc$ok)
      return(res("converged", it, ev, theta, gs, cc))
    if (it == maxit) return(res("maxit", it, ev, theta, gs, cc))
    st <- tr_subproblem(ev$grad, ev$hess, Delta)
    s <- st$s
    pred <- sum(ev$grad * s) + 0.5 * sum(s * drop(ev$hess %*% s))
    evt <- fn(theta + s, 2L)
    if (!evt$ok) {
      rho <- -Inf
    } else if (pred <= eps_f * max(1, abs(ev$value))) {
      gst <- std_grad(evt)
      rho <- if (all(is.finite(gst)) && max(abs(gst)) < max(abs(gs))) 1 else -Inf
    } else {
      rho <- (evt$value - ev$value) / pred
    }
    if (rho > 0.1) {
      theta <- theta + s
      ev <- evt
    }
    if (rho < 0.25) {
      Delta <- Delta / 4
    } else if (rho > 0.75 && st$boundary) {
      Delta <- min(2 * Delta, Delta_max)
    }
    if (Delta < 1e-14 * (1 + sqrt(sum(theta^2)))) {
      gs <- std_grad(ev)
      cc <- curvature_check(ev$hess, ctol)
      st_ <- if (all(is.finite(gs)) && max(abs(gs)) < gtol && cc$ok) "converged" else "stalled"
      return(res(st_, it + 1L, ev, theta, gs, cc))
    }
  }
}

## Newton refinement of an objective candidate: theta <- theta - H^{-1} grad
## with step-halving (at most max_halving) until the standardized gradient
## decreases. Requires a negative definite Hessian at every iterate.
newton_refine <- function(fn, theta0, maxit = 50L, gtol = 1e-8, ctol = 1e-10,
                          max_halving = 30L) {
  theta <- as.numeric(theta0)
  ev <- fn(theta, 2L)
  if (!ev$ok) return(list(status = "start_infeasible", accepted = FALSE, theta = theta, value = -Inf))
  for (it in 0:maxit) {
    gs <- std_grad(ev)
    cc <- curvature_check(ev$hess, ctol)
    if (!all(is.finite(gs))) break
    if (max(abs(gs)) < gtol && cc$ok) break
    if (!cc$ok || it == maxit) break
    dir <- tryCatch(-solve(ev$hess, ev$grad), error = function(e) NULL)
    if (is.null(dir)) break
    done <- FALSE
    for (h in 0:max_halving) {
      evt <- fn(theta + dir / 2^h, 2L)
      if (evt$ok) {
        gst <- std_grad(evt)
        if (all(is.finite(gst)) && max(abs(gst)) < max(abs(gs))) {
          theta <- theta + dir / 2^h
          ev <- evt
          done <- TRUE
          break
        }
      }
    }
    if (!done) break
  }
  gs <- std_grad(ev)
  cc <- curvature_check(ev$hess, ctol)
  ok <- all(is.finite(gs)) && max(abs(gs)) < gtol && cc$ok
  isp <- info_spectrum(ev$iexp)
  list(status = if (ok) "converged" else "newton_failed", accepted = ok,
       theta = theta, value = ev$value, loglik = ev$loglik, penalty = ev$penalty,
       grad = ev$grad, hess = ev$hess, std_grad_max = max(abs(gs)),
       lambda_max = cc$lambda_max, hess_norm = cc$norm2,
       hess_eigenvalues = cc$eigenvalues, info_eigenvalues = isp$values,
       info_cond = isp$cond, iterations = it)
}

## Jacobian check for an adjusted-score root:
## sigma_min(J) > stol * max(1, sigma_max(J)).
jacobian_check <- function(J, stol = 1e-10) {
  sv <- svd(J, nu = 0, nv = 0)$d
  list(ok = min(sv) > stol * max(1, max(sv)), sigma_min = min(sv), sigma_max = max(sv),
       values = sv)
}

as_root_result <- function(status, model, beta, ev, it, gtol, stol) {
  beta <- unname(beta)
  if (!ev$ok) return(list(status = status, accepted = FALSE, beta = beta, iterations = it))
  ej <- if (is.null(ev$J)) bp_adjscore(model, beta, jacobian = TRUE) else ev
  S <- isqrt_sym(ev$iexp)
  rmax <- if (is.null(S)) Inf else max(abs(drop(S %*% ev$U)))
  jc <- if (ej$ok) jacobian_check(ej$J, stol) else
    list(ok = FALSE, sigma_min = NA_real_, sigma_max = NA_real_, values = NA_real_)
  isp <- info_spectrum(ev$iexp)
  list(status = status, accepted = is.finite(rmax) && rmax < gtol && jc$ok,
       beta = beta, U = ev$U, std_resid_max = rmax,
       sigma_min = jc$sigma_min, sigma_max = jc$sigma_max, jac_singular_values = jc$values,
       info_eigenvalues = isp$values, info_cond = isp$cond, iterations = it)
}

## Standardized residual of the adjusted score. merit = U' i^{-1} U is
## computed from the equilibrated Cholesky factor of i; the sup-norm of
## i^{-1/2} U (symmetric inverse square root) is computed only when it can
## decide convergence, since ||r||_inf <= ||r||_2 <= sqrt(p) ||r||_inf.
root_resid <- function(ev, gtol) {
  merit <- chol_quadform(ev$ifac, ev$U)
  if (!is.finite(merit)) return(c(merit = Inf, rmax = Inf))      # never compared as NA
  rmax <- sqrt(merit)
  if (rmax >= gtol && rmax < gtol * sqrt(length(ev$U))) {
    S <- isqrt_sym(ev$iexp)
    rmax <- if (is.null(S)) Inf else max(abs(drop(S %*% ev$U)))
  }
  c(merit = merit, rmax = rmax)
}

## Quasi-Fisher scoring for S + A = 0: beta <- beta + 2^{-h} i^{-1}(S + A),
## h = 0, ..., max_halving, accepting the first step that reduces the squared
## standardized residual U' i^{-1} U.
qfs_solve <- function(model, beta0, maxit = 500L, gtol = 1e-8, stol = 1e-10,
                      max_halving = 30L) {
  beta <- as.numeric(beta0)
  ev <- bp_adjscore(model, beta)
  if (!ev$ok) return(list(status = "start_infeasible", accepted = FALSE, beta = beta, iterations = 0L))
  cur <- root_resid(ev, gtol)
  for (it in 0:maxit) {
    if (cur["rmax"] < gtol) return(as_root_result("converged", model, beta, ev, it, gtol, stol))
    if (it == maxit) break
    dir <- chol_solve(ev$ifac, ev$U)
    done <- FALSE
    for (h in 0:max_halving) {
      evt <- bp_adjscore(model, beta + dir / 2^h)
      if (evt$ok) {
        nr <- root_resid(evt, gtol)
        if (nr["merit"] < cur["merit"]) {
          beta <- beta + dir / 2^h
          ev <- evt
          cur <- nr
          done <- TRUE
          break
        }
      }
    }
    if (!done) return(as_root_result("stalled", model, beta, ev, it, gtol, stol))
  }
  as_root_result("maxit", model, beta, ev, maxit, gtol, stol)
}

## Newton root refinement for S + A = 0 with the analytic Jacobian and
## residual-decreasing step-halving.
newton_root <- function(model, beta0, maxit = 500L, gtol = 1e-8, stol = 1e-10,
                        max_halving = 30L) {
  beta <- as.numeric(beta0)
  ev <- bp_adjscore(model, beta, jacobian = TRUE)
  if (!ev$ok) return(list(status = "start_infeasible", accepted = FALSE, beta = beta, iterations = 0L))
  cur <- root_resid(ev, gtol)
  for (it in 0:maxit) {
    if (cur["rmax"] < gtol) return(as_root_result("converged", model, beta, ev, it, gtol, stol))
    if (it == maxit) break
    dir <- tryCatch(-solve(ev$J, ev$U), error = function(e) NULL)
    if (is.null(dir)) return(as_root_result("singular_jacobian", model, beta, ev, it, gtol, stol))
    done <- FALSE
    for (h in 0:max_halving) {
      evt <- bp_adjscore(model, beta + dir / 2^h, jacobian = TRUE)
      if (evt$ok) {
        nr <- root_resid(evt, gtol)
        if (nr["merit"] < cur["merit"]) {
          beta <- beta + dir / 2^h
          ev <- evt
          cur <- nr
          done <- TRUE
          break
        }
      }
    }
    if (!done) return(as_root_result("stalled", model, beta, ev, it, gtol, stol))
  }
  as_root_result("maxit", model, beta, ev, maxit, gtol, stol)
}


## =============================================================================
## 5. Protocol design, starts, multistart fitting, selection, inference
## =============================================================================
## Starting values, multistart fitting, the deterministic selection rules and
## the rescue sequence of the simulation protocol.

## Fixed design of the protocol: intercept,
## x1 = Phi^{-1}{0.05 + 0.90 (i - 1/2)/n}, and x2 = 1 at the indices
## floor((l - 1/2) n / k_n) + 1, l = 1, ..., k_n, k_n = floor(0.30 n + 1/2).
## Columns are neither centred nor rescaled.
protocol_design <- function(n) {
  x1 <- qnorm(0.05 + 0.90 * (seq_len(n) - 0.5) / n)
  kn <- floor(0.30 * n + 0.5)
  x2 <- numeric(n)
  x2[floor((seq_len(kn) - 0.5) * n / kn) + 1] <- 1
  cbind(intercept = 1, x1 = x1, x2 = x2)
}

## True parameter vectors (alpha, beta1, beta2) of settings S1 and S2.
protocol_truth <- list(S1 = c(-0.5, 0.5, 1), S2 = c(-1.5, 1, 2))

## Routine and rescue starts from a 50 x p matrix Z of independent N(0, 1)
## draws: routine = {0, +2 e_r, -2 e_r, 2 Z[1:10, ]}, rescue = 4 Z[11:50, ].
protocol_starts <- function(p, Z) {
  Z <- as.matrix(Z)
  if (ncol(Z) != p || nrow(Z) < 10L) stop("Z must have p columns and at least 10 rows")
  E <- diag(p)
  list(routine = rbind(numeric(p), 2 * E, -2 * E, 2 * Z[1:10, , drop = FALSE]),
       rescue = if (nrow(Z) >= 50L) 4 * Z[11:50, , drop = FALSE] else NULL)
}

## Two solutions are the same when their relative Euclidean distance
## ||b - b'|| / (1 + max(||b||, ||b'||)) is at most dtol.
same_solution <- function(b1, b2, dtol = 1e-5) {
  sqrt(sum((b1 - b2)^2)) / (1 + max(sqrt(sum(b1^2)), sqrt(sum(b2^2)))) <= dtol
}

## Group candidate vectors into distinct solutions (greedy, in input order).
distinct_solutions <- function(B, dtol = 1e-5) {
  if (!length(B)) return(list())
  groups <- list()
  for (k in seq_along(B)) {
    hit <- FALSE
    for (gi in seq_along(groups)) {
      if (same_solution(B[[k]], B[[groups[[gi]][1]]], dtol)) {
        groups[[gi]] <- c(groups[[gi]], k)
        hit <- TRUE
        break
      }
    }
    if (!hit) groups[[length(groups) + 1L]] <- k
  }
  groups
}

## One entry per distinct solution: its best member (highest criterion for
## maxima, smallest standardized residual for roots) and the number of starts
## that reached it.
summarize_distinct <- function(acc, kind = c("max", "root"), dtol = 1e-5) {
  kind <- match.arg(kind)
  groups <- distinct_solutions(lapply(acc, function(cc) cc$beta), dtol)
  lapply(groups, function(gg) {
    k <- if (kind == "max") gg[which.max(vapply(acc[gg], function(cc) cc$value, 0))] else
      gg[which.min(vapply(acc[gg], function(cc) cc$std_resid_max, 0))]
    out <- acc[[k]][if (kind == "max") c("beta", "value") else c("beta", "std_resid_max")]
    out$n_starts <- length(gg)
    out
  })
}

## Order indices by Euclidean norm, then lexicographically.
order_norm_lex <- function(B) {
  M <- do.call(rbind, B)
  do.call(order, c(list(sqrt(rowSums(M^2))), lapply(seq_len(ncol(M)), function(j) M[, j])))
}

## Selection among accepted maxima: highest criterion; values within
## tie_tol * (1 + |Q_max|) of the maximum are ties, resolved by minimum
## Euclidean norm and then lexicographic order.
select_maximum <- function(cands, tie_tol = 1e-10) {
  acc <- Filter(function(cc) isTRUE(cc$accepted), cands)
  if (!length(acc)) return(NULL)
  v <- vapply(acc, function(cc) cc$value, 0)
  vmax <- max(v)
  tie <- which(v >= vmax - tie_tol * (1 + abs(vmax)))
  pick <- tie[order_norm_lex(lapply(acc[tie], function(cc) cc$beta))[1]]
  acc[[pick]]
}

## Selection among accepted adjusted-score roots: closest to ref (the selected
## l^E estimate, or zero), then lexicographic order.
select_root <- function(cands, ref) {
  acc <- Filter(function(cc) isTRUE(cc$accepted), cands)
  if (!length(acc)) return(NULL)
  d <- vapply(acc, function(cc) sqrt(sum((cc$beta - ref)^2)), 0)
  M <- do.call(rbind, lapply(acc, function(cc) cc$beta))
  o <- do.call(order, c(list(d), lapply(seq_len(ncol(M)), function(j) M[, j])))
  acc[[o[1]]]
}

## Multistart trust-region fit of a criterion. starts: matrix with p columns
## (full-dimensional; for constrained fits only the free coordinates are used,
## i.e. starts are projected onto the nuisance coordinates).
fit_criterion <- function(model, type = "O", a = 0.5, starts, maxit = 500L,
                          gtol = 1e-8, ctol = 1e-10, fixed = NULL, info_scale = NULL) {
  crit <- bp_criterion(model, type, a, fixed, info_scale)
  starts <- as.matrix(starts)
  cands <- vector("list", nrow(starts))
  for (k in seq_len(nrow(starts))) {
    r <- tr_maximize(crit$fn, starts[k, crit$free], maxit, gtol, ctol)
    r$beta <- crit$full(r$theta)
    r$start_index <- k
    r$start <- starts[k, ]
    cands[[k]] <- r
  }
  sel <- select_maximum(cands)
  acc <- Filter(function(cc) isTRUE(cc$accepted), cands)
  dist <- summarize_distinct(acc, "max")
  list(selected = sel, candidates = cands, n_accepted = length(acc),
       n_distinct = length(dist),
       distinct = dist,
       type = type, a = a, fixed = fixed)
}

## Protocol fit of a criterion with the fixed rescue sequence: routine starts
## (ceiling maxit); if no candidate is accepted, the routine starts again with
## ceiling maxit_rescue, then the rescue starts with that ceiling, then Newton
## refinement of every unaccepted candidate whose standardized gradient is
## within 10 * gtol. Tolerances are never relaxed.
## candidates: the set the selection is made from (the last stage reached, with
## Newton-refined candidates in place of their unrefined versions).
## history: every candidate of every stage, each labelled by $stage ("routine",
## "rerun", "rescue", "newton"); refined candidates are added, not substituted.
tag_stage <- function(cands, stage) lapply(cands, function(cc) { cc$stage <- stage; cc })

fit_criterion_protocol <- function(model, type = "O", a = 0.5, routine, rescue = NULL,
                                   maxit = 500L, maxit_rescue = 2000L, gtol = 1e-8,
                                   ctol = 1e-10, fixed = NULL, info_scale = NULL) {
  f1 <- fit_criterion(model, type, a, routine, maxit, gtol, ctol, fixed, info_scale)
  f1$candidates <- tag_stage(f1$candidates, "routine")
  if (!is.null(f1$selected)) {
    f1$selected$stage <- "routine"
    return(c(f1, list(stage = "routine", history = f1$candidates)))
  }
  f2 <- fit_criterion(model, type, a, routine, maxit_rescue, gtol, ctol, fixed, info_scale)
  cands <- tag_stage(f2$candidates, "rerun")
  if (!is.null(rescue)) {
    f3 <- fit_criterion(model, type, a, rescue, maxit_rescue, gtol, ctol, fixed, info_scale)
    cands <- c(cands, tag_stage(f3$candidates, "rescue"))
  }
  history <- c(f1$candidates, cands)
  stage <- "rescue"
  if (!any(vapply(cands, function(cc) isTRUE(cc$accepted), TRUE))) {
    crit <- bp_criterion(model, type, a, fixed, info_scale)
    for (k in seq_along(cands)) {
      cc <- cands[[k]]
      if (is.finite(cc$value) && !is.null(cc$std_grad_max) &&
          is.finite(cc$std_grad_max) && cc$std_grad_max < 10 * gtol) {
        r <- newton_refine(crit$fn, cc$theta, gtol = gtol, ctol = ctol)
        r$beta <- crit$full(r$theta)
        r$start_index <- cc$start_index
        r$start <- cc$start
        r$stage <- "newton"
        r$refined_from <- cc$stage
        cands[[k]] <- r
        history[[length(history) + 1L]] <- r
      }
    }
    stage <- "rescue_newton"
  }
  sel <- select_maximum(cands)
  acc <- Filter(function(cc) isTRUE(cc$accepted), cands)
  dist <- summarize_distinct(acc, "max")
  list(selected = sel, candidates = cands, history = history, n_accepted = length(acc),
       n_distinct = length(dist),
       distinct = dist,
       type = type, a = a, fixed = fixed,
       stage = if (is.null(sel)) "unresolved" else stage)
}

## Multistart solution of S + A = 0 by quasi-Fisher scoring; a start whose
## scoring run from a feasible start is not accepted is continued by Newton root refinement from its
## end point. Roots are selected by select_root(ref).
fit_adjscore <- function(model, starts, ref = NULL, maxit = 500L, gtol = 1e-8,
                         stol = 1e-10) {
  starts <- as.matrix(starts)
  if (is.null(ref)) ref <- numeric(model$p)
  cands <- vector("list", nrow(starts))
  for (k in seq_len(nrow(starts))) {
    r <- qfs_solve(model, starts[k, ], maxit, gtol, stol)
    r$method <- "qfs"
    if (!isTRUE(r$accepted) && r$status != "start_infeasible") {
      r2 <- newton_root(model, r$beta, maxit, gtol, stol)
      r2$method <- "qfs+newton"
      r2$qfs_status <- r$status
      r <- r2
    }
    r$start_index <- k
    r$start <- starts[k, ]
    cands[[k]] <- r
  }
  acc <- Filter(function(cc) isTRUE(cc$accepted), cands)
  dist <- summarize_distinct(acc, "root")
  list(selected = select_root(cands, ref), candidates = cands,
       n_accepted = length(acc), n_distinct = length(dist),
       distinct = dist)
}

## Adjusted-score protocol fit: routine starts; if no root is accepted, the
## routine starts with ceiling maxit_rescue and then the rescue starts.
fit_adjscore_protocol <- function(model, routine, rescue = NULL, ref = NULL,
                                  maxit = 500L, maxit_rescue = 2000L,
                                  gtol = 1e-8, stol = 1e-10) {
  f1 <- fit_adjscore(model, routine, ref, maxit, gtol, stol)
  f1$candidates <- tag_stage(f1$candidates, "routine")
  if (!is.null(f1$selected)) {
    f1$selected$stage <- "routine"
    return(c(f1, list(stage = "routine", history = f1$candidates)))
  }
  f2 <- fit_adjscore(model, routine, ref, maxit_rescue, gtol, stol)
  cands <- tag_stage(f2$candidates, "rerun")
  if (!is.null(rescue))
    cands <- c(cands, tag_stage(fit_adjscore(model, rescue, ref, maxit_rescue, gtol, stol)$candidates, "rescue"))
  if (is.null(ref)) ref <- numeric(model$p)
  sel <- select_root(cands, ref)
  acc <- Filter(function(cc) isTRUE(cc$accepted), cands)
  dist <- summarize_distinct(acc, "root")
  list(selected = sel, candidates = cands, history = c(f1$candidates, cands),
       n_accepted = length(acc), n_distinct = length(dist), distinct = dist,
       stage = if (is.null(sel)) "unresolved" else "rescue")
}

## Penalized-profile statistic for coordinate idx at value t:
## L(t) = 2{Q(beta_hat) - sup_{beta_-idx} Q(t, beta_-idx)}, with the full
## determinant at the constrained parameter. full_fit is the output of
## fit_criterion() or fit_criterion_protocol(); its criterion and multiplier
## are used. The constrained fit starts from the
## projections onto the nuisance coordinates of the supplied full starts and of
## every distinct accepted unrestricted maximum, and uses the rescue sequence
## with the projected rescue starts. A constrained maximum exceeding the
## unrestricted one by more than 1e-8 (1 + |Q|) flags an unresolved
## unrestricted fit; the statistic is never truncated at zero.
## info_scale: optional fixed standardization for the constrained fits (see
## bp_criterion()); NULL uses i at each candidate.
profile_stat <- function(model, full_fit, idx, t, starts, rescue = NULL, maxit = 500L,
                         maxit_rescue = 2000L, gtol = 1e-8, ctol = 1e-10, info_scale = NULL) {
  if (is.null(full_fit$selected))
    return(list(stat = NA_real_, fit = NULL, flag = "unrestricted_missing"))
  fixed <- rep(NA_real_, model$p)
  fixed[idx] <- t
  maxima <- do.call(rbind, lapply(full_fit$distinct, function(s) s$beta))
  fc <- fit_criterion_protocol(model, full_fit$type, full_fit$a,
                               routine = rbind(as.matrix(starts), maxima), rescue = rescue,
                               maxit = maxit, maxit_rescue = maxit_rescue, gtol = gtol,
                               ctol = ctol, fixed = fixed, info_scale = info_scale)
  if (is.null(fc$selected)) return(list(stat = NA_real_, fit = fc, flag = "constrained_unresolved"))
  Q <- full_fit$selected$value
  flag <- if (fc$selected$value > Q + 1e-8 * (1 + abs(Q))) "unrestricted_unresolved" else "ok"
  list(stat = 2 * (Q - fc$selected$value), fit = fc, flag = flag)
}

## Expected-information Wald intervals: beta_r +/- z sqrt{(i(beta)^{-1})_rr}.
wald_ci <- function(model, beta, level = 0.95) {
  I <- bp_info(model, beta)$expected
  se <- sqrt(diag(solve(I)))
  z <- qnorm(1 - (1 - level) / 2)
  cbind(estimate = beta, se = se, lower = beta - z * se, upper = beta + z * se)
}


## Separated samples. Under complete or quasi-complete
## separation a criterion can have two local maxima along the separating
## direction, separated by a shallow dip, the higher one in a narrow basin that
## random starts rarely reach. separation_direction() returns the direction of
## the program that decided the category (lp1 for complete, lp2 for
## quasi-complete separation; original coordinates) and NULL otherwise.
## ray_starts() places starts at beta_hat + t d/|d| for t = +-(0.5, 1, 2, 4, 8)
## and t = +-(0.25, 0.5, 1, 2) |beta_hat|. fit_criterion_search() is
## fit_criterion_protocol() followed, when a direction is given and a maximum
## was accepted, by fits from those starts (stage "ray") and selection over
## all candidates; $ray_starts holds the starts used, and $ray_improved is
## TRUE when the selection is a different maximum from the protocol fit's.
separation_direction <- function(sep) {
  if (is.null(sep) || !isTRUE(sep$category %in% c("complete", "quasi-complete"))) return(NULL)
  res <- if (!is.null(sep$rechecked)) sep$rechecked else sep$scaled
  d <- if (sep$category == "complete") res$lp1$direction else res$lp2$direction
  if (is.null(d) || !all(is.finite(d)) || sum(d^2) == 0) NULL else d
}
ray_starts <- function(beta_hat, d) {
  u <- d / sqrt(sum(d^2)); nb <- sqrt(sum(beta_hat^2))
  tt <- sort(unique(c(0.5, 1, 2, 4, 8, c(0.25, 0.5, 1, 2) * nb)))
  tt <- c(tt, -tt)
  t(vapply(tt, function(t) beta_hat + t * u, numeric(length(beta_hat))))
}
fit_criterion_search <- function(model, type = "O", a = 0.5, routine, rescue = NULL, direction = NULL,
                                 maxit = 500L, maxit_rescue = 2000L, gtol = 1e-8, ctol = 1e-10) {
  f <- fit_criterion_protocol(model, type, a, routine, rescue, maxit, maxit_rescue, gtol, ctol)
  if (is.null(direction) || is.null(f$selected)) return(f)
  rs <- ray_starts(f$selected$beta, direction)
  rc <- tag_stage(fit_criterion(model, type, a, rs, maxit, gtol, ctol)$candidates, "ray")
  cands <- c(f$candidates, rc)
  acc <- Filter(function(cc) isTRUE(cc$accepted), cands)
  dist <- summarize_distinct(acc, "max")
  sel <- select_maximum(cands)
  f$ray_improved <- !same_solution(sel$beta, f$selected$beta)   # a different, higher maximum
  f$selected <- sel; f$candidates <- cands; f$history <- c(f$history, rc)
  f$n_accepted <- length(acc); f$n_distinct <- length(dist); f$distinct <- dist; f$ray_starts <- rs
  if (f$ray_improved) f$stage <- paste0(f$stage, "+ray")
  f
}


## =============================================================================
## 6. First-order biases, control variate, exact intercept-only means
## =============================================================================
## Formal first-order biases and displacement (Proposition 2), the control
## variate of the protocol, and exact finite-sample means by enumeration for
## intercept-only models.

## First-order quantities at the true value beta0 for design X, totals m and
## link. With h_i = m_i omega_i x_i' i^{-1} x_i:
##   A = X'(h l1)/2,  zetaE - A = X'(h tau)/2,  zetaE = A + (zetaE - A),
##   b_ML = -i^{-1} A,
##   bE(a) = i^{-1}(2 a zetaE - A),
##   bO(a) = i^{-1}{2 a (2 zetaE - A) - A},
##   displacement(a) = bO(a) - bE(a) = 2 a i^{-1}(zetaE - A).
## zetaE - A is evaluated from tau, not by subtraction. Also returned: the
## fixed-design quantities computed before generation, namely the
## probabilities, event rate, eigenvalues and condition number of i(beta0), and
## the rank of X.
first_order_bias <- function(X, beta0, link = "probit", m = NULL, a = 0.5) {
  X <- as.matrix(X)
  L <- if (is.character(link)) make_link(link) else link
  if (is.null(m)) m <- rep(1, nrow(X))
  b <- L$base(drop(X %*% beta0))
  kexp <- m * exp(L$lomega(b)$l)
  i <- crossprod(X, kexp * X)
  iinv <- solve(i)
  h <- kexp * rowSums((X %*% iinv) * X)
  A <- 0.5 * drop(crossprod(X, h * b$d1))
  zmA <- 0.5 * drop(crossprod(X, h * L$tau(b)))
  zeta <- A + zmA
  pr <- L$prob(drop(X %*% beta0))
  isp <- info_spectrum(i)
  list(i = i, i_inv = iinv, A = A, zetaE = zeta, zetaE_minus_A = zmA,
       b_ML = -drop(iinv %*% A),
       b_E = drop(iinv %*% (2 * a * zeta - A)),
       b_O = drop(iinv %*% (2 * a * (2 * zeta - A) - A)),
       displacement = 2 * a * drop(iinv %*% zmA),
       a = a, prob = pr, event_rate = sum(m * pr) / sum(m),
       info_eigenvalues = isp$values, info_cond = isp$cond, rank = qr(X)$rank)
}

## Control variate c = i(beta0)^{-1} S(beta0); E(c) = 0, var(c) = i(beta0)^{-1}.
control_variate <- function(model, beta0) {
  inf <- bp_info(model, beta0)
  drop(solve(inf$expected, inf$score))
}

## ---------------------------------------------------------------------------
## Intercept-only models: n Bernoulli observations with common predictor beta
## are equivalent to one binomial observation with total n and t successes.
## Criterion and first two derivatives, vectorized over beta.
intercept_criterion <- function(beta, t, n, link, type = c("E", "O", "ML"),
                                a = 0.5, deriv = 2L) {
  type <- match.arg(type)
  L <- if (is.character(link)) make_link(link) else link
  b <- L$base(beta)
  f <- cmul(t, b$lG) + cmul(n - t, b$lGb)
  pd <- if (type == "O") deriv else 0L
  p1 <- L$psi1(b, pd)
  p0 <- L$psi0(b, pd)
  K <- cmul(t, p1$v) + cmul(n - t, p0$v)              # observed information
  f1 <- cmul(t, b$r1) - cmul(n - t, b$r0)
  f2 <- -K
  if (type == "E") {
    lo <- L$lomega(b, p1$v, p0$v, deriv = 2L)
    f <- f + a * (log(n) + lo$l)
    f1 <- f1 + a * lo$d1
    f2 <- f2 + a * lo$d2
  } else if (type == "O") {
    K1 <- cmul(t, p1$d1) + cmul(n - t, p0$d1)
    K2 <- cmul(t, p1$d2) + cmul(n - t, p0$d2)
    pos <- K > 0
    f <- ifelse(pos, f + a * log(ifelse(pos, K, 1)), -Inf)
    f1 <- f1 + a * K1 / K
    f2 <- f2 + a * (K2 / K - (K1 / K)^2)
  }
  list(f = f, f1 = f1, f2 = f2, info = n * exp(L$lomega(b)$l))
}

## Global maximization in one dimension: grid on [lo, hi] (extended while the
## best grid point lies on the boundary), golden-section refinement of every
## grid-local maximum within its two neighbouring cells, then Newton polishing
## on the derivative. Returns the highest refined maximum and the number of
## distinct local maxima found.
fit_intercept_only <- function(t, n, link, type = "E", a = 0.5, lo = -15, hi = 15,
                               step = 0.02, gtol = 1e-10, max_extend = 6L) {
  L <- if (is.character(link)) make_link(link) else link
  crit <- function(bb, d = 0L) intercept_criterion(bb, t, n, L, type, a, d)
  for (ext in 0:max_extend) {
    grid <- seq(lo, hi, by = step)
    fv <- crit(grid, 0L)$f
    fv[!is.finite(fv)] <- -Inf
    k <- which.max(fv)
    if (k > 1L && k < length(grid)) break
    w <- hi - lo
    if (k == 1L) lo <- lo - w else hi <- hi + w
  }
  ng <- length(grid)
  locmax <- which(fv[2:(ng - 1)] >= fv[1:(ng - 2)] & fv[2:(ng - 1)] >= fv[3:ng] &
                    is.finite(fv[2:(ng - 1)])) + 1L
  gr <- (sqrt(5) - 1) / 2
  sols <- lapply(locmax, function(k) {
    x1 <- grid[k - 1L]; x4 <- grid[k + 1L]
    x2 <- x4 - gr * (x4 - x1); x3 <- x1 + gr * (x4 - x1)
    f2 <- crit(x2)$f; f3 <- crit(x3)$f
    for (it in 1:60) {
      if (f2 >= f3) { x4 <- x3; x3 <- x2; f3 <- f2; x2 <- x4 - gr * (x4 - x1); f2 <- crit(x2)$f }
      else { x1 <- x2; x2 <- x3; f2 <- f3; x3 <- x1 + gr * (x4 - x1); f3 <- crit(x3)$f }
    }
    x <- (x1 + x4) / 2
    cv <- crit(x, 2L)
    for (it in 1:20) {
      if (!isTRUE(cv$f2 < 0)) break
      std <- abs(cv$f1) / sqrt(cv$info)
      if (!is.finite(std) || std < gtol) break
      xn <- x - cv$f1 / cv$f2
      cn <- crit(xn, 2L)
      if (!is.finite(cn$f) || !isTRUE(abs(cn$f1) / sqrt(cn$info) < std)) break
      x <- xn; cv <- cn
    }
    c(beta = x, value = cv$f, std_grad = abs(cv$f1) / sqrt(cv$info), f2 = cv$f2)
  })
  S <- do.call(rbind, sols)
  best <- which.max(S[, "value"])
  nd <- length(distinct_solutions(as.list(S[, "beta"])))
  list(beta = unname(S[best, "beta"]), value = unname(S[best, "value"]),
       std_grad = unname(S[best, "std_grad"]), f2 = unname(S[best, "f2"]),
       n_local = nd, all = S)
}

## Exact mean of the selected maximizer in the intercept-only model with n
## Bernoulli observations and true predictor eta0, by enumeration of the
## response total t = 0, ..., n. Returns the exact bias, the estimates, the
## binomial probabilities and diagnostics.
exact_bias_intercept <- function(link, eta0, n, type = c("E", "O"), a = 0.5, ...) {
  type <- match.arg(type)
  L <- if (is.character(link)) make_link(link) else link
  pr <- L$prob(eta0)
  w <- dbinom(0:n, n, pr)
  fits <- lapply(0:n, function(t) fit_intercept_only(t, n, L, type, a, ...))
  est <- vapply(fits, function(f) f$beta, 0)
  list(bias = sum(w * est) - eta0, estimates = est, weights = w,
       max_std_grad = max(vapply(fits, function(f) f$std_grad, 0)),
       max_f2 = max(vapply(fits, function(f) f$f2, 0)),
       max_local = max(vapply(fits, function(f) f$n_local, 0L)))
}


## =============================================================================
## 7. Separation by linear programming
## =============================================================================
## Detection of complete and quasi-complete separation by linear programming
## (Konis, 2007), with a dense-tableau primal simplex method in base R.

## Maximize c'x subject to A x <= b, x >= 0, with b >= 0 (the slack basis is
## feasible). Bland's rule prevents cycling under the heavy degeneracy of the
## separation programs. Returns the primal solution, the optimal value and the
## dual vector y >= 0 (A'y >= c, b'y = optimal value at optimality).
lp_simplex_max <- function(cvec, A, b, tol = 1e-9, maxit = 100000L) {
  A <- as.matrix(A)
  m <- nrow(A); nv <- ncol(A)
  if (any(b < 0)) stop("b must be non-negative")
  nc <- nv + m + 1L
  Tab <- cbind(A, diag(m), b)
  z <- c(-cvec, numeric(m), 0)
  basis <- nv + seq_len(m)
  status <- "iteration_limit"
  for (it in seq_len(maxit)) {
    enter <- which(z[-nc] < -tol)
    if (!length(enter)) { status <- "optimal"; break }
    j <- enter[1L]
    colj <- Tab[, j]
    pos <- which(colj > tol)
    if (!length(pos)) { status <- "unbounded"; break }
    ratio <- Tab[pos, nc] / colj[pos]
    rmin <- min(ratio)
    cand <- pos[ratio <= rmin + tol * (1 + abs(rmin))]
    i <- cand[which.min(basis[cand])]
    Tab[i, ] <- Tab[i, ] / Tab[i, j]
    f <- Tab[, j]
    f[i] <- 0
    Tab <- Tab - outer(f, Tab[i, ])
    z <- z - z[j] * Tab[i, ]
    basis[i] <- j
  }
  xall <- numeric(nv + m)
  xall[basis] <- Tab[, nc]
  list(status = status, x = xall[seq_len(nv)], value = z[nc],
       dual = z[nv + seq_len(m)], iterations = it)
}

## Rows a_i: x_i for every observed success and -x_i for every observed
## failure (for binary data a_i = (2 y_i - 1) x_i).
separation_rows <- function(X, y, m = rep(1, length(y))) {
  X <- as.matrix(X)
  rbind(X[y > 0, , drop = FALSE], -X[m - y > 0, , drop = FALSE])
}

## LP1: maximize d s.t. a_i'gamma >= d, ||gamma||_1 <= 1   (d* > 0: complete)
## LP2: maximize sum_i a_i'gamma s.t. a_i'gamma >= 0, ||gamma||_1 <= 1
##      (positive optimum with d* = 0: quasi-complete; zero: overlap).
## gamma = gamma+ - gamma-, gamma+-, d >= 0 (d* >= 0 because gamma = 0 is
## feasible). Certificates are rechecked from the returned primal and dual
## vectors: for LP1, lambda >= 0, sum(lambda) >= 1, |sum lambda_i a_i| <= mu;
## for LP2, lambda >= 0, |sum (1 + lambda_i) a_i| <= mu; the duality gap is
## mu minus the primal value.
sep_lp <- function(Arows, tol = 1e-9) {
  N <- nrow(Arows); p <- ncol(Arows)
  A1 <- rbind(cbind(-Arows, Arows, 1), c(rep(1, 2 * p), 0))
  b1 <- c(numeric(N), 1)
  r1 <- lp_simplex_max(c(numeric(2 * p), 1), A1, b1, tol)
  g1 <- r1$x[1:p] - r1$x[p + 1:p]
  d1 <- r1$x[2 * p + 1]
  lam1 <- r1$dual[1:N]; mu1 <- r1$dual[N + 1]
  cert1 <- c(primal_viol = max(0, max(d1 - drop(Arows %*% g1)), sum(abs(g1)) - 1),
             dual_viol = max(0, -min(lam1), -mu1, 1 - sum(lam1),
                             max(abs(drop(crossprod(Arows, lam1)))) - mu1),
             gap = mu1 - d1)
  s <- colSums(Arows)
  A2 <- rbind(cbind(-Arows, Arows), rep(1, 2 * p))
  r2 <- lp_simplex_max(c(s, -s), A2, b1, tol)
  g2 <- r2$x[1:p] - r2$x[p + 1:p]
  v2 <- sum(s * g2)
  lam2 <- r2$dual[1:N]; mu2 <- r2$dual[N + 1]
  cert2 <- c(primal_viol = max(0, -min(drop(Arows %*% g2)), sum(abs(g2)) - 1),
             dual_viol = max(0, -min(lam2), -mu2,
                             max(abs(drop(crossprod(Arows, 1 + lam2)))) - mu2),
             gap = mu2 - v2)
  list(lp1 = list(value = d1, gamma = g1, dual = c(lam1, mu1), status = r1$status,
                  certificate = cert1, iterations = r1$iterations),
       lp2 = list(value = v2, gamma = g2, dual = c(lam2, mu2), status = r2$status,
                  certificate = cert2, iterations = r2$iterations))
}

## Classification: "complete", "quasi-complete", "overlap" or "ambiguous".
## The separating directions are returned in the original coordinates
## (scaled$lp1$direction for complete, scaled$lp2$direction otherwise).
## Values above pos_tol are positive and below zero_tol are zero; values in
## between are rechecked on the unscaled rows (the primary solve uses columns
## scaled to unit maximum absolute value, which leaves separation unchanged).
separation_check <- function(X, y, m = rep(1, length(y)), pos_tol = 1e-8,
                             zero_tol = 1e-10, lp_tol = 1e-9) {
  Ar <- separation_rows(X, y, m)
  sc <- 1 / apply(abs(Ar), 2, max)
  classify <- function(res) {
    v1 <- res$lp1$value; v2 <- res$lp2$value
    if (v1 > pos_tol) return("complete")
    if (v1 >= zero_tol) return("ambiguous")
    if (v2 > pos_tol) return("quasi-complete")
    if (v2 < zero_tol) return("overlap")
    "ambiguous"
  }
  res <- sep_lp(Ar * rep(sc, each = nrow(Ar)), lp_tol)
  res$lp1$direction <- sc * res$lp1$gamma      # in the original coordinates
  res$lp2$direction <- sc * res$lp2$gamma
  cls <- classify(res)
  rechecked <- NULL
  if (cls == "ambiguous") {
    rechecked <- sep_lp(Ar, lp_tol)
    rechecked$lp1$direction <- rechecked$lp1$gamma
    rechecked$lp2$direction <- rechecked$lp2$gamma
    cls <- classify(rechecked)
  }
  cert_ok <- all(c(res$lp1$certificate[1:2], res$lp2$certificate[1:2]) <= lp_tol) &&
    all(abs(c(res$lp1$certificate["gap"], res$lp2$certificate["gap"])) <= lp_tol)
  list(category = cls, scaled = res, rechecked = rechecked, column_scale = sc,
       certificates_ok = cert_ok)
}


## =============================================================================
## 8. Penalized-profile confidence sets
## =============================================================================
## C = {t : L(t) <= q} for one coefficient beta_idx, with
## L(t) = 2 {Q(beta_hat) - max Q(beta) subject to beta_idx = t} and q the
## chi-square(1) quantile at `level`, reconstructed on u = (t - beta_hat_idx)/s,
## s^2 = {i(beta_hat)^{-1}}_{idx,idx}, h = control$grid_h:
##  1. grid of spacing h on [-W0, W0] plus the projections of the other
##     computed unrestricted maxima, evaluated outwards from the estimate, then
##     an inward sweep;
##  2. extensions to each W in control$extensions, each band evaluated outwards
##     and followed by an inward sweep, until on two successive extensions L
##     exceeds q throughout the new bands and increases outwards at their outer
##     points;
##  3. a bounded scalar minimization around every sampled local minimum above q
##     on the current grid, inserting any new below-q minimum before
##     refinement; then for each refinement to h/2, h/4, ...: refinement
##     points, the same search, every crossing refined by bisection to width
##     bisect_width (1 + |t|) with |L - q| < bisect_resid, and stability on two
##     successive grids (same component count, endpoint movement below
##     stab_endpoint (1 + |t|), relative length change below stab_length).
## Starts at a new point: the accepted solutions at the nearest evaluated
## points on both sides and the projections of every distinct accepted
## unrestricted maximum (added by profile_stat()), with, for starts_policy =
## "protocol", the projections of `starts` and the rescue starts. Points
## evaluated outwards have no evaluated outer neighbour: the inward sweep
## refits each from its outer neighbour's solution (when it is one the point
## has not yet been started from) and keeps the selection over both fits, so
## every point receives neighbouring solutions from both directions.
## Refinement, search and bisection points lie between evaluated points and
## receive both at once. An unresolved constrained fit makes the set
## unresolved. Under gradient_scale = "estimate", the nuisance gradient is
## standardized using the free block of expected information at the unrestricted
## estimate. "candidate" evaluates the block at each constrained candidate.
## Possible statuses are "resolved", "tail_check_failed", "unstable",
## "crossing_unresolved", "profile_fit_unresolved", "unrestricted_unresolved",
## "unrestricted_missing".
## Point records: points (u, L, constrained objective Qc, inward-refit gain,
## constrained solution b1, ..., bp) and points_int (resolved, profile_stat
## flag, kind, lower and upper neighbour solution used, inward refits,
## improvement beyond the tie tolerance), with the codes below.
PROFILE_KINDS <- c("grid", "projection", "extension", "refinement", "search", "bisection")
PROFILE_FLAGS <- c("ok", "constrained_unresolved", "unrestricted_unresolved")
code_of <- function(x, levels) if (is.null(x)) 0L else { i <- match(x, levels); if (is.na(i)) 0L else i }

profile_control <- function(...) {
  ctl <- list(maxit = 500L, maxit_rescue = 2000L, grid_h = 0.25, W0 = 8, extensions = c(16, 32, 64),
              n_refine = 3L, bisect_maxit = 200L, bisect_width = 1e-6, bisect_resid = 1e-6,
              stab_endpoint = 1e-4, stab_length = 1e-3, optimize_tol = .Machine$double.eps^0.25,
              gradient_scale = "estimate")
  a <- list(...)
  if (length(a) && (is.null(names(a)) || !all(names(a) %in% names(ctl))))
    stop("unknown profile_control() setting")
  ctl[names(a)] <- a
  ctl
}

profile_set <- function(model, fit, idx, starts, rescue = NULL, level = 0.95,
                        control = profile_control(), starts_policy = c("protocol", "warm")) {
  starts_policy <- match.arg(starts_policy)
  if (is.null(fit$selected)) return(list(status = "unrestricted_missing"))
  num <- control; h <- num$grid_h; cut <- qchisq(level, 1)
  b <- fit$selected$beta
  ib <- bp_info(model, b)$expected
  s <- sqrt(solve(ib)[idx, idx])
  iscale <- if (identical(num$gradient_scale, "estimate")) ib[-idx, -idx, drop = FALSE]
  base <- if (starts_policy == "protocol") as.matrix(starts)
  resc <- if (starts_policy == "protocol") rescue
  nosol <- rep(NA_real_, length(b))
  gu <- gL <- gQ <- ggain <- numeric(0); gok <- gimp <- logical(0)
  gflag <- gkind <- glo <- ghi <- ginw <- integer(0)
  gsol <- used_up <- used_lo <- list()
  flag_unres <- FALSE; n_fits <- 0L; W <- num$W0; tails <- NULL; stab <- list()
  searches <- matrix(numeric(0), 0, 4, dimnames = list(NULL, c("lower", "upper", "u", "L")))
  crossings <- matrix(numeric(0), 0, 3, dimnames = list(NULL, c("t", "iterations", "residual")))
  profile_at <- function(u, st, rs) {
    n_fits <<- n_fits + 1L
    r <- profile_stat(model, fit, idx, b[idx] + s * u, st, rs, maxit = num$maxit, maxit_rescue = num$maxit_rescue,
                      info_scale = iscale)
    if (identical(r$flag, "unrestricted_unresolved")) flag_unres <<- TRUE
    r
  }
  Lfun <- function(u, kind) {
    k <- which(abs(gu - u) < 1e-12)
    if (length(k)) return(gL[k[1]])
    lo <- which(gu < u); hi <- which(gu > u)
    klo <- if (length(lo)) lo[which.max(gu[lo])]; khi <- if (length(hi)) hi[which.min(gu[hi])]
    slo <- if (!is.null(klo) && !anyNA(gsol[[klo]])) gsol[[klo]]
    shi <- if (!is.null(khi) && !anyNA(gsol[[khi]])) gsol[[khi]]
    r <- profile_at(u, rbind(base, slo, shi, b), resc)
    ok <- identical(r$flag, "ok") && is.finite(r$stat)
    sel <- r$fit$selected
    gu <<- c(gu, u); gL <<- c(gL, if (ok) r$stat else NA_real_); gok <<- c(gok, ok)
    gQ <<- c(gQ, if (is.null(sel)) NA_real_ else sel$value)
    gflag <<- c(gflag, code_of(r$flag, PROFILE_FLAGS)); gkind <<- c(gkind, code_of(kind, PROFILE_KINDS))
    glo <<- c(glo, as.integer(!is.null(slo))); ghi <<- c(ghi, as.integer(!is.null(shi)))
    ginw <<- c(ginw, 0L); gimp <<- c(gimp, FALSE); ggain <<- c(ggain, 0)
    i <- length(gu)
    gsol[[i]] <<- if (is.null(sel)) nosol else sel$beta
    used_up[[i]] <<- if (is.null(shi)) nosol else shi
    used_lo[[i]] <<- if (is.null(slo)) nosol else slo
    gL[i]
  }
  inward <- function() {
    for (side in c(1, -1)) {
      ix <- which(side * gu >= 0)
      ix <- ix[order(side * gu[ix], decreasing = TRUE)]          # outermost first
      if (length(ix) < 2L) next
      for (j in 2:length(ix)) {
        k <- ix[j]; so <- gsol[[ix[j - 1L]]]
        if (anyNA(so) || identical(if (side > 0) used_up[[k]] else used_lo[[k]], so)) next
        if (side > 0) { used_up[[k]] <<- so; ghi[k] <<- 1L } else { used_lo[[k]] <<- so; glo[k] <<- 1L }
        ginw[k] <<- ginw[k] + 1L
        r <- profile_at(gu[k], rbind(so), NULL)
        sel <- r$fit$selected
        if (is.null(sel)) next
        cur <- if (anyNA(gsol[[k]])) NULL else list(accepted = TRUE, value = gQ[k], beta = gsol[[k]])
        pick <- select_maximum(list(cur, sel))
        if (!identical(pick$beta, gsol[[k]])) {     # a higher maximum, or a tie decided by norm/lexicographic order
          ok <- identical(r$flag, "ok") && is.finite(r$stat)
          gain <- if (is.null(cur)) NA_real_ else sel$value - gQ[k]
          gimp[k] <<- is.null(cur) || gain > 1e-10 * (1 + abs(gQ[k]))
          ggain[k] <<- gain
          gL[k] <<- if (ok) r$stat else NA_real_; gok[k] <<- ok; gQ[k] <<- sel$value
          gflag[k] <<- code_of(r$flag, PROFILE_FLAGS); gsol[[k]] <<- sel$beta
        }
      }
    }
  }
  check <- function() if (flag_unres) "unrestricted_unresolved" else if (!all(gok)) "profile_fit_unresolved" else NULL
  done <- function(status, comp = NULL) {
    o <- order(gu); sol <- do.call(rbind, gsol)
    colnames(sol) <- paste0("b", seq_len(ncol(sol)))
    pts <- cbind(u = gu, L = gL, Qc = gQ, inward_gain = ggain, sol)[o, , drop = FALSE]
    pts_int <- cbind(ok = as.integer(gok), flag = gflag, kind = gkind, lower_nb = glo, upper_nb = ghi,
                     inward_refits = ginw, improved = as.integer(gimp))[o, , drop = FALSE]
    c(list(status = status, idx = idx, level = level, cutoff = cut, s = s, t_hat = b[idx], W = W,
           tail_checks = tails, points = pts, points_int = pts_int, searches = searches,
           crossings = crossings, stability = stab, n_evaluations = length(gu), n_profile_fits = n_fits),
      if (!is.null(comp)) list(components = comp, length = sum(comp[, 2] - comp[, 1]),
                               n_components = nrow(comp)))
  }
  search_minima <- function() {
    o <- order(gu); U <- gu[o]; LL <- gL[o]; K <- length(U)
    if (K < 3L) return(invisible())
    i <- 2:(K - 1)
    for (k in i[LL[i] <= LL[i - 1] & LL[i] <= LL[i + 1] & LL[i] > cut]) {
      r <- tryCatch(optimize(function(u) {
        v <- Lfun(u, "search")
        if (!is.finite(v)) stop("unresolved profile point")       # stops the search at once
        v
      }, c(U[k - 1], U[k + 1]), tol = num$optimize_tol), error = function(e) NULL)
      searches <<- rbind(searches, c(U[k - 1], U[k + 1], if (is.null(r)) c(NA, NA) else c(r$minimum, r$objective)))
      if (is.null(r)) return(invisible())                          # reported by check()
    }
  }
  crossing <- function(a, c2) {
    fa <- Lfun(a, "bisection") - cut
    fm <- NA_real_
    for (it in seq_len(num$bisect_maxit)) {
      m <- (a + c2) / 2; fm <- Lfun(m, "bisection") - cut
      if (is.na(fm)) break
      tm <- b[idx] + s * m
      if (s * abs(c2 - a) <= num$bisect_width * (1 + abs(tm)) && abs(fm) < num$bisect_resid) {
        crossings <<- rbind(crossings, c(tm, it, fm))
        return(tm)
      }
      if (sign(fm) == sign(fa)) { a <- m; fa <- fm } else c2 <- m
    }
    crossings <<- rbind(crossings, c(NA, it, fm))
    NA_real_
  }
  ## 1. initial grid, evaluated outwards, then inwards; the selected maximum's
  ##    own group in fit$distinct projects to u = 0 and is not evaluated again
  us <- seq(h, W, by = h)
  for (u in c(0, as.vector(rbind(us, -us)))) Lfun(u, "grid")
  for (d in fit$distinct) if (!same_solution(d$beta, b)) Lfun((d$beta[idx] - b[idx]) / s, "projection")
  inward()
  if (!is.null(st <- check())) return(done(st))
  ## 2. tail extensions
  prev <- FALSE; tail_ok <- FALSE
  for (Wn in num$extensions) {
    band <- seq(W + h, Wn, by = h)
    for (u in as.vector(rbind(band, -band))) Lfun(u, "extension")
    inward()
    if (!is.null(st <- check())) return(done(st))
    inb <- abs(gu) > W + 1e-12 & abs(gu) <= Wn + 1e-12
    above <- all(gL[inb] > cut)
    inc_pos <- Lfun(Wn, "extension") > Lfun(Wn - h, "extension")
    inc_neg <- Lfun(-Wn, "extension") > Lfun(-Wn + h, "extension")
    pass <- above && inc_pos && inc_neg
    tails <- rbind(tails, c(W = Wn, above = above, increasing_upper = inc_pos,
                            increasing_lower = inc_neg, pass = pass))
    W <- Wn
    if (pass && prev) { tail_ok <- TRUE; break }
    prev <- pass
  }
  if (!tail_ok) return(done("tail_check_failed"))
  ## 3. minimum search before refinement; refinements, crossings, stability
  search_minima()
  if (!is.null(st <- check())) return(done(st))
  previous <- NULL
  for (j in seq_len(num$n_refine)) {
    hj <- h / 2^j
    for (u in seq(-W + hj, W - hj, by = 2 * hj)) Lfun(u, "refinement")
    if (!is.null(st <- check())) return(done(st))
    search_minima()
    if (!is.null(st <- check())) return(done(st))
    o <- order(gu); U <- gu[o]; below <- gL[o] <= cut
    k <- which(diff(below) != 0)
    ends <- vapply(k, function(i) crossing(U[i], U[i + 1]), 0)
    if (!is.null(st <- check()) || anyNA(ends)) return(done(if (is.null(st)) "crossing_unresolved" else st))
    comp <- matrix(ends, ncol = 2, byrow = TRUE)            # below at both outer ends is excluded by step 2
    stab[[j]] <- comp
    if (!is.null(previous) && nrow(comp) == nrow(previous) &&
        all(abs(comp - previous) <= num$stab_endpoint * (1 + abs(comp))) &&
        abs(sum(comp[, 2] - comp[, 1]) - sum(previous[, 2] - previous[, 1])) /
          (1 + sum(previous[, 2] - previous[, 1])) < num$stab_length)
      return(done("resolved", comp))
    previous <- comp
  }
  done("unstable")
}


## =============================================================================
## 9. Independent re-evaluation by singular value decompositions
## =============================================================================
## Recomputes the quantities that acceptance depends on from singular value
## decompositions of diag(sqrt(k)) X instead of Cholesky factors of
## X' diag(k) X: the objective (log det from the singular values), gradient and
## Hessian of the penalty (leverages from Z = X V D^{-1}, Z Z' = X K^{-1} X'),
## the standardized gradient with i^{-1/2} = V D^{-1} V' and the curvature
## check; for roots (kind "A"), the adjusted score, its standardized residual
## and the singular values of its Jacobian. Only the link-level functions are
## shared with bp_eval() and bp_adjscore(). Use it when i or j is badly
## conditioned (condition number above 1e12).
svd_factor <- function(X, k) {
  if (!all(is.finite(k)) || any(k < 0)) return(NULL)
  s <- svd(sqrt(k) * X)
  if (!(min(s$d) > 0)) return(NULL)
  Z <- X %*% (s$v %*% diag(1 / s$d, length(s$d)))
  list(logdet = 2 * sum(log(s$d)), cond = (max(s$d) / min(s$d))^2,
       isqrt = s$v %*% (t(s$v) / s$d), Z = Z, q = rowSums(Z * Z))
}
stable_reevaluation <- function(model, beta, kind = c("E", "O", "ML", "A"), a = 0.5) {
  kind <- match.arg(kind)
  X <- model$X; y <- model$y; n0 <- model$n0; L <- model$link
  b <- L$base(drop(X %*% beta))
  p1 <- L$psi1(b, 2L); p0 <- L$psi0(b, 2L)
  w <- cmul(y, p1$v) + cmul(n0, p0$v)
  lo <- L$lomega(b, p1$v, p0$v, 2L)
  k <- model$m * exp(lo$l)
  fi <- svd_factor(X, k)
  if (is.null(fi)) return(list(ok = FALSE, reason = "SVD of diag(sqrt(omega)) X failed"))
  fj <- svd_factor(X, w)
  u <- cmul(y, b$r1) - cmul(n0, b$r0)
  if (kind == "A") {
    h <- k * fi$q
    U <- drop(crossprod(X, u)) + 0.5 * drop(crossprod(X, cmul(h, b$d1)))
    k1 <- cmul(k, lo$d1)
    Q <- tcrossprod(fi$Z)
    dh <- (k1 * fi$q) * X - k * ((Q * Q) %*% (k1 * X))
    d1h <- b$d1; d1h[k == 0] <- 0
    J <- -crossprod(X, w * X) + 0.5 * (crossprod(X, cmul(h, b$d2) * X) + crossprod(X, d1h * dh))
    sv <- svd(J, nu = 0, nv = 0)$d
    rmax <- max(abs(fi$isqrt %*% U))
    return(list(ok = TRUE, std_resid_max = rmax, jac_singular_values = sv,
                accept = rmax < 1e-8 && min(sv) > 1e-10 * max(1, max(sv)),
                cond_i = fi$cond, cond_j = if (is.null(fj)) NA_real_ else fj$cond))
  }
  val <- sum(cmul(y, b$lG)) + sum(cmul(n0, b$lGb))
  g <- drop(crossprod(X, u)); H <- -crossprod(X, w * X)
  if (kind %in% c("E", "O")) {
    if (kind == "E") { fK <- fi; c1 <- cmul(k, lo$d1); c2 <- cmul(k, lo$d2 + lo$d1^2) } else {
      if (is.null(fj)) return(list(ok = FALSE, reason = "SVD of diag(sqrt(w)) X failed"))
      fK <- fj; c1 <- cmul(y, p1$d1) + cmul(n0, p0$d1); c2 <- cmul(y, p1$d2) + cmul(n0, p0$d2)
    }
    Q <- tcrossprod(fK$Z)
    val <- val + a * fK$logdet
    g <- g + a * drop(crossprod(X, c1 * fK$q))
    H <- H + a * (crossprod(X, (c2 * fK$q) * X) - crossprod(c1 * X, (Q * Q) %*% (c1 * X)))
  }
  H <- (H + t(H)) / 2
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  gmax <- max(abs(fi$isqrt %*% g))
  list(ok = TRUE, value = val, std_grad_max = gmax, lambda_max = max(ev), hess_norm = max(abs(ev)),
       accept = gmax < 1e-8 && max(ev) < -1e-10 * max(1, max(abs(ev))),
       cond_i = fi$cond, cond_j = if (is.null(fj)) NA_real_ else fj$cond)
}


## =============================================================================
## 10. Data analysis interface
## =============================================================================
## obsinfo_fit(formula, data, link, method) fits a binomial-response GLM by
##   method "O"  l + a log det j   (observed-information penalty; a = 1/2)
##          "E"  l + a log det i   (expected-information penalty; Jeffreys-type)
##          "A"  Firth's adjusted score S + A = 0 (mean bias reduction)
##          "ML" maximum likelihood (overlap data only)
## with the protocol's multistart search, rescue sequence, acceptance checks
## and selection rules. The response is binary (0/1, logical or a two-level
## factor whose second level is the event) or cbind(successes, failures).
## With scale = TRUE (default) every non-constant column of the model matrix is
## divided by its standard deviation before fitting and the estimates are
## transformed back: the four estimators are equivariant under this
## reparametrization (the penalties change by a constant), so the estimates
## are unchanged and only the start set, defined on the scaled coefficients,
## is adapted to the units of the covariates. Starts: 0, +-2 e_r and 2 Z_k,
## k = 1, ..., 10, with Z a 50 x p matrix of N(0, 1) draws from `seed` (the
## random-number state of the session is restored); rescue starts 4 Z_k,
## k = 11, ..., 50. For "O" and "A" the "E" estimate is an additional start,
## and the reference for selecting among adjusted-score roots.
## Standard errors are from the inverse expected information at the estimate.
## confint() gives penalized-profile (or likelihood-ratio, for "ML") sets by
## profile_set(), or Wald intervals; profile sets are not available for "A".
with_seed <- function(seed, expr) {
  genv <- globalenv()
  old_seed <- if (exists(".Random.seed", envir = genv, inherits = FALSE)) get(".Random.seed", envir = genv)
  old_kind <- RNGkind()
  on.exit({
    RNGkind(old_kind[1], old_kind[2], old_kind[3])
    if (is.null(old_seed)) { if (exists(".Random.seed", envir = genv, inherits = FALSE)) rm(".Random.seed", envir = genv) }
    else assign(".Random.seed", old_seed, envir = genv)
  })
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  expr
}

obsinfo_fit <- function(formula, data, link = c("probit", "logit", "cloglog", "loglog", "cauchit"),
                        method = c("O", "E", "A", "ML"), a = 0.5, scale = TRUE, seed = 1L,
                        maxit = 500L, maxit_rescue = 2000L) {
  link <- match.arg(link); method <- match.arg(method)
  cl <- match.call()
  mf <- model.frame(formula, data, na.action = na.omit)
  tt <- attr(mf, "terms")
  X <- model.matrix(tt, mf)
  Y <- model.response(mf)
  if (is.matrix(Y)) {
    if (ncol(Y) != 2L) stop("a matrix response must be cbind(successes, failures)")
    y <- as.numeric(Y[, 1]); m <- as.numeric(Y[, 1] + Y[, 2])
  } else {
    y <- if (is.factor(Y)) {
      if (nlevels(Y) != 2L) stop("a factor response must have two levels")
      as.numeric(Y != levels(Y)[1])
    } else as.numeric(Y)
    if (!all(y %in% c(0, 1))) stop("the response must be binary (0/1, logical or a two-level factor) or cbind(successes, failures)")
    m <- rep(1, length(y))
  }
  keep <- m > 0                                         # rows with no trials contribute nothing
  X <- X[keep, , drop = FALSE]; y <- y[keep]; m <- m[keep]
  p <- ncol(X); cn <- colnames(X)
  if (qr(X)$rank < p) stop("the model matrix does not have full column rank")
  d <- if (scale) apply(X, 2, function(x) { s <- sd(x); if (is.finite(s) && s > 0) s else 1 }) else rep(1, p)
  Xs <- sweep(X, 2, d, "/")
  mod <- bp_model(Xs, y, m, link)
  sep <- separation_check(X, y, m)
  if (method == "ML" && sep$category %in% c("complete", "quasi-complete"))
    stop("the maximum likelihood estimate does not exist: the data show ", sep$category, " separation")
  st <- protocol_starts(p, with_seed(seed, matrix(rnorm(50 * p), 50, p)))
  dir <- separation_direction(sep)                      # direction on the scaled coefficients
  if (!is.null(dir)) dir <- dir * d
  fE <- if (method != "ML") fit_criterion_search(mod, "E", if (method == "E") a else 0.5, st$routine, st$rescue, dir,
                                                 maxit = maxit, maxit_rescue = maxit_rescue)
  refE <- if (!is.null(fE$selected)) fE$selected$beta
  routine <- if (method %in% c("O", "A")) rbind(st$routine, refE) else st$routine
  f <- switch(method, E = fE,
              O = fit_criterion_search(mod, "O", a, routine, st$rescue, dir, maxit = maxit, maxit_rescue = maxit_rescue),
              ML = fit_criterion_protocol(mod, "ML", a, routine, st$rescue, maxit = maxit, maxit_rescue = maxit_rescue),
              A = fit_adjscore_protocol(mod, routine, st$rescue, ref = refE, maxit = maxit, maxit_rescue = maxit_rescue))
  sel <- f$selected
  resolved <- !is.null(sel)
  gam <- if (resolved) sel$beta else rep(NA_real_, p)
  beta <- setNames(gam / d, cn)
  vc <- matrix(NA_real_, p, p, dimnames = list(cn, cn))
  diag_info <- list(stage = f$stage, n_accepted = f$n_accepted, n_distinct = f$n_distinct,
                    values_distinct = if (method != "A") vapply(f$distinct, function(s) s$value, 0))
  ree <- NULL
  if (resolved) {
    b <- mod$link$base(drop(Xs %*% gam))
    fi <- svd_factor(Xs, mod$m * exp(mod$link$lomega(b, NULL, NULL, 0L)$l))
    vc[] <- (fi$isqrt %*% fi$isqrt) / outer(d, d)
    ej <- eigen(bp_info(mod, gam)$observed, symmetric = TRUE, only.values = TRUE)$values
    diag_info$cond_i <- fi$cond
    diag_info$cond_j <- if (min(ej) > 0) max(ej) / min(ej) else Inf
    if (method == "A") {
      diag_info$std_resid_max <- sel$std_resid_max; diag_info$jac_singular_values <- sel$jac_singular_values
    } else {
      diag_info$std_grad_max <- sel$std_grad_max; diag_info$lambda_max <- sel$lambda_max
    }
    if (diag_info$cond_i > 1e12 || diag_info$cond_j > 1e12) {
      ree <- stable_reevaluation(mod, gam, method, a)
      diag_info$reevaluation <- if (!ree$ok) "failed" else if (ree$accept) "agree" else "disagree"
    }
  } else warning("no accepted ", if (method == "A") "root" else "maximum", " was found; see $fit for the candidates")
  shift <- if (method %in% c("E", "O")) 2 * a * sum(log(d)) else 0     # Q(beta) = Q(gamma) + shift
  structure(list(call = cl, terms = tt, xlevels = .getXlevels(tt, mf), link = link, method = method, a = a,
                 coefficients = beta, vcov = vc, se = sqrt(diag(vc)), resolved = resolved,
                 value = if (resolved && method != "A") sel$value + shift else NA_real_,
                 loglik = if (resolved) bp_eval(mod, gam, "ML", a, 0L)$value else NA_real_,
                 separation = sep$category, separation_detail = sep, diagnostics = diag_info, reevaluation = ree,
                 n = length(y), trials = sum(m), y = y, m = m, scale = d, seed = seed,
                 model = mod, fit = f, starts = rbind(routine, f$ray_starts), rescue = st$rescue),
            class = "obsinfo")
}

obsinfo_method_label <- function(x) switch(x$method,
  O = sprintf("penalized likelihood l + %g log det j (observed information)", x$a),
  E = sprintf("penalized likelihood l + %g log det i (expected information)", x$a),
  A = "Firth's adjusted score (mean bias reduction)", ML = "maximum likelihood")

print.obsinfo <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Call:\n"); print(x$call)
  cat(sprintf("\nMethod: %s\nLink: %s\n", obsinfo_method_label(x), x$link))
  if (!x$resolved) { cat("\nNo accepted solution was found.\n"); return(invisible(x)) }
  cat("\nCoefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L, quote = FALSE)
  cat(sprintf("\n%d observations; separation: %s\n", x$n, x$separation))
  invisible(x)
}

summary.obsinfo <- function(object, ...) {
  est <- object$coefficients; se <- object$se; z <- est / se
  structure(list(call = object$call, label = obsinfo_method_label(object), link = object$link,
                 resolved = object$resolved,
                 coefficients = cbind(Estimate = est, `Std. Error` = se, `z value` = z, `Pr(>|z|)` = 2 * pnorm(-abs(z))),
                 value = object$value, loglik = object$loglik, method = object$method, n = object$n,
                 separation = object$separation, diagnostics = object$diagnostics),
            class = "summary.obsinfo")
}

print.summary.obsinfo <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Call:\n"); print(x$call)
  cat(sprintf("\nMethod: %s\nLink: %s\n", x$label, x$link))
  if (!x$resolved) { cat("\nNo accepted solution was found.\n"); return(invisible(x)) }
  cat("\nCoefficients (Wald tests; standard errors from the inverse expected information):\n")
  printCoefmat(x$coefficients, digits = digits)
  dg <- x$diagnostics
  cat(sprintf("\n%d observations; separation: %s\n", x$n, x$separation))
  if (x$method != "A") cat(sprintf("Penalized log-likelihood %.6g; log-likelihood %.6g\n", x$value, x$loglik)) else
    cat(sprintf("Log-likelihood at the root %.6g\n", x$loglik))
  sol <- if (x$method == "A") "roots" else "local maxima"
  cat(sprintf("Distinct accepted %s: %d (%s stage)", sol, dg$n_distinct, dg$stage))
  if (dg$n_distinct > 1 && x$method != "A")
    cat(sprintf("; criterion values %s; the highest is selected", paste(signif(sort(dg$values_distinct, TRUE), 8), collapse = ", ")))
  cat("\n")
  if (x$method == "A") cat(sprintf("Max standardized residual %.1e\n", dg$std_resid_max)) else
    cat(sprintf("Max standardized gradient %.1e; largest Hessian eigenvalue %.3g\n", dg$std_grad_max, dg$lambda_max))
  cat(sprintf("Condition numbers of i and j: %.3g, %.3g", dg$cond_i, dg$cond_j))
  if (!is.null(dg$reevaluation)) cat(sprintf("; stable re-evaluation: %s", dg$reevaluation))
  cat("\n")
  invisible(x)
}

vcov.obsinfo <- function(object, ...) object$vcov

predict.obsinfo <- function(object, newdata = NULL, type = c("link", "response"), ...) {
  type <- match.arg(type)
  X <- if (is.null(newdata)) object$model$X * rep(object$scale, each = nrow(object$model$X)) else {
    tt <- delete.response(object$terms)
    model.matrix(tt, model.frame(tt, newdata, na.action = na.pass, xlev = object$xlevels))
  }
  eta <- drop(X %*% object$coefficients)
  if (type == "link") eta else make_link(object$link)$prob(eta)
}

## Confidence sets for the coefficients `parm` (names or indices): type
## "profile" (penalized-profile sets for "O" and "E",
## likelihood-ratio sets for "ML") or "wald". Returns the interval for each
## coefficient; for a disconnected set its hull, with every component in
## attr(, "components") and a warning. starts = "warm" uses only the
## neighbouring solutions and the unrestricted maxima at each profile point:
## faster than the protocol's full start set, with less protection against
## missed constrained maxima. control: settings of profile_set(), see
## profile_control(). Profile sets take one to a few minutes per coefficient
## for moderate n.
confint.obsinfo <- function(object, parm, level = 0.95, type = c("profile", "wald"),
                            starts = c("protocol", "warm"), control = profile_control(), ...) {
  type <- match.arg(type); starts <- match.arg(starts)
  cf <- object$coefficients; pn <- names(cf)
  if (missing(parm)) parm <- pn else if (is.numeric(parm)) parm <- pn[parm]
  idx <- match(parm, pn)
  if (anyNA(idx)) stop("unknown coefficient names in parm")
  lab <- paste(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2), trim = TRUE, scientific = FALSE, digits = 3), "%")
  ci <- matrix(NA_real_, length(idx), 2L, dimnames = list(parm, lab))
  if (!object$resolved) { warning("no accepted solution: intervals not available"); return(ci) }
  if (type == "wald") {
    z <- qnorm(1 - (1 - level) / 2)
    ci[] <- cbind(cf[idx] - z * object$se[idx], cf[idx] + z * object$se[idx])
    return(ci)
  }
  if (object$method == "A") stop("profile sets need an objective function; use type = \"wald\" for the adjusted score")
  comps <- setNames(vector("list", length(idx)), parm); status <- setNames(character(length(idx)), parm)
  for (j in seq_along(idx)) {
    r <- profile_set(object$model, object$fit, idx[j], object$starts, object$rescue, level, control, starts)
    status[j] <- r$status
    if (identical(r$status, "resolved")) {
      comps[[j]] <- r$components / object$scale[idx[j]]
      ci[j, ] <- range(comps[[j]])
    }
  }
  attr(ci, "status") <- status; attr(ci, "components") <- comps
  if (any(status != "resolved"))
    warning("profile set not resolved for ", paste(sprintf("%s (%s)", parm[status != "resolved"], status[status != "resolved"]), collapse = ", "))
  multi <- vapply(comps, function(cm) !is.null(cm) && nrow(cm) > 1L, TRUE)
  if (any(multi)) warning("disconnected confidence set for ", paste(parm[multi], collapse = ", "),
                          ": the interval shown is its hull; see attr(, \"components\")")
  ci
}
