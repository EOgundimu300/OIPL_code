#!/usr/bin/env Rscript
# Reproduce the endometrial cancer probit illustration. Base R only.
# From the extracted bundle, run: Rscript Endometrial_example.R Endometrial.dat
# Data: Alan Agresti, https://www.alanagresti.com/glm/data/Endometrial.dat
args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args)) args[1L] else "Endometrial.dat"
dat <- read.table(data_file, header = TRUE)
stopifnot(identical(names(dat), c("NV", "PI", "EH", "HG")), nrow(dat) == 79L,
          identical(as.integer(table(dat$NV, dat$HG)), c(49L, 0L, 17L, 13L)))

# Centre and scale only the two continuous predictors. This is an affine
# change of coordinates. The NV coefficient and fitted probabilities agree
# with those in the original units of PI and EH.
ctr <- vapply(dat[c("PI", "EH")], mean, numeric(1))
scl <- vapply(dat[c("PI", "EH")], sd, numeric(1))
X <- cbind(intercept = 1, NV = dat$NV, PI = (dat$PI - ctr[1L]) / scl[1L],
           EH = (dat$EH - ctr[2L]) / scl[2L])
sign_y <- 2 * dat$HG - 1
p <- ncol(X)

# Observed curvature of log(Phi(u)) and its first two derivatives.
curvature <- function(u) {
    r <- exp(dnorm(u, log = TRUE) - pnorm(u, log.p = TRUE))
    w <- r * (u + r)
    dw <- r - (u + 2*r)*w
    list(r = r, w = w, dw = dw, d2w = 2*w*(w-1) - (u+2*r)*dw)
}

terms_at <- function(beta) {
    eta <- drop(X %*% beta)
    pos <- curvature(eta)
    neg <- curvature(-eta)
    obs <- curvature(sign_y * eta)
    logw <- 2*dnorm(eta, log = TRUE) - pnorm(eta, log.p = TRUE) -
            pnorm(eta, lower.tail = FALSE, log.p = TRUE)
    w <- exp(logw)
    dlogw <- -2*eta - pos$r + neg$r
    list(eta = eta, loglik = sum(pnorm(sign_y*eta, log.p = TRUE)),
         score = drop(crossprod(X, sign_y*obs$r)),
         fisher = crossprod(X, w*X), observed = crossprod(X, obs$w*X),
         E = list(w = w, d = w*dlogw,
                  dd = w*(dlogw^2 - 2 + pos$w + neg$w)),
         O = list(w = obs$w, d = sign_y*obs$dw, dd = obs$d2w))
}

# Analytic gradient and Hessian of MINUS the penalised log likelihood.
penalty_fit_terms <- function(beta, method) {
    z <- terms_at(beta)
    w <- z[[method]]
    K <- if (method == "E") z$fisher else z$observed
    upper <- chol(K)
    V <- chol2inv(upper)
    diag_leverage <- rowSums((X %*% V)*X)
    grad <- z$score + drop(crossprod(X, w$d*diag_leverage))/2
    # Hessian of log det K: tr(K^-1 K_rs) - tr(K^-1 K_r K^-1 K_s).
    deriv <- lapply(seq_len(p), function(r) V %*% crossprod(X, (X[,r]*w$d)*X))
    trace_product <- sapply(deriv, function(a)
        vapply(deriv, function(b) sum(a*t(b)), numeric(1)))
    hess <- -z$observed + (crossprod(X, (w$dd*diag_leverage)*X) - trace_product)/2
    list(value = -z$loglik - sum(log(diag(upper))), grad = -grad, hess = -hess)
}

adjusted <- function(beta) {
    z <- terms_at(beta)
    V <- chol2inv(chol(z$fisher))
    hat <- z$E$w*rowSums((X %*% V)*X)
    list(score = z$score - drop(crossprod(X, hat*z$eta))/2, fisher = z$fisher)
}
whitened_norm <- function(score, fisher) sqrt(sum(forwardsolve(t(chol(fisher)), score)^2))

qf_root <- function(start) {
    b <- start
    for (iter in seq_len(500L)) {
        z <- adjusted(b)
        old <- whitened_norm(z$score, z$fisher)
        if (old < 1e-10) return(b)
        direction <- drop(solve(z$fisher, z$score))
        accepted <- FALSE
        for (k in 0:30) {
            candidate <- b + direction/2^k
            trial <- tryCatch(adjusted(candidate), error = function(e) NULL)
            if (!is.null(trial) && is.finite(whitened_norm(trial$score, trial$fisher)) &&
                whitened_norm(trial$score, trial$fisher) < old) {
                b <- candidate
                accepted <- TRUE
                break
            }
        }
        if (!accepted) stop("Adjusted-score step failed")
    }
    stop("Adjusted-score iteration limit")
}

objective_solution <- function(start, method) {
    fn <- function(b) penalty_fit_terms(b, method)
    fit <- nlminb(start, function(b) fn(b)$value,
                  gradient = function(b) fn(b)$grad,
                  hessian = function(b) fn(b)$hess,
                  control = list(iter.max = 500, eval.max = 1000, rel.tol = 1e-14,
                                 x.tol = 1e-12))
    b <- fit$par
    # Refine a stationary candidate while decreasing its gradient norm.
    for (k in 1:4) {
        z <- fn(b)
        if (max(abs(z$grad)) < 1e-11) break
        proposal <- b - drop(solve(z$hess, z$grad))
        if (max(abs(fn(proposal)$grad)) >= max(abs(z$grad))) break
        b <- proposal
    }
    z <- fn(b)
    stopifnot(whitened_norm(z$grad, terms_at(b)$fisher) < 1e-8,
              min(eigen(z$hess, symmetric=TRUE, only.values=TRUE)$values) > 0)
    b
}

# All 25 starts are deterministic. No random simulation is used here.
starts <- rbind(rep(0, p), diag(p), -diag(p),
                as.matrix(expand.grid(rep(list(c(-1,1)), p))))
from_starts <- function(fun, starts) {
    candidates <- lapply(seq_len(nrow(starts)), function(k)
        tryCatch(fun(as.numeric(starts[k,])), error = function(e) NULL))
    ok <- !vapply(candidates, is.null, logical(1))
    cat("Accepted starts:", sum(ok), "/", length(ok), "\n")
    if (!any(ok)) stop("No accepted candidates")
    do.call(rbind, candidates[ok])
}
fits <- list()
for (method in c("E", "O")) {
    candidates <- from_starts(function(b) objective_solution(b, method), starts)
    values <- apply(candidates, 1L, function(b) penalty_fit_terms(b, method)$value)
    fits[[method]] <- candidates[which.min(values),]
    cat(method, "maximum coefficient discrepancy from selected solution:",
        max(abs(sweep(candidates, 2L, fits[[method]], "-"))), "\n")
}
candidates <- from_starts(qf_root, starts)
distances <- rowSums(sweep(candidates, 2L, fits$E, "-")^2)
fits$AS <- candidates[which.min(distances),]
cat("AS maximum coefficient discrepancy from selected root:",
    max(abs(sweep(candidates, 2L, fits$AS, "-"))), "\n")

jacobian <- function(f,b,h=1e-5) {
    sapply(seq_along(b), function(k) {
        e <- numeric(length(b)); e[k] <- h
        (f(b+e)-f(b-e))/(2*h)
    })
}
stopifnot(min(svd(jacobian(function(b) adjusted(b)$score, fits$AS))$d)>1e-8)

# The two within-NV=1 medians form an actually observed patient profile.
profile <- vapply(dat[dat$NV==1,c("PI","EH")], median, numeric(1))
stopifnot(any(dat$NV==1 & dat$PI==profile[1L] & dat$EH==profile[2L]))
x1 <- unname(c(1,1,(profile-ctr)/scl))
x0 <- x1; x0[2L] <- 0
out <- t(vapply(fits, function(b) {
    V <- chol2inv(chol(terms_at(b)$fisher)); eta <- sum(x1*b)
    se <- sqrt(drop(crossprod(x1, V %*% x1))); se_nv <- sqrt(V[2L,2L]); z <- qnorm(.975)
    c(beta_NV = unname(b[2L]), NV_lo = unname(b[2L])-z*se_nv,
      NV_hi = unname(b[2L])+z*se_nv, p0 = pnorm(sum(x0*b)), p1 = pnorm(eta),
      p1_lo = pnorm(eta-z*se), p1_hi = pnorm(eta+z*se))
}, numeric(7)))
print(profile)
print(out, digits = 9)

# Alternative Wald interval for the same observed-information penalized estimate.
nv_estimate <- unname(fits$O[2L])
nv_variance <- chol2inv(chol(terms_at(fits$O)$observed))[2L, 2L]
nv_interval <- nv_estimate + qnorm(c(0.025, 0.975)) * sqrt(nv_variance)
cat("NV Wald interval using observed information at the O estimate:", nv_interval, "\n")
