#!/usr/bin/env Rscript
# Probit likelihood profiles for the neovascularization coefficient.
# Rscript Endometrial_profile.R [data_source [figure.pdf]]
# The default input is Endometrial.dat. A URL or the brglm2 data name is also accepted.
# The figure uses the Cairo PDF device.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop("Usage: Rscript Endometrial_profile.R [data_source [figure.pdf]]")
data_source <- if (length(args)) args[1L] else "Endometrial.dat"
figure_file <- if (length(args) == 2L) args[2L] else "Endometrial_profile.pdf"
if (!capabilities("cairo")) stop("This figure requires an R installation with Cairo support.")
dat <- if (identical(data_source, "brglm2")) {
    if (!requireNamespace("brglm2", quietly = TRUE)) stop("Package brglm2 is not installed.")
    env <- new.env()
    utils::data("endometrial", package = "brglm2", envir = env)
    env$endometrial
} else {
    read.table(data_source, header = TRUE)    # local path or complete URL
}
stopifnot(identical(names(dat), c("NV", "PI", "EH", "HG")), nrow(dat) == 79L,
          all(is.finite(as.matrix(dat))), all(dat$NV %in% 0:1), all(dat$HG %in% 0:1),
          identical(as.integer(table(dat$NV, dat$HG)), c(49L, 0L, 17L, 13L)))

# Scaling PI and EH changes each log determinant by a constant. The NV
# coefficient, fitted probabilities and relative profiles retain their original scale.
continuous <- scale(dat[c("PI", "EH")])
X <- cbind(intercept = 1, NV = dat$NV, continuous)
sign_y <- 2 * dat$HG - 1
nv <- match("NV", colnames(X))
nuisance <- setdiff(seq_len(ncol(X)), nv)
no_nv <- dat$NV == 0
stopifnot(qr(X)$rank == 4L, qr(X[no_nv, nuisance])$rank == 3L)

# For NV=0 a success lies strictly inside a triangle of failures in (PI,EH).
# Its positive barycentric coordinates certify overlap in the three-parameter model.
observation <- function(pi, eh, hg) {
    k <- which(no_nv & dat$PI == pi & dat$EH == eh & dat$HG == hg)
    stopifnot(length(k) > 0L)
    k[1L]
}
triangle <- c(observation(13, 1.64, 0), observation(8, 3.14, 0), observation(5, 2.31, 0))
interior <- observation(8, 2.64, 1)
barycentric <- solve(t(X[triangle, nuisance]), X[interior, nuisance])
stopifnot(min(barycentric) > 0, abs(sum(barycentric) - 1) < 1e-10)

# Evaluate minus the criterion, its gradient and its Hessian. For penalized
# fits the full design is retained here; restrictions are applied only below.
criterion <- function(beta, method, design = X, signs = sign_y) {
    eta <- drop(design %*% beta)
    u <- signs * eta
    mills <- function(z) exp(dnorm(z, log = TRUE) - pnorm(z, log.p = TRUE))
    r <- mills(u)
    psi <- r * (u + r)
    observed <- crossprod(design, psi * design)
    value <- -sum(pnorm(u, log.p = TRUE))
    grad <- -drop(crossprod(design, signs * r))
    if (method == "L") return(list(value = value, grad = grad, hess = observed))
    if (method == "E") {
        r_plus <- mills(eta)
        r_minus <- mills(-eta)
        w <- exp(2 * dnorm(eta, log = TRUE) - pnorm(eta, log.p = TRUE) -
                 pnorm(-eta, log.p = TRUE))
        dlogw <- -2 * eta - r_plus + r_minus
        dw <- w * dlogw
        ddw <- w * (dlogw^2 - 2 + r_plus * (eta + r_plus) + r_minus * (-eta + r_minus))
    } else {
        stopifnot(method == "O")
        w <- psi
        dpsi <- r - (u + 2 * r) * psi
        dw <- signs * dpsi
        ddw <- 2 * psi * (psi - 1) - (u + 2 * r) * dpsi
    }
    factor <- chol(crossprod(design, w * design))
    inverse <- chol2inv(factor)
    leverage <- rowSums((design %*% inverse) * design)
    value <- value - sum(log(diag(factor)))        # minus one half log det
    grad <- grad - drop(crossprod(design, dw * leverage)) / 2
    # d_rs log det K = tr(K^-1 K_rs) - tr(K^-1 K_r K^-1 K_s).
    products <- lapply(seq_len(ncol(design)), function(k)
        inverse %*% crossprod(design, (dw * design[, k]) * design))
    traces <- vapply(products, function(a)
        vapply(products, function(b) sum(a * t(b)), 0.0), numeric(ncol(design)))
    hess <- observed - (crossprod(design, (ddw * leverage) * design) - traces) / 2
    list(value = value, grad = grad, hess = hess)
}

# Constrained and unrestricted fits use analytic derivatives and the same acceptance criteria.
fit_from <- function(terms, initial) {
    last <- NULL
    cached <- function(b) {
        if (is.null(last) || !identical(as.numeric(b), last$b))
            last <<- list(b = as.numeric(b), z = terms(b))
        last$z
    }
    fit <- nlminb(as.numeric(initial), function(b) cached(b)$value,
                  gradient = function(b) cached(b)$grad,
                  hessian = function(b) cached(b)$hess,
                  control = list(iter.max = 400L, eval.max = 800L,
                                 rel.tol = 1e-14, x.tol = 1e-12))
    b <- fit$par
    z <- terms(b)
    for (k in seq_len(5L)) {
        if (max(abs(z$grad)) < 1e-10) break
        candidate <- b - drop(solve(z$hess, z$grad))
        trial <- terms(candidate)
        tol <- 32 * .Machine$double.eps * (1 + abs(z$value))
        if (trial$value > z$value + tol || max(abs(trial$grad)) >= max(abs(z$grad))) break
        b <- candidate
        z <- trial
    }
    score_norm <- max(abs(z$grad))
    curvature <- min(eigen(z$hess, symmetric = TRUE, only.values = TRUE)$values)
    if (!is.finite(z$value) || !(score_norm < 1e-6 && curvature > 0))
        stop(sprintf("Unverified optimum: nlminb code %d, score %.3g, curvature %.3g",
                     fit$convergence, score_norm, curvature))
    list(par = b, value = -z$value, score = score_norm, curvature = curvature)
}
best_of <- function(terms, starts) {
    fits <- lapply(seq_len(nrow(starts)), function(k) fit_from(terms, starts[k, ]))
    values <- vapply(fits, function(fit) fit$value, 0.0)
    best <- fits[[which.max(values)]]
    best$spread <- diff(range(values))
    best
}
start4 <- rbind(numeric(4), diag(4), -diag(4), as.matrix(expand.grid(rep(list(c(-1, 1)), 4))))
full <- lapply(c(E = "E", O = "O"), function(m)
    best_of(function(b) criterion(b, m), start4))
stopifnot(max(vapply(full, function(f) f$spread, 0.0)) < 1e-7)
# For every finite nuisance vector the 13 NV=1 contributions increase strictly
# to zero. Maximizing on the overlapping NV=0 subset gives the exact reference.
limit <- fit_from(function(b) criterion(b, "L", X[no_nv, nuisance], sign_y[no_nv]), numeric(3))
tops <- c(L = limit$value, E = full$E$value, O = full$O$value)

profile_fit <- function(t, method, starts) {
    restricted <- function(theta) {
        b <- numeric(4); b[nv] <- t; b[nuisance] <- theta
        z <- criterion(b, method)
        list(value = z$value, grad = z$grad[nuisance], hess = z$hess[nuisance, nuisance])
    }
    best_of(restricted, unique(as.matrix(starts)))
}
grid <- seq(-1, 7, by = 0.025)
cutoff <- qchisq(0.95, df = 1) / 2
trace_profile <- function(method) {
    anchor <- if (method == "L") limit$par else full[[method]]$par[nuisance]
    centre <- if (method == "L") 0 else full[[method]]$par[nv]
    pivot <- which.min(abs(grid - centre))
    out <- vector("list", length(grid))
    # Start at the preceding fit, zero and the anchor; discard duplicate starts.
    for (indices in list(pivot:length(grid), rev(seq_len(pivot - 1L)))) {
        previous <- anchor
        for (k in indices) {
            out[[k]] <- profile_fit(grid[k], method, rbind(previous, numeric(3), anchor))
            previous <- out[[k]]$par
        }
    }
    out
}
methods <- c(L = "L", E = "E", O = "O")
profiles <- lapply(methods, trace_profile)
relative <- sapply(methods, function(m)
    vapply(profiles[[m]], function(f) f$value, 0.0) - tops[[m]])
stopifnot(all(diff(relative[, "L"]) >= -1e-10), max(relative) <= 1e-8,
          abs(tail(relative[, "L"], 1L)) < 1e-7)

# Independent, deterministic checks from all eight sign combinations at spaced
# grid points. Stationarity and positive curvature are local checks; multistart
# agreement does not prove global uniqueness of a penalized maximum.
corners <- 1.5 * as.matrix(expand.grid(rep(list(c(-1, 1)), 3)))
gain <- 0
for (m in methods) for (k in seq(1L, length(grid), by = 32L)) {
    audited <- profile_fit(grid[k], m, corners)
    gain <- max(gain, audited$value - profiles[[m]][[k]]$value)
}
stopifnot(gain < 1e-8)
at <- function(t, m) {
    j <- which.min(abs(grid - t))
    profile_fit(t, m, rbind(profiles[[m]][[j]]$par, numeric(3)))$value - tops[[m]]
}
# Locate every threshold crossing on the plotting grid and solve the re-fitted
# profile equation in each bracket. The ML upper limit is infinite by separation.
interval <- function(m) {
    above <- relative[, m] >= -cutoff
    changes <- which(above[-1L] != head(above, -1L))
    stopifnot(length(changes) == if (m == "L") 1L else 2L)
    roots <- vapply(changes, function(k) {
        f <- function(t) at(t, m) + cutoff
        root <- uniroot(f, grid[k + 0:1], tol = 1e-9)$root
        stopifnot(abs(f(root)) < 1e-7)
        root
    }, 0.0)
    if (m == "L") c(roots, Inf) else roots
}
ci <- t(vapply(methods, interval, numeric(2)))
colnames(ci) <- c("lower", "upper")
bhat <- vapply(full, function(f) f$par[nv], 0.0)
peak_error <- vapply(names(full), function(m) abs(at(bhat[[m]], m)), 0.0)
x_patient <- X[which(dat$NV == 1 & dat$PI == 21 & dat$EH == 0.98)[1L], ]
phat <- vapply(full, function(f) pnorm(sum(x_patient * f$par)), 0.0)
stopifnot(max(peak_error) < 1e-8)
allfits <- unlist(profiles, recursive = FALSE)
print(with(dat, table(NV, HG)))
print(cbind(beta_NV = bhat, probability = phat, ci[c("E", "O"), ]), digits = 8)
cat(sprintf("Unpenalized supremum %.9f; nominal 95%% LR set [%.6f, Inf); LR at zero %.6f.\n",
            tops["L"], ci["L", "lower"], -2 * at(0, "L")))
cat(sprintf("Grid: %d points. Maximum nuisance score %.3g; minimum eigenvalue of negative nuisance Hessian %.4g.\n",
            length(grid), max(vapply(allfits, function(f) f$score, 0.0)),
            min(vapply(allfits, function(f) f$curvature, 0.0))))
cat(sprintf("Largest objective spread across 25 full starts %.3g; peak error %.3g; restart gain %.3g.\n",
            max(vapply(full, function(f) f$spread, 0.0)), max(peak_error), gain))

# Vector curves; line patterns distinguish the almost coincident penalized fits in print.
colour <- c(L = "#222222", E = "#00689D", O = "#B94D12")
pattern <- c(L = "13", E = "solid", O = "52")
draw_profiles <- function() {
    old <- par(mar = c(3.75, 4.15, 0.45, 0.55), mgp = c(2.45, 0.55, 0), tcl = -0.25,
               las = 1, xaxs = "i", yaxs = "i", lend = "round")
    on.exit(par(old))
    plot(NA, type = "n", xlim = range(grid), ylim = c(-6.05, 0.35), axes = FALSE,
         xlab = expression("Neovascularization coefficient, " * beta[v]),
         ylab = "Relative profile log-likelihood")
    axis(1, at = -1:7, cex.axis = 0.92)
    axis(2, at = seq(-6, 0, by = 1), cex.axis = 0.92)
    box(bty = "l", col = "grey40")
    abline(h = 0, col = "grey85", lwd = 0.65)
    abline(h = -cutoff, col = "grey45", lty = "83", lwd = 0.8)
    for (m in methods) lines(grid, relative[, m], col = colour[m], lty = pattern[m],
                             lwd = if (m == "E") 2.4 else 1.8)
    legend("bottom", inset = 0.04, bty = "n", cex = 0.9, seg.len = 2.5, x.intersp = 0.65, y.intersp = 1.03,
           legend = c("Unpenalized", "Expected information", "Observed information"),
           col = colour, lty = pattern, lwd = c(1.8, 2.4, 1.8))
}

# Save the figure.
cairo_pdf(figure_file, width = 5.9, height = 3.8, family = "sans", pointsize = 10)
draw_profiles()
invisible(dev.off())
cat("Figure:", figure_file, "\n")

# Also draw it on screen when the script is run in an interactive R session.
if (interactive()) draw_profiles()

