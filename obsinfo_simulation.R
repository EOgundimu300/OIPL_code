#!/usr/bin/env Rscript
## Simulation and fixed-design calculations for the binomial-regression study.
## Run from the directory containing both obsinfo R files. For available stages:
##   Rscript obsinfo_simulation.R --help
## In an R session, source("obsinfo_simulation.R") to load the definitions.
## Outputs go to cfg$out_dir. Completed replicate blocks can be resumed with
## the same two scripts. Scenario order and stream allocation are fixed.

if (!exists(".obsinfo_worker")) source("obsinfo_core.R", local = TRUE)
if (!exists("profile_set") || !exists("tag_stage"))
  stop("this obsinfo_simulation.R needs the obsinfo_core.R delivered with it")


## =============================================================================
## Part 1. Settings
## =============================================================================
## Workers: the CPUs allocated by SLURM (SLURM_CPUS_PER_TASK) when set, otherwise
## all but one of the machine's cores. parallel::detectCores() counts every core
## of a cluster node, not the CPUs allocated to the job.
default_cores <- function() {
  s <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
  if (!is.na(s) && s >= 1L) return(s)
  available <- parallel::detectCores()
  if (is.na(available)) 1L else max(1L, available - 1L)
}
cfg <- list(
  out_dir = "sim_out",
  code_revision = 0L,                 # archive key for the standard production streams
  cores = default_cores(),
  seed = 20261101L,
  n_pilot = 2000L,                    # pilot samples per scenario
  n_profile_coverage = 10000L,        # production indices with profile coverage at the truth
  n_profile_sets = 2000L,             # production indices with set reconstruction
  audit_every = 100L,                 # production audit frequency
  chunk_size = 200L,                  # replicates per archive file
  recon_starts = "protocol",          # use all specified starts for profile fits
  keep_candidates = TRUE,             # archive every candidate of every start and stage
  store_starts = TRUE,                # archive the 50 x 3 auxiliary draws and the audit draws
  num = list(
    maxit = 500L,
    maxit_rescue = 2000L,
    grid_h = 0.25,                    # coarse profile grid spacing on the u scale
    W0 = 8,                           # initial grid [-W0, W0]
    extensions = c(16, 32, 64),       # tail extensions
    n_refine = 3L,                    # refinements to grid_h/2, grid_h/4, grid_h/8
    bisect_maxit = 200L,
    bisect_width = 1e-6,              # bracket width <= bisect_width (1 + |t|)
    bisect_resid = 1e-6,              # |L - cutoff| < bisect_resid
    stab_endpoint = 1e-4,             # endpoint movement < stab_endpoint (1 + |t|)
    stab_length = 1e-3,               # relative length change < stab_length
    optimize_tol = .Machine$double.eps^0.25,
    gradient_scale = "estimate",      # constrained fits: i at the unrestricted estimate
    reeval_cond = 1e12))              # condition number triggering re-evaluation
CUT <- 3.841458820694124              # chi-square(1) 0.95 quantile
Z975 <- 1.959963984540054             # Wilson intervals
MC_Z975 <- 1.96                       # Monte Carlo intervals for means

## Integer codes of the candidate matrices (also in archive_notes.txt); the
## reconstruction codes PROFILE_KINDS and PROFILE_FLAGS are in the core.
STAGES <- c("routine", "rerun", "rescue", "newton", "audit", "ray")
STATUSES <- c("converged", "maxit", "stalled", "start_infeasible", "newton_failed", "singular_jacobian")
METHODS <- c("qfs", "qfs+newton")


## =============================================================================
## Part 2. Scenarios and random-number streams
## =============================================================================
## The probability difference compares x2 = 1 with x2 = 0 at x1 = 0.
d_contrast <- function(beta, link) {
  G <- make_link(link)$prob
  G(beta[1] + beta[3]) - G(beta[1])
}
make_scenarios <- function() {
  g <- expand.grid(n = c(25L, 50L, 100L, 200L), setting = c("S1", "S2"),
                   link = c("probit", "cloglog"), stringsAsFactors = FALSE)
  lapply(seq_len(nrow(g)), function(k) {
    n <- g$n[k]; X <- protocol_design(n); b0 <- protocol_truth[[g$setting[k]]]
    ref <- first_order_bias(X, b0, g$link[k])   # b_E, b_O, event rate, eigenvalues, rank
    d0 <- d_contrast(b0, g$link[k])
    list(id = k, link = g$link[k], setting = g$setting[k], n = n, X = X, beta0 = b0,
         prob = ref$prob, ref = ref, d0 = d0)
  })
}

## Eight streams per scenario, allocated in scenario order by nextRNGStream from
## the master seed; replicate b uses the (b - 1)-th substream of each stream.
## The two reserved restart streams retain their positions in the allocation.
## Restart audits use the production audit stream, independently of their fits.
stream_names <- c("pilot_resp", "pilot_start", "pilot_audit", "prod_resp", "prod_start",
                  "prod_audit", "restart_resp", "restart_start")
phase_streams <- list(pilot = c("pilot_resp", "pilot_start", "pilot_audit"),
                      production = c("prod_resp", "prod_start", "prod_audit"),
                      restart = c("restart_resp", "restart_start", "prod_audit"))
make_streams <- function(seed) {
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(seed)
  s <- .Random.seed
  lapply(1:16, function(k) {
    st <- list()
    for (nm in stream_names) { s <<- parallel::nextRNGStream(s); st[[nm]] <- s }
    st
  })
}
substreams <- function(s, R) {
  out <- vector("list", R)
  out[[1]] <- s
  if (R > 1) for (b in 2:R) out[[b]] <- parallel::nextRNGSubStream(out[[b - 1]])
  out
}
use_stream <- function(s) assign(".Random.seed", s, envir = globalenv())


## =============================================================================
## Part 3. One replicate: records, re-evaluation, audits, profile inference
## =============================================================================
## ---- Archive records of fits ----------------------------------------------------
sel_beta <- function(f) if (is.null(f) || is.null(f$selected)) rep(NA_real_, 3) else f$selected$beta
nz <- function(x) if (is.null(x) || !length(x)) NA_real_ else x

## One row per candidate: stage, start index (within the start matrix of its
## stage), status, accepted, objective, coefficients, standardized gradient,
## curvature (lambda_max, ||H||_2), iterations, condition number of i.
cand_matrix_max <- function(cands) {
  if (!length(cands)) return(NULL)
  t(vapply(cands, function(cc) c(
    stage = code_of(cc$stage, STAGES), start = nz(cc$start_index), status = code_of(cc$status, STATUSES),
    accepted = as.numeric(isTRUE(cc$accepted)), value = nz(cc$value),
    beta = if (length(cc$beta) == 3) cc$beta else rep(NA_real_, 3),
    std_grad_max = nz(cc$std_grad_max), lambda_max = nz(cc$lambda_max), hess_norm = nz(cc$hess_norm),
    iterations = nz(cc$iterations), info_cond = nz(cc$info_cond)), numeric(13)))
}
## Roots: as above with the method (qfs, qfs+newton), the scoring status before
## Newton continuation, the standardized residual and the Jacobian singular values.
cand_matrix_root <- function(cands) {
  if (!length(cands)) return(NULL)
  t(vapply(cands, function(cc) c(
    stage = code_of(cc$stage, STAGES), start = nz(cc$start_index), status = code_of(cc$status, STATUSES),
    accepted = as.numeric(isTRUE(cc$accepted)), method = code_of(cc$method, METHODS),
    qfs_status = code_of(cc$qfs_status, STATUSES),
    beta = if (length(cc$beta) == 3) cc$beta else rep(NA_real_, 3), std_resid_max = nz(cc$std_resid_max),
    jac_sv = if (length(cc$jac_singular_values) == 3) cc$jac_singular_values else rep(NA_real_, 3),
    iterations = nz(cc$iterations), info_cond = nz(cc$info_cond)), numeric(15)))
}
obs_info_spectrum <- function(mod, beta) {
  ev <- eigen(bp_info(mod, beta)$observed, symmetric = TRUE, only.values = TRUE)$values
  list(values = ev, cond = if (min(ev) > 0) max(ev) / min(ev) else Inf)
}
record_fit <- function(f, mod, cfg) {
  if (is.null(f)) return(NULL)
  s <- f$selected
  out <- list(stage = f$stage, n_accepted = f$n_accepted, n_distinct = f$n_distinct,
              distinct = t(vapply(f$distinct, function(d) c(d$beta, d$value, d$n_starts), numeric(5))),
              status = table(vapply(f$history, function(cc) cc$status, "")))
  if (!is.null(s)) {
    out[c("beta", "value", "grad", "std_grad_max", "lambda_max", "hess_norm", "hess_eigenvalues",
          "info_eigenvalues", "info_cond", "iterations")] <- s[c("beta", "value", "grad", "std_grad_max",
          "lambda_max", "hess_norm", "hess_eigenvalues", "info_eigenvalues", "info_cond", "iterations")]
    out$selected_stage <- s$stage; out$selected_start <- s$start_index
    oj <- obs_info_spectrum(mod, s$beta)
    out$obs_info_eigenvalues <- oj$values; out$obs_info_cond <- oj$cond
  }
  if (cfg$keep_candidates) out$candidates <- cand_matrix_max(f$history)
  out
}
record_root <- function(f, mod, cfg) {
  if (is.null(f)) return(NULL)
  s <- f$selected
  out <- list(stage = f$stage, n_accepted = f$n_accepted, n_distinct = f$n_distinct,
              distinct = t(vapply(f$distinct, function(d) c(d$beta, d$std_resid_max, d$n_starts), numeric(5))),
              status = table(vapply(f$history, function(cc) cc$status, "")))
  if (!is.null(s)) {
    out[c("beta", "U", "std_resid_max", "jac_singular_values", "info_eigenvalues", "info_cond",
          "iterations", "method")] <- s[c("beta", "U", "std_resid_max", "jac_singular_values",
          "info_eigenvalues", "info_cond", "iterations", "method")]
    out$selected_stage <- s$stage; out$selected_start <- s$start_index
    oj <- obs_info_spectrum(mod, s$beta)
    out$obs_info_eigenvalues <- oj$values; out$obs_info_cond <- oj$cond
  }
  if (cfg$keep_candidates) out$candidates <- cand_matrix_root(f$history)
  out
}

## ---- Re-evaluation of badly conditioned fits ---------------------
## When the condition number of i or j at a selected estimate exceeds
## cfg$num$reeval_cond, stable_reevaluation() recomputes the
## acceptance quantities from singular value decompositions. The re-evaluation
## never replaces or deletes an estimate: "disagree" (acceptance fails, or the
## objective moves by more than 1e-8 (1 + |Q|)) and "failed" are reported as
## unresolved numerical checks.
reevaluate <- function(mod, fr, kind, cfg) {
  if (is.null(fr) || is.null(fr$beta)) return(NULL)
  ci <- fr$info_cond; cj <- fr$obs_info_cond; thr <- cfg$num$reeval_cond
  if (!(ci > thr || cj > thr)) return(NULL)                  # Inf triggers; NaN cannot occur
  r <- stable_reevaluation(mod, fr$beta, kind)
  if (!r$ok) return(list(status = "failed", reason = r$reason, cond_i = ci, cond_j = cj))
  if (kind == "A") {
    diffs <- c(std_resid_max = r$std_resid_max - fr$std_resid_max,
               sigma_min = min(r$jac_singular_values) - min(fr$jac_singular_values))
    agree <- r$accept
  } else {
    diffs <- c(value = r$value - fr$value, std_grad_max = r$std_grad_max - fr$std_grad_max,
               lambda_max = r$lambda_max - fr$lambda_max, hess_norm = r$hess_norm - fr$hess_norm)
    agree <- r$accept && abs(diffs[["value"]]) <= 1e-8 * (1 + abs(fr$value))
  }
  list(status = if (agree) "agree" else "disagree", cond_i = ci, cond_j = cj,
       svd = r[setdiff(names(r), c("ok", "accept"))], differences = diffs)
}

## ---- Audits -------------------------------------------------------
## 40 additional N(0, 16 I) starts with the routine 500-iteration ceiling and no
## rescue. Records the audit candidates, improved objectives, extra maxima or
## roots, and the selection over the routine and audit candidates together
## (select_maximum, select_root), including ties within 1e-10 (1 + |Q|) that
## the norm and lexicographic rules decide; whether that selection differs
## from the routine selection, and by how much. Diagnostic only: audit fits
## never replace routine estimates. ML is not audited: for probit and cloglog
## the log-likelihood is strictly concave for a full-rank design, so an
## overlap maximum is unique.
selection_change <- function(sel, union) {
  a <- sel_beta(list(selected = sel)); u <- sel_beta(list(selected = union))
  if (anyNA(a) || anyNA(u))
    return(list(selection_changed = anyNA(a) != anyNA(u), change = rep(NA_real_, 3), change_rel = NA_real_))
  list(selection_changed = !same_solution(a, u), change = u - a,
       change_rel = sqrt(sum((u - a)^2)) / (1 + max(sqrt(sum(a^2)), sqrt(sum(u^2)))))
}
run_audit <- function(mod, fits, A_st, ref, cfg) {
  out <- list()
  for (k in c("E", "O")) {
    f <- fits[[k]]
    au <- fit_criterion(mod, k, 0.5, A_st, maxit = cfg$num$maxit)
    au$candidates <- tag_stage(au$candidates, "audit")
    pool <- c(f$candidates, au$candidates)
    union <- select_maximum(pool)
    acc <- Filter(function(cc) isTRUE(cc$accepted), pool)
    Q <- if (is.null(f$selected)) NA_real_ else f$selected$value
    Qa <- if (is.null(au$selected)) NA_real_ else au$selected$value
    out[[k]] <- c(list(n_distinct_union = length(distinct_solutions(lapply(acc, `[[`, "beta"))),
                       Q = Q, Q_audit = Qa, Q_union = if (is.null(union)) NA_real_ else union$value,
                       beta_audit = sel_beta(au), beta_union = sel_beta(list(selected = union))),
                  selection_change(f$selected, union),
                  list(material = isTRUE(Qa - Q > 1e-8 * (1 + abs(Q))) || (is.na(Q) && !is.na(Qa)),
                       candidates = if (cfg$keep_candidates) cand_matrix_max(au$candidates)))
  }
  au <- fit_adjscore(mod, A_st, ref = ref, maxit = cfg$num$maxit)
  au$candidates <- tag_stage(au$candidates, "audit")
  pool <- c(fits$A$candidates, au$candidates)
  union <- select_root(pool, if (is.null(ref)) numeric(mod$p) else ref)
  acc <- Filter(function(cc) isTRUE(cc$accepted), pool)
  out$A <- c(list(n_distinct_union = length(distinct_solutions(lapply(acc, `[[`, "beta"))),
                  beta_audit = sel_beta(au), beta_union = sel_beta(list(selected = union))),
             selection_change(fits$A$selected, union),
             list(candidates = if (cfg$keep_candidates) cand_matrix_root(au$candidates)))
  out
}

## ---- Profile inference -----------------------------------------------
## Penalized-profile statistic for beta2 at its true value; the constrained
## fit uses the gradient scale of the reconstruction (cfg$num$gradient_scale).
pl_at_truth <- function(mod, f, t, starts, rescue, cfg) {
  if (is.null(f$selected)) return(list(stat = NA_real_, flag = "unrestricted_missing", covered = NA))
  isc <- if (identical(cfg$num$gradient_scale, "estimate")) bp_info(mod, f$selected$beta)$expected[-3, -3, drop = FALSE]
  r <- profile_stat(mod, f, 3, t, starts, rescue, maxit = cfg$num$maxit, maxit_rescue = cfg$num$maxit_rescue,
                    info_scale = isc)
  ok <- identical(r$flag, "ok") && is.finite(r$stat)
  sel <- r$fit$selected
  list(stat = r$stat, flag = r$flag, covered = if (ok) r$stat <= CUT else NA, Q = f$selected$value,
       Q_constrained = if (is.null(sel)) NA_real_ else sel$value,
       beta_constrained = if (is.null(sel)) NULL else sel$beta,
       constrained_stage = if (is.null(r$fit)) NA_character_ else r$fit$stage)
}

## Reconstruction of C = {t : L(t) <= 3.8415} for beta2 by profile_set() with the settings cfg$num and the start policy cfg$recon_starts.
reconstruct_set <- function(mod, f, starts, rescue, cfg)
  profile_set(mod, f, 3, starts, rescue, level = 0.95, control = cfg$num, starts_policy = cfg$recon_starts)

## ---- One replicate --------------------------------------------------------------
## The response vector, the RNG states and the auxiliary starts are stored
## before any fitting. Each component (separation, control variate, each
## primary fit and its record and re-evaluation, each Wald interval, the audit,
## each profile coverage check and each reconstruction) has its own error
## capture; an R error leaves that component missing and is logged in
## rec$errors, never discarding the replicate.
record_sep <- function(sep) {
  one <- function(res) list(d = res$lp1$value, q = res$lp2$value, status = c(res$lp1$status, res$lp2$status),
                            iterations = c(res$lp1$iterations, res$lp2$iterations),
                            direction_complete = res$lp1$direction, direction_quasi = res$lp2$direction,
                            certificate_lp1 = res$lp1$certificate, certificate_lp2 = res$lp2$certificate,
                            dual_lp1 = res$lp1$dual, dual_lp2 = res$lp2$dual)
  dec <- if (is.null(sep$rechecked)) sep$scaled else sep$rechecked
  list(category = sep$category, d = dec$lp1$value, q = dec$lp2$value,
       decided_by = if (is.null(sep$rechecked)) "column-scaled programs" else "unscaled recheck",
       column_scale = sep$column_scale, certificates_ok_scaled = sep$certificates_ok,
       scaled = one(sep$scaled), rechecked = if (!is.null(sep$rechecked)) one(sep$rechecked))
}

analyse_replicate <- function(sc, b, seeds, phase, audit, prof_cov, prof_set, cfg) {
  t0 <- proc.time()[["elapsed"]]
  tm <- c(fits = NA, audit = 0, coverage = 0, sets = 0)
  errs <- character(0)
  part <- function(what, expr) tryCatch(expr, error = function(e) { errs[[what]] <<- conditionMessage(e); NULL })
  use_stream(seeds$resp)
  y <- as.numeric(runif(sc$n) < sc$prob)
  use_stream(seeds$start)
  Zs <- matrix(rnorm(150), 50, 3)                  # filled by column: start k uses draws k, 50 + k, 100 + k
  st <- protocol_starts(3, Zs)
  rec <- list(scenario = sc$id, b = b, phase = phase, seeds = seeds, y = as.integer(y),
              starts = if (cfg$store_starts) Zs,
              scheduled = c(audit = audit, prof_cov = prof_cov, prof_set = prof_set))
  mod <- bp_model(sc$X, y, link = sc$link)
  num <- cfg$num
  sep <- part("separation", separation_check(sc$X, y))
  cvar <- part("control_variate", control_variate(mod, sc$beta0))
  dir <- separation_direction(sep)                  # additional starts along a separating direction
  fE <- part("fit_E", fit_criterion_search(mod, "E", 0.5, st$routine, st$rescue, dir,
                                           maxit = num$maxit, maxit_rescue = num$maxit_rescue))
  refE <- if (!is.null(fE$selected)) fE$selected$beta
  st2 <- rbind(st$routine, refE)
  fO <- part("fit_O", fit_criterion_search(mod, "O", 0.5, st2, st$rescue, dir,
                                           maxit = num$maxit, maxit_rescue = num$maxit_rescue))
  stE <- rbind(st$routine, fE$ray_starts); stO <- rbind(st2, fO$ray_starts)    # full-fit starts
  fA <- part("fit_A", fit_adjscore_protocol(mod, st2, st$rescue, ref = refE,
                                            maxit = num$maxit, maxit_rescue = num$maxit_rescue))
  fML <- if (identical(sep$category, "overlap"))
    part("fit_ML", fit_criterion_protocol(mod, "ML", 0.5, st$routine, st$rescue,
                                          maxit = num$maxit, maxit_rescue = num$maxit_rescue))
  fits <- list(E = fE, O = fO, A = fA, ML = fML)
  est <- lapply(fits, sel_beta)
  cvec <- if (is.null(cvar)) rep(NA_real_, 3) else cvar
  err <- lapply(est, function(e) e - sc$beta0)                 # NA rows for unresolved fits
  dhat <- c(E = d_contrast(est$E, sc$link), O = d_contrast(est$O, sc$link))
  rec$sep <- if (is.null(sep)) list(category = "error") else record_sep(sep)
  rec$est <- est
  rec$c <- cvec
  rec$derived <- list(err = err, z = lapply(err[c("E", "O", "A")], function(e) e - cvec),
                      D = est$O - est$E, H = err$O - 2 * err$E, ZH = err$O - 2 * err$E + cvec,
                      d0 = sc$d0, d_hat = dhat, T = dhat - sc$d0, W = unname(dhat["O"] - dhat["E"]))
  rec$fits <- list()
  for (k in c("E", "O", "A", "ML")) {
    fr <- part(paste0("record_", k), if (k == "A") record_root(fits[[k]], mod, cfg) else record_fit(fits[[k]], mod, cfg))
    if (!is.null(fr$beta)) fr$reevaluation <- part(paste0("reevaluation_", k), reevaluate(mod, fr, k, cfg))
    rec$fits[k] <- list(fr)
  }
  rec$wald <- list()
  for (k in names(est)) rec$wald[k] <- list(if (!anyNA(est[[k]])) part(paste0("wald_", k), wald_ci(mod, est[[k]])))
  tm["fits"] <- proc.time()[["elapsed"]] - t0
  if (audit) {
    t1 <- proc.time()[["elapsed"]]
    use_stream(seeds$audit)
    Za <- matrix(rnorm(120), 40, 3)
    rec$audit <- part("audit", run_audit(mod, fits, 4 * Za, refE, cfg))
    if (!is.null(rec$audit) && cfg$store_starts) rec$audit$starts <- Za
    tm["audit"] <- proc.time()[["elapsed"]] - t1
  }
  if (prof_cov) {
    t1 <- proc.time()[["elapsed"]]
    rec$pl_cov <- list(E = part("coverage_E", pl_at_truth(mod, fE, sc$beta0[3], stE, st$rescue, cfg)),
                       O = part("coverage_O", pl_at_truth(mod, fO, sc$beta0[3], stO, st$rescue, cfg)))
    tm["coverage"] <- proc.time()[["elapsed"]] - t1
  }
  if (prof_set) {
    t1 <- proc.time()[["elapsed"]]
    rec$pl_set <- list(E = part("set_E", reconstruct_set(mod, fE, stE, st$rescue, cfg)),
                       O = part("set_O", reconstruct_set(mod, fO, stO, st$rescue, cfg)))
    tm["sets"] <- proc.time()[["elapsed"]] - t1
  }
  rec$errors <- errs
  rec$time <- tm
  rec
}

## Last-resort capture: an error outside the component captures (a defect in
## this file) still yields a record with the seeds and the regenerated
## response vector; it counts as unresolved for every component and blocks
## release.
run_job <- function(j) tryCatch(
  analyse_replicate(j$sc, j$b, j$seeds, j$phase, j$audit, j$prof_cov, j$prof_set, j$cfg),
  error = function(e) {
    use_stream(j$seeds$resp)
    list(scenario = j$sc$id, b = j$b, phase = j$phase, seeds = j$seeds,
         y = as.integer(runif(j$sc$n) < j$sc$prob),
         scheduled = c(audit = j$audit, prof_cov = j$prof_cov, prof_set = j$prof_set),
         fatal_error = conditionMessage(e))
  })


## =============================================================================
## Part 4. Running a phase, archive, precision plan, session records
## =============================================================================
## ---- Running a phase ------------------------------------------------------------
code_files <- c("obsinfo_core.R", "obsinfo_simulation.R")
code_md5 <- function() setNames(unname(tools::md5sum(code_files)), code_files)
code_tag <- function(md5) paste0(substr(md5[[1]], 1, 6), substr(md5[[2]], 1, 6))

## Socket (PSOCK) workers on Windows, where forking is unavailable; each loads
## both files by their absolute paths. Set psock = TRUE to use sockets on other platforms.
start_workers <- function(cfg, psock = .Platform$OS.type == "windows") {
  if (cfg$cores <= 1L || !psock) return(NULL)
  files <- normalizePath(code_files, mustWork = TRUE)
  cl <- parallel::makeCluster(cfg$cores)
  loaded <- FALSE
  on.exit(if (!loaded) parallel::stopCluster(cl))          # no orphaned workers if loading fails
  parallel::clusterCall(cl, function(dir, files) {
    setwd(dir); assign(".obsinfo_worker", TRUE, envir = globalenv())
    source(files[1]); source(files[2]); NULL }, getwd(), files)
  loaded <- TRUE
  cl
}
## Forked workers except on Windows. mclapply runs a single job in the calling
## process, so a lone job (e.g. a rerun) is forked explicitly: a worker killed
## by the system then returns NULL instead of ending the run. With cores = 1
## jobs run in this process, without that isolation.
par_run <- function(jobs, cfg, cl) {
  if (!is.null(cl)) return(parallel::parLapplyLB(cl, jobs, run_job))
  if (cfg$cores <= 1L) return(lapply(jobs, run_job))
  if (length(jobs) == 1L) return(unname(parallel::mccollect(parallel::mcparallel(run_job(jobs[[1]])), wait = TRUE)))
  parallel::mclapply(jobs, run_job, mc.cores = cfg$cores, mc.preschedule = FALSE)
}
valid_record <- function(r, b) is.list(r) && identical(as.integer(r$b), as.integer(b)) && length(r$y) > 0

## Chunks of a phase and scenario, with their index ranges and code tags.
list_chunks <- function(cfg, phase, k) {
  f <- list.files(cfg$out_dir, sprintf("^%s_sc%02d_[0-9]{6}-[0-9]{6}_[0-9a-f]{12}\\.rds$", phase, k))
  m <- regmatches(f, regexec("_([0-9]{6})-([0-9]{6})_([0-9a-f]{12})\\.rds$", f))
  data.frame(file = file.path(cfg$out_dir, f), from = as.integer(vapply(m, `[`, "", 2)),
             to = as.integer(vapply(m, `[`, "", 3)), tag = vapply(m, `[`, "", 4), stringsAsFactors = FALSE)
}
## Archived indices; stops on chunks from other code or overlapping ranges.
archived_indices <- function(cfg, phase, k, md5) {
  ch <- list_chunks(cfg, phase, k)
  if (any(ch$tag != code_tag(md5)))
    stop(sprintf("%s scenario %d: chunks written by other code are present; use a new out_dir", phase, k))
  idx <- unlist(Map(seq, ch$from, ch$to))
  if (anyDuplicated(idx)) stop(sprintf("%s scenario %d: overlapping chunks", phase, k))
  idx
}

run_phase <- function(phase, ids, R_of, cfg, scen, streams, cl, md5) {
  nm <- phase_streams[[phase]]
  for (k in ids) {
    R <- R_of(k)
    if (is.na(R)) { cat(sprintf("%s scenario %2d: not run (pilot status not ok)\n", phase, k)); next }
    sc <- scen[[k]]
    todo <- setdiff(seq_len(R), archived_indices(cfg, phase, k, md5))
    if (!length(todo)) next
    S <- lapply(nm, function(x) substreams(streams[[k]][[x]], max(todo)))
    runs <- split(todo, cumsum(c(1, diff(todo) != 1)))
    blocks <- unlist(lapply(runs, function(r) split(r, ceiling(seq_along(r) / cfg$chunk_size))), recursive = FALSE)
    for (blk in blocks) {
      jobs <- lapply(blk, function(b) list(
        sc = sc, b = b, seeds = list(resp = S[[1]][[b]], start = S[[2]][[b]], audit = S[[3]][[b]]),
        phase = phase, audit = phase == "pilot" || b %% cfg$audit_every == 0,
        prof_cov = phase != "pilot" && b <= cfg$n_profile_coverage,
        prof_set = phase != "pilot" && b <= cfg$n_profile_sets, cfg = cfg))
      t0 <- proc.time()[["elapsed"]]
      res <- par_run(jobs, cfg, cl)
      bad <- which(!vapply(seq_along(blk), function(i) valid_record(res[[i]], blk[i]), TRUE))
      if (length(bad)) {                    # e.g. a forked process killed for memory
        cat(sprintf("%s scenario %d: rerunning replicates %s\n", phase, k, paste(blk[bad], collapse = ", ")))
        res[bad] <- par_run(jobs[bad], cfg, cl)
        bad <- which(!vapply(seq_along(blk), function(i) valid_record(res[[i]], blk[i]), TRUE))
        if (length(bad)) stop(sprintf("%s scenario %d: replicates %s returned no record; chunk not written",
                                      phase, k, paste(blk[bad], collapse = ", ")))
      }
      attr(res, "meta") <- list(phase = phase, scenario = k, from = min(blk), to = max(blk), R = R,
                                md5 = md5, code_revision = cfg$code_revision, recon_starts = cfg$recon_starts,
                                num = cfg$num, settings = cfg[c("n_pilot", "n_profile_coverage",
                                "n_profile_sets", "audit_every", "keep_candidates", "store_starts")],
                                written = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), R_version = R.version.string)
      file <- file.path(cfg$out_dir, sprintf("%s_sc%02d_%06d-%06d_%s.rds", phase, k, min(blk), max(blk), code_tag(md5)))
      saveRDS(res, tmp <- paste0(file, ".tmp")); file.rename(tmp, file)
      cat(sprintf("%s scenario %2d (%s %s n = %3d): replicates %d-%d of %d done in %.0f s\n", phase, k,
                  sc$link, sc$setting, sc$n, min(blk), max(blk), R, proc.time()[["elapsed"]] - t0))
    }
  }
}

## Records of a phase and scenario, optionally reduced chunk by chunk by
## `extract`. Checks that all chunks come from the same code and, if R is
## given, whether the indices are exactly 1, ..., R.
read_phase <- function(cfg, phase, k, R = NULL, extract = NULL) {
  ch <- list_chunks(cfg, phase, k)
  if (length(unique(ch$tag)) > 1) stop(sprintf("%s scenario %d: chunks from different code versions", phase, k))
  recs <- vector("list", nrow(ch)); meta <- NULL
  for (i in seq_len(nrow(ch))) {
    x <- readRDS(ch$file[i])
    m <- attr(x, "meta")
    if (is.null(meta)) meta <- m else if (!identical(m$md5, meta$md5))
      stop(sprintf("%s scenario %d: chunks with different code hashes", phase, k))
    recs[[i]] <- if (is.null(extract)) x else lapply(x, extract)
  }
  recs <- unlist(recs, recursive = FALSE)
  b <- vapply(recs, function(r) as.numeric(r$b), 0)
  if (anyDuplicated(b)) stop(sprintf("%s scenario %d: duplicated replicate indices", phase, k))
  recs <- recs[order(b)]
  missing <- if (is.null(R)) integer(0) else setdiff(seq_len(R), b)
  list(records = recs, meta = meta, R = R, missing = missing,
       complete = !is.null(R) && !length(missing) && length(b) == R)
}

## ---- Precision rule -------------------------------------------------
## An incomplete pilot, unresolved fits, errors or a material change in the
## selected maximum on audit preclude a production size for that scenario.
pilot_extract <- function(r) list(
  b = r$b, n_err = as.integer(!is.null(r$fatal_error)) + length(r$errors),
  eE = if (is.null(r$est)) rep(NA_real_, 3) else r$est$E, eO = if (is.null(r$est)) rep(NA_real_, 3) else r$est$O,
  c = if (is.null(r$c)) rep(NA_real_, 3) else r$c,
  mE = isTRUE(r$audit$E$material), mO = isTRUE(r$audit$O$material))
plan_precision <- function(cfg, scen, ids, md5) {
  do.call(rbind, lapply(ids, function(k) {
    sc <- scen[[k]]
    rd <- read_phase(cfg, "pilot", k, R = cfg$n_pilot, extract = pilot_extract)
    rs <- rd$records
    g <- function(f, v) vapply(rs, f, v)
    eE <- matrix(unlist(lapply(rs, `[[`, "eE")), ncol = 3, byrow = TRUE) - rep(sc$beta0, each = length(rs))
    eO <- matrix(unlist(lapply(rs, `[[`, "eO")), ncol = 3, byrow = TRUE) - rep(sc$beta0, each = length(rs))
    C <- matrix(unlist(lapply(rs, `[[`, "c")), ncol = 3, byrow = TRUE)
    n_err <- sum(g(function(r) r$n_err, 0L) > 0)
    uE <- sum(!complete.cases(eE)); uO <- sum(!complete.cases(eO)); uC <- sum(!complete.cases(C))
    mE <- sum(g(function(r) r$mE, TRUE)); mO <- sum(g(function(r) r$mO, TRUE))
    why <- c(if (!rd$complete) "incomplete", if (n_err) "R errors", if (uE + uO + uC) "unresolved fits",
             if (mE + mO) "material audit")
    ok <- !length(why)
    V <- if (ok) c(var(eE[, 3] - C[, 3]), var(eO[, 3] - C[, 3]), var(eO[, 3] - 2 * eE[, 3] + C[, 3])) else rep(NA_real_, 3)
    R_unc <- 1.25 * sc$n^2 * max(V) / 0.05^2
    data.frame(scenario = k, link = sc$link, setting = sc$setting, n = sc$n,
               status = if (ok) "ok" else paste(why, collapse = "; "), records = length(rs),
               replicates_with_errors = n_err, unresolved_E = uE, unresolved_O = uO, material_E = mE,
               material_O = mO, V_ZE2 = V[1], V_ZO2 = V[2], V_ZH2 = V[3], R_uncapped = R_unc,
               R = if (ok) min(1e5, 1000 * ceiling(max(1e4, R_unc) / 1000)) else NA_real_,
               cap_binds = R_unc > 1e5, md5_core = md5[[1]], md5_simulation = md5[[2]],
               code_revision = cfg$code_revision, row.names = NULL)
  }))
}
## Files that concurrent jobs (e.g. a SLURM array) may write: a temporary file
## renamed into place, so that a reader never sees a partial file.
save_atomic <- function(object, file, writer = saveRDS) {
  tmp <- tempfile(pattern = ".tmp_", tmpdir = dirname(file))
  writer(object, tmp)
  if (!file.rename(tmp, file)) {
    unlink(file); if (!file.rename(tmp, file)) { unlink(tmp); stop("cannot write ", file) }
  }
  invisible(file)
}
write_csv_atomic <- function(x, file) save_atomic(x, file, function(o, f) write.csv(o, f, row.names = FALSE))
## The precision plan is kept as one file per scenario (pilot_plan_scXX.rds),
## so that pilot jobs for different scenarios can run at the same time;
## pilot_plan.csv collects the rows present, for reading.
plan_file <- function(cfg, k) file.path(cfg$out_dir, sprintf("pilot_plan_sc%02d.rds", k))
read_plan_files <- function(cfg) {
  f <- sort(list.files(cfg$out_dir, "^pilot_plan_sc[0-9]{2}\\.rds$", full.names = TRUE))
  if (!length(f)) return(NULL)
  do.call(rbind, lapply(f, readRDS))
}
read_plan <- function(cfg, md5) {
  plan <- read_plan_files(cfg)
  if (is.null(plan)) stop("no pilot plan in ", cfg$out_dir, ": run the pilot first")
  if (!all(plan$md5_core == md5[[1]] & plan$md5_simulation == md5[[2]] & plan$code_revision == cfg$code_revision))
    stop("the pilot plan was written by other code: repeat the pilot with this code")
  plan
}

## ---- Session records and archive notes -----------------------------------------
check_out_dir <- function(cfg, md5) {
  f <- file.path(cfg$out_dir, "code_version.rds")
  cur <- list(md5 = md5, code_revision = cfg$code_revision)
  if (!file.exists(f)) return(invisible(save_atomic(c(cur, created = format(Sys.time())), f)))
  old <- readRDS(f)
  if (!identical(old$md5, md5) || !identical(old$code_revision, cfg$code_revision))
    stop(cfg$out_dir, " was written by other code or another code revision; use a new out_dir")
}
log_session <- function(cfg, stage, md5) {
  cat(c(sprintf("==== %s  stage %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), stage),
        paste(names(md5), md5),
        sprintf("code revision %d; recon_starts %s; cores %d", cfg$code_revision, cfg$recon_starts, cfg$cores),
        if (nzchar(Sys.getenv("SLURM_JOB_ID"))) sprintf("SLURM job %s, array task %s, CPUs per task %s",
          Sys.getenv("SLURM_JOB_ID"), Sys.getenv("SLURM_ARRAY_TASK_ID", "-"), Sys.getenv("SLURM_CPUS_PER_TASK", "-")),
        paste(names(Sys.info()), Sys.info(), sep = ": "), capture.output(sessionInfo()), ""),
      file = file.path(cfg$out_dir, "session_log.txt"), sep = "\n", append = TRUE)
}
archive_notes <- c(
  "Chunks contain indexed replicates, per-part status, candidate histories and stream states.",
  "Stage codes: routine, rerun, rescue, Newton, audit, separating direction.",
  "Candidate and profile-point codes follow STAGES, STATUSES, METHODS, PROFILE_KINDS and PROFILE_FLAGS.",
  "Profiles use the full determinant under restriction and a fixed nuisance gradient scale.",
  "Separation records include certificates on the column-scaled design and any unscaled recheck.")

## =============================================================================
## Part 5. Summary statistics and figures
## =============================================================================
## Fields used by the summaries, extracted chunk by chunk and packed into
## plain vectors: a full record takes about 55 KB in memory (more with a
## reconstruction), a reduced one about 3 KB. Coefficient estimates are packed
## as 12 values (blocks E, O, A, ML) and Wald limits as 24 (per method: lower
## then upper, three coefficients each; NA when unavailable).
METH <- c("E", "O", "A", "ML")
slim_record <- function(r) {
  wl <- function(k) { m <- r$wald[[k]]
    if (is.null(m) || anyNA(m[, c("lower", "upper")])) rep(NA_real_, 6) else unname(as.vector(m[, c("lower", "upper")])) }
  au <- r$audit
  list(b = r$b, fatal = !is.null(r$fatal_error), errors = names(r$errors),
       scheduled = c(audit = isTRUE(r$scheduled[["audit"]]), prof_cov = isTRUE(r$scheduled[["prof_cov"]]),
                     prof_set = isTRUE(r$scheduled[["prof_set"]])),
       sep = if (is.null(r$sep)) "error" else r$sep$category,
       est = if (is.null(r$est)) rep(NA_real_, 12) else unname(unlist(r$est[METH])),
       c = if (is.null(r$c)) rep(NA_real_, 3) else unname(r$c),
       d_hat = if (is.null(r$derived)) c(NA_real_, NA_real_) else unname(r$derived$d_hat),
       nd_A = if (is.null(r$fits$A)) NA_integer_ else r$fits$A$n_distinct,
       reeval = vapply(METH, function(k) {
         f <- r$fits[[k]]
         if (is.null(f) || is.null(f$beta)) "no fit" else if (is.null(f$reevaluation)) "not triggered" else
           f$reevaluation$status }, "", USE.NAMES = FALSE),
       audit = if (!is.null(au)) c(nE = au$E$n_distinct_union, nO = au$O$n_distinct_union,
                                   nA = au$A$n_distinct_union, mE = isTRUE(au$E$material),
                                   mO = isTRUE(au$O$material), cE = isTRUE(au$E$selection_changed),
                                   cO = isTRUE(au$O$selection_changed), cA = isTRUE(au$A$selection_changed)),
       wald = unlist(lapply(METH, wl)),
       pl_cov = lapply(r$pl_cov, function(x) if (is.null(x)) NULL else x[c("stat", "flag", "covered")]),
       pl_set = lapply(r$pl_set, function(z) if (is.null(z)) NULL else
         z[intersect(names(z), c("status", "length", "n_components", "components"))]),
       time = if (is.null(r$time)) rep(NA_real_, 4) else unname(r$time))
}
est_of <- function(r, k) r$est[(match(k, METH) - 1L) * 3L + 1:3]
wald_of <- function(r, k) {
  m <- r$wald[(match(k, METH) - 1L) * 6L + 1:6]
  if (anyNA(m)) NULL else matrix(m, 3, 2, dimnames = list(NULL, c("lower", "upper")))
}

wilson <- function(x, M, z = Z975) {
  p <- x / M; d <- 1 + z^2 / M; h <- z * sqrt(p * (1 - p) / M + z^2 / (4 * M^2))
  c((p + z^2 / (2 * M) - h) / d, (p + z^2 / (2 * M) + h) / d)
}
q_interval <- function(W, q) {
  W <- sort(W); M <- length(W); l <- qbinom(0.025, M, q); u <- qbinom(0.975, M, q) + 1
  c(estimate = unname(quantile(W, q)), lower = if (l < 1) -Inf else W[l], upper = if (u > M) Inf else W[u])
}
prop <- function(item, x, M) {
  if (M == 0) return(data.frame(item = item, x = x, M = 0, pct = NA, mcse = NA, lower = NA, upper = NA))
  p <- x / M; w <- wilson(x, M)
  data.frame(item = item, x = x, M = M, pct = 100 * p, mcse = 100 * sqrt(p * (1 - p) / M),
             lower = 100 * w[1], upper = 100 * w[2])
}
precision_status <- function(se) if (is.na(se)) "not assessable" else if (se <= 0.05) "met" else "unmet"

## Separation and numerical diagnostics: numerators, denominators, MCSE and
## Wilson limits. Audit items use the audited samples; an audit that stopped
## with an R error is counted under audit_errors. Re-evaluation items use the
## resolved fits of the method.
diagnostics_table <- function(recs) {
  R <- length(recs)
  cnt <- function(f, rs) sum(vapply(rs, function(r) isTRUE(f(r)), TRUE))
  sepc <- vapply(recs, function(r) r$sep, "")
  acc <- function(k) vapply(recs, function(r) !anyNA(est_of(r, k)), TRUE)
  sched <- vapply(recs, function(r) isTRUE(r$scheduled[["audit"]]), TRUE)
  aud <- Filter(function(r) !is.null(r$audit), recs[sched])
  ree <- function(k, what) cnt(function(r) r$reeval[match(k, METH)] %in% what, recs)
  fitted <- function(k) cnt(function(r) r$reeval[match(k, METH)] != "no fit", recs)
  rbind(
    prop("complete", sum(sepc == "complete"), R),
    prop("quasi_complete", sum(sepc == "quasi-complete"), R),
    prop("ambiguous", sum(sepc == "ambiguous"), R),
    prop("separation_error", sum(sepc == "error"), R),
    prop("accepted_E", sum(acc("E")), R), prop("accepted_O", sum(acc("O")), R), prop("accepted_A", sum(acc("A")), R),
    prop("accepted_ML_overlap", sum(acc("ML") & sepc == "overlap"), sum(sepc == "overlap")),
    prop("multiple_roots_A", cnt(function(r) r$nd_A > 1, recs), R),
    prop("multiple_E", cnt(function(r) r$audit[["nE"]] > 1, aud), length(aud)),
    prop("multiple_O", cnt(function(r) r$audit[["nO"]] > 1, aud), length(aud)),
    prop("audit_extra_roots_A", cnt(function(r) r$audit[["nA"]] > r$nd_A, aud), length(aud)),
    prop("material_E", cnt(function(r) r$audit[["mE"]] == 1, aud), length(aud)),
    prop("material_O", cnt(function(r) r$audit[["mO"]] == 1, aud), length(aud)),
    prop("audit_selection_changed_E", cnt(function(r) r$audit[["cE"]] == 1, aud), length(aud)),
    prop("audit_selection_changed_O", cnt(function(r) r$audit[["cO"]] == 1, aud), length(aud)),
    prop("audit_selection_changed_A", cnt(function(r) r$audit[["cA"]] == 1, aud), length(aud)),
    prop("audit_errors", sum(sched) - length(aud), sum(sched)),
    do.call(rbind, lapply(c("E", "O", "A", "ML"), function(k) rbind(
      prop(paste0("reevaluation_triggered_", k), ree(k, c("agree", "disagree", "failed")), fitted(k)),
      prop(paste0("reevaluation_unresolved_", k), ree(k, c("disagree", "failed")), fitted(k))))),
    prop("replicates_with_R_errors", cnt(function(r) r$fatal || length(r$errors) > 0, recs), R))
}

## Production or restart summaries for one scenario. recs: reduced records
## with indices exactly 1, ..., R. A component that stopped with an R error is
## unresolved for that component; the replicate stays in every denominator.
summarise_scenario <- function(recs, sc) {
  n <- sc$n; b0 <- sc$beta0; R <- length(recs); cn <- c("alpha", "beta1", "beta2")
  qs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  w <- rep(1, R)
  wm <- function(M, rows = rep(TRUE, nrow(M))) colSums(M[rows, , drop = FALSE] * w[rows]) / sum(w[rows])
  getm <- function(f) matrix(as.numeric(unlist(lapply(recs, f))), ncol = 3, byrow = TRUE)
  C <- getm(function(r) r$c); c_ok <- all(complete.cases(C))
  E <- lapply(c(E = "E", O = "O", A = "A", ML = "ML"), function(k) getm(function(r) est_of(r, k)) - rep(b0, each = R))
  ok <- lapply(E, complete.cases)
  sepc <- vapply(recs, function(r) r$sep, "")
  sched <- function(nm) which(vapply(recs, function(r) isTRUE(r$scheduled[[nm]]), TRUE))
  out <- list(diag = diagnostics_table(recs))
  ## Bias and RMSE, all coefficients. The controlled mean is reported only when
  ## every fit of the method and every c are available; otherwise the raw
  ## conditional mean and the accepted-sample mean of c (a selection
  ## diagnostic) are reported.
  out$bias <- do.call(rbind, lapply(c("E", "O", "A"), function(k) {
    e <- E[[k]][ok[[k]], , drop = FALSE]; M <- nrow(e); all_ok <- M == R && c_ok
    V <- e^2; rmse <- sqrt(wm(E[[k]]^2, ok[[k]])); z <- if (all_ok) E[[k]] - C
    top1 <- apply(V, 2, function(v) sum(sort(v, decreasing = TRUE)[seq_len(ceiling(0.01 * length(v)))]) / sum(v))
    data.frame(method = k, coef = cn, generated = R, M = M, all_resolved = all_ok,
               raw = wm(E[[k]], ok[[k]]), raw_se = apply(e, 2, sd) / sqrt(M), var_raw = apply(e, 2, var),
               cv = if (all_ok) wm(z) else NA, cv_se = if (all_ok) apply(z, 2, sd) / sqrt(M) else NA,
               var_cv = if (all_ok) apply(z, 2, var) else NA,
               cv_identity = if (all_ok) wm(E[[k]]) - wm(z) - wm(C) else NA,   # zero up to rounding
               cbar_accepted = wm(C, ok[[k]]),
               rmse = rmse, rmse_se = ifelse(rmse > 0, apply(V, 2, sd) / (2 * sqrt(M) * rmse), NA),
               max_abs_err = apply(abs(e), 2, max), top1pct_mse_share = top1, row.names = NULL)
  }))
  ## Paired displacement D and H = E_O - 2 E_E on common resolved pairs; the
  ## controlled n Z_H when every pair and every c are available.
  pr <- ok$E & ok$O; Mp <- sum(pr); all_pairs <- Mp == R; all_cv <- all_pairs && c_ok
  Dm <- E$O - E$E; Hm <- E$O - 2 * E$E
  D <- Dm[pr, , drop = FALSE]; H <- Hm[pr, , drop = FALSE]
  ZH <- if (all_cv) Hm + C
  dsq <- (E$O^2 - E$E^2)[pr, , drop = FALSE]
  sdi <- sqrt(diag(sc$ref$i_inv))
  out$paired <- data.frame(coef = cn, generated = R, M = Mp, all_resolved = all_pairs, all_resolved_cv = all_cv,
    nD = n * wm(Dm, pr), nD_se = n * apply(D, 2, sd) / sqrt(Mp), nbE = n * sc$ref$b_E,
    D_mean = wm(Dm, pr), D_mean_over_sd = wm(Dm, pr) / sdi, bE_over_sd = sc$ref$b_E / sdi,
    nH = n * wm(Hm, pr), nH_se = n * apply(H, 2, sd) / sqrt(Mp),
    nZH = if (all_cv) n * wm(ZH) else NA, nZH_se = if (all_cv) n * apply(ZH, 2, sd) / sqrt(R) else NA,
    sqerr_diff = colMeans(dsq), sqerr_diff_se = apply(dsq, 2, sd) / sqrt(Mp),
    sqerr_diff_lower = colMeans(dsq) - MC_Z975 * apply(dsq, 2, sd) / sqrt(Mp),
    sqerr_diff_upper = colMeans(dsq) + MC_Z975 * apply(dsq, 2, sd) / sqrt(Mp), row.names = NULL)
  out$quantiles <- if (Mp) do.call(rbind, lapply(1:3, function(r) rbind(
    data.frame(coef = cn[r], quantity = "D", q = qs, all_resolved = all_pairs,
               t(vapply(qs, function(p) q_interval(D[, r], p), numeric(3)))),
    data.frame(coef = cn[r], quantity = "sqrt(n)(nD - nbE)", q = qs, all_resolved = all_pairs,
               t(vapply(qs, function(p) q_interval(sqrt(n) * (n * D[, r] - n * sc$ref$b_E[r]), p), numeric(3)))))))
  bz <- out$bias[out$bias$coef == "beta2", ]
  se3 <- c(n * bz$cv_se[bz$method == "E"], n * bz$cv_se[bz$method == "O"], out$paired$nZH_se[3])
  out$precision <- data.frame(se_nZE2 = se3[1], se_nZO2 = se3[2], se_nZH2 = se3[3],
                              met_E = se3[1] <= 0.05, met_O = se3[2] <= 0.05, met_H = se3[3] <= 0.05,
                              status_E = precision_status(se3[1]), status_O = precision_status(se3[2]),
                              status_H = precision_status(se3[3]))
  ## Probability difference d at the design centre, percentage points.
  dh <- matrix(as.numeric(unlist(lapply(recs, `[[`, "d_hat"))), ncol = 2, byrow = TRUE,
               dimnames = list(NULL, c("E", "O")))
  dsum <- function(label, x, paired) {
    x <- x[is.finite(x)]; M <- length(x); ms <- mean(x^2)
    data.frame(quantity = label, generated = R, M = M, all_resolved = M == R, d0_pp = 100 * sc$d0,
               mean_pp = 100 * mean(x), mean_se_pp = 100 * sd(x) / sqrt(M),
               rmse_pp = if (paired) NA else 100 * sqrt(ms),
               rmse_se_pp = if (paired) NA else 100 * sd(x^2) / (2 * sqrt(M * ms)))
  }
  out$prob <- rbind(dsum("bias E", dh[, "E"] - sc$d0, FALSE), dsum("bias O", dh[, "O"] - sc$d0, FALSE),
                    dsum("paired O-E", dh[, "O"] - dh[, "E"], TRUE))
  ## Wald intervals: all samples, and the profile index sets for Wald-profile
  ## comparisons. eligible: samples of the index set (overlap samples for ML);
  ## unresolved: eligible samples without an interval (unresolved fit or R error).
  all_i <- seq_len(R); cov_i <- sched("prof_cov"); set_i <- sched("prof_set")
  wald_one <- function(k, idx, label) {
    el <- if (k == "ML") sum(sepc[idx] == "overlap") else length(idx)
    wl <- Filter(Negate(is.null), lapply(recs[idx], wald_of, k = k))
    M <- length(wl)
    if (!M) return(NULL)
    cv <- t(vapply(wl, function(m) as.numeric(m[, "lower"] <= b0 & b0 <= m[, "upper"]), numeric(3)))
    ln <- t(vapply(wl, function(m) m[, "upper"] - m[, "lower"], numeric(3)))
    p <- colMeans(cv); wi <- t(vapply(colSums(cv), wilson, numeric(2), M = M))
    data.frame(method = k, coef = cn, index_set = label, eligible = el, M = M, unresolved = el - M,
               coverage = 100 * p, coverage_mcse = 100 * sqrt(p * (1 - p) / M), lower = 100 * wi[, 1],
               upper = 100 * wi[, 2], length = colMeans(ln), length_se = apply(ln, 2, sd) / sqrt(M), row.names = NULL)
  }
  out$wald <- rbind(do.call(rbind, lapply(c("E", "O", "A", "ML"), wald_one, idx = all_i, label = "all")),
                    do.call(rbind, lapply(c("E", "O"), wald_one, idx = cov_i, label = "coverage indices")),
                    do.call(rbind, lapply(c("E", "O"), wald_one, idx = set_i, label = "reconstruction indices")))
  ## Penalized-profile coverage at the truth (C covered, F unresolved among the
  ## M prescribed samples) and set reconstruction, beta2.
  tabstr <- function(x) if (length(x)) paste(names(table(x)), table(x), sep = " ", collapse = "; ") else ""
  cov_of <- function(r, k) { x <- r$pl_cov[[k]]; if (is.null(x)) list(flag = "error", covered = NA) else x }
  set_of <- function(r, k) { z <- r$pl_set[[k]]; if (is.null(z)) list(status = "error") else z }
  out$profile <- do.call(rbind, lapply(c("E", "O"), function(k) {
    xs <- lapply(recs[cov_i], cov_of, k = k)
    v <- vapply(xs, function(x) as.logical(x$covered), NA)
    M <- length(v); Cn <- sum(v %in% TRUE); Fn <- sum(is.na(v)); pc <- Cn / (M - Fn)
    s <- lapply(recs[set_i], set_of, k = k)
    oks <- Filter(function(z) identical(z$status, "resolved"), s); len <- vapply(oks, function(z) z$length, 0)
    data.frame(method = k, M = M, covered = Cn, unresolved = Fn, conditional = 100 * pc,
               conditional_mcse = 100 * sqrt(pc * (1 - pc) / (M - Fn)),
               range_low = 100 * Cn / M, range_high = 100 * (Cn + Fn) / M,
               envelope_low = 100 * wilson(Cn, M)[1], envelope_high = 100 * wilson(Cn + Fn, M)[2],
               coverage_flags = tabstr(vapply(xs, function(x) x$flag, "")),
               sets = length(s), sets_resolved = length(oks), length = mean(len), length_se = sd(len) / sqrt(length(len)),
               disconnected = sum(vapply(oks, function(z) z$n_components > 1, TRUE)),
               reconstruction_status = tabstr(vapply(s, function(z) z$status, "")))
  }))
  ## Paired interval contrasts on common resolved samples: lO - lE
  ## for each procedure, and Wald - profile on the same prespecified indices.
  wcov <- function(k, idx) vapply(recs[idx], function(r) { m <- wald_of(r, k)
    if (is.null(m)) NA_real_ else as.numeric(m[3, "lower"] <= b0[3] && b0[3] <= m[3, "upper"]) }, 0)
  wlen <- function(k, idx) vapply(recs[idx], function(r) { m <- wald_of(r, k)
    if (is.null(m)) NA_real_ else m[3, "upper"] - m[3, "lower"] }, 0)
  pcov <- function(k, idx) vapply(recs[idx], function(r) as.numeric(cov_of(r, k)$covered), 0)
  plen <- function(k, idx) vapply(recs[idx], function(r) { z <- set_of(r, k)
    if (identical(z$status, "resolved")) z$length else NA_real_ }, 0)
  pdiff <- function(a, b, contrast, procedure, measure, index_set) {
    okp <- is.finite(a) & is.finite(b); d <- a[okp] - b[okp]; M <- length(d)
    s <- if (measure == "coverage") 100 else 1; m <- mean(d); se <- sd(d) / sqrt(M)
    data.frame(contrast = contrast, procedure = procedure, measure = measure, index_set = index_set,
               eligible = length(a), M = M, mean = s * m, se = s * se, lower = s * (m - MC_Z975 * se), upper = s * (m + MC_Z975 * se))
  }
  out$paired_intervals <- rbind(
    pdiff(wcov("O", all_i), wcov("E", all_i), "O-E", "Wald", "coverage", "all"),
    pdiff(wlen("O", all_i), wlen("E", all_i), "O-E", "Wald", "length", "all"),
    pdiff(pcov("O", cov_i), pcov("E", cov_i), "O-E", "profile", "coverage", "coverage indices"),
    pdiff(plen("O", set_i), plen("E", set_i), "O-E", "profile", "length", "reconstruction indices"),
    pdiff(wcov("E", cov_i), pcov("E", cov_i), "Wald-profile", "lE", "coverage", "coverage indices"),
    pdiff(wcov("O", cov_i), pcov("O", cov_i), "Wald-profile", "lO", "coverage", "coverage indices"),
    pdiff(wlen("E", set_i), plen("E", set_i), "Wald-profile", "lE", "length", "reconstruction indices"),
    pdiff(wlen("O", set_i), plen("O", set_i), "Wald-profile", "lO", "length", "reconstruction indices"))
  ## Separation strata: raw conditional means and RMSE of every coefficient
  ## with their Monte Carlo standard errors, stratum denominators.
  out$strata <- do.call(rbind, lapply(c("complete", "quasi-complete", "overlap"), function(g)
    do.call(rbind, lapply(c("E", "O", "A"), function(k) {
      e <- E[[k]][sepc == g & ok[[k]], , drop = FALSE]; M <- nrow(e); V <- e^2; rmse <- sqrt(colMeans(V))
      data.frame(stratum = g, method = k, coef = cn, generated = sum(sepc == g), M = M,
                 bias = if (M) colMeans(e) else NA, bias_se = if (M > 1) apply(e, 2, sd) / sqrt(M) else NA,
                 rmse = if (M) rmse else NA,
                 rmse_se = if (M > 1) ifelse(rmse > 0, apply(V, 2, sd) / (2 * sqrt(M) * rmse), NA) else NA,
                 row.names = NULL)
    }))))
  ## ML panel: common accepted overlap samples for every method.
  ov <- sepc == "overlap" & ok$E & ok$O & ok$A & ok$ML
  out$overlap <- do.call(rbind, lapply(c("E", "O", "A", "ML"), function(k) {
    e <- E[[k]][ov, , drop = FALSE]; M <- nrow(e); V <- e^2; rmse <- sqrt(colMeans(V))
    data.frame(method = k, coef = cn, overlap_generated = sum(sepc == "overlap"), M = M, bias = colMeans(e),
               bias_se = apply(e, 2, sd) / sqrt(M), rmse = rmse, rmse_se = apply(V, 2, sd) / (2 * sqrt(M) * rmse),
               row.names = NULL)
  }))
  tms <- matrix(as.numeric(unlist(lapply(recs, `[[`, "time"))), ncol = 4, byrow = TRUE,
                dimnames = list(NULL, c("fits", "audit", "coverage", "sets")))
  out$time <- c(colMeans(tms, na.rm = TRUE),
                audit_per_sample = mean(tms[which(tms[, "audit"] > 0), "audit"]),
                coverage_per_sample = mean(tms[which(tms[, "coverage"] > 0), "coverage"]),
                sets_per_sample = mean(tms[which(tms[, "sets"] > 0), "sets"]))
  out
}

## Sensitivity analysis: fits whose stable re-evaluation
## disagreed or failed are treated as unresolved (lE, lO, S + A).
mask_reeval <- function(r) {
  for (i in 1:3) if (r$reeval[i] %in% c("disagree", "failed")) r$est[(i - 1L) * 3L + 1:3] <- NA_real_
  r
}

## Completeness, numerical outcomes and agreement of the pilot and production archives.
release_row <- function(s, rd, plan_row, phase) {
  dg <- s$diag; x <- function(item) dg$x[dg$item == item]
  md5 <- unname(rd$meta$md5)
  why <- c(
    if (!is.null(plan_row) && plan_row$status != "ok") paste("pilot status:", plan_row$status),
    if (!is.null(plan_row) && plan_row$material_E + plan_row$material_O > 0) "material pilot audit",
    if (!rd$complete) "archive incomplete",
    if (!is.null(plan_row) && !identical(c(plan_row$md5_core, plan_row$md5_simulation), md5))
      "pilot and production code differ",
    if (x("replicates_with_R_errors") > 0) "R errors",
    if (x("material_E") + x("material_O") > 0) "material audit")
  ree <- sum(dg$x[grepl("^reevaluation_unresolved_", dg$item)])
  data.frame(phase = phase, generated = if (is.null(rd$R)) NA else rd$R, archived = length(rd$records),
             md5_core = md5[1], md5_simulation = md5[2],
             unresolved_E = length(rd$records) - x("accepted_E"), unresolved_O = length(rd$records) - x("accepted_O"),
             unresolved_A = length(rd$records) - x("accepted_A"), reevaluation_unresolved = ree,
             material_audits = x("material_E") + x("material_O"),
             pilot_material_audits = if (is.null(plan_row)) NA else plan_row$material_E + plan_row$material_O,
             replicates_with_R_errors = x("replicates_with_R_errors"),
             blocked = length(why) > 0, reasons = paste(why, collapse = "; "))
}

## ---- Tables and figures -------------------------------------------------------------
summarise_all <- function(cfg, scen, ids, phase = "production", R_of = NULL) {
  sdir <- file.path(cfg$out_dir, if (phase == "production") "summary" else paste0("summary_", phase))
  dir.create(sdir, showWarnings = FALSE, recursive = TRUE)
  plan <- if (is.null(R_of)) read_plan_files(cfg)
  if (is.null(R_of) && is.null(plan)) stop("No pilot plan found in ", cfg$out_dir)
  S <- list(); done_ids <- integer(0)
  for (k in ids) {
    prow <- if (!is.null(plan)) plan[plan$scenario == k, ]
    R <- if (is.null(R_of)) prow$R else R_of(k)
    if (!length(R) || is.na(R)) { cat(sprintf("scenario %d: no production size (pilot status not ok)\n", k)); next }
    rd <- read_phase(cfg, phase, k, R = R, extract = slim_record)
    if (!rd$complete) stop(sprintf("%s scenario %d: archive incomplete (%d of %d indices)", phase, k,
                                   length(rd$records), R))
    s <- summarise_scenario(rd$records, scen[[k]])
    s$release <- release_row(s, rd, prow, phase)
    flagged <- vapply(rd$records, function(r) any(r$reeval[1:3] %in% c("disagree", "failed")), TRUE)
    if (any(flagged)) {                      # sensitivity analysis for unresolved re-evaluations
      s2 <- summarise_scenario(lapply(rd$records, mask_reeval), scen[[k]])
      s$reeval_sensitivity <- list(bias = s2$bias, paired = s2$paired)
    }
    S[[length(S) + 1L]] <- s; done_ids <- c(done_ids, k)
  }
  if (!length(S)) { cat("no scenario to summarise\n"); return(invisible(NULL)) }
  lab <- function(k) data.frame(scenario = k, link = scen[[k]]$link, setting = scen[[k]]$setting, n = scen[[k]]$n)
  bind <- function(f) do.call(rbind, Map(function(s, k) {
    x <- f(s); if (is.null(x) || !nrow(x)) NULL else cbind(lab(k)[rep(1, nrow(x)), ], x, row.names = NULL) }, S, done_ids))
  tabs <- list(diagnostics = bind(function(s) s$diag), bias = bind(function(s) s$bias),
               paired = bind(function(s) s$paired), quantiles = bind(function(s) s$quantiles),
               probability = bind(function(s) s$prob), wald = bind(function(s) s$wald),
               profile = bind(function(s) s$profile), paired_intervals = bind(function(s) s$paired_intervals),
               strata = bind(function(s) s$strata), overlap_ml = bind(function(s) s$overlap),
               precision = bind(function(s) s$precision), release = bind(function(s) s$release),
               time = bind(function(s) as.data.frame(t(s$time))),
               bias_reeval_sensitivity = bind(function(s) s$reeval_sensitivity$bias),
               paired_reeval_sensitivity = bind(function(s) s$reeval_sensitivity$paired))
  tabs <- tabs[!vapply(tabs, is.null, TRUE)]
  for (nm in names(tabs)) write.csv(tabs[[nm]], file.path(sdir, paste0(nm, ".csv")), row.names = FALSE)
  plot_figure1(tabs$bias, tabs$paired, tabs$quantiles, file.path(sdir, "figure_bias.pdf"),
               file.path(sdir, "figure_bias_companion.pdf"))
  invisible(tabs)
}

## Pilot outputs: planning variances (in the plan), failure counts,
## audits and timing only. No pilot contrasts are summarised.
summarise_pilot_all <- function(cfg, scen, ids, plan) {
  sdir <- file.path(cfg$out_dir, "summary_pilot"); dir.create(sdir, showWarnings = FALSE, recursive = TRUE)
  lab <- function(k) data.frame(scenario = k, link = scen[[k]]$link, setting = scen[[k]]$setting, n = scen[[k]]$n)
  for (k in ids) {
    rs <- read_phase(cfg, "pilot", k, R = cfg$n_pilot, extract = slim_record)$records
    if (!length(rs)) next
    d <- diagnostics_table(rs)
    write_csv_atomic(cbind(lab(k)[rep(1, nrow(d)), ], d, row.names = NULL),
                     file.path(sdir, sprintf("pilot_diagnostics_sc%02d.csv", k)))
  }
  f <- sort(list.files(sdir, "^pilot_diagnostics_sc[0-9]{2}\\.csv$", full.names = TRUE))
  dg <- do.call(rbind, lapply(f, read.csv))
  write_csv_atomic(dg, file.path(sdir, "pilot_diagnostics.csv"))
  write_csv_atomic(read_plan_files(cfg), file.path(sdir, "pilot_plan.csv"))
  save_atomic("Pilot planning variances, fit counts and audits.",
              file.path(sdir, "README.txt"), writeLines)
  invisible(dg)
}

## Figures 1 and S1. The inputs are the bias, paired and quantile summaries.
plot_figure1 <- function(bias, paired, quant, file, file2) {
  bias <- subset(bias, coef == "beta2")
  paired <- subset(paired, coef == "beta2")
  quant <- if (NROW(quant)) subset(quant, coef == "beta2" & quantity == "sqrt(n)(nD - nbE)")
  panels <- unique(paired[c("link", "setting")])
  label <- c(probit = "Probit", cloglog = "Complementary log-log")
  colours <- c(E = "#246596", O = "#B84D1C", D = "#24704C")
  take <- function(x, k) {
    if (!NROW(x)) return(x)
    rows <- x$link == panels$link[k] & x$setting == panels$setting[k]
    x <- x[rows, , drop = FALSE]
    x[order(x$n), , drop = FALSE]
  }
  panel <- function(k, values, ylabel) {
    yr <- range(c(0, values), finite = TRUE)
    plot(c(log2(25) - 0.23, log2(200) + 0.23), yr, type = "n", xaxt = "n",
         xlab = "n", ylab = ylabel,
         main = paste(label[panels$link[k]], panels$setting[k], sep = "\n"))
    sizes <- c(25, 50, 100, 200)
    axis(1, at = log2(sizes), labels = sizes)
    abline(h = 0, col = "grey75", lty = 3)
  }
  bounds <- function(y, se) c(y - MC_Z975 * se, y + MC_Z975 * se)
  estimate <- function(n, y, se, offset, symbol, colour) {
    x <- log2(n) + offset
    lines(x, y, col = colour, lwd = 0.7)
    lo <- y - MC_Z975 * se
    hi <- y + MC_Z975 * se
    segments(x, lo, x, hi, col = colour)
    segments(x - 0.035, lo, x + 0.035, lo, col = colour)
    segments(x - 0.035, hi, x + 0.035, hi, col = colour)
    points(x, y, pch = symbol, col = colour)
  }
  draw_bias <- function() {
    pdf(file, width = 8.2, height = 7.5)
    on.exit(dev.off(), add = TRUE)
    par(mfrow = c(2, 2), mar = c(3.6, 4.1, 3.3, 0.8), oma = c(2.7, 0, 0, 0))
    for (k in seq_len(nrow(panels))) {
      b <- take(bias, k)
      e <- subset(b, method == "E" & all_resolved)
      o <- subset(b, method == "O" & all_resolved)
      p <- take(paired, k)
      d <- subset(p, all_resolved)
      panel(k, c(bounds(e$n * e$cv, e$n * e$cv_se), bounds(o$n * o$cv, o$n * o$cv_se),
                 bounds(d$nD, d$nD_se), p$nbE, 2 * p$nbE), expression(n %*% "bias of " * beta[2]))
      lines(log2(p$n), p$nbE, lty = 2)
      lines(log2(p$n), 2 * p$nbE, lty = 3)
      estimate(e$n, e$n * e$cv, e$n * e$cv_se, 0, 16, colours["E"])
      estimate(o$n, o$n * o$cv, o$n * o$cv_se, 0, 17, colours["O"])
      estimate(d$n, d$nD, d$nD_se, 0, 15, colours["D"])
    }
    par(fig = c(0, 1, 0, 0.075), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), new = TRUE)
    plot.new()
    legend("center", c("Expected", "Observed", "Mean paired difference", "n b(E)", "2 n b(E)"),
           pch = c(16, 17, 15, NA, NA), lty = c(1, 1, 1, 2, 3),
           col = c(colours, "black", "black"), ncol = 3, bty = "n", cex = 0.75)
  }
  draw_bias()
  pdf(file2, width = 11.5, height = 6.5)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(2, nrow(panels)), mar = c(3.5, 3.8, 3.3, 0.8))
  for (k in seq_len(nrow(panels))) {
    p <- take(paired, k)
    h <- subset(p, all_resolved)
    z <- subset(p, all_resolved_cv)
    panel(k, c(bounds(h$nH, h$nH_se), bounds(z$nZH, z$nZH_se)), expression(n %*% "mean of " * H))
    estimate(h$n, h$nH, h$nH_se, 0.07, 1, "black")
    estimate(z$n, z$nZH, z$nZH_se, -0.07, 16, colours["E"])
  }
  for (k in seq_len(nrow(panels))) {
    q <- take(quant, k)
    if (NROW(q)) q <- subset(q, all_resolved)
    panel(k, q$estimate, expression(sqrt(n) * (n * D - n * b^E)))
    for (size in unique(q$n)) {
      sample_q <- q[q$n == size, ]
      values <- sample_q$estimate[match(c(0.025, 0.25, 0.5, 0.75, 0.975), sample_q$q)]
      x <- log2(size)
      segments(x, values[1], x, values[5])
      segments(x - 0.07, values[c(1, 5)], x + 0.07, values[c(1, 5)])
      rect(x - 0.15, values[2], x + 0.15, values[4], col = "grey92")
      points(x, values[3], pch = 4)
    }
  }
  invisible(NULL)
}

## First-order bias of the probability difference at x1 = 0, in percentage points.
probability_bias_terms <- function(sc) {
  V <- rbind(c(1, 0, 0), c(1, 0, 1))
  link <- make_link(sc$link)$base(drop(V %*% sc$beta0))
  g <- exp(link$lg)
  gp <- g * link$d1
  variance <- rowSums((V %*% sc$ref$i_inv) * V)
  C <- diff(gp * variance) / 2
  L <- diff(g * drop(V %*% sc$ref$b_E))
  data.frame(link = sc$link, setting = sc$setting, n = sc$n,
             C = 100 * C, L = 100 * L, E = 100 * (C + L), O = 100 * (C + 2 * L))
}

paper_tables <- function() {
  old <- options(width = max(getOption("width"), 120))
  on.exit(options(old), add = TRUE)
  sample_sizes <- c(25, 100, 400, 1600)

  fmt_tab <- function(M, digits = 4)
    noquote(matrix(formatC(M, digits = digits, format = "g"), nrow(M), dimnames = dimnames(M)))
  round_df <- function(d, digits = 4) { k <- sapply(d, is.numeric); d[k] <- lapply(d[k], round, digits); d }
  link_label <- c(logit = "Logit", probit = "Probit", cloglog = "Cloglog", loglog = "Loglog", cauchit = "Cauchit")

  ## Starts for the two-parameter examples: 0, +-2 e_r, 10 N(0, 4 I) and 40
  ## N(0, 16 I) draws (the protocol's routine and rescue starts), used together.
  set.seed(20261101)
  st2 <- protocol_starts(2, matrix(rnorm(100), 50, 2))
  starts2 <- rbind(st2$routine, st2$rescue)

  ## Tables 2 and S2
  xc <- c(-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2); yc <- as.numeric(xc > 0)
  xq <- c(-2, -1.5, -1, 0, 0, 1, 1.5, 2);       yq <- c(0, 0, 0, 0, 1, 1, 1, 1)
  cat("=== Tables 2 and S2: penalized estimates for the separated data sets ===\n")
  cat("Separation (linear programming): complete data:", separation_check(cbind(1, xc), yc)$category,
      "| quasi data:", separation_check(cbind(1, xq), yq)$category, "\n\n")
  est_names <- c("lE(1/2)", "lO(1/2)", "S+A=0", "lR_i N=n", "lR_j N=n", "lR_i N=n^.5", "lR_j N=n^.5")
  fit_separated_fits <- function(x, y, link) {
    mod <- bp_model(cbind(1, x), y, link = link)
    n <- mod$n
    fE <- fit_criterion_protocol(mod, "E", 0.5, starts2)
    st <- rbind(starts2, fE$selected$beta)
    fits <- list(fE,
                 fit_criterion_protocol(mod, "O", 0.5, st),
                 fit_adjscore_protocol(mod, st, ref = fE$selected$beta),
                 fit_criterion_protocol(mod, "Ri", 1 / n, st),
                 fit_criterion_protocol(mod, "Rj", 1 / n, st),
                 fit_criterion_protocol(mod, "Ri", 1 / sqrt(n), st),
                 fit_criterion_protocol(mod, "Rj", 1 / sqrt(n), st))
    est <- sapply(fits, function(f) f$selected$beta)
    dimnames(est) <- list(c("Intercept", "Slope"), est_names)
    list(model = mod, fits = fits, est = est, n_distinct = sapply(fits, `[[`, "n_distinct"))
  }
  separated <- list()
  for (lk in c("logit", "probit", "cloglog", "loglog", "cauchit"))
    for (ds in c("Complete", "Quasi"))
      separated[[paste(lk, ds)]] <- if (ds == "Complete") fit_separated_fits(xc, yc, lk) else fit_separated_fits(xq, yq, lk)
  rows <- list()
  for (lk in c("logit", "probit", "cloglog", "cauchit")) for (ds in c("Complete", "Quasi")) {
    e <- separated[[paste(lk, ds)]]$est
    if (lk == "cloglog") rows[[paste(link_label[lk], ds, "Intercept")]] <- e["Intercept", ]
    rows[[paste(link_label[lk], ds, "Slope")]] <- e["Slope", ]
  }
  separated_fits <- do.call(rbind, rows)
  print(round(separated_fits, 4))

  ## Table S3
  cat("\n=== Table S3: log-likelihood and penalty (1/2) log det j along b +/- t (0, 1), complete data ===\n")
  ray <- function(res, tt = c(0, 2, 5, 10, 20), gam = c(0, 1)) {
    b <- res$fits[[2]]$selected$beta                       # computed maximum of l^O(.; 1/2)
    f <- function(bb) { ev <- bp_eval(res$model, bb, "O", 0.5, 0L); c(ev$loglik, ev$penalty) }
    out <- t(sapply(tt, function(t) c(t, f(b + t * gam), f(b - t * gam))))
    colnames(out) <- c("t", "l(b+tg)", "pen(b+tg)", "l(b-tg)", "pen(b-tg)")
    list(b = b, table = out)
  }
  for (lk in c("probit", "cloglog")) {
    r3 <- ray(separated[[paste(lk, "Complete")]])
    cat(sprintf("%s: b = (%.4f, %.4f)\n", link_label[lk], r3$b[1], r3$b[2]))
    print(fmt_tab(r3$table), right = TRUE)
  }

  ## Supplementary Section D.1
  cat("\n=== Supplementary Section D.1: l - tr(j^{-1} e)/2 along beta = (0, t), probit, complete data ===\n")
  tt <- c(1, 2, 4, 8)
  trace_criterion <- sapply(tt, function(t) bp_eval(separated[["probit Complete"]]$model, c(0, t), "Ri", 0, 0L)$value)
  print(fmt_tab(cbind(t = tt, value = trace_criterion)), right = TRUE)

  ## Table S1
  cat("\n=== Table S1: Probit slope vs multiplier a on (xc, yc) ===\n")
  modp <- separated[["probit Complete"]]$model
  a_grid <- c(0.5, 1, 2, 4, 6, 8, 16, 64)
  slope_sweep <- t(sapply(a_grid, function(a) {
    fE <- fit_criterion_protocol(modp, "E", a, starts2)
    fO <- fit_criterion_protocol(modp, "O", a, rbind(starts2, fE$selected$beta))
    c(a = a, slope_lO = fO$selected$beta[2], slope_lE = fE$selected$beta[2])
  }))
  print(round(slope_sweep, 4))
  g_l <- bp_eval(modp, c(0, 0), "ML", level = 1L)$grad
  g_ld <- bp_eval(modp, c(0, 0), "O", a = 1, level = 1L)$grad - g_l
  cat(sprintf("Origin stationary for l^O(.; a) at a = %.2f; at a = 1/2 the slopes differ by %.1f%%\n",
              -g_l[2] / g_ld[2], 100 * (slope_sweep[1, 2] / slope_sweep[1, 3] - 1)))

  ## Table 3
  cat("\n=== Table 3: intercept-only models; n x first-order bias (FO) and n x exact bias by enumeration ===\n")
  exact_biases <- NULL
  for (cf in list(list("probit", 0.5), list("cloglog", -0.5), list("cloglog", 0.5))) {
    fo <- first_order_bias(matrix(1, 1, 1), cf[[2]], cf[[1]], m = 1)   # n x first-order bias, free of n
    for (n in sample_sizes) {
      bE <- exact_bias_intercept(cf[[1]], cf[[2]], n, "E")$bias
      bO <- exact_bias_intercept(cf[[1]], cf[[2]], n, "O")$bias
      exact_biases <- rbind(exact_biases, data.frame(Link = link_label[cf[[1]]], eta0 = cf[[2]], FO_ML = fo$b_ML, FO_lE = fo$b_E,
                                     FO_lO = fo$b_O, n = n, nbias_lE = n * bE, nbias_lO = n * bO, Ratio = bO / bE))
    }
  }
  print(round_df(exact_biases), row.names = FALSE)

  ## Table S4
  cat("\n=== Table S4: first-order paired difference 10^3 x i^{-1}(zetaE - A); x equally spaced on [-1.6, 1.6] ===\n")
  displacements <- NULL
  for (cf in list(list("probit", c(0.4, 0.8)), list("cloglog", c(-0.4, 0.8)))) for (n in c(30, 60, 120, 240)) {
    d <- 1e3 * first_order_bias(cbind(1, seq(-1.6, 1.6, length.out = n)), cf[[2]], cf[[1]], a = 0.5)$displacement
    displacements <- rbind(displacements, data.frame(Link = link_label[cf[[1]]], n = n, Intercept = d[1], Slope = d[2]))
  }
  print(round_df(displacements, 3), row.names = FALSE)

  ## Table S6 and Supplementary Section D.3
  cat("\n=== Table S6: strict local maxima of the probit penalized log-likelihoods (data of Table S5) ===\n")
  d7 <- data.frame(y = c(0, 0, 1, 0, 1, 0, 0), x1 = c(0.15, 0.37, -0.01, -0.87, 0.35, 0.13, -1.68),
                   x2 = c(1, 1, 0, 1, 1, 1, 1))
  d12 <- data.frame(y = c(0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0),
                    x1 = c(-1.10, -1.31, -0.39, -0.57, 2.56, 0.08, -0.85, -0.72, -0.21, -0.23, 0.40, 1.30),
                    x2 = c(0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1))
  Z3 <- matrix(rnorm(350 * 3), 350, 3)
  st3 <- protocol_starts(3, Z3[1:50, ])
  starts3 <- rbind(st3$routine, st3$rescue, 3 * Z3[51:350, ])     # 357 starts
  maxima <- function(d, type, a, label) {
    mod <- bp_model(cbind(1, d$x1, d$x2), d$y, link = "probit")
    f <- fit_criterion_protocol(mod, type, a, starts3)
    out <- do.call(rbind, lapply(f$distinct, function(s) data.frame(
      Data = label, Objective = type, alpha = s$beta[1], beta1 = s$beta[2], beta2 = s$beta[3], Value = s$value,
      lambda_max = max(eigen(bp_eval(mod, s$beta, type, a, 2L)$hess, symmetric = TRUE, only.values = TRUE)$values),
      starts = s$n_starts)))
    out[order(-out$Value), ]
  }
  local_maxima <- rbind(maxima(d7, "E", 0.5, "n = 7"), maxima(d7, "O", 0.5, "n = 7"),
                maxima(d12, "E", 0.5, "n = 12"), maxima(d12, "O", 0.5, "n = 12"))
  local_maxima$Objective <- c(E = "lE(1/2)", O = "lO(1/2)")[local_maxima$Objective]
  print(round_df(local_maxima), row.names = FALSE)
  cat("starts = number of the 357 starts that reached each maximum\n")

  cat("\n=== Supplementary Section D.3: local maxima of l^R_i with N = n, n = 12 data ===\n")
  empirical_maxima <- maxima(d12, "Ri", 1 / 12, "n = 12")
  empirical_maxima$Objective <- "lR_i N=n"
  print(round_df(empirical_maxima, 5), row.names = FALSE)
  cat(sprintf("Difference between the two near-tied values: %.1e\n",
              abs(diff(sort(empirical_maxima$Value, decreasing = TRUE)[1:2]))))

  cat("\nTable S7: first-order probability-difference biases (percentage points)\n")
  print(do.call(rbind, lapply(make_scenarios(), probability_bias_terms)), row.names = FALSE)
  invisible(NULL)
}

main <- function(stage, ids = 1:16, phase = "production") {
  if (length(stage) != 1L || !stage %in% c("pilot", "production", "restart", "summary", "tables"))
    stop("Unknown stage. See --help.")
  if (!all(file.exists(code_files))) stop("Run from the directory containing both obsinfo R files.")
  if (!length(ids) || anyNA(ids) || any(!ids %in% 1:16)) stop("Scenario IDs must be in 1:16.")
  if (!phase %in% c("production", "restart")) stop("Summary phase must be production or restart.")
  if (stage == "tables") return(invisible(paper_tables()))
  scen <- make_scenarios()
  streams <- if (stage != "summary") make_streams(cfg$seed)
  rc <- cfg
  if (stage == "summary" && !dir.exists(rc$out_dir)) stop("Simulation archive directory missing: ", rc$out_dir)
  if (stage == "production" && rc$code_revision != 0L) stop("Use the standard production configuration.")
  if (stage == "restart" && rc$code_revision != 1L) stop("Reserved restart streams require code_revision = 1.")
  if (stage != "summary" && rc$recon_starts != "protocol") stop("Simulation profile fits require all specified starts.")
  md5 <- code_md5()
  dir.create(rc$out_dir, showWarnings = FALSE, recursive = TRUE)
  if (stage != "summary") {
    check_out_dir(rc, md5)
    save_atomic(archive_notes, file.path(rc$out_dir, "archive_notes.txt"), writeLines)
    save_atomic(lapply(scen, function(s) s[c("id", "link", "setting", "n", "X", "beta0", "prob", "ref", "d0")]),
                file.path(rc$out_dir, "design_record.rds"))
  }
  log_session(rc, stage, md5)
  cl <- if (stage != "summary") start_workers(rc)
  if (!is.null(cl)) on.exit(parallel::stopCluster(cl))
  if (stage == "pilot") {
    run_phase("pilot", ids, function(k) rc$n_pilot, rc, scen, streams, cl, md5)
    plan <- plan_precision(rc, scen, ids, md5)
    print(plan[, c("scenario", "link", "setting", "n", "status", "V_ZE2", "V_ZO2", "V_ZH2", "R_uncapped", "R")])
    for (i in seq_len(nrow(plan))) save_atomic(plan[i, ], plan_file(rc, plan$scenario[i]))
    write_csv_atomic(read_plan_files(rc), file.path(rc$out_dir, "pilot_plan.csv"))
    summarise_pilot_all(rc, scen, ids, plan)
  }
  if (stage %in% c("production", "restart")) {
    plan <- read_plan(rc, md5)
    R_of <- function(k) {
      row <- plan[plan$scenario == k, ]
      if (nrow(row) == 1L && row$status == "ok") row$R else NA_real_
    }
    run_phase(stage, ids, R_of, rc, scen, streams, cl, md5)
  }
  if (stage == "summary") return(invisible(summarise_all(rc, scen, ids, phase = phase)))
  invisible(NULL)
}

scenario_ids <- function(spec) {
  if (length(spec) != 1L || !grepl("^[0-9]+(:[0-9]+)?(,[0-9]+(:[0-9]+)?)*$", spec))
    stop("Use scenario numbers or ranges, for example 1,3:6.")
  decode <- function(token) {
    x <- suppressWarnings(as.integer(strsplit(token, ":", fixed = TRUE)[[1L]]))
    if (anyNA(x) || any(x < 1L | x > 16L)) stop("Scenario IDs must be in 1:16.")
    seq.int(x[1L], x[length(x)])
  }
  unique(unlist(lapply(strsplit(spec, ",", fixed = TRUE)[[1L]], decode), use.names = FALSE))
}

if (sys.nframe() == 0L && !exists(".obsinfo_worker")) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || args[1L] %in% c("-h", "--help", "help")) {
    cat("Run from the directory containing both obsinfo R files.\n",
        "Usage: Rscript obsinfo_simulation.R STAGE [SCENARIOS [PHASE]]\n",
        "  tables       Fixed-design examples and first-order calculations\n",
        "  pilot        Independent precision pilot for the selected scenarios\n",
        "  production   Main simulation using the pilot replication counts\n",
        "  summary      Summaries from saved replicate blocks\n",
        "  restart      Reserved streams for a separately configured restart\n",
        "SCENARIOS defaults to 1:16. PHASE is production (default) or restart.\n", sep = "")
  } else {
    if (length(args) > 3L) stop("Too many arguments. See --help.")
    ids <- if (length(args) >= 2L) scenario_ids(args[2L]) else seq_len(16L)
    main(args[1L], ids, if (length(args) >= 3L) args[3L] else "production")
  }
}
