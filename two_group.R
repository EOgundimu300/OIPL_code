#!/usr/bin/env Rscript
# Reproduce the two-group enumeration in the supplementary material.
# From this directory run Rscript two_group.R. Only base R is required.
source("obsinfo_core.R")
group_sizes <- list(c(17L, 8L), c(35L, 15L))

# A fitted intercept for a given response total does not depend on its true value.
intercept_fits <- function(size, link, method) {
    vapply(0:size, function(total) {
        fit <- fit_intercept_only(total, size, link, type = method, a = 0.5)
        stopifnot(is.finite(fit$beta), fit$f2 < 0, fit$std_grad < 1e-7)
        fit$beta
    }, numeric(1))
}

results <- list()
for (link_name in c("probit", "cloglog")) {
    link <- make_link(link_name)
    for (sizes in group_sizes) {
        paired_fits <- lapply(sizes, function(size) {
            intercept_fits(size, link, "O") - intercept_fits(size, link, "E")
        })
        for (setting in names(protocol_truth)) {
            truth <- protocol_truth[[setting]]
            eta <- c(truth[1], truth[1] + truth[3])
            group_means <- vapply(seq_along(sizes), function(g) {
                size <- sizes[g]
                prob <- dbinom(0:size, size, link$prob(eta[g]))
                sum(prob * paired_fits[[g]])
            }, numeric(1))
            first_order <- vapply(seq_along(sizes), function(g) {
                # Use one binomial total, the same likelihood as repeated Bernoulli rows.
                first_order_bias(matrix(1, 1, 1), eta[g], link,
                                 m = sizes[g], a = 0.5)$displacement[1]
            }, numeric(1))
            mean_difference <- diff(group_means)
            approximation <- diff(first_order)
            results[[length(results) + 1L]] <- data.frame(
                link = link_name, setting = setting, n0 = sizes[1], n1 = sizes[2],
                mean = mean_difference, first_order = approximation,
                ratio = mean_difference / approximation)
        }
    }
}
print(do.call(rbind, results), digits = 7, row.names = FALSE)
