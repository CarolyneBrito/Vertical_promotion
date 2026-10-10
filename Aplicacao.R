# =============================================================================
# Vertical promotion-time cure model for the melanoma competing-risks data
# Estimation by the EM algorithm
#
# Data coding (boot::melanoma):
#   status = 1: death from melanoma
#   status = 2: censored/alive
#   status = 3: death from other causes
#
# Model:
#   N_i | x_i ~ Poisson(theta_i), log(theta_i) = x_i' gamma
#   T_i = min(Y_i1,...,Y_iN_i), with N_i = 0 defining cure
#   Y_ij has a piecewise-exponential distribution WITHOUT covariates
#   P(D_i = 3 | T_i=t, event, q_i) is logistic (vertical component)
#
# The script needs only the recommended packages boot and survival. All
# numerical differentiation is performed with base R.
# =============================================================================

required <- c("boot", "survival", "dplyr", "xtable", "tidyverse")
missing_packages <- required[!vapply(required, requireNamespace,
                                     quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages)) {
  stop("Install the following packages first: ",
       paste(missing_packages, collapse = ", "))
}

pacman::p_load("boot", "survival", "dplyr", "xtable", "tidyverse")

# ---------- Utilities for a piecewise-exponential activation distribution ----

make_intervals <- function(time, event, K = 4L, cuts = NULL) {
  stopifnot(length(time) == length(event), all(time > 0), K >= 1L)
  if (is.null(cuts)) {
    probs <- seq(0, 1, length.out = K + 1L)
    inner <- as.numeric(stats::quantile(time[event == 1L], probs[-c(1, K + 1L)],
                                        names = FALSE, type = 2))
    cuts <- c(0, unique(inner), Inf)
  } else {
    cuts <- sort(unique(c(0, cuts[cuts > 0 & is.finite(cuts)], Inf)))
  }
  if (length(cuts) < 2L) stop("The interval partition is invalid.")
  cuts
}

exposure_matrix <- function(time, cuts) {
  K <- length(cuts) - 1L
  E <- matrix(0, nrow = length(time), ncol = K)
  for (k in seq_len(K)) {
    E[, k] <- pmax(0, pmin(time, cuts[k + 1L]) - cuts[k])
  }
  E
}

interval_index <- function(time, cuts) {
  findInterval(time, cuts, rightmost.closed = TRUE, all.inside = TRUE)
}

activation_terms <- function(time, lambda, cuts, E = NULL) {
  if (is.null(E)) E <- exposure_matrix(time, cuts)
  k <- interval_index(time, cuts)
  H0 <- as.vector(E %*% lambda)
  S0 <- exp(-H0)
  F0 <- -expm1(-H0)                 # numerically stable 1-exp(-H0)
  f0 <- lambda[k] * S0
  list(H0 = H0, S0 = S0, F0 = F0, f0 = f0, interval = k, E = E)
}

# -------------------------- Observed log-likelihood --------------------------

loglik_total_time <- function(gamma, lambda, X, time, event, cuts, E = NULL) {
  if (any(!is.finite(gamma)) || any(!is.finite(lambda)) || any(lambda <= 0)) {
    return(-Inf)
  }
  eta <- as.vector(X %*% gamma)
  if (any(eta > 700)) return(-Inf)
  theta <- exp(eta)
  a <- activation_terms(time, lambda, cuts, E)
  if (any(event == 1L & a$f0 <= 0)) return(-Inf)
  sum(event * (eta + log(a$f0)) - theta * a$F0)
}

loglik_cause <- function(beta, Z, event, cause_other) {
  ind <- event == 1L
  eta <- as.vector(Z[ind, , drop = FALSE] %*% beta)
  # y*eta-log(1+exp(eta)), evaluated stably
  sum(cause_other[ind] * eta - pmax(eta, 0) - log1p(exp(-abs(eta))))
}

loglik_joint_unconstrained <- function(par, X, Z, time, event, cause_other,
                                       cuts, E = NULL) {
  p <- ncol(X); K <- length(cuts) - 1L; q <- ncol(Z)
  gamma <- par[seq_len(p)]
  log_lambda <- par[p + seq_len(K)]
  beta <- par[p + K + seq_len(q)]
  loglik_total_time(gamma, exp(log_lambda), X, time, event, cuts, E) +
    loglik_cause(beta, Z, event, cause_other)
}

# ------------------------------- EM algorithm --------------------------------

fit_vertical_pt_em <- function(data, K = 4L, cuts = NULL,
                               count_formula = ~ sexo_reg + IDADE + ano + catend_reg  + 
                                 ecgrup_2_3_reg + ecgrup_4_reg + m_reg + 
                                 RADIO + QUIMIO + CIRURGIA,
                               cause_formula = ~ log_time + sexo_reg + IDADE + ano + catend_reg  + 
                                 ecgrup_2_3_reg + ecgrup_4_reg + m_reg + 
                                 RADIO + QUIMIO + CIRURGIA,
                               start = NULL, tol = 1e-6, maxit = 1000L,
                               trace = FALSE, compute_se = TRUE, scaling = NULL) {
  vars <- c("time", "status", "sexo_reg", "IDADE", "ano", "catend_reg", 
            "ecgrup_2_3_reg", "ecgrup_4_reg", "m_reg", 
            "RADIO", "QUIMIO", "CIRURGIA")
  if (!all(vars %in% names(data))) stop("Missing variables: ",
                                        paste(setdiff(vars, names(data)), collapse = ", "))
  dat <- data[, vars]
  dat <- dat[stats::complete.cases(dat), ]
  if (!all(dat$status %in% 1:3)) stop("status must take values 1, 2 or 3.")
  if (any(dat$time <= 0)) stop("All observed times must be positive.")
  
  raw_dat <- dat                       # retained in days for bootstrap resampling
  dat$time <- dat$time / 365.25        # all modeling is performed in years
  # Scaling improves the conditioning of the Hessian. Coefficients refer to
  # one standard deviation for continuous covariates.
  scale_info <- scaling
  if (is.null(scale_info)) {
    scale_info <- lapply(c("IDADE", "ano"), function(v) {
      c(center = mean(dat[[v]]), scale = stats::sd(dat[[v]]))
    })
    names(scale_info) <- c("IDADE", "ano")
  }
  for (v in names(scale_info)) {
    dat[[v]] <- (dat[[v]] - scale_info[[v]]["center"]) / scale_info[[v]]["scale"]
  }
  dat$log_time <- log(dat$time)
  event <- as.integer(dat$status != 2L)
  cause_other <- as.integer(dat$status == 3L)
  X <- stats::model.matrix(count_formula, dat)
  Z <- stats::model.matrix(cause_formula, dat)
  if (qr(X)$rank < ncol(X)) stop("The count-component design matrix is not full rank.")
  
  cuts <- make_intervals(dat$time, event, K, cuts)
  K <- length(cuts) - 1L
  E <- exposure_matrix(dat$time, cuts)
  ik <- interval_index(dat$time, cuts)
  d <- tabulate(ik[event == 1L], nbins = K)
  if (any(d == 0L)) {
    stop("At least one activation-time interval has no failures. Reduce K or supply other cuts.")
  }
  
  if (is.null(start)) {
    gamma <- rep(0, ncol(X))
    gamma[1] <- log(max(-log(mean(event == 0L)), 0.2))
    lambda <- d / pmax(colSums(E), .Machine$double.eps)
  } else {
    gamma <- start$gamma
    lambda <- start$lambda
  }
  if (length(gamma) != ncol(X) || length(lambda) != K || any(lambda <= 0)) {
    stop("Invalid starting values.")
  }
  
  ll_history <- numeric(maxit + 1L)
  ll_history[1] <- loglik_total_time(gamma, lambda, X, dat$time, event, cuts, E)
  converged <- FALSE
  
  for (iter in seq_len(maxit)) {
    theta <- exp(as.vector(X %*% gamma))
    a <- activation_terms(dat$time, lambda, cuts, E)
    
    # E-step: E(N_i | observed data)
    zhat <- event + theta * a$S0
    
    # M-step for gamma: Poisson pseudo-likelihood with fractional response.
    gamma_fit <- stats::glm.fit(x = X, y = zhat,
                                family = stats::poisson(link = "log"),
                                start = gamma, control = stats::glm.control(maxit = 100L))
    gamma_new <- as.vector(gamma_fit$coefficients)
    if (any(!is.finite(gamma_new))) stop("The gamma M-step failed; inspect scaling/separation.")
    
    # M-step for the common, covariate-free activation hazard.
    lambda_new <- d / pmax(colSums(E * zhat), .Machine$double.eps)
    ll_new <- loglik_total_time(gamma_new, lambda_new, X, dat$time, event, cuts, E)
    ll_history[iter + 1L] <- ll_new
    
    rel_change <- abs(ll_new - ll_history[iter]) / (1 + abs(ll_history[iter]))
    if (trace) message(sprintf("EM %4d: logLik = %.8f; relative change = %.3e",
                               iter, ll_new, rel_change))
    gamma <- gamma_new; lambda <- lambda_new
    if (rel_change < tol) { converged <- TRUE; break }
  }
  ll_history <- ll_history[seq_len(iter + 1L)]
  
  # Vertical cause component: among failures, model cause 3 versus cause 1.
  ind <- event == 1L
  if (length(unique(cause_other[ind])) < 2L) stop("Both failure causes must be observed.")
  cause_fit <- stats::glm.fit(x = Z[ind, , drop = FALSE], y = cause_other[ind],
                              family = stats::binomial(link = "logit"),
                              control = stats::glm.control(maxit = 100L))
  beta <- as.vector(cause_fit$coefficients)
  if (any(!is.finite(beta)) || !cause_fit$converged) {
    stop("The cause model did not converge; possible complete/quasi-complete separation.")
  }
  
  theta <- exp(as.vector(X %*% gamma))
  cure <- exp(-theta)
  a <- activation_terms(dat$time, lambda, cuts, E)
  zhat <- event + theta * a$S0
  
  fit <- list(call = match.call(), data = dat, raw_data = raw_dat,
              scale_info = scale_info, X = X, Z = Z,
              count_formula = count_formula, cause_formula = cause_formula,
              gamma = gamma, lambda = lambda, beta = beta, cuts = cuts,
              event = event, cause_other = cause_other, zhat = zhat,
              theta = theta, cure = cure, logLik = tail(ll_history, 1) +
                loglik_cause(beta, Z, event, cause_other),
              logLik_history = ll_history, iterations = iter,
              converged = converged, d_interval = d,
              cause_converged = cause_fit$converged, se = NULL)
  class(fit) <- "vertical_pt_em"
  
  if (compute_se) fit$se <- observed_information(fit)
  fit
}

# ---------------- Standard errors from observed information ------------------

observed_information <- function(object, condition_limit = 1e10) {
  stopifnot(inherits(object, "vertical_pt_em"))
  p <- length(object$gamma); K <- length(object$lambda); q <- length(object$beta)
  par_hat <- c(object$gamma, log(object$lambda), object$beta)
  E <- exposure_matrix(object$data$time, object$cuts)
  nll <- function(par) -loglik_joint_unconstrained(
    par, object$X, object$Z, object$data$time, object$event,
    object$cause_other, object$cuts, E)
  H <- stats::optimHess(par_hat, nll)
  H <- (H + t(H)) / 2
  eig <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  positive_definite <- all(is.finite(eig)) && min(eig) > 0
  condition_number <- if (positive_definite) max(eig) / min(eig) else Inf
  reliable <- positive_definite && condition_number < condition_limit
  V_u <- if (positive_definite) solve(H) else matrix(NA_real_, length(par_hat), length(par_hat))
  
  # Delta transformation from log(lambda) to lambda.
  J <- diag(length(par_hat)); J[p + seq_len(K), p + seq_len(K)] <- object$lambda
  V <- J %*% V_u %*% J
  estimates <- c(object$gamma, object$lambda, object$beta)
  nm <- c(paste0("gamma:", colnames(object$X)),
          paste0("lambda[", seq_len(K), "]"),
          paste0("beta:", colnames(object$Z)))
  names(estimates) <- nm; dimnames(V) <- list(nm, nm)
  se <- sqrt(pmax(diag(V), 0)); names(se) <- nm
  table <- data.frame(parameter = nm, estimate = estimates, std.error = se,
                      lower.95 = estimates - 1.96 * se,
                      upper.95 = estimates + 1.96 * se,
                      row.names = NULL)
  
  # Cure probability and delta-method SE for every individual:
  # dp_i/dgamma = -p_i theta_i x_i.
  Vgamma <- V[seq_len(p), seq_len(p), drop = FALSE]
  Gcure <- -object$cure * object$theta * object$X
  se_cure <- sqrt(pmax(rowSums((Gcure %*% Vgamma) * Gcure), 0))
  cure_table <- data.frame(cure = object$cure, std.error = se_cure,
                           lower.95 = pmax(0, object$cure - 1.96 * se_cure),
                           upper.95 = pmin(1, object$cure + 1.96 * se_cure))
  list(vcov = V, vcov_unconstrained = V_u, information = H,
       eigenvalues = eig, condition_number = condition_number,
       positive_definite = positive_definite, reliable = reliable,
       coefficients = table, cure = cure_table)
}

# ---------------------------- Bootstrap standard errors ----------------------

bootstrap_em <- function(object, B = 500L, seed = 2026L, trace = TRUE) {
  stopifnot(inherits(object, "vertical_pt_em"), B >= 2L)
  set.seed(seed)
  n <- nrow(object$data)
  original_names <- c(paste0("gamma:", colnames(object$X)),
                      paste0("lambda[", seq_along(object$lambda), "]"),
                      paste0("beta:", colnames(object$Z)))
  draws <- matrix(NA_real_, B, length(original_names), dimnames = list(NULL, original_names))
  failures <- character(B)
  
  for (b in seq_len(B)) {
    id <- sample.int(n, n, replace = TRUE)
    ans <- try(fit_vertical_pt_em(
      object$raw_data[id, ], cuts = object$cuts,
      count_formula = object$count_formula, cause_formula = object$cause_formula,
      tol = 1e-8, maxit = 1000L, compute_se = FALSE,
      scaling = object$scale_info), silent = TRUE)
    if (!inherits(ans, "try-error") && ans$converged) {
      draws[b, ] <- c(ans$gamma, ans$lambda, ans$beta)
    } else {
      failures[b] <- if (inherits(ans, "try-error")) as.character(ans) else "EM did not converge"
    }
    if (trace && (b %% 25L == 0L || b == B)) message("Bootstrap: ", b, "/", B)
  }
  ok <- stats::complete.cases(draws)
  if (sum(ok) < max(30L, ceiling(0.8 * B))) {
    warning("Few valid bootstrap replications: ", sum(ok), " of ", B, ".")
  }
  list(B = B, successful = sum(ok), failure_rate = mean(!ok), draws = draws[ok, , drop = FALSE],
       std.error = apply(draws[ok, , drop = FALSE], 2, stats::sd),
       percentile.95 = t(apply(draws[ok, , drop = FALSE], 2, stats::quantile,
                               probs = c(.025, .975), na.rm = TRUE)),
       failures = failures[!ok])
}

# ------------------------- Predictions and summaries -------------------------

predict.vertical_pt_em <- function(object, times, newdata = NULL, type = c("survival", "cif"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    X <- object$X; Zbase <- object$data
  } else {
    Zbase <- newdata
    for (v in names(object$scale_info)) {
      Zbase[[v]] <- (Zbase[[v]] - object$scale_info[[v]]["center"]) /
        object$scale_info[[v]]["scale"]
    }
    Zbase$log_time <- 0 # placeholder; updated over the CIF integration grid
    X <- stats::model.matrix(object$count_formula, Zbase)
  }
  times <- sort(unique(times[times >= 0]))
  theta <- exp(as.vector(X %*% object$gamma))
  A <- activation_terms(times, object$lambda, object$cuts)
  S <- outer(theta, A$F0, function(th, F) exp(-th * F))
  if (type == "survival") return(S)
  
  # Numerical integration of theta*f0(t)*pi_j(t) for the two causes.
  grid <- sort(unique(c(0, times, object$cuts[is.finite(object$cuts)])))
  if (length(grid) < 2L) stop("At least one positive prediction time is required.")
  # Individual CIFs are returned by a direct trapezoidal integration grid.
  fine <- sort(unique(c(grid, seq(0, max(times), length.out = 1001L))))
  fine[1] <- max(fine[1], 1e-8)
  af <- activation_terms(fine, object$lambda, object$cuts)
  out1 <- out3 <- matrix(0, nrow(X), length(times))
  for (i in seq_len(nrow(X))) {
    zd <- Zbase[rep(i, length(fine)), , drop = FALSE]
    zd$log_time <- log(fine)
    Zi <- stats::model.matrix(object$cause_formula, zd)
    pi3 <- stats::plogis(as.vector(Zi %*% object$beta)); pi1 <- 1 - pi3
    dens <- theta[i] * af$f0 * exp(-theta[i] * af$F0)
    cum1 <- c(0, cumsum(diff(fine) * (head(dens * pi1, -1) + tail(dens * pi1, -1)) / 2))
    cum3 <- c(0, cumsum(diff(fine) * (head(dens * pi3, -1) + tail(dens * pi3, -1)) / 2))
    out1[i, ] <- stats::approx(fine, cum1, xout = times, rule = 2)$y
    out3[i, ] <- stats::approx(fine, cum3, xout = times, rule = 2)$y
  }
  list(cause1_melanoma = out1, cause3_other = out3, times = times)
}

# Quantities used in Figures 2--4 of Nicolaie et al. (2019), adapted to the
# promotion-time model. For N_i ~ Poisson(theta_i), the susceptible population
# is the subpopulation N_i > 0 and has survival
# [S_pop(t)-exp(-theta_i)]/[1-exp(-theta_i)].
predict_profile_quantities <- function(object, profile, times) {
  times <- sort(unique(times[times >= 0]))
  S_pop <- as.numeric(predict(object, times, newdata = profile, type = "survival"))
  cif <- predict(object, times, newdata = profile, type = "cif")
  
  prof_scaled <- profile
  for (v in names(object$scale_info)) {
    prof_scaled[[v]] <- (prof_scaled[[v]] - object$scale_info[[v]]["center"]) /
      object$scale_info[[v]]["scale"]
  }
  prof_scaled$log_time <- 0
  Xp <- stats::model.matrix(object$count_formula, prof_scaled)
  theta <- exp(as.numeric(Xp %*% object$gamma))
  p_cure <- exp(-theta)
  S_sus <- pmax(0, (S_pop - p_cure) / (1 - p_cure))
  F1_pop <- as.numeric(cif$cause1_melanoma)
  F3_pop <- as.numeric(cif$cause3_other)
  F1_sus <- F1_pop / (1 - p_cure)
  F3_sus <- F3_pop / (1 - p_cure)
  
  zd <- profile[rep(1, length(times)), , drop = FALSE]
  for (v in names(object$scale_info)) {
    zd[[v]] <- (zd[[v]] - object$scale_info[[v]]["center"]) /
      object$scale_info[[v]]["scale"]
  }
  zd$log_time <- log(pmax(times, 1e-8))
  Zp <- stats::model.matrix(object$cause_formula, zd)
  pi_other <- stats::plogis(as.vector(Zp %*% object$beta))
  
  data.frame(time = times, survival_population = S_pop,
             survival_susceptible = S_sus,
             cif_melanoma_population = F1_pop,
             cif_other_population = F3_pop,
             cif_melanoma_susceptible = F1_sus,
             cif_other_susceptible = F3_sus,
             pi_melanoma = 1 - pi_other, pi_other = pi_other,
             cure_probability = p_cure)
}

# Bootstrap for both model parameters and profile-specific curves. The same
# resamples are used throughout, so the uncertainty summaries are coherent.
bootstrap_full <- function(object, profile, times, B = 500L, seed = 2026L,
                           trace = TRUE, min_success = 0.80) {
  stopifnot(inherits(object, "vertical_pt_em"), B >= 2L)
  set.seed(seed)
  n <- nrow(object$raw_data)
  par_names <- c(paste0("gamma:", colnames(object$X)),
                 paste0("lambda[", seq_along(object$lambda), "]"),
                 paste0("beta:", colnames(object$Z)))
  curve_names <- setdiff(names(predict_profile_quantities(object, profile, times)), "time")
  par_draws <- matrix(NA_real_, B, length(par_names), dimnames = list(NULL, par_names))
  curve_draws <- array(NA_real_, c(B, length(times), length(curve_names)),
                       dimnames = list(NULL, NULL, curve_names))
  failures <- character(B)
  
  for (b in seq_len(B)) {
    id <- sample.int(n, n, replace = TRUE)
    ans <- try(fit_vertical_pt_em(
      object$raw_data[id, ], cuts = object$cuts,
      count_formula = object$count_formula, cause_formula = object$cause_formula,
      tol = 1e-8, maxit = 1000L, compute_se = FALSE,
      scaling = object$scale_info), silent = TRUE)
    if (!inherits(ans, "try-error") && ans$converged) {
      pq <- try(predict_profile_quantities(ans, profile, times), silent = TRUE)
      if (!inherits(pq, "try-error")) {
        par_draws[b, ] <- c(ans$gamma, ans$lambda, ans$beta)
        curve_draws[b, , ] <- as.matrix(pq[, curve_names, drop = FALSE])
      } else failures[b] <- as.character(pq)
    } else {
      failures[b] <- if (inherits(ans, "try-error")) as.character(ans) else "EM did not converge"
    }
    if (trace && (b %% 25L == 0L || b == B)) message("Full bootstrap: ", b, "/", B)
  }
  ok <- stats::complete.cases(par_draws)
  if (mean(ok) < min_success) {
    warning(sprintf("Only %.1f%% of bootstrap fits were successful.", 100 * mean(ok)))
  }
  pd <- par_draws[ok, , drop = FALSE]
  parameter_table <- data.frame(
    parameter = par_names,
    estimate = c(object$gamma, object$lambda, object$beta),
    bootstrap.se = apply(pd, 2, stats::sd),
    lower.95 = apply(pd, 2, stats::quantile, probs = .025, na.rm = TRUE),
    upper.95 = apply(pd, 2, stats::quantile, probs = .975, na.rm = TRUE),
    row.names = NULL)
  
  estimate <- predict_profile_quantities(object, profile, times)
  bands <- list()
  for (j in seq_along(curve_names)) {
    arr <- curve_draws[ok, , j, drop = FALSE][, , 1]
    bands[[curve_names[j]]] <- data.frame(
      time = times, estimate = estimate[[curve_names[j]]],
      lower = apply(arr, 2, stats::quantile, probs = .025, na.rm = TRUE),
      upper = apply(arr, 2, stats::quantile, probs = .975, na.rm = TRUE))
  }
  list(B = B, successful = sum(ok), failure_rate = mean(!ok),
       parameter_draws = pd, parameter_table = parameter_table,
       estimate = estimate, bands = bands, failures = failures[!ok])
}

add_band <- function(band, col, alpha = 0.18) {
  graphics::polygon(c(band$time, rev(band$time)),
                    c(band$lower, rev(band$upper)),
                    col = grDevices::adjustcolor(col, alpha.f = alpha), border = NA)
}

plot_nicolaie_figures <- function(object, boot_object, file = NULL,
                                  width = 10, height = 8) {
  if (!is.null(file)) grDevices::pdf(file, width = width, height = height)
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit({graphics::par(oldpar); if (!is.null(file)) grDevices::dev.off()}, add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.6, 1), las = 1)
  cols <- c(melanoma = "#0072B2", other = "#D55E00", population = "#009E73",
            susceptible = "#CC79A7")
  
  # Figure 1 analogue: overall Kaplan--Meier curve.
  km <- survival::survfit(survival::Surv(object$data$time, object$event) ~ 1)
  graphics::plot(km, conf.int = TRUE, mark.time = TRUE, xlab = "Years since diagnosis",
                 ylab = "Surviving fraction", main = "(a) Overall survival (Kaplan--Meier)",
                 xlim = range(boot_object$estimate$time), ylim = c(0, 1),
                 col = "black", lwd = 2)
  
  # Figure 2 analogue: population cumulative incidences.
  b1 <- boot_object$bands$cif_melanoma_population
  b2 <- boot_object$bands$cif_other_population
  graphics::plot(b1$time, b1$estimate, type = "n", ylim = c(0, 1),
                 xlab = "Years since diagnosis", ylab = "Population cumulative incidence",
                 main = "(b) Population cumulative incidences")
  add_band(b1, cols["melanoma"]); add_band(b2, cols["other"])
  graphics::lines(b1$time, b1$estimate, col = cols["melanoma"], lwd = 2)
  graphics::lines(b2$time, b2$estimate, col = cols["other"], lwd = 2)
  graphics::legend("topleft", c("Primary bone cancer", "Other causes"),
                   col = cols[c("melanoma", "other")], lwd = 2, bty = "n")
  
  # Figure 3 analogue: population versus susceptible survival.
  bs1 <- boot_object$bands$survival_population
  bs2 <- boot_object$bands$survival_susceptible
  graphics::plot(bs1$time, bs1$estimate, type = "n", ylim = c(0, 1),
                 xlab = "Years since diagnosis", ylab = "Surviving fraction",
                 main = "(c) Promotion-time survival curves")
  add_band(bs1, cols["population"]); add_band(bs2, cols["susceptible"])
  graphics::lines(bs1$time, bs1$estimate, col = cols["population"], lwd = 2)
  graphics::lines(bs2$time, bs2$estimate, col = cols["susceptible"], lwd = 2, lty = 2)
  graphics::abline(h = boot_object$estimate$cure_probability[1], lty = 3, col = "gray35")
  graphics::legend("topright", c("Population", "Susceptible (N>0)", "Cure fraction"),
                   col = c(cols["population"], cols["susceptible"], "gray35"),
                   lwd = c(2, 2, 1), lty = c(1, 2, 3), bty = "n")
  
  # Figure 4 analogue: cumulative incidences among susceptible patients.
  bc1 <- boot_object$bands$cif_melanoma_susceptible
  bc2 <- boot_object$bands$cif_other_susceptible
  graphics::plot(bc1$time, bc1$estimate, type = "n", ylim = c(0, 1),
                 xlab = "Years since diagnosis", ylab = "Conditional cumulative incidence",
                 main = "(d) Cumulative incidences given N>0")
  add_band(bc1, cols["melanoma"]); add_band(bc2, cols["other"])
  graphics::lines(bc1$time, bc1$estimate, col = cols["melanoma"], lwd = 2)
  graphics::lines(bc2$time, bc2$estimate, col = cols["other"], lwd = 2)
  graphics::legend("topleft", c("Primary bone cancer", "Other causes"),
                   col = cols[c("melanoma", "other")], lwd = 2, bty = "n")
  invisible(boot_object)
}

print.vertical_pt_em <- function(x, ...) {
  cat("Vertical promotion-time cure model fitted by EM\n")
  cat("n =", nrow(x$data), "; failures =", sum(x$event),
      "; EM iterations =", x$iterations, "; converged =", x$converged, "\n")
  cat("log-likelihood =", format(x$logLik, digits = 8), "\n")
  if (!is.null(x$se)) {
    cat("Observed information positive definite =", x$se$positive_definite,
        "; condition number =", format(x$se$condition_number, digits = 4), "\n")
  }
  invisible(x)
}

# ------------------------------- Data analysis -------------------------------

mel <- dados_org %>% select(time, status, sexo_reg, IDADE, ano, catend_reg,  
                            ecgrup_2_3_reg, ecgrup_4_reg, m_reg, RADIO, QUIMIO, CIRURGIA) %>% 
  mutate(status = case_when(status == 0 ~ 2,
                            status == 1 ~ 1,
                            status == 2 ~ 3))

# Four intervals (in years) determined by quartiles of observed failure times. If an
# interval is sparse, use K=3 or manually specify clinically meaningful cuts.
# Ajuste dos modelos com K = 1,...,6
modelos <- lapply(1:10, function(k) {
  cat("\n-----------------------------------\n")
  cat("Iniciando modelo com K =", k, "\n")
  cat("Horário de início:", format(Sys.time(), "%H:%M:%S"), "\n")
  inicio <- Sys.time()
  
  fit <- fit_vertical_pt_em(mel, K = k, trace = FALSE, compute_se = FALSE)
  tempo <- as.numeric(difftime(Sys.time(), inicio, units = "secs"))
  
  cat("Modelo K =", k, "finalizado!\n")
  cat("Convergiu:", fit$converged, "\n")
  cat("Iterações do EM:", fit$iterations, "\n")
  cat("Tempo gasto:", round(tempo, 2), "segundos\n")
  cat("-----------------------------------\n")
  
  fit$tempo <- tempo
  
  return(fit)
})

names(modelos) <- paste0("K=", 1:10)

# Comparação dos modelos
comparacao <- do.call(rbind, lapply(1:10, function(k) {
  fit <- modelos[[k]]
  npar <- length(fit$gamma) + length(fit$lambda) + length(fit$beta)
  n <- nrow(fit$data)
  data.frame(K = k, logLik = fit$logLik,
             AIC = -2 * fit$logLik + 2 * npar,
             BIC = -2 * fit$logLik + log(n) * npar,
             convergiu = fit$converged)
}))

print(comparacao)

# Quantidade de eventos em cada intervalo
lapply(modelos, function(fit) fit$d_interval)
# Estimativas das hazards de ativação
lapply(modelos, function(fit) round(fit$lambda, 4))
# Pontos de corte, em anos
lapply(modelos, function(fit) round(fit$cuts, 2))

data.frame(K = 1:length(modelos),
           Cura_media = sapply(modelos, function(fit) mean(fit$cure)),
           Eventos_min = sapply(modelos, function(fit) min(fit$d_interval)),
           Eventos_max = sapply(modelos, function(fit) max(fit$d_interval)))

# Seleção do melhor modelo pelo BIC
melhor_k <- comparacao$K[which.min(ifelse(comparacao$convergiu,
                                          comparacao$BIC, Inf))]

fit_melanoma <- modelos[[melhor_k]]

cat("Melhor modelo pelo BIC: K =", melhor_k, "\n")
print(fit_melanoma)

# Profile used by Nicolaie et al. (2019): male, ulceration present and mean
# values of all continuous covariates (on their original scales).
profile_nicolaie <- data.frame(
  sexo_reg = 1,
  IDADE = mean(mel$IDADE),
  ano = mean(mel$ano),
  catend_reg = 1,
  ecgrup_2_3_reg = 0,
  ecgrup_4_reg =1,
  m_reg = 1,  
  RADIO =0,
  QUIMIO=0,
  CIRURGIA=1
)

# Use 500 replicates for standard errors and at least 1000 for final percentile
# intervals. The environment variable B_BOOT can be set for a quick pilot, e.g.
# Sys.setenv(B_BOOT = 100), before sourcing this file.
B_BOOT <- as.integer(Sys.getenv("B_BOOT", "500"))
time_grid <- seq(0, max(mel$time) / 365.25, length.out = 301)
boot_melanoma <- bootstrap_full(
  fit_melanoma, profile = profile_nicolaie, times = time_grid,
  B = B_BOOT, seed = 2026, trace = TRUE
)

# Main inferential table: nonparametric-bootstrap standard errors and percentile
# 95% confidence intervals. Unlike Wald intervals, lambda limits remain positive.
print(boot_melanoma$parameter_table, digits = 4, row.names = FALSE)
cat("Bootstrap successful fits:", boot_melanoma$successful, "of", B_BOOT, "\n")

# Profile-specific cure probability and its bootstrap percentile interval.
cure_band <- boot_melanoma$bands$cure_probability
cure_summary <- data.frame(
  profile = "Female, SUS, EC group 4, metastase, with only surgery, continuous covariates at their means",
  estimate = cure_band$estimate[1], lower.95 = cure_band$lower[1],
  upper.95 = cure_band$upper[1])
print(cure_summary, digits = 4, row.names = FALSE)

# Four-panel PDF analogous to Figures 1--4 in Nicolaie et al. (2019), adapted
# to the promotion-time susceptible stratum N>0.
plot_nicolaie_figures(
  fit_melanoma, boot_melanoma,
  file = "melanoma_vertical_promotion_time_figures.pdf"
)

# The same figure on the active graphics device.
plot_nicolaie_figures(fit_melanoma, boot_melanoma)
# 9 6

# Basic diagnostics required before reporting bootstrap inference.
diagnostics <- list(
  em_converged = fit_melanoma$converged,
  loglik_monotone = all(diff(fit_melanoma$logLik_history) >= -1e-7),
  events_per_interval = fit_melanoma$d_interval,
  events_by_cause = table(factor(mel$status[mel$status != 2], levels = c(1, 3))),
  bootstrap_successful = boot_melanoma$successful,
  bootstrap_failure_rate = boot_melanoma$failure_rate,
  bootstrap_reliable = boot_melanoma$failure_rate <= 0.20
)
print(diagnostics)

# Conditions for interpreting the bootstrap standard errors and intervals:
# 1. independent, non-informative right censoring;
# 2. adequate follow-up to identify the cure plateau;
# 3. full-rank count and cause design matrices;
# 4. no complete/quasi-complete separation in the cause model;
# 5. failures and positive exposure in every activation interval;
# 6. an interior solution: finite gamma/beta and strictly positive lambdas;
# 7. at least 80%--90% successful bootstrap refits;
# 8. convergence to the same maximum from several starting values;
# 9. B >= 500 for SEs and preferably B >= 1000 for percentile intervals.
