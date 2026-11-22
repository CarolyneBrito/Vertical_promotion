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
  J <- diag(length(par_hat))
  idx_lambda <- p + seq_len(K)
  J[cbind(idx_lambda, idx_lambda)] <- object$lambda
  V <- J %*% V_u %*% t(J)
  #J[p + seq_len(K), p + seq_len(K)] <- object$lambda
  #V <- J %*% V_u %*% J
  
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

# 3. Matriz de Informação Observada Analítica via Fórmula de Louis (1982)
compute_louis_se <- function(gamma, lambda, beta, X, Z, event, cause_other, time, E, ik, cuts) {
  p <- length(gamma)
  K <- length(lambda)
  q <- length(beta)
  n <- nrow(X)
  
  theta <- exp(as.vector(X %*% gamma))
  a     <- activation_terms(time, lambda, cuts, E)
  vhat  <- theta * a$S0  # Var[N_i | dados_i]
  
  # Informação Completa (I_comp)
  I_comp <- matrix(0, nrow = p + K, ncol = p + K)
  I_comp[1:p, 1:p] <- t(X) %*% (theta * X)
  for (k in seq_len(K)) {
    events_k <- sum(event == 1L & ik == k)
    I_comp[p + k, p + k] <- events_k / (lambda[k]^2)
  }
  
  # Informação Ausente (I_miss)
  I_miss <- matrix(0, nrow = p + K, ncol = p + K)
  for (i in seq_len(n)) {
    A_i <- c(X[i, ], -E[i, ])
    I_miss <- I_miss + vhat[i] * outer(A_i, A_i)
  }
  
  # Informação Observada
  I_obs_gl <- I_comp - I_miss
  vcov_gl  <- tryCatch(solve(I_obs_gl), error = function(e) matrix(NA, p + K, p + K))
  se_gl    <- sqrt(diag(vcov_gl))
  
  # Componente de causas (Regressão Logística)
  ind <- event == 1L
  Z_ev <- Z[ind, , drop = FALSE]
  eta_z <- as.vector(Z_ev %*% beta)
  p_cause <- 1 / (1 + exp(-eta_z))
  w <- p_cause * (1 - p_cause)
  I_obs_beta <- t(Z_ev) %*% (w * Z_ev)
  
  vcov_beta <- tryCatch(solve(I_obs_beta), error = function(e) matrix(NA, q, q))
  se_beta   <- sqrt(diag(vcov_beta))
  
  se_all <- c(se_gl[1:p], se_gl[p + seq_len(K)], se_beta)
  names(se_all) <- c(paste0("gamma_", 0:(p-1)), paste0("lambda_", 1:K), paste0("beta_", 0:(q-1)))
  
  list(se = se_all, vcov_gl = vcov_gl, vcov_beta = vcov_beta)
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

# ============================== DATA ANALYSIS ================================
setwd("~/Estudos/Doutorado/Materias/Sobrevivência/Projeto Vertical Promotion/Resultados_codigo")

mel <- dados_org %>%
  dplyr::select(time, status, sexo_reg, IDADE, ano, catend_reg,
                ecgrup_2_3_reg, ecgrup_4_reg, m_reg,
                RADIO, QUIMIO, CIRURGIA) %>%
  dplyr::mutate(status = dplyr::case_when(
    status == 0 ~ 2L, status == 1 ~ 1L, status == 2 ~ 3L,
    TRUE ~ NA_integer_))

# Numero de reamostragens.
B_BOOT <- 500L

tabela_regressao <- function(fit, boot) {
  tab <- fit$se$coefficients
  nm <- tab$parameter
  ep_louis <- unname(fit$louis$se)
  stopifnot(length(ep_louis) == nrow(tab))
  if (!is.null(boot$parameter_table)) {  # resultado de bootstrap_full()
    b <- boot$parameter_table[match(nm, boot$parameter_table$parameter), ]
    ep_b <- b$bootstrap.se
    li_b <- b$lower.95
    ls_b <- b$upper.95
  } else {                                # resultado de bootstrap_em()
    ep_b <- unname(boot$std.error[nm])
    li_b <- boot$percentile.95[nm, 1]
    ls_b <- boot$percentile.95[nm, 2]
  }
  
  p_valor <- 2 * pnorm(-abs(tab$estimate / tab$std.error))
  p_valor[grepl("^lambda", nm)] <- NA_real_
  
  data.frame(
    componente = ifelse(grepl("^gamma:", nm), "Cura",
                        ifelse(grepl("^lambda", nm), "Ativacao", "Causa")),
    parametro = nm,
    estimativa = tab$estimate,
    EP_Hessiana = tab$std.error,
    LI95_Hessiana = tab$lower.95,
    LS95_Hessiana = tab$upper.95,
    EP_Louis = ep_louis,
    LI95_Louis = tab$estimate - 1.96 * ep_louis,
    LS95_Louis = tab$estimate + 1.96 * ep_louis,
    EP_bootstrap = ep_b,
    LI95_bootstrap = li_b,
    LS95_bootstrap = ls_b,
    p_valor = p_valor,
    row.names = NULL
  )
}

# 1. Comparacao de K=1,...,10 
valores_K <- 1:10
modelos <- lapply(valores_K, function(k) {
  cat('\n--------------------------------------------\n')
  cat('Iniciando K =', k, 'as', format(Sys.time(), '%H:%M:%S'), '\n')
  inicio <- Sys.time()
  fit <- fit_vertical_pt_em(mel, K = k, trace = FALSE, compute_se = FALSE)
  fit$tempo <- as.numeric(difftime(Sys.time(), inicio, units = 'secs'))
  cat('Convergiu:', fit$converged,
      '| Iteracoes:', fit$iterations,
      '| Tempo (s):', round(fit$tempo, 2), '\n')
  fit
})
names(modelos) <- paste0('K=', valores_K)

comparacao <- do.call(rbind, lapply(valores_K, function(k) {
  fit <- modelos[[paste0('K=', k)]]
  p <- length(fit$gamma) + length(fit$lambda) + length(fit$beta)
  n <- nrow(fit$data)
  aic <- -2 * fit$logLik + 2 * p
  data.frame(K = k, logLik = fit$logLik, npar = p, AIC = aic,
             AICc = if (n > p + 1) aic + 2*p*(p+1)/(n-p-1) else Inf,
             BIC = -2*fit$logLik + log(n)*p,
             convergiu = fit$converged, tempo_seg = fit$tempo)
}))

print(comparacao, digits = 6, row.names = FALSE)
write.csv(comparacao, 'tabela_6_comparacao_K1_a_K10.csv', row.names = FALSE)
saveRDS(modelos, 'modelos_K1_a_K10.rds')
cat('Menor BIC:', comparacao$K[which.min(
  ifelse(comparacao$convergiu, comparacao$BIC, Inf))], '\n')

# 2. Modelo final com k=5
fit_k5 <- modelos[['K=5']]
stopifnot(fit_k5$converged)
print(fit_k5)

cat('\nInformacao observada do modelo completo K=5...\n')
fit_k5$se <- observed_information(fit_k5)
fit_k5$louis <- compute_louis_se(
  gamma = fit_k5$gamma, lambda = fit_k5$lambda, beta = fit_k5$beta,
  X = fit_k5$X, Z = fit_k5$Z,
  event = fit_k5$event, cause_other = fit_k5$cause_other,
  time = fit_k5$data$time,
  E = exposure_matrix(fit_k5$data$time, fit_k5$cuts),
  ik = interval_index(fit_k5$data$time, fit_k5$cuts),
  cuts = fit_k5$cuts
)

cat('\nBootstrap do modelo completo K=5: B =', B_BOOT, '\n')
boot_completo <- bootstrap_em(fit_k5, B = B_BOOT, seed = 2026, trace = TRUE)
saveRDS(boot_completo, 'bootstrap_modelo_completo_K5.rds')

regressao_k5 <- tabela_regressao(fit_k5, boot_completo)
print(regressao_k5, digits = 5, row.names = FALSE)
write.csv(regressao_k5, 'tabela_7_regressao_completa_K5.csv', row.names = FALSE)

cat('Bootstrap completo:', boot_completo$successful, '/', B_BOOT,
    '| falha:', boot_completo$failure_rate, '\n')

# 3. Investigacao do estadio IV e metastase 
frequencias_m <- with(fit_k5$data[fit_k5$event == 1, ],
                      table(ecgrup_4_reg, m_reg, status))
print(frequencias_m)
tabela_frequencias <- as.data.frame(frequencias_m,
                                    responseName = 'frequencia')
print(tabela_frequencias, row.names = FALSE)
write.csv(tabela_frequencias, 'tabela_8_frequencias_estadio_metastase.csv',
          row.names = FALSE)

# Instabilidade dos dois coeficientes entre reamostragens bootstrap.
coef_diagnostico <- c('beta:ecgrup_4_reg', 'beta:m_reg')
draws <- boot_completo$draws
cor_boot <- cor(draws[, coef_diagnostico[1]],
                draws[, coef_diagnostico[2]])
cat('Correlacao bootstrap estadio IV x metastase:', cor_boot, '\n')

quantis_boot <- t(apply(draws[, coef_diagnostico, drop = FALSE], 2,
                        quantile, probs = c(.025, .5, .975)))
tabela_instabilidade <- data.frame(
  parametro = coef_diagnostico,
  EP_Hessiana = regressao_k5$EP_Hessiana[match(
    coef_diagnostico, regressao_k5$parametro)],
  EP_Louis = regressao_k5$EP_Louis[match(
    coef_diagnostico, regressao_k5$parametro)],
  LI95_Louis = regressao_k5$LI95_Louis[match(
    coef_diagnostico, regressao_k5$parametro)],
  LS95_Louis = regressao_k5$LS95_Louis[match(
    coef_diagnostico, regressao_k5$parametro)],
  EP_bootstrap = unname(boot_completo$std.error[coef_diagnostico]),
  LI_bootstrap = quantis_boot[, 1],
  Mediana_bootstrap = quantis_boot[, 2],
  LS_bootstrap = quantis_boot[, 3],
  correlacao_entre_coeficientes = cor_boot,
  row.names = NULL)
print(tabela_instabilidade, digits = 5, row.names = FALSE)
write.csv(tabela_instabilidade, 'diagnostico_bootstrap_estadio_metastase.csv',
          row.names = FALSE)

# 4. Modelo alternativo K=5: retirar m_reg APENAS da componente da causa 
# A componente de cura continua com m_reg, sem alterar gamma nem a cura.
form_causa_sem_m <- update(fit_k5$cause_formula, ~ . - m_reg)
print(form_causa_sem_m)

fit_sem_m <- fit_vertical_pt_em(
  fit_k5$raw_data,
  K = 5,
  cuts = fit_k5$cuts,
  count_formula = fit_k5$count_formula,
  cause_formula = form_causa_sem_m,
  scaling = fit_k5$scale_info,
  compute_se = FALSE)
stopifnot(fit_sem_m$converged)

criterios <- function(fit) {
  p <- length(fit$gamma) + length(fit$lambda) + length(fit$beta)
  n <- nrow(fit$data)
  aic <- -2 * fit$logLik + 2 * p
  c(logLik = fit$logLik, parametros = p,
    AIC = aic, AICc = aic + 2*p*(p+1)/(n-p-1),
    BIC = -2 * fit$logLik + log(n) * p)
}

comparacao_m <- as.data.frame(rbind(
  Completo = criterios(fit_k5),
  Sem_m_na_causa = criterios(fit_sem_m)))
comparacao_m$Modelo <- rownames(comparacao_m)
comparacao_m <- comparacao_m[, c('Modelo', 'parametros', 'logLik',
                                 'AIC', 'AICc', 'BIC')]

LR <- 2 * (fit_k5$logLik - fit_sem_m$logLik)
p_LR <- pchisq(max(LR, 0), df = 1, lower.tail = FALSE)

comparacao_m$LR <- c(NA_real_, LR)
comparacao_m$p_LR <- c(NA_real_, p_LR)
print(comparacao_m, digits = 6, row.names = FALSE)
write.csv(comparacao_m, 'tabela_9_completo_versus_sem_m.csv', row.names = FALSE)

fit_sem_m$se <- observed_information(fit_sem_m)
fit_sem_m$louis <- compute_louis_se(
  gamma = fit_sem_m$gamma, lambda = fit_sem_m$lambda, beta = fit_sem_m$beta,
  X = fit_sem_m$X, Z = fit_sem_m$Z,
  event = fit_sem_m$event, cause_other = fit_sem_m$cause_other,
  time = fit_sem_m$data$time,
  E = exposure_matrix(fit_sem_m$data$time, fit_sem_m$cuts),
  ik = interval_index(fit_sem_m$data$time, fit_sem_m$cuts),
  cuts = fit_sem_m$cuts
)
diag_modelos <- data.frame(
  Modelo = c('Completo', 'Sem_m_na_causa'),
  Hessiana_positiva = c(fit_k5$se$positive_definite,
                        fit_sem_m$se$positive_definite),
  EP_confiaveis = c(fit_k5$se$reliable, fit_sem_m$se$reliable),
  Numero_condicao = c(fit_k5$se$condition_number,
                      fit_sem_m$se$condition_number),
  EP_beta_estadio_IV = c(
    fit_k5$se$coefficients$std.error[match('beta:ecgrup_4_reg',
                                           fit_k5$se$coefficients$parameter)],
    fit_sem_m$se$coefficients$std.error[match('beta:ecgrup_4_reg',
                                              fit_sem_m$se$coefficients$parameter)]),
  EP_Louis_beta_estadio_IV = c(
    regressao_k5$EP_Louis[match('beta:ecgrup_4_reg', regressao_k5$parametro)],
    unname(fit_sem_m$louis$se)[match('beta:ecgrup_4_reg',
                                     fit_sem_m$se$coefficients$parameter)]))
print(diag_modelos, digits = 6, row.names = FALSE)
write.csv(diag_modelos, 'diagnosticos_completo_versus_sem_m.csv',
          row.names = FALSE)

# 5. Modelo FINAL: K=5 com m_reg apenas na cura 
fit_final <- fit_sem_m
fit_melanoma <- fit_final  
mel_analise <- fit_final$raw_data

vars_perfil <- c('sexo_reg', 'IDADE', 'ano', 'catend_reg',
                 'ecgrup_2_3_reg', 'ecgrup_4_reg', 'm_reg',
                 'RADIO', 'QUIMIO', 'CIRURGIA')
escolher_perfil <- function(ec23, ec4, metastase) {
  idx <- which(mel_analise$ecgrup_2_3_reg == ec23 &
                 mel_analise$ecgrup_4_reg == ec4 &
                 mel_analise$m_reg == metastase)
  if (!length(idx)) stop('Nenhum registro observado para este estrato.')
  d <- mel_analise[idx, , drop = FALSE]
  escala_idade <- max(stats::sd(mel_analise$IDADE), 1)
  escala_ano <- max(stats::sd(mel_analise$ano), 1)
  distancia <- abs(d$IDADE - median(d$IDADE))/escala_idade +
    abs(d$ano - median(d$ano))/escala_ano
  ans <- d[which.min(distancia), vars_perfil, drop = FALSE]
  rownames(ans) <- NULL
  ans
}
perfis <- list('Estadio I, sem metastase' = escolher_perfil(0, 0, 0),
               'Estadio IV, com metastase' = escolher_perfil(0, 1, 1))
time_grid <- seq(0, max(mel_analise$time) / 365.25, length.out = 301)

# 6. Bootstrap do modelo FINAL 
#Primeiro perfil
cat('\nBootstrap final K=5 - perfil I | B =', B_BOOT, '\n')
boot_final_1 <- bootstrap_full(fit_final, profile = perfis[[1]], times = time_grid,
                               B = B_BOOT, seed = 2026, trace = TRUE)
saveRDS(boot_final_1, 'bootstrap_final_perfil_I_K5.rds')

regressao_final <- tabela_regressao(fit_final, boot_final_1)
print(regressao_final, digits = 5, row.names = FALSE)
write.csv(regressao_final, 'tabela_10_regressao_final_K5.csv', row.names = FALSE)

print(subset(regressao_final, parametro == 'beta:ecgrup_4_reg'),
      digits = 6, row.names = FALSE)

#Segundo perfil
cat('\nBootstrap final K=5 - perfil IV | B =', B_BOOT, '\n')
boot_final_2 <- bootstrap_full(
  fit_final, profile = perfis[[2]], times = time_grid,
  B = B_BOOT, seed = 2026, trace = TRUE)
saveRDS(boot_final_2, 'bootstrap_final_perfil_IV_K5.rds')

boots_perfis <- list('Estadio I, sem metastase' = boot_final_1,
                     'Estadio IV, com metastase' = boot_final_2)

tabela_bootstrap_validos <- data.frame(
  modelo = c("Completo K5","Final K5 / perfil I","Final K5 / perfil IV"),
  tentativas = B_BOOT,
  validas = c(boot_completo$successful, boot_final_1$successful, boot_final_2$successful),
  falhas = B_BOOT - c(boot_completo$successful, boot_final_1$successful, boot_final_2$successful))

print(tabela_bootstrap_validos, row.names = FALSE)
write.csv(tabela_bootstrap_validos, "diagnostico_bootstrap_validos.csv",
          row.names = FALSE)

# 7. Tabelas dos perfis e probabilidades de cura 
tabela_perfis <- do.call(rbind, lapply(names(perfis), function(nm) {
  data.frame(perfil = nm, perfis[[nm]], row.names = NULL)
}))
print(tabela_perfis, row.names = FALSE)
write.csv(tabela_perfis, 'tabela_11_caracteristicas_perfis.csv', row.names = FALSE)

cura_perfis <- do.call(rbind, lapply(names(boots_perfis), function(nm) {
  b <- boots_perfis[[nm]]$bands$cure_probability
  data.frame(perfil = nm, cura = b$estimate[1],
             LI95 = b$lower[1], LS95 = b$upper[1],
             bootstrap_sucesso = boots_perfis[[nm]]$successful)
}))
print(cura_perfis, digits = 5, row.names = FALSE)
write.csv(cura_perfis, 'tabela_12_cura_dois_perfis_K5_final.csv',
          row.names = FALSE)

# 8. Curvas marginais para TODOS os pacientes - MODELO FINAL 
cat('\nCalculando as curvas marginais para todos os pacientes...\n')
inicio_geral <- Sys.time()
S_ind <- predict(fit_final, times = time_grid, type = 'survival')
CIF_ind <- predict(fit_final, times = time_grid, type = 'cif')
p0_medio <- mean(fit_final$cure)
p_susc <- 1 - p0_medio
if (p_susc <= 0) stop('Nao ha massa suscetivel estimada.')

S_geral <- colMeans(S_ind)
F1_geral <- colMeans(CIF_ind$cause1_melanoma)
F3_geral <- colMeans(CIF_ind$cause3_other)

# CIFs do modelo completo com K=5
CIF_completo_K5 <- predict(fit_k5,times = time_grid,type = "cif")

# Comparacao com o modelo final reduzido
comparacao_CIF_m <- data.frame(
  curva = c("Sobrevivencia","CIF cancer osseo","CIF outras causas"),
  diferenca_maxima = c(max(abs(predict(fit_k5, time_grid, type = "survival") -  S_ind)),
                       max(abs(colMeans(CIF_completo_K5$cause1_melanoma) -F1_geral)),
                       max(abs(colMeans(CIF_completo_K5$cause3_other) - F3_geral))))

print(comparacao_CIF_m, digits = 8, row.names = FALSE)
write.csv(
  comparacao_CIF_m,
  "comparacao_CIF_completo_sem_m_K5.csv",
  row.names = FALSE)


pred_geral <- data.frame(
  time = time_grid,
  survival_population = S_geral,
  survival_susceptible = pmax(0, (S_geral - p0_medio) / p_susc),
  cif_melanoma_population = F1_geral,
  cif_other_population = F3_geral,
  cif_melanoma_susceptible = F1_geral / p_susc,
  cif_other_susceptible = F3_geral / p_susc,
  cure_probability = rep(p0_medio, length(time_grid)))
cat('Tempo CIF geral (s):',
    round(as.numeric(difftime(Sys.time(), inicio_geral, units = 'secs')), 1), '\n')
erro_CIF <- max(abs(with(pred_geral, survival_population +
                           cif_melanoma_population + cif_other_population - 1)))
cat('Erro maximo S+F1+F3-1:', erro_CIF, '\n')
write.csv(pred_geral, 'curvas_marginais_K5_final.csv', row.names = FALSE)


# 9. Graficos finais
plot_nicolaie_adaptado <- function(estimativa, titulo,
                                   bandas = NULL, objeto = NULL,
                                   geral = FALSE, arquivo = NULL,
                                   grupo_estadio = NULL) {
  if (is.null(objeto)) stop('Informe objeto = fit_final para calcular o painel A.')
  if (!geral && (is.null(grupo_estadio) || !grupo_estadio %in% c('I', 'IV'))) {
    stop("Nos graficos por perfil, informe grupo_estadio = 'I' ou 'IV'.")
  }
  if (!is.null(arquivo)) grDevices::pdf(arquivo, width = 11, height = 8.5)
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(oldpar)
    if (!is.null(arquivo)) grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.4, 3.3, 1),
                oma = c(0, 0, 2, 0), las = 1)
  cores <- c(cancer = '#0072B2', outras = '#D55E00',
             pop = '#009E73', susc = '#CC79A7')
  t <- estimativa$time
  
  # (a) KM e sobrevivencia prevista para exatamente o mesmo conjunto de dados
  if (geral) {
    idx <- rep(TRUE, nrow(objeto$data))
    S_ajustada <- estimativa$survival_population
  } else {
    idx <- if (grupo_estadio == 'I') {
      objeto$data$ecgrup_2_3_reg == 0 & objeto$data$ecgrup_4_reg == 0
    } else {
      objeto$data$ecgrup_4_reg == 1
    }
    if (!any(idx)) stop('Nao ha pacientes no grupo de estadio solicitado.')
    # predict(..., newdata=NULL) preserva a ordem das linhas de objeto$data
    S_todos <- predict(objeto, times = t, type = 'survival')
    S_ajustada <- colMeans(S_todos[idx, , drop = FALSE])
  }
  km <- survival::survfit(
    survival::Surv(objeto$data$time[idx], objeto$event[idx]) ~ 1
  )
  graphics::plot(km, conf.int = FALSE, mark.time = FALSE,
                 col = 'gray30', lwd = 2, xlim = range(t), ylim = c(0, 1),
                 xlab = 'Tempo (anos)', ylab = 'Sobrevivencia',
                 main = '(a) Kaplan-Meier e modelo ajustado')
  graphics::lines(t, S_ajustada, col = cores['pop'], lwd = 2)
  graphics::legend('topright', c('Kaplan-Meier', 'Media prevista (K=5)'),
                   col = c('gray30', cores['pop']), lty = 1, lwd = 2, bty = 'n')
  
  # Funcao auxiliar: desenha banda percentil quando fornecida
  desenhar_banda <- function(campo, cor) {
    if (!is.null(bandas) && !is.null(bandas[[campo]]))
      add_band(bandas[[campo]], cor)
  }
  graphics::plot(t, estimativa$cif_melanoma_population, type = 'n',
                 ylim = c(0,1), xlab = 'Tempo (anos)', ylab = 'Incidencia acumulada',
                 main = '(b) CIF populacional')
  desenhar_banda('cif_melanoma_population', cores['cancer'])
  desenhar_banda('cif_other_population', cores['outras'])
  graphics::lines(t, estimativa$cif_melanoma_population, col = cores['cancer'], lwd = 2)
  graphics::lines(t, estimativa$cif_other_population, col = cores['outras'], lwd = 2)
  graphics::legend('topleft', c('Cancer osseo','Outras causas'),
                   col = cores[c('cancer','outras')], lty = 1, lwd = 2, bty = 'n')
  
  graphics::plot(t, estimativa$survival_population, type = 'n',
                 ylim = c(0,1), xlab = 'Tempo (anos)', ylab = 'Sobrevivencia',
                 main = '(c) Populacao e suscetiveis')
  desenhar_banda('survival_population', cores['pop'])
  desenhar_banda('survival_susceptible', cores['susc'])
  graphics::lines(t, estimativa$survival_population, col = cores['pop'], lwd = 2)
  graphics::lines(t, estimativa$survival_susceptible, col = cores['susc'], lwd = 2, lty = 2)
  graphics::abline(h = estimativa$cure_probability[1], col = 'gray40', lty = 3)
  graphics::legend('topright', c('Populacao','Suscetiveis','Fracao de cura'),
                   col = c(cores['pop'],cores['susc'],'gray40'),
                   lty = c(1,2,3), lwd = c(2,2,1), bty = 'n')
  
  graphics::plot(t, estimativa$cif_melanoma_susceptible, type = 'n',
                 ylim = c(0,1), xlab = 'Tempo (anos)', ylab = 'Incidencia acumulada',
                 main = '(d) CIF entre suscetiveis')
  desenhar_banda('cif_melanoma_susceptible', cores['cancer'])
  desenhar_banda('cif_other_susceptible', cores['outras'])
  graphics::lines(t, estimativa$cif_melanoma_susceptible,
                  col = cores['cancer'], lwd = 2)
  graphics::lines(t, estimativa$cif_other_susceptible,
                  col = cores['outras'], lwd = 2)
  graphics::legend('topleft', c('Cancer osseo','Outras causas'),
                   col = cores[c('cancer','outras')], lty = 1, lwd = 2, bty = 'n')
  graphics::mtext(titulo, outer = TRUE, line = 0.4, cex = 1.1)
  invisible(estimativa)
}

# Figuras finais K=5: dois perfis com bandas e grupo geral pontual.
plot_nicolaie_adaptado(boot_final_1$estimate, 'Estadio I',
                       bandas = boot_final_1$bands, objeto = fit_final,
                       grupo_estadio = 'I')
plot_nicolaie_adaptado(boot_final_2$estimate, 'Estadio IV',
                       bandas = boot_final_2$bands, objeto = fit_final,
                       grupo_estadio = 'IV')
plot_nicolaie_adaptado(pred_geral, 'Geral', objeto = fit_final, geral = TRUE)

# Salvar 
plot_nicolaie_adaptado(
  estimativa = pred_geral,
  titulo = "K=5 - Todos os pacientes",
  objeto = fit_final,
  geral = TRUE,
  arquivo = "figuras_gerais_K5_final.pdf")

plot_nicolaie_adaptado(
  estimativa = boot_final_1$estimate,
  titulo = "K=5 - Estágio I, sem metástase",
  bandas = boot_final_1$bands,
  objeto = fit_final,
  geral = FALSE,
  grupo_estadio = "I",
  arquivo = "figuras_estadio_I_K5_final.pdf")

plot_nicolaie_adaptado(
  estimativa = boot_final_2$estimate,
  titulo = "K=5 - Estágio IV, com metástase",
  bandas = boot_final_2$bands,
  objeto = fit_final,
  geral = FALSE,
  grupo_estadio = "IV",
  arquivo = "figuras_estadio_IV_K5_final.pdf")


# 10. Diagnosticos finais 
print(list(
  convergencia_K5_completo = fit_k5$converged,
  convergencia_K5_final = fit_final$converged,
  logLik_monotona_K5 = all(diff(fit_final$logLik_history) >= -1e-7),
  eventos_por_intervalo_K5 = fit_final$d_interval,
  hessiana_positiva_final = fit_final$se$positive_definite,
  EP_confiaveis_final = fit_final$se$reliable,
  numero_condicao_final = fit_final$se$condition_number,
  simetria_covariancia_final = max(abs(fit_final$se$vcov -
                                         t(fit_final$se$vcov))),
  erro_maximo_CIF = erro_CIF,
  bootstrap_completo_sucesso = boot_completo$successful,
  bootstrap_perfis_sucesso = sapply(boots_perfis, function(b) b$successful),
  bootstrap_perfis_falha = sapply(boots_perfis, function(b) b$failure_rate)
))

saveRDS(fit_final, 'modelo_final_K5_sem_m_na_causa.rds')