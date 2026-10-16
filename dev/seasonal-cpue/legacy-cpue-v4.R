# Frozen CPUE implementation for the v4/v5 sensitivity analysis only.
# Copied verbatim, apart from its name, from R/likelihoods.R at
# 2dd8e3de2bbfb4f00a9489c1ae456d8dd2ce5e14 (main / opal 0.0.4).
# Source: https://github.com/N-DucharmeBarth-NOAA/opal/blob/2dd8e3de2bbfb4f00a9489c1ae456d8dd2ce5e14/R/likelihoods.R
# GPL-3 or later, as in the source package. Never used by the package model.
legacy_cpue_like_v4 <- function(cpue_data, parameters, number_ysa, sel_fya, weight_fya, cpue_switch = 1L) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  data_n_index <- NULL
  if (is.list(cpue_data) && !is.data.frame(cpue_data) && !is.null(cpue_data$cpue_data)) {
    if (!is.null(cpue_data$cpue_switch)) cpue_switch <- cpue_data$cpue_switch
    data_n_index <- cpue_data$n_index
    cpue_data <- cpue_data$cpue_data
  }
  log_cpue_q <- parameters$log_cpue_q
  log_cpue_tau <- parameters$log_cpue_tau
  log_cpue_omega <- parameters$log_cpue_omega
  cpue_creep <- parameters$cpue_creep
  n_cpue <- nrow(cpue_data)
  if (!("index" %in% names(cpue_data))) cpue_data$index <- rep(1L, n_cpue)
  n_index <- max(
    cpue_data$index,
    length(log_cpue_q),
    length(log_cpue_tau),
    length(log_cpue_omega),
    if (is.null(cpue_creep)) 1L else length(cpue_creep),
    if (is.null(data_n_index)) 1L else data_n_index
  )
  if (is.null(cpue_creep)) cpue_creep <- rep(0, n_index)
  if (length(log_cpue_q) == 1L && n_index > 1L) log_cpue_q <- rep(log_cpue_q, n_index)
  if (length(log_cpue_tau) == 1L && n_index > 1L) log_cpue_tau <- rep(log_cpue_tau, n_index)
  if (length(log_cpue_omega) == 1L && n_index > 1L) log_cpue_omega <- rep(log_cpue_omega, n_index)
  if (length(cpue_creep) == 1L && n_index > 1L) cpue_creep <- rep(cpue_creep, n_index)
  cpue_log_pred <- cpue_adjust <- cpue_sigma <- lp <- numeric(n_cpue)
  for (idx in seq_len(n_index)) {
    rows <- which(cpue_data$index == idx)
    if (length(rows) == 0) next
    cpue_adjust[rows[1]] <- 1
    if (length(rows) > 1) {
      for (j in 2:length(rows)) {
        cpue_adjust[rows[j]] <- cpue_adjust[rows[j - 1]] + cpue_creep[idx]
      }
    }
  }
  for (i in seq_len(n_cpue)) {
    y <- cpue_data$ts[i]
    f <- cpue_data$fishery[i]
    idx <- cpue_data$index[i]
    cpue_n <- number_ysa[y, 1, ] * sel_fya[f, y, ]
    if (cpue_data$units[i] == 1) cpue_n <- cpue_n * weight_fya[f, y,] # 1=weight, 2=numbers
    sum_n <- sum(cpue_n) + 1e-6
    cpue_log_pred[i] <- log(cpue_adjust[i]) + exp(log_cpue_omega[idx]) * log(sum_n)
  }
  for (idx in seq_len(n_index)) {
    rows <- which(cpue_data$index == idx)
    if (length(rows) == 0) next
    centre <- log(mean(exp(cpue_log_pred[rows])))
    cpue_log_pred[rows] <- cpue_log_pred[rows] - centre + log_cpue_q[idx]
    tau_idx <- exp(log_cpue_tau[idx])
    cpue_sigma[rows] <- sqrt(cpue_data$se[rows]^2 + tau_idx^2)
  }
  cpue_log_obs <- log(cpue_data$value)
  cpue_log_obs <- OBS(cpue_log_obs)
  if (cpue_switch > 0) {
    lp[] <- -dnorm(x = cpue_log_obs, mean = cpue_log_pred, sd = cpue_sigma, log = TRUE)
  }
  cpue_pred <- exp(cpue_log_pred)
  REPORT(cpue_pred)
  REPORT(cpue_sigma)
  return(lp = lp)
}
