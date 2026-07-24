#' Natural mortality at age
#'
#' Constructs a vector of M-at-age values using a declining early-age curve and late-age increase.
#'
#' @param min_age minimum model age.
#' @param max_age maximum model age.
#' @param age_increase_M age at which M begins to increase again.
#' @param m0 M at age-1 and 2.
#' @param m4 M at age-4 (controls slope of decline).
#' @param m10 M at age-10 (base for flat zone).
#' @param m30 M at age-30 (terminal age M).
#' @return vector of M at age.
#' @export
#'
get_M <- function(min_age, max_age, age_increase_M, m0, m4, m10, m30) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  n_age <- max_age - min_age + 1
  n_increase_M <- age_increase_M - min_age + 1
  M_slope <- log((m0 - m4) / (m0 - m10)) / log(0.33333)
  M_inc <- (m30 - m10) / (max_age - age_increase_M)
  M_age <- numeric(n_age)
  M_age[1:2] <- m0
  for (i in 3:10) {
    M_age[i] <- m0 - (m0 - m10) * ((i - 2) / 9)^M_slope
  }
  if (n_increase_M >= 10) {
    M_age[11:n_increase_M] <- m10
  }
  if (n_age > n_increase_M) {
    for (i in (n_increase_M + 1):(n_age - 1)) {
      M_age[i] <- M_age[i - 1] + M_inc
    }
  }
  M_age[n_age] <- m30
  return(M_age)
}

#' Natural mortality at length
#'
#' Constructs a vector of M-at-age values using a declining early-age curve and late-age increase.
#'
#' @param min_age minimum model age.
#' @param max_age maximum model age.
#' @param age_increase_M age at which M begins to increase again.
#' @param m0 M at age-1 and 2.
#' @param m30 M at age-30 (terminal age M).
#' @param length_mu_ysa M at age-30 (terminal age M).
#' @return vector of M at age.
#' @export
#'
get_M_length <- function(min_age, max_age, age_increase_M, m0, m30, length_mu_ysa) {
  "[<-" <- ADoverload("[<-")
  n_age <- max_age - min_age + 1
  M_age <- numeric(n_age)
  L0 <- length_mu_ysa[1, 1, 1]
  Mc <- -1.5
  Ma <- log(m0) - Mc * log(L0)
  for (a in seq_len(n_age)) {
    M_age[a] <- exp(Ma + Mc * log(length_mu_ysa[1, 1, a]))
  }
  return(M_age)
}

#' Plot natural mortality
#' 
#' Plot natural mortality (M) by age.
#' 
#' @param data A model data list passed to \code{MakeADFun}.
#' @param object The AD object created using \code{MakeADFun}.
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @export
#' 
plot_natural_mortality <- function(data, object) {
  mortality <- data.frame(
    age = data$min_age:data$max_age,
    value = object$report()$M_a
  )
  
  ggplot(mortality, aes(x = .data$age, y = .data$value)) +
    geom_line(linetype = "dashed") +
    labs(x = "Age", y = "Natural mortality") +
    scale_x_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    )
}
