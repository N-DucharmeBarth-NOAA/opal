
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
