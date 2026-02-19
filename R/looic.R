#' Obtain the leave-one-out information criterion (LOO IC)
#' 
#' This function extracts point-wise log-likelihood values and then calculates LOO statistics.
#' 
#' @param data a \code{list} containing the data that was passed to \code{MakeADFun}.
#' @param object A \code{list} specifying the AD object created using \code{MakeADFun}.
#' @param posterior An \code{snutsfit} object created using the \code{sample_snuts} function.
#' @return a \code{psis_loo} object.
#' @import dplyr
#' @importFrom loo relative_eff loo
#' @export
#' 
get_loo <- function(data, object, posterior) {
  pars <- c("lp_af", "lp_lf", "lp_cpue", "lp_aerial", "lp_troll", "lp_tags", "lp_pop", "lp_hsp", "lp_gt")
  pars <- pars[c(TRUE, TRUE, data$cpue_switch > 0, data$aerial_switch > 0, data$troll_switch > 0, data$tag_switch > 0, data$pop_switch > 0, data$hsp_switch > 0, data$gt_switch > 0)]
  nll_list <- get_posterior(object = object, posterior = posterior, pars = pars, option = 2, type = "list")
  n_chains <- length(nll_list)
  for (j in 1:n_chains) {
    nll_list[[j]]$lp_af <- nll_list[[j]]$lp_af[,data$af_n > 0]
    nll_list[[j]]$lp_lf <- nll_list[[j]]$lp_lf[,data$lf_n > 0]
  }
  df <- melt(nll_list) %>%
    select(parameter = L2, index, chain = L1, iteration, value) %>%
    tibble() %>%
    mutate(names = paste(parameter, index, sep = "_")) %>% 
    select(chain, iteration, names, value) %>%
    pivot_wider(names_from = names, values_from = value)
  x <- df %>% select(-chain, -iteration) %>% as.matrix()
  r_eff <- relative_eff(x = exp(-x), chain_id = df$chain) # Convert from negative log-likelihood to likelihood
  loo_1 <- loo(x = -x, r_eff = r_eff) # Convert from negative log-likelihood to log-likelihood
  return(loo_1)
}


#' Plot PSIS LOO
#' 
#' This is the same as \code{plot(loo1)} except it colours/groups the different 
#' types of data.
#' 
#' @param x A \code{psis_loo} object.
#' @param exclude Any variables to exclude from the plot (e.g., aerial, af, lf, cpue).
#' @return a \code{ggplot2}.
#' @importFrom rlang .data
#' @importFrom tidyr separate
#' @import dplyr
#' @import ggplot2
#' @export
#' 
plot_loo <- function(x, exclude = NULL) {
  df <- x$pointwise %>%
    data.frame() %>%
    # filter(!is.infinite(influence_pareto_k)) %>%
    mutate(name = row.names(.)) %>%
    mutate(name = gsub("lp_", "", .data$name)) %>%
    separate(name, into = c("Variable", "point")) %>%
    mutate(Variable = factor(.data$Variable))
  if (!is.null(exclude)) df <- df %>% filter(!.data$Variable %in% exclude)
  df <- df %>% mutate(x = 1:n())
  p <- ggplot(data = df, aes(x = .data$x, y = .data$influence_pareto_k, colour = .data$Variable)) +
    geom_point(aes(shape = .data$Variable)) +
    scale_shape_manual(values = 1:nlevels(df$Variable)) +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    labs(x = "Data point", y = "Pareto shape k", title = "PSIS diagnostic plot")
  return(p)
}
