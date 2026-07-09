# Tests for the length-engine equilibrium initialisation
library(RTMB)

# opakapaka-like fixtures
n_age <- 30L
A1 <- 1L; A2 <- 30L
L1 <- 20.0; L2 <- 90.0
log_k <- log(0.15)
len_lower <- seq(5, by = 5, length.out = 20)
len_upper <- len_lower + 5
len_mid   <- (len_lower + len_upper) / 2
n_len <- length(len_lower)
B0 <- 5000; h <- 0.8
M_a <- rep(0.2, n_age)

# maturity x fecundity at length
mat_l <- 1 / (1 + exp(-0.3 * (len_mid - 50)))
fec_l <- len_mid^3 / 1e5
sp_l  <- mat_l * fec_l

make_bits <- function(log_CV1, log_CV2, pg = TRUE) {
  mu_a <- get_growth(n_age, A1, A2, L1, L2, log_k)
  sd_a <- get_sd_at_age(mu_a, L1, L2, log_CV1, log_CV2)
  pla  <- get_pla(len_lower, len_upper, mu_a, sd_a)
  G    <- get_growth_matrix(len_lower, len_upper, len_mid, sd_a, L1, L2, log_k, A1, A2,
                            plus_group_growth = pg)
  rec  <- get_recruit_length_dist(len_lower, len_upper, mu_a[1], sd_a[1])
  list(mu_a = mu_a, sd_a = sd_a, pla = pla, G = G, rec = rec)
}

test_that("unfished equilibrium SSB0 equals B0 (M at age)", {
  b <- make_bits(log(0.1), log(0.1))
  init <- get_initial_numbers_length(B0, h, M_a, sp_l, b$pla, b$G, b$rec)
  SSB0 <- sum(colSums(init$Ninit0_al) * sp_l)
  expect_equal(SSB0, B0, tolerance = 1e-6)
})

test_that("no-fishing fished equilibrium equals unfished (M at age)", {
  b <- make_bits(log(0.1), log(0.1))
  init <- get_initial_numbers_length(B0, h, M_a, sp_l, b$pla, b$G, b$rec,
                                     init_F_f = 0, sel_fl = rep(0, n_len))
  expect_equal(init$Ninit_al, init$Ninit0_al, tolerance = 1e-8)
})

test_that("length equilibrium (M at age) reduces exactly to age-only", {
  b <- make_bits(log(0.12), log(0.08))     # realistic, non-degenerate SD
  sp_a <- as.vector(t(b$pla) %*% sp_l)
  age  <- get_initial_numbers(B0, h, M_a, sp_a)
  len  <- get_initial_numbers_length(B0, h, M_a, sp_l, b$pla, b$G, b$rec)
  expect_equal(len$R0,    age$R0,    tolerance = 1e-10)
  expect_equal(len$alpha, age$alpha, tolerance = 1e-10)
  expect_equal(len$beta,  age$beta,  tolerance = 1e-10)
  # age marginal of the length-structured initial numbers matches age-only exactly
  expect_equal(rowSums(len$Ninit0_al), age$Ninit0, tolerance = 1e-10)
})

test_that("initial recruitment deviations scale the age rows", {
  b <- make_bits(log(0.1), log(0.1))
  dev <- rep(0, n_age); dev[3] <- 0.5
  base <- get_initial_numbers_length(B0, h, M_a, sp_l, b$pla, b$G, b$rec)
  pert <- get_initial_numbers_length(B0, h, M_a, sp_l, b$pla, b$G, b$rec,
                                     init_rdev_a = dev, sigma_r = 0.6)
  ratio <- rowSums(pert$Ninit0_al) / rowSums(base$Ninit0_al)
  expect_equal(ratio[3], exp(0.5), tolerance = 1e-8)
  expect_equal(ratio[-3], rep(1, n_age - 1), tolerance = 1e-8)
})

test_that("M-at-length path builds a valid equilibrium (SSB0 = B0)", {
  b <- make_bits(log(0.1), log(0.1))
  M_l <- rep(0.2, n_len)
  init <- get_initial_numbers_length(B0, h, M_l, sp_l, b$pla, b$G, b$rec,
                                     plus_group_growth = TRUE, M_basis = "length")
  SSB0 <- sum(colSums(init$Ninit0_al) * sp_l)
  expect_equal(SSB0, B0, tolerance = 1e-6)
  expect_true(all(init$Ninit0_al >= -1e-10))
})

test_that("get_initial_numbers_length differentiates through B0 (both M bases)", {
  b_num <- make_bits(log(0.1), log(0.1))
  for (basis in c("age", "length")) {
    Mvec <- if (basis == "age") M_a else rep(0.2, n_len)
    f <- function(p) {
      mu_a <- get_growth(n_age, A1, A2, L1, L2, log_k)
      sd_a <- get_sd_at_age(mu_a, L1, L2, log(0.1), log(0.1))
      pla  <- get_pla(len_lower, len_upper, mu_a, sd_a)
      G    <- get_growth_matrix(len_lower, len_upper, len_mid, sd_a, L1, L2, log_k, A1, A2)
      rec  <- get_recruit_length_dist(len_lower, len_upper, mu_a[1], sd_a[1])
      init <- get_initial_numbers_length(p[1], h, Mvec, sp_l, pla, G, rec, M_basis = basis)
      sum(init$Ninit0_al)
    }
    tap <- RTMB::MakeTape(f, c(B0))
    g <- tap$jacobian(c(B0))
    expect_true(all(is.finite(g)), info = basis)
  }
})
