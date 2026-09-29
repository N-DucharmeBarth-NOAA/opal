test_that("posterior summaries reproduce draw-wise reports and preserve source identity", {
  x <- opal_fit(small_opal_object())
  x <- opal_attach_mcmc(x, mock_opal_sampler(opal_rtmb(x), num_samples=120), check=FALSE)
  x <- opal_posterior(x, quantities=c("B0","spawning_biomass_y"), draws=c(1,3,121,123))
  result <- opal_derived(x,"posterior")
  expect_equal(result$draws$chain,c(1L,1L,2L,2L))
  selected <- c(x$mcmc$samples[c(11,13),1,1],x$mcmc$samples[c(11,13),2,1])
  expect_equal(result$reports$mean[result$reports$quantity=="B0"],mean(exp(selected)))
  expect_equal(result$parameters$mean,mean(selected))
  expect_equal(result$dimensions$spawning_biomass_y,3L)
  path <- tempfile(fileext=".rds"); on.exit(unlink(path))
  opal_save(x,path)
  expect_identical(opal_derived(opal_read(path),"posterior"),result)
  expect_error(opal_posterior(x,draws=c(1,1)),"distinct")
  expect_error(opal_posterior(x,quantities="missing"),"numeric model reports")
})

test_that("profiles re-optimise a fixed parameter on its stored scale", {
  x <- opal_fit(small_opal_object())
  values <- x$fit$opt$par + c(-.1,0,.1)
  x <- opal_profile(x,"log_B0",as.numeric(values))
  result <- opal_derived(x,"profile")
  expect_true(all(is.na(result$table$error)))
  expect_equal(result$table$delta[2],0,tolerance=1e-7)
  expect_true(all(result$table$delta >= -1e-7))
  expect_equal(result$table$objective, vapply(result$fits,function(z) opal_rtmb(z)$fn(numeric()),numeric(1)))
  expect_equal(as.numeric(rowsum(result$components$objective,result$components$value)),result$table$objective)
  expect_s3_class(plot_opal_profile(x),"ggplot")
})

test_that("grid checkpoints match configurations and balanced draws retain model weights", {
  x <- small_opal_object()
  directory <- tempfile(); dir.create(directory); on.exit(unlink(directory,recursive=TRUE))
  grid <- opal_grid(x,list(base=list(), high=list(priors=list(log_B0=list(type="normal",par1=15,par2=1,index=1)))),directory)
  expect_true(all(grid$summary$passes))
  expect_false(any(grid$summary$reused))
  again <- opal_grid(x,grid$settings$scenarios,directory)
  expect_true(all(again$summary$reused))
  grid <- opal_grid_mcmc(grid,sampler=mock_opal_sampler,seed=32,chains=2)
  expect_true(all(grid$summary$mcmc_passes))
  ids <- opal_grid_draws(grid,10)
  expect_equal(as.integer(table(ids$model)),c(10L,10L))
  expect_equal(sum(ids$weight),1)
  expect_true(all(table(ids$model,ids$chain)==5))
})

test_that("equilibrium MSY agrees with unfished quantities and a dense harvest grid", {
  x <- opal_fit(small_opal_object())
  r <- opal_report(x)
  eq <- function(u) opal:::.opal_equilibrium(u,c(.4,.6),
    matrix(r$sel_fya[,2,],2,5),matrix(r$weight_fya_mod[,2,],2,5),
    r$M_a,r$spawning_potential_a,r$alpha,r$beta,1)
  expect_equal(eq(0)$spawning,r$B0)
  expect_equal(eq(0)$recruitment,r$R0)
  expect_equal(eq(0)$yield,0)
  x <- opal_msy(x,c(.4,.6))
  result <- opal_derived(x,"msy")$draws
  expect_true(result$resolved)
  brute <- max(vapply(seq(0,1,length.out=10001),function(u)eq(u)$yield,numeric(1)))
  expect_equal(result$msy,brute,tolerance=1e-6)
  expect_equal(eq(1)$recruitment,0)
  x <- opal_attach_mcmc(x,mock_opal_sampler(opal_rtmb(x),num_samples=10),check=FALSE)
  x <- opal_msy(x,c(.4,.6),uncertainty="mcmc",draws=c(1,11))
  expect_equal(nrow(opal_derived(x,"msy")$draws),2L)
})

test_that("simulation examples cover both seasons and all composition distributions", {
  for(family in c("multinomial","Dirichlet","Dirichlet-multinomial")) {
    set.seed(14); rng <- .Random.seed
    inputs <- opal_example_inputs(family)
    expect_identical(.Random.seed,rng)
    expect_equal(range(inputs$data$cpue_data$ts),c(1L,24L))
    x <- opal_fit(opal_obj(inputs$data,inputs$parameters,inputs$map))
    expect_true(x$validation$fit$passes)
    x <- opal_osa(x)
    expect_equal(nrow(opal_derived(x,"osa")$summary),6L)
    expect_true(all(opal_derived(x,"osa")$summary$failures==0))
  }
})

test_that("posterior prior densities retain parameter scales and reject ambiguous mappings", {
  inputs <- opal_example_inputs()
  x <- opal_fit(opal_obj(inputs$data,inputs$parameters,inputs$map))
  x <- opal_attach_mcmc(x,mock_opal_sampler(opal_rtmb(x),num_samples=30),check=FALSE)
  p <- plot_prior_posterior(x)
  curve <- p$layers[[2]]$data
  expect_equal(curve$density,dnorm(curve$value,log(1e5),.3))
  expect_equal(unique(curve$parameter),"log_B0")
  expect_error(plot_prior_posterior(x,parameters="missing"),"No unambiguous")
})

test_that("MSY seasonal equilibrium matches independent repeated population transitions", {
  sel <- rbind(c(.1,.5,1),c(.2,.6,.8)); weight <- rbind(c(1,2,3),c(1,2,4))
  M <- c(.2,.3,.4); spawning <- c(0,1,2)
  initial <- get_unfished_init(1000,.75,M,spawning)
  eq <- opal:::.opal_equilibrium(.2,c(.4,.6),sel,weight,M,spawning,initial$alpha,initial$beta,2)
  N <- initial$rel_N * initial$R0
  harvest <- sweep(sel,1,c(.08,.12),"*")
  for(year in 1:1000) {
    catch <- 0
    for(season in 1:2) {
      catch <- catch + sum(sweep(harvest*weight,2,N,"*"))
      N <- N*(1-colSums(harvest))*exp(-M/2)
    }
    N <- c(0,N[1],N[2]+N[3])
    B <- sum(N*spawning)
    N[1] <- initial$alpha*B/(initial$beta+B)
  }
  expect_equal(N,eq$numbers,tolerance=1e-9)
  expect_equal(B,eq$spawning,tolerance=1e-9)
  expect_equal(catch,eq$yield,tolerance=1e-9)
})

test_that("shared-map profiles fix the whole group and retain remaining bounds", {
  inputs <- opal_example_inputs()
  map <- inputs$map
  map$log_cpue_q <- factor(c(1,1))
  x <- opal_fit(opal_obj(inputs$data,inputs$parameters,map))
  value <- x$fit$parameters$log_cpue_q[1]
  x <- opal_profile(x,"log_cpue_q",value,element=2)
  result <- opal_derived(x,"profile")
  expect_true(result$table$passes)
  expect_equal(result$mapped_elements,1:2)
  expect_equal(as.numeric(result$fits[[1]]$fit$parameters$log_cpue_q),rep(value,2))
  expect_identical(names(result$fits[[1]]$bounds$lower),"log_B0")
})

test_that("multiple priors on one block are not misrepresented as a single density", {
  inputs <- opal_example_inputs()
  inputs$data$priors$extra <- inputs$data$priors$log_B0
  x <- opal_fit(opal_obj(inputs$data,inputs$parameters,inputs$map))
  x <- opal_attach_mcmc(x,mock_opal_sampler(opal_rtmb(x),num_samples=10),check=FALSE)
  expect_warning(expect_error(plot_prior_posterior(x),"No unambiguous"),"Ambiguous prior blocks")
})

test_that("growth uncertainty reaches posterior biology and equilibrium reference points", {
  inputs <- opal_example_inputs()
  inputs$data$weight <- c(1,2,4)
  inputs$data$M <- c(.6,.3,.2)
  inputs$map$log_L1 <- NULL
  x <- opal_build(opal_obj(inputs$data,inputs$parameters,inputs$map))
  o <- opal_rtmb(x)
  draws <- rbind(o$par,o$par)
  colnames(draws) <- opal:::.opal_expand_parameter_names(names(o$par))
  draws[2,"log_L1"] <- draws[2,"log_L1"] + .15
  x <- opal_attach_mcmc(x,draws,check=FALSE)
  x <- opal_posterior(x,quantities=c("M_a","weight_fya_mod"))
  reports <- opal_derived(x,"posterior")$reports
  expect_true(any(reports$sd[reports$quantity=="M_a"] > .001))
  expect_true(any(reports$sd[reports$quantity=="weight_fya_mod"] > .001))
  x <- opal_msy(x,c(.5,.5),uncertainty="mcmc")
  points <- opal_derived(x,"msy")$draws
  expect_true(all(points$resolved))
  expect_gt(abs(diff(points$msy)),1)
})
