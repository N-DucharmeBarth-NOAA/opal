test_that("plot_catch handles named array dimensions", {
  observed <- array(
    c(10, 20, 30, 40),
    dim = c(2, 1, 2),
    dimnames = list(
      year = c("1", "2"),
      season = "1",
      fishery = c("1", "2")
    )
  )
  predicted <- observed + 1
  data <- list(
    years = 2001:2002,
    n_year = 2L,
    n_fishery = 2L,
    catch_obs_ysf = observed
  )
  obj <- list(report = function() list(catch_pred_ysf = predicted))

  expect_message(
    plot <- plot_catch(data, obj),
    "maximum absolute catch difference"
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(sort(unique(plot$data$year)), 2001:2002)
  expect_equal(unique(as.character(plot$data$season)), "Season: 1")
  expect_equal(
    levels(plot$data$fishery),
    c("Fishery: 1 (weight units)", "Fishery: 2 (weight units)")
  )
  expect_equal(plot$data$resid, rep(-1, 4))
})

test_that("plot_catch derives years when data$years is absent", {
  observed <- array(c(10, 20), dim = c(2, 1, 1))
  data <- list(
    first_yr = 1990L,
    n_year = 2L,
    n_fishery = 1L,
    catch_obs_ysf = observed
  )
  obj <- list(report = function() list(catch_pred_ysf = observed))

  plot <- suppressMessages(plot_catch(data, obj, plot_resid = TRUE))

  expect_s3_class(plot, "ggplot")
  expect_equal(plot$data$year, 1990:1991)
  expect_equal(plot$data$resid, c(0, 0))
})

test_that("catch panels and spawning plots respect model units and display scaling", {
  data <- list(first_yr=2000,last_yr=2001,n_year=2,n_fishery=2,
    catch_obs_ysf=array(10,c(2,1,2)),catch_units_f=c(1L,2L))
  object <- list(report=function()list(catch_pred_ysf=data$catch_obs_ysf,
    B0=1000,spawning_biomass_y=c(1000,800,600)))
  p <- suppressMessages(plot_catch(data,object,weight_units="kg"))
  expect_equal(levels(p$data$fishery),c("Fishery: 1 (kg)","Fishery: 2 (numbers)"))
  p <- plot_biomass_spawning(list(data),list(object),relative=FALSE,units="eggs",scale=1000)
  expect_equal(p$data$value,c(1,.8,.6))
  expect_equal(p$labels$y,"Spawning output (eggs / 1000)")
  p <- plot_biomass_spawning(list(data),list(object),relative=FALSE)
  expect_equal(p$data$value,c(1000,800,600))
  expect_equal(p$labels$y,"Spawning output (model units)")
  expect_error(plot_biomass_spawning(list(data),list(object),scale=0),"positive")
})
