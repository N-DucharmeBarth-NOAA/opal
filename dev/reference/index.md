# Package index

## Assessment objects

Create, fit, sample, inspect, and save one portable assessment.

- [`opal_obj()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_obj.md)
  : Create a portable Opal assessment object
- [`opal_build()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_build.md)
  [`opal_rtmb()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_build.md)
  : Build or access an Opal runtime objective
- [`opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_fit.md)
  : Fit an Opal assessment
- [`opal_update()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_update.md)
  : Update an Opal model configuration
- [`opal_attach_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_fit.md)
  : Attach an externally optimised Opal fit
- [`opal_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_mcmc.md)
  : Run MCMC for an Opal object
- [`opal_attach_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_mcmc.md)
  : Attach posterior draws to an Opal object
- [`opal_check()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_check.md)
  : Check fitted or sampled Opal results
- [`opal_report()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_report.md)
  : Report a configured or fitted Opal model
- [`opal_save()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_io.md)
  [`opal_read()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_io.md)
  : Save and read staged Opal objects
- [`summary(`*`<opal_obj>`*`)`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/summary.opal_obj.md)
  [`print(`*`<summary.opal_obj>`*`)`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/summary.opal_obj.md)
  [`print(`*`<opal_obj>`*`)`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/summary.opal_obj.md)
  : Summarise a staged Opal assessment
- [`validate_opal_obj()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/validate_opal_obj.md)
  : Validate a staged Opal object
- [`opal_from_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_from_fit.md)
  : Convert a legacy Opal fit to the staged object workflow
- [`opal_as_tmbfit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_as_tmbfit.md)
  : Convert normalized opal posterior draws to a tmbfit

## Assessment diagnostics and summaries

- [`opal_diagnose()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_diagnose.md)
  : Diagnose biological feasibility of an Opal model
- [`opal_osa()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_osa.md)
  : Calculate one-step-ahead observation residuals
- [`plot_osa_sdnr()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_osa_sdnr.md)
  : Compare OSA residual dispersion across datasets
- [`plot_osa_residuals()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_osa_residuals.md)
  : Plot stored OSA residuals
- [`plot_composition()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_composition.md)
  : Plot observed and fitted compositions
- [`plot_prior_posterior()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_prior_posterior.md)
  : Compare parameter priors and posterior distributions
- [`opal_posterior()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_posterior.md)
  : Summarise posterior parameters and derived model quantities
- [`opal_derived()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_derived.md)
  : Retrieve a stored assessment analysis
- [`opal_example_inputs()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_example_inputs.md)
  : Simulated inputs for the assessment-tools tutorial

## Sensitivities and reference points

- [`opal_profile()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_profile.md)
  : Profile an assessment parameter
- [`plot_opal_profile()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_opal_profile.md)
  : Plot an objective profile
- [`opal_grid()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_grid.md)
  : Fit a reproducible grid of assessment scenarios
- [`opal_grid_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_grid_mcmc.md)
  : Sample accepted members of an assessment grid
- [`opal_grid_draws()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_grid_draws.md)
  : Select balanced posterior draws from a model grid
- [`opal_msy()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_msy.md)
  : Calculate deterministic equilibrium reference points

## Projections

Specify future assumptions and retain projections with their source
assessment.

- [`opal_project()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_project.md)
  : Project an Opal object and retain the result
- [`project_dynamics()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/project_dynamics.md)
  : Project dynamics
- [`project_rec_devs()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/project_rec_devs.md)
  : Project recruitment deviates
- [`project_selectivity()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/project_selectivity.md)
  : Project selectivity

## Plots

Pass an assessment object directly to the model plots.

- [`plot_biomass_spawning()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_biomass_spawning.md)
  : Plot spawning biomass
- [`plot_catch()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_catch.md)
  : Plot catch
- [`plot_cpue()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_cpue.md)
  : Plot CPUE
- [`plot_hrate()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_hrate.md)
  : Plot harvest rate
- [`plot_initial_numbers()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_initial_numbers.md)
  : Plot initial numbers
- [`plot_recruitment()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_recruitment.md)
  : Plot recruitment

## Model inputs and diagnostics

Prepare data, configure parameters, and inspect the underlying RTMB
model.

- [`opaka_quickstart_inputs()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opaka_quickstart_inputs.md)
  : Opakapaka quickstart model inputs
- [`get_data()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_data.md)
  : Get bundled model data
- [`get_parameters()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_parameters.md)
  : Get bundled initial parameter values
- [`get_map()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_map.md)
  : Get default parameter mapping
- [`get_bounds()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_bounds.md)
  : Get default parameter bounds
- [`get_priors()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_priors.md)
  : Get priors
- [`prep_lf_data()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/prep_lf_data.md)
  : Prepare length composition data for model input
- [`prep_wf_data()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/prep_wf_data.md)
  : Prepare weight composition data for model input
- [`check_bounds()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/check_bounds.md)
  : Check if parameters are up against the bounds
- [`check_estimability()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/check_estimability.md)
  : Check for identifiability of fixed effects
- [`get_cor_pairs()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_cor_pairs.md)
  : Return strongly-correlated parameter pairs
- [`get_par_table()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_par_table.md)
  : Summarise model parameters in a table
- [`extract_fixed()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/extract_fixed.md)
  : Extract fixed effects

## Model components

Lower-level functions used by the assessment model and advanced
workflows. The natural-mortality plot takes a data list and an RTMB
runtime.

- [`plot_natural_mortality()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_natural_mortality.md)
  : Plot natural mortality
- [`cmb()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/cmb.md)
  : Helper to make closure
- [`convert_rtmb_selex_to_ss3()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/convert_rtmb_selex_to_ss3.md)
  : Convert RTMB real-line selectivity parameters back to SS3 natural
  scale
- [`convert_ss3_selex_to_rtmb()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/convert_ss3_selex_to_rtmb.md)
  : Convert SS3 selectivity parameters to RTMB real-line
  parameterization
- [`do_dynamics()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/do_dynamics.md)
  : Population dynamics
- [`double_richards_natural()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/double_richards_natural.md)
  : Convert double Richards parameters to natural scale
- [`evaluate_priors()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/evaluate_priors.md)
  : Evaluate priors
- [`get_bias_adj_vector()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_bias_adj_vector.md)
  : Calculate Recruitment Bias Adjustment Ramp
- [`get_cpue_like()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_cpue_like.md)
  : CPUE index likelihood (multi-index)
- [`get_growth()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_growth.md)
  : Compute mean length-at-age using the Schnute parameterization of VB
  growth
- [`get_harvest_rate()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_harvest_rate.md)
  : Harvest rate calculation
- [`get_initial_numbers()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_initial_numbers.md)
  : Initial numbers and Beverton-Holt parameters
- [`get_length_like()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_length_like.md)
  : Length Composition Likelihood
- [`get_maturity_at_age()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_maturity_at_age.md)
  : Convert maturity-at-length to maturity-at-age
- [`get_pla()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_pla.md)
  : Probability of length at age matrix (age-length key)
- [`get_recruitment()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_recruitment.md)
  : Calculate recruitment
- [`get_recruitment_prior()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_recruitment_prior.md)
  : Recruitment prior
- [`get_sd_at_age()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_sd_at_age.md)
  : Compute SD of length-at-age from CV1 and CV2
- [`get_selectivity()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_selectivity.md)
  : Compute selectivity-at-age from length-based selectivity curves
- [`get_unfished_init()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_unfished_init.md)
  : Compute unfished equilibrium quantities from natural mortality only
- [`get_weight_at_length()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_weight_at_length.md)
  : Compute weight at each length bin midpoint
- [`get_weight_like()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_weight_like.md)
  : Weight Composition Likelihood
- [`opal_globals()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_globals.md)
  : The opal globals
- [`opal_model()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_model.md)
  : The opal model
- [`posfun()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/posfun.md)
  : Positive Constraint Penalty Function
- [`rebin_counts()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/rebin_counts.md)
  : Linear area rebinning of frequency data
- [`rebin_matrix()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/rebin_matrix.md)
  : Compute rebinning weight matrix
- [`resolve_bio_vector()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/resolve_bio_vector.md)
  : Resolve a biology vector to age-basis
- [`sel_double_normal()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/sel_double_normal.md)
  : Double-normal selectivity as a function of length (SS3 pattern 24,
  full form)
- [`sel_double_richards()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/sel_double_richards.md)
  : Double Richards selectivity as a function of length
- [`sel_length()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/sel_length.md)
  : Selectivity-at-length for a single fishery, dispatched on type code
- [`sel_logistic()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/sel_logistic.md)
  : Logistic selectivity as a function of length

## Bundled data

Example inputs and numerical baseline datasets.

- [`opaka_data`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opaka_data.md)
  : 'Opakapaka Stock Assessment Data

- [`opaka_lf`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opaka_lf.md)
  : 'Opakapaka Length Frequency Data

- [`opaka_parameters`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opaka_parameters.md)
  : 'Opakapaka Stock Assessment Parameters

- [`opaka_truth`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opaka_truth.md)
  : Opakapaka SS3 OM/EM truth and EM output (extracted)

- [`opal_baseline`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_baseline.md)
  :

  `opal` model regression baseline

- [`opal_baseline_data`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_baseline_data.md)
  : Opal Baseline Data

- [`opal_baseline_map`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_baseline_map.md)
  : Opal Baseline Map

- [`opal_baseline_parameters`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_baseline_parameters.md)
  : Opal Baseline Parameters

- [`wcpo_bet_data`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/wcpo_bet_data.md)
  : West Central Pacific Ocean Bigeye Tuna (WCPO BET) Assessment Data

- [`wcpo_bet_lf`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/wcpo_bet_lf.md)
  : West Central Pacific Ocean Bigeye Tuna (WCPO BET) Length Frequency
  Data

- [`wcpo_bet_parameters`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/wcpo_bet_parameters.md)
  : West Central Pacific Ocean Bigeye Tuna (WCPO BET) Assessment
  Parameters

- [`wcpo_bet_wf`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/wcpo_bet_wf.md)
  : West Central Pacific Ocean Bigeye Tuna (WCPO BET) Weight Frequency
  Data

## Legacy fit compatibility

Compatibility helpers for older opal_fit files; use opal_obj for new
assessments.

- [`opal_fit_compatibility()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_fit_compatibility.md)
  : Check compatibility of a saved opal fit
- [`save_opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_fit_io.md)
  [`read_opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_fit_io.md)
  : Save and read portable opal fits
- [`opal_fit_object()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_fit_object.md)
  : Access the runtime objective for an opal fit
- [`opal_fit_report()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_fit_report.md)
  : Recreate the fitted model report
- [`rebuild_opal_object()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/rebuild_opal_object.md)
  : Rebuild the RTMB objective for a portable opal fit
- [`update_opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/update_opal_fit.md)
  : Update portable results attached to an opal fit
- [`validate_opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/validate_opal_fit.md)
  : Validate a portable opal fit
