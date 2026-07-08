# CLAUDE.md — opal

RTMB-based, age-structured (soon optionally length-structured) fisheries stock-assessment
R package. Working dir: `/home/philipp/dragonfly/Misc/opal`.

## Model architecture

**Entry point** `opal_model(parameters, data)` (`R/model.R`) — the RTMB NLL. Flow:
1. `getAll(data, parameters)`; guard optional switches (`cpue_switch`, `lf_switch`, `wf_switch`, ...).
2. **Growth module**: `mu_a = get_growth(n_age, A1, A2, L1, L2, log_k)` (Schnute VB, `R/growth.R`);
   `sd_a = get_sd_at_age(...)` (linear CV in length, SS3 CV_Growth_Pattern=2);
   `pla = get_pla(len_lower, len_upper, mu_a, sd_a)` ([n_len,n_age], cols sum to 1, `R/selectivity.R`).
3. **Biology to age via PLA**: `weight_a = t(pla) %*% wt_at_len`; `maturity`/`M`/`fecundity`/
   `sex_ratio` via `resolve_bio_vector()` (age vec passes through; length vec => `t(pla) %*% v`).
   `spawning_potential_a = sex_ratio_a*maturity_a*fecundity_a` (or supplied).
4. **Selectivity**: `sel_fya = get_selectivity(data, par_sel, pla, len_mid)` — length curve
   (`sel_logistic`/`sel_double_normal`, real-line pars) collapsed to age via PLA; time-invariant.
   `sel_fa_external` overrides.
5. `get_initial_numbers()` (`R/dynamics.R`) — unfished survivorship => R0, BH alpha/beta; fished
   equilibrium; optional `init_rdev_a`.
6. `do_dynamics()` (`R/dynamics.R`) — the population loop.
7. `comp_pred_fya` = `catch_pred_fya`, but for **zero-catch (survey) fisheries** overwritten with
   `number_ysa[y,1,]*sel_fya[f,y,]` (composition from vulnerable numbers, not catch).
8. Priors (`evaluate_priors`), recruitment prior, CPUE/LF/WF likelihoods; `nll` summed.

**Dynamics** `do_dynamics()` — state `number_ysa` [n_year+1, n_season, n_age]. **Pope-style**
harvest rate (NOT Baranov): `F_f = catch_obs / sum(N*sel*weight)`, `h_rate_fa = F_f*sel_fya`,
capped by `posfun(1 - sum_F)`. Seasonal survival `S_a = exp(-M_a/n_season)`. Ageing = age shift
with plus-group accumulation at `n_age`. BH recruitment to age-1 (`get_recruitment`, `R/recruitment.R`).
`SSB = sum(number_ysa[y,1,]*spawning_potential_a)`. Biology passed as explicit args so AD
gradients reach growth params. ADREPORT: `spawning_biomass_y`, `spawning_biomass0_y`, `dynamic_depletion_y`.

**Length observations** (`R/likelihoods.R`): predicted length comp = `catch_pred_fya %*% t(pla)`;
weight comp = `catch_pred_fya %*% t(wf_rebin_matrix %*% pla)` (`wf_rebin_matrix` from
`prep_wf_data`, a proportional-overlap length->weight rebin, `R/rebin.R`). Likelihood types:
1=multinomial, 2=Dirichlet, 3=Dirichlet-multinomial. CPUE (`get_cpue_like`): log-linear with
per-index q/tau/omega/creep, vuln = `number_ysa[y,1,]*sel_fya[f,y,]` (xweight if units==1).

## Data object (via `get_data(model)`; models: opal_baseline, opakapaka, wcpo_bet)

Dims: `n_age, n_len, n_year, n_season, n_fishery, n_index`. Time: `first_yr, last_yr, years,
first_yr_catch`. Catch: `catch_obs_ysf` [n_year,1,n_fishery], `catch_units_f` (1=weight,2=numbers).
Length bins: `len_bin_start, len_bin_width, n_len` (=> `len_lower/len_upper/len_mid` derived in
`prep_lf_data`). Growth/biology: `A1, A2, lw_a, lw_b`, `maturity` & `fecundity` at **length**
(n_len), `M` at **age** (n_age). Selectivity: `sel_type_f` (1=logistic, 2/24=double-normal).
Comp data: `lf_*` / `wf_*` (flat/ints/prop, per-fishery split lists, min/max bin, switches),
`wt_bin_start/width/n_wt` (BET). `priors` list. Optional: `weight`, `spawning_potential`,
`sel_fa_external`, `sex_ratio`, `init_rdev_a`, `bias_adj_y`.
opakapaka: n_age=44, n_len=17, n_year=75, n_fishery=3. wcpo_bet: n_age=40, n_len=95, n_year=268,
n_fishery=15.

## Parameters (via `get_parameters(model)`)

`log_B0, log_h, log_sigma_r, log_cpue_q[n_index], cpue_creep, log_cpue_tau, log_cpue_omega,
rdev_y[n_year], par_sel[n_fishery,6]`, growth block `log_L1, log_L2, log_k, log_CV1, log_CV2`.
Optional: `log_init_F_f, init_rdev_a, log_lf_tau, log_wf_tau`. (`get_map`/`get_bounds` in
`R/parameters.R` reference some legacy names — treat as stale.)

## Projections (`R/projections.R`)

`project_dynamics(data, object, mcmc=NULL, n_proj, n_iter, rdev_y, sel_fya, catch_ysf,
return_hist=FALSE)`. `mcmc=NULL` => MVN draws from Hessian vcov (Cholesky); else MCMC posterior.
Perturbed terminal state per iter: `rep$number_ysa[n_year+1,1,]` as `init_number_a`. **MVN rdev
handling** uses the conditional expectation `E[rdev|scalar_draw]` to preserve the
log_B0<->rdev correlation and avoid Jensen inflation (see below). `project_rec_devs`,
`project_selectivity` support this.

### Jensen's-inequality warning (historical uncertainty)
Do NOT draw all rdev_y simultaneously from the full-parameter MVN to show historical SBY ribbons:
SE(rdev_y)~sigma_r~0.57, and E[exp(X)]>exp(E[X]) inflates the MVN mean SBY ~76%. Use
`sdreport()` + `ADREPORT(spawning_biomass_y)` (delta method) for the MVN historical ribbon; the
MCMC historical ribbon (`proj_mcmc$hist_sbio`) is fine.

## MCMC

`SparseNUTS::sample_snuts(obj, globals=bet_globals(), num_samples, num_warmup, chains, cores,
metric="dense", seed, control=list(adapt_delta, max_treedepth))`. `bet_globals()` (`R/model.R`)
exports all helper fns to parallel workers — **keep it in sync when adding functions**. Workers
call `library(opal)`, so `R CMD INSTALL` opal before sampling. opakapaka MCMC cache:
`inst/extdata/opaka_mcmc_fit.rds` (regenerate with tuned settings if the model changes —
survey composition / likelihood changes invalidate it).

## Data prep

`wcpo_bet_data` needs manual field additions before MakeADFun (len_lower/upper/mid,
`removal_switch_f`) unless `prep_lf_data()` is run. `prep_lf_data` / `prep_wf_data` are the
single source of truth for bin geometry and the LF/WF split lists; `prep_wf_data` also builds
`wf_rebin_matrix` (length->weight).

## Environment / rendering

- **renv**: set `RENV_CONFIG_AUTOLOADER_ENABLED=FALSE` for all R/render invocations to use the
  global R 4.6.1 library (avoids the project renv isolation).
- Vignettes are **Quarto** (`.qmd`, `quarto::html`); projection vignettes use `eval: true` with a
  cached results bundle. Plots: no titles/subtitles — required text goes in `fig-cap`
  (skyblue/orange/seagreen/purple palette, faded gridlines, patchwork multi-panel).
- The documentation paper `opal-documentation/assets/opal-sc22-wp.qmd` renders to PDF via
  pdflatex; its cache lives at `assets/cached/opakapaka/` (`doc_dir` already resolves to `assets/`
  — do not re-append `assets`).

## Length-dynamics extension (in progress)

An opt-in joint age×length engine (growth transition matrix, mirroring `lbm/`) is being added:
`data$length_dynamics=1` dispatches to `do_dynamics_length()` (state `number_ysal`
[y,s,a,l]) with `get_growth_matrix()` (conditional-SD normal-CDF transition), length-basis
selectivity/biology, and LF/WF likelihoods on the length marginal (via identity PLA). Age-only
path is unchanged when the flag is absent. Growth `G` is built in-model from existing Schnute
params (`Linf = L1 + (L2-L1)/(1-exp(-k*(A2-A1)))`, `rho = exp(-exp(log_k))`). The transition is a
mean-anchored binned normal (same discretisation as `get_pla`) so it reduces to the static PLA
under no fishing — do NOT enforce no-shrinkage (it drifts length-at-age upward and inflates SSB).
Showcased in `vignette("length_dynamics")` (`vignettes/length_dynamics.qmd`, eval:true). See
`.claude/plans/` for signatures, the growth-matrix formula, and the validation checklist.

## Recent model changes (Darcy, on `dev`)

`comp_pred_fya` for zero-catch survey fisheries; optional `init_rdev_a`/`init_bias_adj_a`;
`sel_fa_external` (fixed selectivity-at-age); corrected harvest-fraction penalty (logspace_add);
precomputed composition splits; modernized `get_data()`/`get_parameters()` accessors with
deterministic fishery ordering.
