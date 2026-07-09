# State-space numbers-at-age development plan

## Objective

Evaluate whether opal can estimate annual numbers-at-age as latent states with
process error while retaining the existing seasonal catch, survival,
stock-recruitment, and observation components. The first prototype should be a
reversible extension of the deterministic model, not a replacement for it.

## Prototype formulation

The initial state remains the equilibrium state produced by
`get_initial_numbers()`. For transition year `y` and age `a`, the existing
dynamics produce a deterministic prediction `N_pred[y, a]`. The latent state is
parameterised as `log_number_state_ya[y, a]` and follows

```text
log N[y + 1, a] ~ Normal(
  log N_pred[y, a] - 0.5 * sigma_state[a]^2,
  sigma_state[a]
)
```

The bias correction is configurable. A scalar process standard deviation is
the first estimation target; the implementation also accepts an age-specific
vector for structured experiments. Within-year seasonal dynamics remain
deterministic. The latent state replaces the predicted state at the next year
boundary, so process deviations propagate into catches, indices, compositions,
and subsequent transitions.

Use `state_space_switch = 1L` in the data and add parameters with
`initialize_state_space_parameters()`. In fitted RTMB models, pass
`random = "log_number_state_ya"` to use the Laplace approximation. Keep the
switch absent or zero for the current deterministic behavior.

## Work phases

### 1. Establish deterministic parity — initiated here

- Preserve the current objective and trajectories when state-space mode is
  disabled.
- Initialise latent states from a deterministic report.
- Confirm that matching latent and predicted states reproduce the deterministic
  trajectory when lognormal bias correction is disabled.
- Report transition predictions, standardised residuals, and the process
  likelihood separately.

Exit criterion: unit tests cover switch-off parity, zero-residual transitions,
state propagation, dimensions, and registration of states as RTMB random
effects.

### 2. Fit the smallest identifiable model

- Start with one shared `sigma_state` fixed at a short grid of plausible values
  (for example 0.02, 0.05, 0.10, and 0.20).
- Fix `rdev_y = 0` in the first experiments. Estimating both recruitment
  deviations and an age-1 state innovation would otherwise allocate the same
  variation twice.
- Estimate latent states with the Laplace approximation, then estimate the
  shared process standard deviation only after fixed-sigma fits are stable.
- Use a reduced synthetic case before the full BET and opakapaka data sets.

Exit criterion: stable convergence, positive-definite fixed-effect Hessian,
small gradients, finite random-effect Hessian, and recovery of simulated
states/process variance without material deterministic-limit bias.

### 3. Compare process-error structures

Evaluate only after the scalar model is stable:

1. one shared standard deviation across all ages;
2. recruitment versus survival standard deviations;
3. smooth age-varying standard deviations with few hyperparameters;
4. cohort-correlated or age-correlated innovations; and
5. process error on survival/mortality rates rather than unconstrained state
   innovations.

Prefer the lowest-dimensional structure supported by simulation recovery and
predictive diagnostics. A free standard deviation for every age is available
for exploration but is not the default scientific recommendation.

### 4. Observation and catch treatment

- Test sensitivity to treating catch as exact, as in the current dynamics.
- If catch error is needed, move catch to an observation likelihood before
  interpreting process variance; otherwise process error may absorb catch
  error.
- Check whether CPUE and composition data identify cohort trajectories rather
  than only aggregate biomass.
- Verify that the plus group does not accumulate systematic residuals caused by
  misspecified ageing or mortality.

Exit criterion: process residuals show no strong year, age, cohort, or fishery
pattern and variance is not acting as a substitute for a known observation
error source.

### 5. Diagnostics and model comparison

- Add plots for observed versus predicted indices/compositions, latent versus
  transition-predicted states, standardised innovations, cohort traces, and
  process variance profiles.
- Run simulation-estimation tests across data-rich and data-poor scenarios.
- Compare deterministic and state-space fits using marginal likelihood-based
  criteria only when their treatment of random effects is consistent.
- Measure runtime and memory growth with years, ages, and seasons.

Exit criterion: a written go/no-go decision records recovery, diagnostics,
computational cost, and consequences for management quantities.

### 6. Production integration if supported

- Add a public model-building helper that consistently sets the switch,
  parameter map, random-effect declaration, bounds, and initial states.
- Decide how projections propagate process uncertainty.
- Add end-to-end examples and migration notes.
- Keep the deterministic mode as a tested special case.

## Immediate experiment matrix

| Case | Data | State sigma | Recruitment deviations | Purpose |
|---|---|---:|---|---|
| D0 | synthetic | off | fixed zero | deterministic reference |
| S1 | synthetic | 0.05 fixed | fixed zero | latent-state recovery |
| S2 | synthetic | 0.10 fixed | fixed zero | variance sensitivity |
| S3 | synthetic | estimated scalar | fixed zero | variance estimability |
| S4 | synthetic | estimated scalar | estimated | confounding stress test |
| B1 | BET subset | best fixed value | fixed zero | scaling and diagnostics |
| B2 | BET full | estimated scalar | fixed zero | full prototype evaluation |

## Known limitations of the initiated code

- Process error is applied at annual boundaries, not every season.
- The initial equilibrium state is conditioned on rather than estimated.
- The same state process currently covers recruitment and survival ages.
- Catch is still treated as exact and converted to harvest rates internally.
- No correlation structure is yet imposed across ages or cohorts.
- A convenience wrapper for constructing the full RTMB object is deferred until
  the simplest formulation demonstrates stable estimation.
