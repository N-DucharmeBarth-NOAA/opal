# Indian Ocean yellowfin tuna (simulated) SS3 and truth reference values

Quantities from the SS3 model and from the simulation operating model,
for comparison with opal fits to
[io_yft_data](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/io_yft_data.md).

## Usage

``` r
io_yft_reference
```

## Format

A named list:

- fleet_names:

  Names of the 20 opal fisheries.

- ss3:

  SS3 results: `version`; unfished spawning biomass `ssb0` (t);
  `trajectory` (`ts`, `ssb`, `static_dep`); `recdevs` (`ts`, `recdev`);
  `cpue` (`ts`, `fishery`, `obs`, `expected`, `se`), each quarterly
  series scaled to mean one; `lf` (`fishery`, `bin`, `expected`),
  expected length compositions aggregated per fishery weighted by sample
  size; `sel_age` (`fishery`, `age`, `sel`); `sel_pattern` (`fleet`,
  `name`, `pattern`, `special`), SS3's `age_selex` pattern per fleet (20
  = double normal, 12 = logistic, 15 = mirror of fleet `special`);
  `sel_par` (`fleet`, `pattern`, `par`, `value`, `phase`, `lo`, `hi`,
  `label`), SS3's estimated age selectivity parameters P1–P6, with
  `phase > 0` for the ones SS3 estimated (see
  [`convert_ss3_selex_to_rtmb()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/convert_ss3_selex_to_rtmb.md)
  with `basis = "age"`); and `length_at_age` (`age`, `length`, start of
  the time step).

- truth:

  The operating model (SPM, 159 cells on a 5 degree grid), summed over
  the whole domain: `b0`, unfished equilibrium spawning biomass (t), and
  `trajectory` (`ts`, `ssb`, `static_dep`).

- schnute_length_at_age:

  Length-at-age (`age`, `length`) of the single-k Schnute curve used
  before opal supported age-specific K.

- history:

  B0 and terminal depletion of three opal fits: Schnute growth with
  SS3's selectivity, SS3 growth with SS3's selectivity, and SS3 growth
  with opal's estimated selectivity (C2).

- source:

  Provenance of the export: repository, script, commit, date and package
  versions.

## Source

SS3 Report.sso and data.ss_new expected values; operating model outputs
from
<https://github.com/aaronmberger-nwfsc/Spatial-Assessment-Modeling-Workshop>.
Exported by `code/rtmb/scripts/04-export-io-yft-for-opal.r` in
<https://github.com/N-DucharmeBarth/wcpfc-bet-single-region>.
