# Illness-Death Competing Risks with Two Causes of Death

Simulates data from an illness-death model with two causes of death from
the illness state. Covariate effects can be introduced via relative risk
terms.

## Usage

``` r
sim_multistateII(
  cumhaz,
  death.cumhaz,
  death.cumhaz2,
  n = NULL,
  rr = NULL,
  rd = NULL,
  rd2 = NULL,
  gamma23 = 0,
  gamma24 = 0,
  early2 = 10000,
  gap.time = FALSE,
  max.recurrent = 100,
  cens = NULL,
  rrc = NULL,
  extend = TRUE,
  ...
)
```

## Arguments

- cumhaz:

  Cumulative hazard from state 1 to 2.

- death.cumhaz:

  Cumulative hazard of death from state 1.

- death.cumhaz2:

  Cumulative hazard of death from state 2.

- n:

  Number of simulations.

- rr:

  Relative risks.

- rd:

  Relative risks for death from state 1.

- rd2:

  Relative risks for death from state 2.

- gamma23:

  Early effect parameters for death causes.

- gamma24:

  Early effect parameters for death causes.

- early2:

  Time threshold for early effect.

- gap.time:

  Gap time indicator.

- max.recurrent:

  Maximum recurrent events.

- cens:

  Censoring specification.

- rrc:

  Censoring relative risks.

- extend:

  Extend hazards.

- ...:

  Additional arguments.

## Value

Data frame with multi-state event history.

## Author

Thomas Scheike
