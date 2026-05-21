# Simulation Helper Functions

Internal simulation functions used for generating data from various
survival, competing risks, frailty, and twin/family models. These
functions are primarily intended for use in examples and testing.

## Survival and Competing Risks

- `simrchaz`:

  Simulate from cumulative hazard via inverse CDF.

- `simul_cifs`:

  Simulate competing risks data from cumulative incidence functions.

- `simlogitSurvd`:

  Simulate survival data using logistic model.

- `kumarsim`:

  Simulate competing risks data (Kumaraswamy-type).

- `kumarsimRCT`:

  Simulate RCT competing risks data (Kumaraswamy-type).

## Clayton-Oakes and Frailty Models

- `sim_ClaytonOakes_family_ace`:

  Simulate family data from Clayton-Oakes ACE model.

- `sim_ClaytonOakes_twin_ace`:

  Simulate twin data from Clayton-Oakes ACE model.

- `sim_Compete_twin_ace`:

  Simulate twin competing risks data with ACE frailty.

- `sim_Compete_simple`:

  Simulate competing risks data with shared frailty.

- `sim_Frailty_simple`:

  Simulate survival data with shared frailty.

- `sim_SurvFam`:

  Simulate family survival data with shared frailty.

## Binomial Twin/Family Models

- `sim_BinPlack`:

  Simulate paired binary data (Plackett model).

- `sim_BinFam`:

  Simulate binary family data.

- `sim_BinFam2`:

  Simulate binary family data (alternative parameterization).

- `sim_binClaytonOakes_family_ace`:

  Simulate binary family data (Clayton-Oakes ACE).

- `sim_binClaytonOakes_twin_ace`:

  Simulate binary twin data (Clayton-Oakes ACE).

- `sim_binClaytonOakes_pairs`:

  Simulate binary paired data (Clayton-Oakes).

- `sim_bptwin`:

  Simulate from a biprobit twin model.

## Nordic Twin Studies

- `sim_nordictwin`:

  Simulate Nordic twin registry data.

- `sim_nordic_random`:

  Simulate Nordic twin data with random effects.

- `corsim_prostate_random`:

  Simulate correlated prostate cancer data with random effects.
