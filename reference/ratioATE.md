# Ratio of Average Treatment Effects

Computes the ratio of two Average Treatment Effects (ATEs), typically
comparing the ATE for a specific cause (e.g., RMTL due to cause 1) to
the ATE for the total RMTL.

## Usage

``` r
ratioATE(rmtl, rmtl1, h = NULL, null = 1)
```

## Arguments

- rmtl:

  Object containing the ATE for the total RMTL (from `resmeanATE`).

- rmtl1:

  Object containing the ATE for the specific cause RMTL (from
  `resmeanATE`).

- h:

  Transformation function (e.g., `log`) applied to the ratio. Default is
  identity.

- null:

  Value under the null hypothesis for the ratio (default 1).

## Value

A list containing:

- ratioG:

  Ratio based on the simple IPCW estimator.

- ratioDR:

  Ratio based on the double robust estimator.

## Details

The function transforms the estimates (e.g., using log) to compute the
ratio and its standard error using the delta method and the joint
influence functions.

## See also

[`resmeanATE`](http://kkholst.github.io/mets/reference/resmeanATE.md)

## Author

Thomas Scheike

## Examples

``` r
data(bmt); bmt$event <- bmt$cause!=0; dfactor(bmt) <- tcell~tcell
out <- resmeanATE(Event(time,event)~tcell+platelet, data=bmt, time=40, outcome="rmtl")
out1 <- resmeanATE(Event(time,cause)~tcell+platelet, data=bmt, cause=1, time=40, outcome="rmtl")
ratioATE(out, out1, h=log)
```
