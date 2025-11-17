# Simulate observations from a Weibull distribution

Simulate observations from the model with cumulative hazard given by
\$\$\Lambda(t) = \lambda\cdot t^s\$\$ where \\\lambda\\ is the *rate
parameter* and \\s\\ is the *shape parameter*.

## Usage

``` r
rweibullcox(n, rate, shape)
```

## Arguments

- n:

  (integer) number of observations

- rate:

  (numeric) rate parameter (can be a vector of size n)

- shape:

  (numeric) shape parameter (can be a vector of size n)

## Details

\[stats::rweibull()\] uses a different parametrization with cumulative
hazard given by \$\$H(t) = (t/b)^a,\$\$ i.e., the shape is the same
\\a:=s\\ but the scale paramter \\b\\ is related to rate paramter \\r\\
by \$\$r := b^{-a}\$\$

## See also

\[stats::rweibull()\]
