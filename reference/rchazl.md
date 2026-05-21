# Multiple Cause Piecewise Constant Hazard Simulation

Simulates data from multiple piecewise constant baseline hazards for
competing risks. Takes the minimum of all cause-specific event times and
assigns the corresponding cause.

## Usage

``` r
rchazl(cumhaz, rr, causes = NULL, ...)
```

## Arguments

- cumhaz:

  List of cumulative hazard matrices, one for each cause.

- rr:

  Matrix of relative risks (rows = subjects, columns = causes).

- causes:

  Vector of cause codes to assign (default NULL, uses 1,2,...).

- ...:

  Additional arguments passed to `rchaz`.

## Value

A data frame with event times and status indicating the cause of the
first event.

## See also

[`rchaz`](http://kkholst.github.io/mets/reference/rchaz.md)

## Author

Thomas Scheike
