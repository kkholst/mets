# Dermal ridges data (monozygotic twins)

Data on dermal ridge counts in left and right hand in (nuclear) families

## Format

Data on dermal ridge counts (left and right hand) in 18 monozygotic twin
pairs.

## Source

Sarah B. Holt (1952). Genetics of dermal ridges: bilateral asymmetry in
finger ridge-counts. Annals of Eugenics 17 (1), pp.211–231. DOI:
10.1111/j.1469-1809.1952.tb02513.x

## Examples

``` r
data(dermalridgesMZ)
fast.reshape(dermalridgesMZ,id="id",varying=c("left","right"))
#>       sex left1 right1              group id left2 right2
#> 1  female    95     83 psychotic&neurotic  1    90     85
#> 3  female    93     90 psychotic&neurotic  2    80     80
#> 5    male    99     94 psychotic&neurotic  3    94     99
#> 7  female    55     61 psychotic&neurotic  4    48     54
#> 9  female    41     24 psychotic&neurotic  5    26     42
#> 11   male    39     38 psychotic&neurotic  6    49     38
#> 13 female    69     81 psychotic&neurotic  7    81     71
#> 15   male    78     72 psychotic&neurotic  8    77     77
#> 17   male    64     64 psychotic&neurotic  9    69     55
#> 19   male    74     83 psychotic&neurotic 10    78     87
#> 21 female    37     34 psychotic&neurotic 11    27     46
#> 23 female    39     45 psychotic&neurotic 12    45     53
#> 25 female    62     62             normal 13    69     66
#> 27   male    92     95             normal 14    90     89
#> 29 female    97     90             normal 15    92     99
#> 31 female    28     44             normal 16    44     48
#> 33 female    53     59             normal 17    51     53
#> 35   male    62     60             normal 18    56     51
```
