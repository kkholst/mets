# Transform that allows condition

Defines new variables under condition for data frame

## Usage

``` r
dtransform(data, ...)
```

## Arguments

- data:

  is data frame

- ...:

  new variable definitions including possible if condition

## Examples

``` r
data(mena)

xx <- dtransform(mena,ll=log(agemena)+twinnum)

xx <- dtransform(mena,ll=log(agemena)+twinnum,agemena<15)
xx <- dtransform(xx  ,ll=100+agemena,ll2=1000,agemena>15)
dsummary(xx,ll+ll2~I(agemena>15))
#> I(agemena > 15): FALSE
#>        ll             ll2      
#>  Min.   :3.227   Min.   : NA   
#>  1st Qu.:3.557   1st Qu.: NA   
#>  Median :3.702   Median : NA   
#>  Mean   :4.048   Mean   :NaN   
#>  3rd Qu.:4.557   3rd Qu.: NA   
#>  Max.   :4.707   Max.   : NA   
#>                  NA's   :1859  
#> ------------------------------------------------------------ 
#> I(agemena > 15): TRUE
#>        ll             ll2      
#>  Min.   :115.0   Min.   :1000  
#>  1st Qu.:115.3   1st Qu.:1000  
#>  Median :115.5   Median :1000  
#>  Mean   :115.7   Mean   :1000  
#>  3rd Qu.:115.9   3rd Qu.:1000  
#>  Max.   :118.2   Max.   :1000  
```
