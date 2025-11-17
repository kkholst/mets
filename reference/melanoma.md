# The Melanoma Survival Data

The melanoma data frame has 205 rows and 7 columns. It contains data
relating to survival of patients after operation for malignant melanoma
collected at Odense University Hospital by K.T. Drzewiecki.

## Format

This data frame contains the following columns:

- no:

  a numeric vector. Patient code.

- status:

  a numeric vector code. Survival status. 1: dead from melanoma, 2:
  alive, 3: dead from other cause.

- days:

  a numeric vector. Survival time.

- ulc:

  a numeric vector code. Ulceration, 1: present, 0: absent.

- thick:

  a numeric vector. Tumour thickness (1/100 mm).

- sex:

  a numeric vector code. 0: female, 1: male.

## Source

Andersen, P.K., Borgan O, Gill R.D., Keiding N. (1993), *Statistical
Models Based on Counting Processes*, Springer-Verlag.

Drzewiecki, K.T., Ladefoged, C., and Christensen, H.E. (1980), Biopsy
and prognosis for cutaneous malignant melanoma in clinical stage I.
Scand. J. Plast. Reconstru. Surg. 14, 141-144.

## Examples

``` r
data(melanoma)
names(melanoma)
#> [1] "no"     "status" "days"   "ulc"    "thick"  "sex"   
```
