# Dermal ridges data (families)

Data on dermal ridge counts in left and right hand in (nuclear) families

## Format

Data on 50 families with ridge counts in left and right hand for moter,
father and each child. Family id in 'family' and gender and child number
in 'sex' and 'child'.

## Source

Sarah B. Holt (1952). Genetics of dermal ridges: bilateral asymmetry in
finger ridge-counts. Annals of Eugenics 17 (1), pp.211–231. DOI:
10.1111/j.1469-1809.1952.tb02513.x

## Examples

``` r
data(dermalridges)
fast.reshape(dermalridges,id="family",varying=c("child.left","child.right","sex"))
#>     father.left father.right mother.left mother.right   sex1 child.left1
#> 1            51           63           5           18   male          66
#> 2             7            8          78           82 female          30
#> 3            30           35          95           91   male          54
#> 4            39           17          68           67   male          74
#> 5            76           78          72           76 female          69
#> 6           105          108          97          104   male          88
#> 7            35           42          64           91   male          54
#> 8           102           97          34           36 female          71
#> 9            93          113          62           47 female          67
#> 10           76           78          64           71 female          75
#> 11           61           70          13           15   male          76
#> 12           91           98          74           83 female          84
#> 13           69           84          73           71   male          61
#> 14           59           57          87          110   male          94
#> 15           84           81          84           86 female          78
#> 16           77           82          30           30 female          72
#> 17           94          113          49           60   male          82
#> 18           42           58          61           81   male          45
#> 20           65           43          96           92 female          69
#> 22           25           40          67           53 female          44
#> 24           70           76          50           60   male          71
#> 26           93           95          29           31 female          72
#> 28           77           88          79           83   male          53
#> 30          114          112          94           87 female         103
#> 32           32           23          20           33 female          66
#> 34           61           37          59           62   male          64
#> 36           51           47          84           73   male          41
#> 38            6           23          67           72   male          38
#> 40           61           75          79           83   male          95
#> 42          113          116          91           91   male         106
#> 44           46           39          82           93   male          80
#> 46           82           90          17           25   male          52
#> 48           61           75          59           72 female          59
#> 50           69           68          46           69 female          33
#> 52           95          100          43           44   male         100
#> 54           93           96          46           57   male          76
#> 56           20           26          53           50 female          42
#> 59           28           38          37           37 female          12
#> 62           49           59          62           73 female          31
#> 65           32           34          15           20   male          35
#> 68           57           69          92           86   male          93
#> 71           95          106         115          115 female         116
#> 74           44           57          79           77   male          68
#> 77           52           50          77           98   male          76
#> 80           46           67          38           50   male          81
#> 83           87           85          11            7 female          19
#> 87          110          114          97           84 female         103
#> 91           90           98          20           24   male          87
#> 96           90           87          68           64   male         106
#> 101          66           62          15           35 female          42
#>     child.right1 family child note child.left2 child.right2   sex2 child.left3
#> 1             67      1     1 <NA>          NA           NA   <NA>          NA
#> 2             26      2     1 <NA>          NA           NA   <NA>          NA
#> 3             74      3     1 <NA>          NA           NA   <NA>          NA
#> 4             77      4     1 <NA>          NA           NA   <NA>          NA
#> 5             62      5     1 <NA>          NA           NA   <NA>          NA
#> 6            101      6     1 <NA>          NA           NA   <NA>          NA
#> 7             59      7     1 <NA>          NA           NA   <NA>          NA
#> 8             78      8     1 <NA>          NA           NA   <NA>          NA
#> 9             59      9     1 <NA>          NA           NA   <NA>          NA
#> 10            76     10     1 <NA>          NA           NA   <NA>          NA
#> 11            89     11     1 <NA>          NA           NA   <NA>          NA
#> 12            97     12     1 <NA>          NA           NA   <NA>          NA
#> 13            75     13     1 <NA>          NA           NA   <NA>          NA
#> 14            88     14     1 <NA>          NA           NA   <NA>          NA
#> 15            69     15     1 <NA>          NA           NA   <NA>          NA
#> 16            79     16     1 <NA>          NA           NA   <NA>          NA
#> 17            89     17     1 <NA>          NA           NA   <NA>          NA
#> 18            44     18     1 <NA>          49           53 female          NA
#> 20            69     19     1 <NA>          75           76   male          NA
#> 22            45     20     1 <NA>          68           77 female          NA
#> 24            77     21     1 <NA>          41           54   male          NA
#> 26            75     22     1 <NA>          62           63 female          NA
#> 28            66     23     1 <NA>          87           83   male          NA
#> 30           101     24     1 <NA>         100           93 female          NA
#> 32            63     25     1 <NA>          40           26 female          NA
#> 34            81     26     1   dz          66           65   male          NA
#> 36            37     27     1 <NA>          61           57 female          NA
#> 38            39     28     1 <NA>          12           12   male          NA
#> 40            88     29     1 <NA>          62           57   male          NA
#> 42           106     30     1 <NA>          97           90   male          NA
#> 44            81     31     1 <NA>          82           82 female          NA
#> 46            51     32     1 <NA>          53           73   male          NA
#> 48            78     33     1 <NA>          47           48 female          NA
#> 50            43     34     1 <NA>          81           85   male          NA
#> 52           108     35     1 <NA>          66           45   male          NA
#> 54            80     36     1 <NA>          72           71 female          NA
#> 56            39     37     1 <NA>          28           34 female          49
#> 59            12     38     1 <NA>           5            4 female          17
#> 62            51     39     1 <NA>          46           50 female          58
#> 65            44     40     1 <NA>          26           38   male          16
#> 68           100     41     1 <NA>          99          108 female          53
#> 71           124     42     1 <NA>          89          101 female          99
#> 74            74     43     1 <NA>          57           67 female          52
#> 77            97     44     1 <NA>          89           77   male          64
#> 80            65     45     1 <NA>          59           42   male          57
#> 83            45     46     1 <NA>          76           75   male           4
#> 87           104     47     1 <NA>          99          104 female          99
#> 91            90     48     1 <NA>          43           50 female          15
#> 96            93     49     1 <NA>          96          103   male          76
#> 101           65     50     1 <NA>           0            0 female           0
#>     child.right3   sex3 child.left4 child.right4   sex4 child.left5
#> 1             NA   <NA>          NA           NA   <NA>          NA
#> 2             NA   <NA>          NA           NA   <NA>          NA
#> 3             NA   <NA>          NA           NA   <NA>          NA
#> 4             NA   <NA>          NA           NA   <NA>          NA
#> 5             NA   <NA>          NA           NA   <NA>          NA
#> 6             NA   <NA>          NA           NA   <NA>          NA
#> 7             NA   <NA>          NA           NA   <NA>          NA
#> 8             NA   <NA>          NA           NA   <NA>          NA
#> 9             NA   <NA>          NA           NA   <NA>          NA
#> 10            NA   <NA>          NA           NA   <NA>          NA
#> 11            NA   <NA>          NA           NA   <NA>          NA
#> 12            NA   <NA>          NA           NA   <NA>          NA
#> 13            NA   <NA>          NA           NA   <NA>          NA
#> 14            NA   <NA>          NA           NA   <NA>          NA
#> 15            NA   <NA>          NA           NA   <NA>          NA
#> 16            NA   <NA>          NA           NA   <NA>          NA
#> 17            NA   <NA>          NA           NA   <NA>          NA
#> 18            NA   <NA>          NA           NA   <NA>          NA
#> 20            NA   <NA>          NA           NA   <NA>          NA
#> 22            NA   <NA>          NA           NA   <NA>          NA
#> 24            NA   <NA>          NA           NA   <NA>          NA
#> 26            NA   <NA>          NA           NA   <NA>          NA
#> 28            NA   <NA>          NA           NA   <NA>          NA
#> 30            NA   <NA>          NA           NA   <NA>          NA
#> 32            NA   <NA>          NA           NA   <NA>          NA
#> 34            NA   <NA>          NA           NA   <NA>          NA
#> 36            NA   <NA>          NA           NA   <NA>          NA
#> 38            NA   <NA>          NA           NA   <NA>          NA
#> 40            NA   <NA>          NA           NA   <NA>          NA
#> 42            NA   <NA>          NA           NA   <NA>          NA
#> 44            NA   <NA>          NA           NA   <NA>          NA
#> 46            NA   <NA>          NA           NA   <NA>          NA
#> 48            NA   <NA>          NA           NA   <NA>          NA
#> 50            NA   <NA>          NA           NA   <NA>          NA
#> 52            NA   <NA>          NA           NA   <NA>          NA
#> 54            NA   <NA>          NA           NA   <NA>          NA
#> 56            40   male          NA           NA   <NA>          NA
#> 59            11   male          NA           NA   <NA>          NA
#> 62            74 female          NA           NA   <NA>          NA
#> 65            18   male          NA           NA   <NA>          NA
#> 68            79 female          NA           NA   <NA>          NA
#> 71            98   male          NA           NA   <NA>          NA
#> 74            54 female          NA           NA   <NA>          NA
#> 77            72 female          NA           NA   <NA>          NA
#> 80            50   male          NA           NA   <NA>          NA
#> 83             0 female          35           31 female          NA
#> 87           102 female         101           93   male          NA
#> 91             8   male          54           42 female          36
#> 96            82   male          83           85   male          76
#> 101            8 female          76           77 female          30
#>     child.right5   sex5 child.left6 child.right6   sex6
#> 1             NA   <NA>          NA           NA   <NA>
#> 2             NA   <NA>          NA           NA   <NA>
#> 3             NA   <NA>          NA           NA   <NA>
#> 4             NA   <NA>          NA           NA   <NA>
#> 5             NA   <NA>          NA           NA   <NA>
#> 6             NA   <NA>          NA           NA   <NA>
#> 7             NA   <NA>          NA           NA   <NA>
#> 8             NA   <NA>          NA           NA   <NA>
#> 9             NA   <NA>          NA           NA   <NA>
#> 10            NA   <NA>          NA           NA   <NA>
#> 11            NA   <NA>          NA           NA   <NA>
#> 12            NA   <NA>          NA           NA   <NA>
#> 13            NA   <NA>          NA           NA   <NA>
#> 14            NA   <NA>          NA           NA   <NA>
#> 15            NA   <NA>          NA           NA   <NA>
#> 16            NA   <NA>          NA           NA   <NA>
#> 17            NA   <NA>          NA           NA   <NA>
#> 18            NA   <NA>          NA           NA   <NA>
#> 20            NA   <NA>          NA           NA   <NA>
#> 22            NA   <NA>          NA           NA   <NA>
#> 24            NA   <NA>          NA           NA   <NA>
#> 26            NA   <NA>          NA           NA   <NA>
#> 28            NA   <NA>          NA           NA   <NA>
#> 30            NA   <NA>          NA           NA   <NA>
#> 32            NA   <NA>          NA           NA   <NA>
#> 34            NA   <NA>          NA           NA   <NA>
#> 36            NA   <NA>          NA           NA   <NA>
#> 38            NA   <NA>          NA           NA   <NA>
#> 40            NA   <NA>          NA           NA   <NA>
#> 42            NA   <NA>          NA           NA   <NA>
#> 44            NA   <NA>          NA           NA   <NA>
#> 46            NA   <NA>          NA           NA   <NA>
#> 48            NA   <NA>          NA           NA   <NA>
#> 50            NA   <NA>          NA           NA   <NA>
#> 52            NA   <NA>          NA           NA   <NA>
#> 54            NA   <NA>          NA           NA   <NA>
#> 56            NA   <NA>          NA           NA   <NA>
#> 59            NA   <NA>          NA           NA   <NA>
#> 62            NA   <NA>          NA           NA   <NA>
#> 65            NA   <NA>          NA           NA   <NA>
#> 68            NA   <NA>          NA           NA   <NA>
#> 71            NA   <NA>          NA           NA   <NA>
#> 74            NA   <NA>          NA           NA   <NA>
#> 77            NA   <NA>          NA           NA   <NA>
#> 80            NA   <NA>          NA           NA   <NA>
#> 83            NA   <NA>          NA           NA   <NA>
#> 87            NA   <NA>          NA           NA   <NA>
#> 91            32   male          NA           NA   <NA>
#> 96            74 female          NA           NA   <NA>
#> 101           26 female          59           61 female
```
