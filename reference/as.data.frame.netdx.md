# Extract Timed Edgelists for netdx Objects

This function extracts timed edgelists for objects of class `netdx` into
a data frame using the generic `as.data.frame` function.

## Usage

``` r
# S3 method for class 'netdx'
as.data.frame(x, row.names = NULL, optional = FALSE, sim = NULL, ...)
```

## Arguments

- x:

  An `EpiModel` object of class `netdx`.

- row.names:

  See
  [`as.data.frame.default()`](https://rdrr.io/r/base/as.data.frame.html).

- optional:

  See
  [`as.data.frame.default()`](https://rdrr.io/r/base/as.data.frame.html).

- sim:

  The simulation number to output. If `NULL`, then data from all
  simulations will be output.

- ...:

  See
  [`as.data.frame.default()`](https://rdrr.io/r/base/as.data.frame.html).

## Value

A data frame containing the data from `x`.

## Examples

``` r
# \donttest{
# Initialize and parameterize the network model
nw <- network_initialize(n = 100)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)

# Model estimation
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#> Starting simulated annealing (SAN)
#> Iteration 1 of at most 4
#> Finished simulated annealing
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.

# Simulate the network with netdx
dx <- netdx(est, nsims = 3, nsteps = 10, keep.tedgelist = TRUE,
            verbose = FALSE)

# Extract data from the first simulation
as.data.frame(dx, sim = 1)
#>    onset terminus tail head onset.censored terminus.censored duration edge.id
#> 1      0        7    1   18          FALSE             FALSE        7       1
#> 2      0        5    1   62          FALSE             FALSE        5       2
#> 3      0        2    2   58          FALSE             FALSE        2       3
#> 4      0        1    2   88          FALSE             FALSE        1       4
#> 5      0        6    3   24          FALSE             FALSE        6       5
#> 6      0       11    3   90          FALSE              TRUE       11       6
#> 7      0       11    4   98          FALSE              TRUE       11       7
#> 8      0       11    5   13          FALSE              TRUE       11       8
#> 9      0       10    7   76          FALSE             FALSE       10       9
#> 10     0        4    8   51          FALSE             FALSE        4      10
#> 11     0        4   10   69          FALSE             FALSE        4      11
#> 12     0        6   10  100          FALSE             FALSE        6      12
#> 13     0       11   11   29          FALSE              TRUE       11      13
#> 14     0        3   11   66          FALSE             FALSE        3      14
#> 15     0       11   12   67          FALSE              TRUE       11      15
#> 16     0       11   16   44          FALSE              TRUE       11      16
#> 17     0       11   17   63          FALSE              TRUE       11      17
#> 18     0        4   18   48          FALSE             FALSE        4      18
#> 19     0       11   18   90          FALSE              TRUE       11      19
#> 20     0       11   18   94          FALSE              TRUE       11      20
#> 21     0       11   20   96          FALSE              TRUE       11      21
#> 22     0       11   21   63          FALSE              TRUE       11      22
#> 23     0       11   21   84          FALSE              TRUE       11      23
#> 24     0       11   22   23          FALSE              TRUE       11      24
#> 25     0       11   24   32          FALSE              TRUE       11      25
#> 26     0       11   27   89          FALSE              TRUE       11      26
#> 27     0       11   28   45          FALSE              TRUE       11      27
#> 28     0        1   29   68          FALSE             FALSE        1      28
#> 29     0       11   30   57          FALSE              TRUE       11      29
#> 30     0       11   30   73          FALSE              TRUE       11      30
#> 31     0       11   32   41          FALSE              TRUE       11      31
#> 32     0       11   32   46          FALSE              TRUE       11      32
#> 33     0       11   32   53          FALSE              TRUE       11      33
#> 34     0       11   32   66          FALSE              TRUE       11      34
#> 35     0        9   35   40          FALSE             FALSE        9      35
#> 36     0       11   35   88          FALSE              TRUE       11      36
#> 37     0       11   37   38          FALSE              TRUE       11      37
#> 38     0       11   37   62          FALSE              TRUE       11      38
#> 39     0       10   37   93          FALSE             FALSE       10      39
#> 40     0       11   38   41          FALSE              TRUE       11      40
#> 41     0       11   38   49          FALSE              TRUE       11      41
#> 42     0        2   38   62          FALSE             FALSE        2      42
#> 43     0       11   41   99          FALSE              TRUE       11      43
#> 44     0        6   42   78          FALSE             FALSE        6      44
#> 45     0       11   43   93          FALSE              TRUE       11      45
#> 46     0       11   44   69          FALSE              TRUE       11      46
#> 47     0       11   49   88          FALSE              TRUE       11      47
#> 48     0       11   52   78          FALSE              TRUE       11      48
#> 49     0       11   57   83          FALSE              TRUE       11      49
#> 50     0       10   68   97          FALSE             FALSE       10      50
#> 51     0        7   74   83          FALSE             FALSE        7      51
#> 52     0       11   74  100          FALSE              TRUE       11      52
#> 53     0        1   89   93          FALSE             FALSE        1      53
#> 54     0       11   90   91          FALSE              TRUE       11      54
#> 55     0        2   91   99          FALSE             FALSE        2      55
#> 56     0        5   93   97          FALSE             FALSE        5      56
#> 57     1        6   41   93          FALSE             FALSE        5      57
#> 58     1       10    9   86          FALSE             FALSE        9      58
#> 59     2       11   38   99          FALSE              TRUE        9      59
#> 60     2       11   39   99          FALSE              TRUE        9      60
#> 61     2       10   29   36          FALSE             FALSE        8      61
#> 62     3       11   38   94          FALSE              TRUE        8      62
#> 63     4       11   21   82          FALSE              TRUE        7      63
#> 64     4       11   75   90          FALSE              TRUE        7      64
#> 65     4       11   42   89          FALSE              TRUE        7      65
#> 66     4       11   12   28          FALSE              TRUE        7      66
#> 67     5        8   60   73          FALSE             FALSE        3      67
#> 68     6       11   12   83          FALSE              TRUE        5      68
#> 69     6       11   26   81          FALSE              TRUE        5      69
#> 70     7        8   58   62          FALSE             FALSE        1      70
#> 71     7       11    1   61          FALSE              TRUE        4      71
#> 72     7       11   32   69          FALSE              TRUE        4      72
#> 73     7       11   57   97          FALSE              TRUE        4      73
#> 74     7       11   54   75          FALSE              TRUE        4      74
#> 75     8       11   74   96          FALSE              TRUE        3      75
#> 76     8       11   33   97          FALSE              TRUE        3      76
#> 77     8       11   10   21          FALSE              TRUE        3      77
#> 78     8       11    6   79          FALSE              TRUE        3      78
#> 79     9       11   88   94          FALSE              TRUE        2      79
#> 80    10       11    5   92          FALSE              TRUE        1      80

# Extract data from all simulations
as.data.frame(dx)
#>     sim onset terminus tail head onset.censored terminus.censored duration
#> 1     1     0        7    1   18          FALSE             FALSE        7
#> 2     1     0        5    1   62          FALSE             FALSE        5
#> 3     1     0        2    2   58          FALSE             FALSE        2
#> 4     1     0        1    2   88          FALSE             FALSE        1
#> 5     1     0        6    3   24          FALSE             FALSE        6
#> 6     1     0       11    3   90          FALSE              TRUE       11
#> 7     1     0       11    4   98          FALSE              TRUE       11
#> 8     1     0       11    5   13          FALSE              TRUE       11
#> 9     1     0       10    7   76          FALSE             FALSE       10
#> 10    1     0        4    8   51          FALSE             FALSE        4
#> 11    1     0        4   10   69          FALSE             FALSE        4
#> 12    1     0        6   10  100          FALSE             FALSE        6
#> 13    1     0       11   11   29          FALSE              TRUE       11
#> 14    1     0        3   11   66          FALSE             FALSE        3
#> 15    1     0       11   12   67          FALSE              TRUE       11
#> 16    1     0       11   16   44          FALSE              TRUE       11
#> 17    1     0       11   17   63          FALSE              TRUE       11
#> 18    1     0        4   18   48          FALSE             FALSE        4
#> 19    1     0       11   18   90          FALSE              TRUE       11
#> 20    1     0       11   18   94          FALSE              TRUE       11
#> 21    1     0       11   20   96          FALSE              TRUE       11
#> 22    1     0       11   21   63          FALSE              TRUE       11
#> 23    1     0       11   21   84          FALSE              TRUE       11
#> 24    1     0       11   22   23          FALSE              TRUE       11
#> 25    1     0       11   24   32          FALSE              TRUE       11
#> 26    1     0       11   27   89          FALSE              TRUE       11
#> 27    1     0       11   28   45          FALSE              TRUE       11
#> 28    1     0        1   29   68          FALSE             FALSE        1
#> 29    1     0       11   30   57          FALSE              TRUE       11
#> 30    1     0       11   30   73          FALSE              TRUE       11
#> 31    1     0       11   32   41          FALSE              TRUE       11
#> 32    1     0       11   32   46          FALSE              TRUE       11
#> 33    1     0       11   32   53          FALSE              TRUE       11
#> 34    1     0       11   32   66          FALSE              TRUE       11
#> 35    1     0        9   35   40          FALSE             FALSE        9
#> 36    1     0       11   35   88          FALSE              TRUE       11
#> 37    1     0       11   37   38          FALSE              TRUE       11
#> 38    1     0       11   37   62          FALSE              TRUE       11
#> 39    1     0       10   37   93          FALSE             FALSE       10
#> 40    1     0       11   38   41          FALSE              TRUE       11
#> 41    1     0       11   38   49          FALSE              TRUE       11
#> 42    1     0        2   38   62          FALSE             FALSE        2
#> 43    1     0       11   41   99          FALSE              TRUE       11
#> 44    1     0        6   42   78          FALSE             FALSE        6
#> 45    1     0       11   43   93          FALSE              TRUE       11
#> 46    1     0       11   44   69          FALSE              TRUE       11
#> 47    1     0       11   49   88          FALSE              TRUE       11
#> 48    1     0       11   52   78          FALSE              TRUE       11
#> 49    1     0       11   57   83          FALSE              TRUE       11
#> 50    1     0       10   68   97          FALSE             FALSE       10
#> 51    1     0        7   74   83          FALSE             FALSE        7
#> 52    1     0       11   74  100          FALSE              TRUE       11
#> 53    1     0        1   89   93          FALSE             FALSE        1
#> 54    1     0       11   90   91          FALSE              TRUE       11
#> 55    1     0        2   91   99          FALSE             FALSE        2
#> 56    1     0        5   93   97          FALSE             FALSE        5
#> 57    1     1        6   41   93          FALSE             FALSE        5
#> 58    1     1       10    9   86          FALSE             FALSE        9
#> 59    1     2       11   38   99          FALSE              TRUE        9
#> 60    1     2       11   39   99          FALSE              TRUE        9
#> 61    1     2       10   29   36          FALSE             FALSE        8
#> 62    1     3       11   38   94          FALSE              TRUE        8
#> 63    1     4       11   21   82          FALSE              TRUE        7
#> 64    1     4       11   75   90          FALSE              TRUE        7
#> 65    1     4       11   42   89          FALSE              TRUE        7
#> 66    1     4       11   12   28          FALSE              TRUE        7
#> 67    1     5        8   60   73          FALSE             FALSE        3
#> 68    1     6       11   12   83          FALSE              TRUE        5
#> 69    1     6       11   26   81          FALSE              TRUE        5
#> 70    1     7        8   58   62          FALSE             FALSE        1
#> 71    1     7       11    1   61          FALSE              TRUE        4
#> 72    1     7       11   32   69          FALSE              TRUE        4
#> 73    1     7       11   57   97          FALSE              TRUE        4
#> 74    1     7       11   54   75          FALSE              TRUE        4
#> 75    1     8       11   74   96          FALSE              TRUE        3
#> 76    1     8       11   33   97          FALSE              TRUE        3
#> 77    1     8       11   10   21          FALSE              TRUE        3
#> 78    1     8       11    6   79          FALSE              TRUE        3
#> 79    1     9       11   88   94          FALSE              TRUE        2
#> 80    1    10       11    5   92          FALSE              TRUE        1
#> 81    2     0       11    1   10          FALSE              TRUE       11
#> 82    2     0        2    2   97          FALSE             FALSE        2
#> 83    2     0        5    4   18          FALSE             FALSE        5
#> 84    2     0       11    4   41          FALSE              TRUE       11
#> 85    2     0        2    5   36          FALSE             FALSE        2
#> 86    2     0       11    5   45          FALSE              TRUE       11
#> 87    2     0        3    6   13          FALSE             FALSE        3
#> 88    2     0       11    6   63          FALSE              TRUE       11
#> 89    2     0        3    7   45          FALSE             FALSE        3
#> 90    2     0        4    8   29          FALSE             FALSE        4
#> 91    2     0        2    9   38          FALSE             FALSE        2
#> 92    2     0       11   11   75          FALSE              TRUE       11
#> 93    2     0        3   13   41          FALSE             FALSE        3
#> 94    2     0       11   16   88          FALSE              TRUE       11
#> 95    2     0       11   16   91          FALSE              TRUE       11
#> 96    2     0       11   22   55          FALSE              TRUE       11
#> 97    2     0       11   23   49          FALSE              TRUE       11
#> 98    2     0        2   25   51          FALSE             FALSE        2
#> 99    2     0        8   27   57          FALSE             FALSE        8
#> 100   2     0       11   28   77          FALSE              TRUE       11
#> 101   2     0        6   32   57          FALSE             FALSE        6
#> 102   2     0        7   34   52          FALSE             FALSE        7
#> 103   2     0       11   36   55          FALSE              TRUE       11
#> 104   2     0       11   37   87          FALSE              TRUE       11
#> 105   2     0       11   38   52          FALSE              TRUE       11
#> 106   2     0       11   39   87          FALSE              TRUE       11
#> 107   2     0       11   44   49          FALSE              TRUE       11
#> 108   2     0       11   47   93          FALSE              TRUE       11
#> 109   2     0        3   51   99          FALSE             FALSE        3
#> 110   2     0       11   52   69          FALSE              TRUE       11
#> 111   2     0       11   54   97          FALSE              TRUE       11
#> 112   2     0       11   63   96          FALSE              TRUE       11
#> 113   2     0       11   64   75          FALSE              TRUE       11
#> 114   2     0        2   67   68          FALSE             FALSE        2
#> 115   2     0       11   73   89          FALSE              TRUE       11
#> 116   2     0       11   73  100          FALSE              TRUE       11
#> 117   2     0       11   74   91          FALSE              TRUE       11
#> 118   2     0        6   76   77          FALSE             FALSE        6
#> 119   2     0        3   79   99          FALSE             FALSE        3
#> 120   2     0       10   81  100          FALSE             FALSE       10
#> 121   2     0       11   83  100          FALSE              TRUE       11
#> 122   2     0       11   87   98          FALSE              TRUE       11
#> 123   2     0       11   92   96          FALSE              TRUE       11
#> 124   2     0        3   92  100          FALSE             FALSE        3
#> 125   2     0       11   97   98          FALSE              TRUE       11
#> 126   2     1       11    1  100          FALSE              TRUE       10
#> 127   2     1       11   62   93          FALSE              TRUE       10
#> 128   2     1       11   58   62          FALSE              TRUE       10
#> 129   2     2       10   61   97          FALSE             FALSE        8
#> 130   2     2        7   52   65          FALSE             FALSE        5
#> 131   2     2       11   48   95          FALSE              TRUE        9
#> 132   2     3       11    8   53          FALSE              TRUE        8
#> 133   2     3       11   63   99          FALSE              TRUE        8
#> 134   2     3        7   54   61          FALSE             FALSE        4
#> 135   2     3        5    1   31          FALSE             FALSE        2
#> 136   2     4       11   15   17          FALSE              TRUE        7
#> 137   2     4        9   57   99          FALSE             FALSE        5
#> 138   2     5       11   24   53          FALSE              TRUE        6
#> 139   2     5       11   19   72          FALSE              TRUE        6
#> 140   2     6        7   24   85          FALSE             FALSE        1
#> 141   2     6       11   65   97          FALSE              TRUE        5
#> 142   2     6       11   44   65          FALSE              TRUE        5
#> 143   2     6       11    3    5          FALSE              TRUE        5
#> 144   2     7       11   26   66          FALSE              TRUE        4
#> 145   2     8       11   16   83          FALSE              TRUE        3
#> 146   2     8       11   22   73          FALSE              TRUE        3
#> 147   2     8       11   14   79          FALSE              TRUE        3
#> 148   2     9       10   27   50          FALSE             FALSE        1
#> 149   2    10       11   55   71          FALSE              TRUE        1
#> 150   3     0        8    4   12          FALSE             FALSE        8
#> 151   3     0       10    5   51          FALSE             FALSE       10
#> 152   3     0       11    6   49          FALSE              TRUE       11
#> 153   3     0        8    7   94          FALSE             FALSE        8
#> 154   3     0       11    8   18          FALSE              TRUE       11
#> 155   3     0        2    8   39          FALSE             FALSE        2
#> 156   3     0       11    9   86          FALSE              TRUE       11
#> 157   3     0       11   10   71          FALSE              TRUE       11
#> 158   3     0        8   10   78          FALSE             FALSE        8
#> 159   3     0        5   11   55          FALSE             FALSE        5
#> 160   3     0       11   11   88          FALSE              TRUE       11
#> 161   3     0        1   13   66          FALSE             FALSE        1
#> 162   3     0       11   14   41          FALSE              TRUE       11
#> 163   3     0       11   15   38          FALSE              TRUE       11
#> 164   3     0        8   16   99          FALSE             FALSE        8
#> 165   3     0        8   17   80          FALSE             FALSE        8
#> 166   3     0       11   18   79          FALSE              TRUE       11
#> 167   3     0       11   20   32          FALSE              TRUE       11
#> 168   3     0       10   21   92          FALSE             FALSE       10
#> 169   3     0        1   22   27          FALSE             FALSE        1
#> 170   3     0        9   22   37          FALSE             FALSE        9
#> 171   3     0       11   23   33          FALSE              TRUE       11
#> 172   3     0        1   24   43          FALSE             FALSE        1
#> 173   3     0        5   24   56          FALSE             FALSE        5
#> 174   3     0       11   25   42          FALSE              TRUE       11
#> 175   3     0       11   25   78          FALSE              TRUE       11
#> 176   3     0       11   27   59          FALSE              TRUE       11
#> 177   3     0       11   27   67          FALSE              TRUE       11
#> 178   3     0       11   32   73          FALSE              TRUE       11
#> 179   3     0       11   33   69          FALSE              TRUE       11
#> 180   3     0        7   34   85          FALSE             FALSE        7
#> 181   3     0       11   36   56          FALSE              TRUE       11
#> 182   3     0        1   38   39          FALSE             FALSE        1
#> 183   3     0        5   44   54          FALSE             FALSE        5
#> 184   3     0        2   45   48          FALSE             FALSE        2
#> 185   3     0       11   45   69          FALSE              TRUE       11
#> 186   3     0       10   47   74          FALSE             FALSE       10
#> 187   3     0       11   49   58          FALSE              TRUE       11
#> 188   3     0       11   49   88          FALSE              TRUE       11
#> 189   3     0       11   50   81          FALSE              TRUE       11
#> 190   3     0       11   53   80          FALSE              TRUE       11
#> 191   3     0        8   57   64          FALSE             FALSE        8
#> 192   3     0       11   59   75          FALSE              TRUE       11
#> 193   3     0       11   63   83          FALSE              TRUE       11
#> 194   3     0        4   65   68          FALSE             FALSE        4
#> 195   3     0        6   69   75          FALSE             FALSE        6
#> 196   3     0       11   70   88          FALSE              TRUE       11
#> 197   3     0       11   76   86          FALSE              TRUE       11
#> 198   3     0       11   83   95          FALSE              TRUE       11
#> 199   3     0       11   90   99          FALSE              TRUE       11
#> 200   3     1       10    3   16          FALSE             FALSE        9
#> 201   3     2       11   59   67          FALSE              TRUE        9
#> 202   3     2        7   38   97          FALSE             FALSE        5
#> 203   3     2       11   39   86          FALSE              TRUE        9
#> 204   3     2       11   70   91          FALSE              TRUE        9
#> 205   3     4       11   48   71          FALSE              TRUE        7
#> 206   3     4        7    7   75          FALSE             FALSE        3
#> 207   3     5       11   27   43          FALSE              TRUE        6
#> 208   3     5       11    5   55          FALSE              TRUE        6
#> 209   3     5       11    6   44          FALSE              TRUE        6
#> 210   3     5       11    2   16          FALSE              TRUE        6
#> 211   3     5       11    6   16          FALSE              TRUE        6
#> 212   3     6       11   23   91          FALSE              TRUE        5
#> 213   3     6       11   26   41          FALSE              TRUE        5
#> 214   3     6       11   45   81          FALSE              TRUE        5
#> 215   3     7        9    6   77          FALSE             FALSE        2
#> 216   3     7        9   69   87          FALSE             FALSE        2
#> 217   3     7       11   39   65          FALSE              TRUE        4
#> 218   3     7       11    9   27          FALSE              TRUE        4
#> 219   3     8       11   14   18          FALSE              TRUE        3
#> 220   3     8       11    2   25          FALSE              TRUE        3
#> 221   3    10       11   23   49          FALSE              TRUE        1
#> 222   3    10       11   44   87          FALSE              TRUE        1
#>     edge.id
#> 1         1
#> 2         2
#> 3         3
#> 4         4
#> 5         5
#> 6         6
#> 7         7
#> 8         8
#> 9         9
#> 10       10
#> 11       11
#> 12       12
#> 13       13
#> 14       14
#> 15       15
#> 16       16
#> 17       17
#> 18       18
#> 19       19
#> 20       20
#> 21       21
#> 22       22
#> 23       23
#> 24       24
#> 25       25
#> 26       26
#> 27       27
#> 28       28
#> 29       29
#> 30       30
#> 31       31
#> 32       32
#> 33       33
#> 34       34
#> 35       35
#> 36       36
#> 37       37
#> 38       38
#> 39       39
#> 40       40
#> 41       41
#> 42       42
#> 43       43
#> 44       44
#> 45       45
#> 46       46
#> 47       47
#> 48       48
#> 49       49
#> 50       50
#> 51       51
#> 52       52
#> 53       53
#> 54       54
#> 55       55
#> 56       56
#> 57       57
#> 58       58
#> 59       59
#> 60       60
#> 61       61
#> 62       62
#> 63       63
#> 64       64
#> 65       65
#> 66       66
#> 67       67
#> 68       68
#> 69       69
#> 70       70
#> 71       71
#> 72       72
#> 73       73
#> 74       74
#> 75       75
#> 76       76
#> 77       77
#> 78       78
#> 79       79
#> 80       80
#> 81        1
#> 82        2
#> 83        3
#> 84        4
#> 85        5
#> 86        6
#> 87        7
#> 88        8
#> 89        9
#> 90       10
#> 91       11
#> 92       12
#> 93       13
#> 94       14
#> 95       15
#> 96       16
#> 97       17
#> 98       18
#> 99       19
#> 100      20
#> 101      21
#> 102      22
#> 103      23
#> 104      24
#> 105      25
#> 106      26
#> 107      27
#> 108      28
#> 109      29
#> 110      30
#> 111      31
#> 112      32
#> 113      33
#> 114      34
#> 115      35
#> 116      36
#> 117      37
#> 118      38
#> 119      39
#> 120      40
#> 121      41
#> 122      42
#> 123      43
#> 124      44
#> 125      45
#> 126      46
#> 127      47
#> 128      48
#> 129      49
#> 130      50
#> 131      51
#> 132      52
#> 133      53
#> 134      54
#> 135      55
#> 136      56
#> 137      57
#> 138      58
#> 139      59
#> 140      60
#> 141      61
#> 142      62
#> 143      63
#> 144      64
#> 145      65
#> 146      66
#> 147      67
#> 148      68
#> 149      69
#> 150       1
#> 151       2
#> 152       3
#> 153       4
#> 154       5
#> 155       6
#> 156       7
#> 157       8
#> 158       9
#> 159      10
#> 160      11
#> 161      12
#> 162      13
#> 163      14
#> 164      15
#> 165      16
#> 166      17
#> 167      18
#> 168      19
#> 169      20
#> 170      21
#> 171      22
#> 172      23
#> 173      24
#> 174      25
#> 175      26
#> 176      27
#> 177      28
#> 178      29
#> 179      30
#> 180      31
#> 181      32
#> 182      33
#> 183      34
#> 184      35
#> 185      36
#> 186      37
#> 187      38
#> 188      39
#> 189      40
#> 190      41
#> 191      42
#> 192      43
#> 193      44
#> 194      45
#> 195      46
#> 196      47
#> 197      48
#> 198      49
#> 199      50
#> 200      51
#> 201      52
#> 202      53
#> 203      54
#> 204      55
#> 205      56
#> 206      57
#> 207      58
#> 208      59
#> 209      60
#> 210      61
#> 211      62
#> 212      63
#> 213      64
#> 214      65
#> 215      66
#> 216      67
#> 217      68
#> 218      69
#> 219      70
#> 220      71
#> 221      72
#> 222      73
# }
```
