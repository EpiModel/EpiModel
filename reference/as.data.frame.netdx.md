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
#> 1      0       11    2   34          FALSE              TRUE       11       1
#> 2      0        7    3   23          FALSE             FALSE        7       2
#> 3      0       11    5   18          FALSE              TRUE       11       3
#> 4      0        8    5   52          FALSE             FALSE        8       4
#> 5      0       11    5   67          FALSE              TRUE       11       5
#> 6      0       11    6   47          FALSE              TRUE       11       6
#> 7      0        5    6   64          FALSE             FALSE        5       7
#> 8      0       11    6   67          FALSE              TRUE       11       8
#> 9      0       11    7   62          FALSE              TRUE       11       9
#> 10     0       11    7   96          FALSE              TRUE       11      10
#> 11     0        5    8   64          FALSE             FALSE        5      11
#> 12     0        4   10   38          FALSE             FALSE        4      12
#> 13     0       11   12   58          FALSE              TRUE       11      13
#> 14     0       10   13   58          FALSE             FALSE       10      14
#> 15     0        9   15   36          FALSE             FALSE        9      15
#> 16     0       11   16   39          FALSE              TRUE       11      16
#> 17     0        1   16   68          FALSE             FALSE        1      17
#> 18     0        3   16   97          FALSE             FALSE        3      18
#> 19     0       11   17   83          FALSE              TRUE       11      19
#> 20     0       11   18   27          FALSE              TRUE       11      20
#> 21     0        7   19   69          FALSE             FALSE        7      21
#> 22     0       11   19   92          FALSE              TRUE       11      22
#> 23     0       11   20   59          FALSE              TRUE       11      23
#> 24     0        4   21   67          FALSE             FALSE        4      24
#> 25     0        1   23   58          FALSE             FALSE        1      25
#> 26     0       11   26   65          FALSE              TRUE       11      26
#> 27     0       11   27   34          FALSE              TRUE       11      27
#> 28     0       11   27   72          FALSE              TRUE       11      28
#> 29     0       11   28   70          FALSE              TRUE       11      29
#> 30     0       11   30   36          FALSE              TRUE       11      30
#> 31     0       11   30   41          FALSE              TRUE       11      31
#> 32     0       11   31   33          FALSE              TRUE       11      32
#> 33     0       11   31   51          FALSE              TRUE       11      33
#> 34     0        8   33   81          FALSE             FALSE        8      34
#> 35     0       11   34   38          FALSE              TRUE       11      35
#> 36     0       10   35   67          FALSE             FALSE       10      36
#> 37     0       11   36   55          FALSE              TRUE       11      37
#> 38     0       11   38   96          FALSE              TRUE       11      38
#> 39     0       11   41   91          FALSE              TRUE       11      39
#> 40     0        7   42   80          FALSE             FALSE        7      40
#> 41     0        7   45   79          FALSE             FALSE        7      41
#> 42     0        1   46   97          FALSE             FALSE        1      42
#> 43     0        8   51   74          FALSE             FALSE        8      43
#> 44     0       11   62   67          FALSE              TRUE       11      44
#> 45     0        4   63   66          FALSE             FALSE        4      45
#> 46     0       11   64   94          FALSE              TRUE       11      46
#> 47     0        2   65   75          FALSE             FALSE        2      47
#> 48     0       11   65   89          FALSE              TRUE       11      48
#> 49     0       11   67   79          FALSE              TRUE       11      49
#> 50     0       11   68   81          FALSE              TRUE       11      50
#> 51     0        3   70  100          FALSE             FALSE        3      51
#> 52     0       11   72   81          FALSE              TRUE       11      52
#> 53     0       10   75   99          FALSE             FALSE       10      53
#> 54     0       11   78   96          FALSE              TRUE       11      54
#> 55     0       10   81   84          FALSE             FALSE       10      55
#> 56     0       11   89   98          FALSE              TRUE       11      56
#> 57     1        6   26   55          FALSE             FALSE        5      57
#> 58     1       11   25   47          FALSE              TRUE       10      58
#> 59     1        5   60   92          FALSE             FALSE        4      59
#> 60     1        2    5   82          FALSE             FALSE        1      60
#> 61     2        8   46   63          FALSE             FALSE        6      61
#> 62     2       11   19   70          FALSE              TRUE        9      62
#> 63     2       11   66   89          FALSE              TRUE        9      63
#> 64     2        9   50   84          FALSE             FALSE        7      64
#> 65     3       11   16   34          FALSE              TRUE        8      65
#> 66     3       11   38   65          FALSE              TRUE        8      66
#> 67     3       11   54   91          FALSE              TRUE        8      67
#> 68     3        7   53   74          FALSE             FALSE        4      68
#> 69     3       11   69   89          FALSE              TRUE        8      69
#> 70     4       11   38   77          FALSE              TRUE        7      70
#> 71     4        7   24   32          FALSE             FALSE        3      71
#> 72     4       11   45   92          FALSE              TRUE        7      72
#> 73     5       11   82   98          FALSE              TRUE        6      73
#> 74     6       11   44   64          FALSE              TRUE        5      74
#> 75     7       11   40   92          FALSE              TRUE        4      75
#> 76     7       11   26   33          FALSE              TRUE        4      76
#> 77     7       11    7   65          FALSE              TRUE        4      77
#> 78     8       11   43   50          FALSE              TRUE        3      78
#> 79     8       11   73   92          FALSE              TRUE        3      79
#> 80     9       11   44   48          FALSE              TRUE        2      80
#> 81     9       11   76   91          FALSE              TRUE        2      81
#> 82     9       11   52   83          FALSE              TRUE        2      82
#> 83     9       11   45   50          FALSE              TRUE        2      83
#> 84     9       11   30   76          FALSE              TRUE        2      84
#> 85    10       11   13   66          FALSE              TRUE        1      85

# Extract data from all simulations
as.data.frame(dx)
#>     sim onset terminus tail head onset.censored terminus.censored duration
#> 1     1     0       11    2   34          FALSE              TRUE       11
#> 2     1     0        7    3   23          FALSE             FALSE        7
#> 3     1     0       11    5   18          FALSE              TRUE       11
#> 4     1     0        8    5   52          FALSE             FALSE        8
#> 5     1     0       11    5   67          FALSE              TRUE       11
#> 6     1     0       11    6   47          FALSE              TRUE       11
#> 7     1     0        5    6   64          FALSE             FALSE        5
#> 8     1     0       11    6   67          FALSE              TRUE       11
#> 9     1     0       11    7   62          FALSE              TRUE       11
#> 10    1     0       11    7   96          FALSE              TRUE       11
#> 11    1     0        5    8   64          FALSE             FALSE        5
#> 12    1     0        4   10   38          FALSE             FALSE        4
#> 13    1     0       11   12   58          FALSE              TRUE       11
#> 14    1     0       10   13   58          FALSE             FALSE       10
#> 15    1     0        9   15   36          FALSE             FALSE        9
#> 16    1     0       11   16   39          FALSE              TRUE       11
#> 17    1     0        1   16   68          FALSE             FALSE        1
#> 18    1     0        3   16   97          FALSE             FALSE        3
#> 19    1     0       11   17   83          FALSE              TRUE       11
#> 20    1     0       11   18   27          FALSE              TRUE       11
#> 21    1     0        7   19   69          FALSE             FALSE        7
#> 22    1     0       11   19   92          FALSE              TRUE       11
#> 23    1     0       11   20   59          FALSE              TRUE       11
#> 24    1     0        4   21   67          FALSE             FALSE        4
#> 25    1     0        1   23   58          FALSE             FALSE        1
#> 26    1     0       11   26   65          FALSE              TRUE       11
#> 27    1     0       11   27   34          FALSE              TRUE       11
#> 28    1     0       11   27   72          FALSE              TRUE       11
#> 29    1     0       11   28   70          FALSE              TRUE       11
#> 30    1     0       11   30   36          FALSE              TRUE       11
#> 31    1     0       11   30   41          FALSE              TRUE       11
#> 32    1     0       11   31   33          FALSE              TRUE       11
#> 33    1     0       11   31   51          FALSE              TRUE       11
#> 34    1     0        8   33   81          FALSE             FALSE        8
#> 35    1     0       11   34   38          FALSE              TRUE       11
#> 36    1     0       10   35   67          FALSE             FALSE       10
#> 37    1     0       11   36   55          FALSE              TRUE       11
#> 38    1     0       11   38   96          FALSE              TRUE       11
#> 39    1     0       11   41   91          FALSE              TRUE       11
#> 40    1     0        7   42   80          FALSE             FALSE        7
#> 41    1     0        7   45   79          FALSE             FALSE        7
#> 42    1     0        1   46   97          FALSE             FALSE        1
#> 43    1     0        8   51   74          FALSE             FALSE        8
#> 44    1     0       11   62   67          FALSE              TRUE       11
#> 45    1     0        4   63   66          FALSE             FALSE        4
#> 46    1     0       11   64   94          FALSE              TRUE       11
#> 47    1     0        2   65   75          FALSE             FALSE        2
#> 48    1     0       11   65   89          FALSE              TRUE       11
#> 49    1     0       11   67   79          FALSE              TRUE       11
#> 50    1     0       11   68   81          FALSE              TRUE       11
#> 51    1     0        3   70  100          FALSE             FALSE        3
#> 52    1     0       11   72   81          FALSE              TRUE       11
#> 53    1     0       10   75   99          FALSE             FALSE       10
#> 54    1     0       11   78   96          FALSE              TRUE       11
#> 55    1     0       10   81   84          FALSE             FALSE       10
#> 56    1     0       11   89   98          FALSE              TRUE       11
#> 57    1     1        6   26   55          FALSE             FALSE        5
#> 58    1     1       11   25   47          FALSE              TRUE       10
#> 59    1     1        5   60   92          FALSE             FALSE        4
#> 60    1     1        2    5   82          FALSE             FALSE        1
#> 61    1     2        8   46   63          FALSE             FALSE        6
#> 62    1     2       11   19   70          FALSE              TRUE        9
#> 63    1     2       11   66   89          FALSE              TRUE        9
#> 64    1     2        9   50   84          FALSE             FALSE        7
#> 65    1     3       11   16   34          FALSE              TRUE        8
#> 66    1     3       11   38   65          FALSE              TRUE        8
#> 67    1     3       11   54   91          FALSE              TRUE        8
#> 68    1     3        7   53   74          FALSE             FALSE        4
#> 69    1     3       11   69   89          FALSE              TRUE        8
#> 70    1     4       11   38   77          FALSE              TRUE        7
#> 71    1     4        7   24   32          FALSE             FALSE        3
#> 72    1     4       11   45   92          FALSE              TRUE        7
#> 73    1     5       11   82   98          FALSE              TRUE        6
#> 74    1     6       11   44   64          FALSE              TRUE        5
#> 75    1     7       11   40   92          FALSE              TRUE        4
#> 76    1     7       11   26   33          FALSE              TRUE        4
#> 77    1     7       11    7   65          FALSE              TRUE        4
#> 78    1     8       11   43   50          FALSE              TRUE        3
#> 79    1     8       11   73   92          FALSE              TRUE        3
#> 80    1     9       11   44   48          FALSE              TRUE        2
#> 81    1     9       11   76   91          FALSE              TRUE        2
#> 82    1     9       11   52   83          FALSE              TRUE        2
#> 83    1     9       11   45   50          FALSE              TRUE        2
#> 84    1     9       11   30   76          FALSE              TRUE        2
#> 85    1    10       11   13   66          FALSE              TRUE        1
#> 86    2     0       11    1   86          FALSE              TRUE       11
#> 87    2     0       11    3   62          FALSE              TRUE       11
#> 88    2     0       11    4   65          FALSE              TRUE       11
#> 89    2     0       11    5   46          FALSE              TRUE       11
#> 90    2     0       11    6   97          FALSE              TRUE       11
#> 91    2     0       11    9   52          FALSE              TRUE       11
#> 92    2     0        8   10   24          FALSE             FALSE        8
#> 93    2     0       11   10   60          FALSE              TRUE       11
#> 94    2     0       11   11   69          FALSE              TRUE       11
#> 95    2     0       11   14   21          FALSE              TRUE       11
#> 96    2     0       11   16   73          FALSE              TRUE       11
#> 97    2     0        4   21   23          FALSE             FALSE        4
#> 98    2     0       11   22   65          FALSE              TRUE       11
#> 99    2     0        9   23   89          FALSE             FALSE        9
#> 100   2     0        2   24   86          FALSE             FALSE        2
#> 101   2     0        3   27   34          FALSE             FALSE        3
#> 102   2     0        5   27   81          FALSE             FALSE        5
#> 103   2     0       11   29   88          FALSE              TRUE       11
#> 104   2     0       11   30   68          FALSE              TRUE       11
#> 105   2     0       11   32   71          FALSE              TRUE       11
#> 106   2     0        1   33   55          FALSE             FALSE        1
#> 107   2     0        2   34   45          FALSE             FALSE        2
#> 108   2     0       11   35   57          FALSE              TRUE       11
#> 109   2     0        8   37   43          FALSE             FALSE        8
#> 110   2     0        9   37   49          FALSE             FALSE        9
#> 111   2     0        5   38   95          FALSE             FALSE        5
#> 112   2     0        5   40   75          FALSE             FALSE        5
#> 113   2     0       11   43   52          FALSE              TRUE       11
#> 114   2     0       11   45   58          FALSE              TRUE       11
#> 115   2     0       11   45  100          FALSE              TRUE       11
#> 116   2     0        6   48   53          FALSE             FALSE        6
#> 117   2     0        6   48   82          FALSE             FALSE        6
#> 118   2     0       11   50   82          FALSE              TRUE       11
#> 119   2     0        9   51   69          FALSE             FALSE        9
#> 120   2     0       11   52   64          FALSE              TRUE       11
#> 121   2     0       11   55   83          FALSE              TRUE       11
#> 122   2     0        1   56   89          FALSE             FALSE        1
#> 123   2     0       11   63   85          FALSE              TRUE       11
#> 124   2     0        2   68   70          FALSE             FALSE        2
#> 125   2     0        7   74   94          FALSE             FALSE        7
#> 126   2     0       11   78   85          FALSE              TRUE       11
#> 127   2     0       11   78   96          FALSE              TRUE       11
#> 128   2     0        9   93   97          FALSE             FALSE        9
#> 129   2     1       11   29   79          FALSE              TRUE       10
#> 130   2     1       11   42   59          FALSE              TRUE       10
#> 131   2     1       11   18   81          FALSE              TRUE       10
#> 132   2     1       11   17   90          FALSE              TRUE       10
#> 133   2     2        8   41   83          FALSE             FALSE        6
#> 134   2     2       11   23   88          FALSE              TRUE        9
#> 135   2     2       11   15   73          FALSE              TRUE        9
#> 136   2     3       11    4   29          FALSE              TRUE        8
#> 137   2     3       11   39   74          FALSE              TRUE        8
#> 138   2     4       11   14   41          FALSE              TRUE        7
#> 139   2     4       11   17   56          FALSE              TRUE        7
#> 140   2     4        7    4   35          FALSE             FALSE        3
#> 141   2     5       11   51   63          FALSE              TRUE        6
#> 142   2     6       11   28   63          FALSE              TRUE        5
#> 143   2     7       11   57   87          FALSE              TRUE        4
#> 144   2     7       11   65   90          FALSE              TRUE        4
#> 145   2     7        9    6   19          FALSE             FALSE        2
#> 146   2     8       11    1   79          FALSE              TRUE        3
#> 147   2     8       11   30   50          FALSE              TRUE        3
#> 148   2     8       11    9   89          FALSE              TRUE        3
#> 149   2     9       11   32   59          FALSE              TRUE        2
#> 150   2    10       11   27   95          FALSE              TRUE        1
#> 151   2    10       11   97   98          FALSE              TRUE        1
#> 152   2    10       11   45   51          FALSE              TRUE        1
#> 153   2    10       11   65   75          FALSE              TRUE        1
#> 154   2    10       11    1   78          FALSE              TRUE        1
#> 155   3     0       11    1   20          FALSE              TRUE       11
#> 156   3     0       11    1   58          FALSE              TRUE       11
#> 157   3     0        9    2   12          FALSE             FALSE        9
#> 158   3     0        3    2   72          FALSE             FALSE        3
#> 159   3     0       10    3   53          FALSE             FALSE       10
#> 160   3     0        6    6   26          FALSE             FALSE        6
#> 161   3     0        5    7   22          FALSE             FALSE        5
#> 162   3     0       11   10   15          FALSE              TRUE       11
#> 163   3     0       11   13   57          FALSE              TRUE       11
#> 164   3     0        3   15   36          FALSE             FALSE        3
#> 165   3     0       11   16   28          FALSE              TRUE       11
#> 166   3     0        3   16   79          FALSE             FALSE        3
#> 167   3     0       11   17   35          FALSE              TRUE       11
#> 168   3     0        8   19   47          FALSE             FALSE        8
#> 169   3     0       11   19   59          FALSE              TRUE       11
#> 170   3     0       11   20   46          FALSE              TRUE       11
#> 171   3     0        6   21   38          FALSE             FALSE        6
#> 172   3     0        2   22   68          FALSE             FALSE        2
#> 173   3     0        2   24   67          FALSE             FALSE        2
#> 174   3     0       11   25   26          FALSE              TRUE       11
#> 175   3     0        6   25   29          FALSE             FALSE        6
#> 176   3     0       11   25   99          FALSE              TRUE       11
#> 177   3     0        6   26   67          FALSE             FALSE        6
#> 178   3     0       11   26   79          FALSE              TRUE       11
#> 179   3     0        9   27   50          FALSE             FALSE        9
#> 180   3     0        3   27   58          FALSE             FALSE        3
#> 181   3     0       11   29   32          FALSE              TRUE       11
#> 182   3     0        1   30   55          FALSE             FALSE        1
#> 183   3     0       11   35   67          FALSE              TRUE       11
#> 184   3     0       11   39   62          FALSE              TRUE       11
#> 185   3     0       10   44   69          FALSE             FALSE       10
#> 186   3     0       11   47   75          FALSE              TRUE       11
#> 187   3     0        1   47   97          FALSE             FALSE        1
#> 188   3     0        5   49   77          FALSE             FALSE        5
#> 189   3     0       11   49   92          FALSE              TRUE       11
#> 190   3     0       11   50   78          FALSE              TRUE       11
#> 191   3     0        5   53   54          FALSE             FALSE        5
#> 192   3     0        3   57   71          FALSE             FALSE        3
#> 193   3     0       11   58   68          FALSE              TRUE       11
#> 194   3     0       11   63   75          FALSE              TRUE       11
#> 195   3     0        2   69   96          FALSE             FALSE        2
#> 196   3     0        3   71   84          FALSE             FALSE        3
#> 197   3     0       11   86   90          FALSE              TRUE       11
#> 198   3     0        5   88   99          FALSE             FALSE        5
#> 199   3     1       11   36   81          FALSE              TRUE       10
#> 200   3     1       11   20   80          FALSE              TRUE       10
#> 201   3     2       11   79   98          FALSE              TRUE        9
#> 202   3     2       11   73   83          FALSE              TRUE        9
#> 203   3     2       11   21   43          FALSE              TRUE        9
#> 204   3     3       11   40   66          FALSE              TRUE        8
#> 205   3     3       11   21   34          FALSE              TRUE        8
#> 206   3     3       11   51   68          FALSE              TRUE        8
#> 207   3     5        6   52   60          FALSE             FALSE        1
#> 208   3     5       11   31   81          FALSE              TRUE        6
#> 209   3     5       11   23   80          FALSE              TRUE        6
#> 210   3     5       11   33   67          FALSE              TRUE        6
#> 211   3     6       11   53   55          FALSE              TRUE        5
#> 212   3     7       11   28   84          FALSE              TRUE        4
#> 213   3     8       11   25   61          FALSE              TRUE        3
#> 214   3     8       11   54   96          FALSE              TRUE        3
#> 215   3     8       11   78  100          FALSE              TRUE        3
#> 216   3     8       11   39   47          FALSE              TRUE        3
#> 217   3     9       11   62   75          FALSE              TRUE        2
#> 218   3    10       11   32   79          FALSE              TRUE        1
#> 219   3    10       11   45   93          FALSE              TRUE        1
#> 220   3    10       11   52   66          FALSE              TRUE        1
#> 221   3    10       11   24   70          FALSE              TRUE        1
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
#> 81       81
#> 82       82
#> 83       83
#> 84       84
#> 85       85
#> 86        1
#> 87        2
#> 88        3
#> 89        4
#> 90        5
#> 91        6
#> 92        7
#> 93        8
#> 94        9
#> 95       10
#> 96       11
#> 97       12
#> 98       13
#> 99       14
#> 100      15
#> 101      16
#> 102      17
#> 103      18
#> 104      19
#> 105      20
#> 106      21
#> 107      22
#> 108      23
#> 109      24
#> 110      25
#> 111      26
#> 112      27
#> 113      28
#> 114      29
#> 115      30
#> 116      31
#> 117      32
#> 118      33
#> 119      34
#> 120      35
#> 121      36
#> 122      37
#> 123      38
#> 124      39
#> 125      40
#> 126      41
#> 127      42
#> 128      43
#> 129      44
#> 130      45
#> 131      46
#> 132      47
#> 133      48
#> 134      49
#> 135      50
#> 136      51
#> 137      52
#> 138      53
#> 139      54
#> 140      55
#> 141      56
#> 142      57
#> 143      58
#> 144      59
#> 145      60
#> 146      61
#> 147      62
#> 148      63
#> 149      64
#> 150      65
#> 151      66
#> 152      67
#> 153      68
#> 154      69
#> 155       1
#> 156       2
#> 157       3
#> 158       4
#> 159       5
#> 160       6
#> 161       7
#> 162       8
#> 163       9
#> 164      10
#> 165      11
#> 166      12
#> 167      13
#> 168      14
#> 169      15
#> 170      16
#> 171      17
#> 172      18
#> 173      19
#> 174      20
#> 175      21
#> 176      22
#> 177      23
#> 178      24
#> 179      25
#> 180      26
#> 181      27
#> 182      28
#> 183      29
#> 184      30
#> 185      31
#> 186      32
#> 187      33
#> 188      34
#> 189      35
#> 190      36
#> 191      37
#> 192      38
#> 193      39
#> 194      40
#> 195      41
#> 196      42
#> 197      43
#> 198      44
#> 199      45
#> 200      46
#> 201      47
#> 202      48
#> 203      49
#> 204      50
#> 205      51
#> 206      52
#> 207      53
#> 208      54
#> 209      55
#> 210      56
#> 211      57
#> 212      58
#> 213      59
#> 214      60
#> 215      61
#> 216      62
#> 217      63
#> 218      64
#> 219      65
#> 220      66
#> 221      67
# }
```
