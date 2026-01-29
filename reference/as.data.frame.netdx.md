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
# Initialize and parameterize the network model
nw <- network_initialize(n = 100)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)

# Model estimation
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#> Starting maximum pseudolikelihood estimation (MPLE):
#> Obtaining the responsible dyads.
#> Evaluating the predictor and response matrix.
#> Maximizing the pseudolikelihood.
#> Finished MPLE.

# Simulate the network with netdx
dx <- netdx(est, nsims = 3, nsteps = 10, keep.tedgelist = TRUE,
            verbose = FALSE)
#> Warning: NAs introduced by coercion to integer range
#> Warning: NAs introduced by coercion to integer range
#> Warning: NAs introduced by coercion to integer range

# Extract data from the first simulation
as.data.frame(dx, sim = 1)
#>    onset terminus tail head onset.censored terminus.censored duration edge.id
#> 1      0        6    1   46          FALSE             FALSE        6       1
#> 2      0        7    1   86          FALSE             FALSE        7       2
#> 3      0       10    2   97          FALSE             FALSE       10       3
#> 4      0       11    3    8          FALSE              TRUE       11       4
#> 5      0        3    6   79          FALSE             FALSE        3       5
#> 6      0       11    7   86          FALSE              TRUE       11       6
#> 7      0        1    8   67          FALSE             FALSE        1       7
#> 8      0        7    9   14          FALSE             FALSE        7       8
#> 9      0       11    9   64          FALSE              TRUE       11       9
#> 10     0       11    9   82          FALSE              TRUE       11      10
#> 11     0        6   10   25          FALSE             FALSE        6      11
#> 12     0       11   11   15          FALSE              TRUE       11      12
#> 13     0       11   11   90          FALSE              TRUE       11      13
#> 14     0       11   12   81          FALSE              TRUE       11      14
#> 15     0       11   13   44          FALSE              TRUE       11      15
#> 16     0       11   14   16          FALSE              TRUE       11      16
#> 17     0       11   15   43          FALSE              TRUE       11      17
#> 18     0       11   15   76          FALSE              TRUE       11      18
#> 19     0       11   18   79          FALSE              TRUE       11      19
#> 20     0       11   19   98          FALSE              TRUE       11      20
#> 21     0       11   20   55          FALSE              TRUE       11      21
#> 22     0       11   22   73          FALSE              TRUE       11      22
#> 23     0        8   22   86          FALSE             FALSE        8      23
#> 24     0        9   23   54          FALSE             FALSE        9      24
#> 25     0        3   23   59          FALSE             FALSE        3      25
#> 26     0        2   24   57          FALSE             FALSE        2      26
#> 27     0        4   24   79          FALSE             FALSE        4      27
#> 28     0       11   25   40          FALSE              TRUE       11      28
#> 29     0       11   26   60          FALSE              TRUE       11      29
#> 30     0       11   30   50          FALSE              TRUE       11      30
#> 31     0        9   30   75          FALSE             FALSE        9      31
#> 32     0       11   31   67          FALSE              TRUE       11      32
#> 33     0        5   31   95          FALSE             FALSE        5      33
#> 34     0       11   34   54          FALSE              TRUE       11      34
#> 35     0       11   34   80          FALSE              TRUE       11      35
#> 36     0       11   38   39          FALSE              TRUE       11      36
#> 37     0        5   38   56          FALSE             FALSE        5      37
#> 38     0       11   39   85          FALSE              TRUE       11      38
#> 39     0       10   41   69          FALSE             FALSE       10      39
#> 40     0       11   42   76          FALSE              TRUE       11      40
#> 41     0       11   43   97          FALSE              TRUE       11      41
#> 42     0       11   45   50          FALSE              TRUE       11      42
#> 43     0        7   46   82          FALSE             FALSE        7      43
#> 44     0        1   51   67          FALSE             FALSE        1      44
#> 45     0        4   53   75          FALSE             FALSE        4      45
#> 46     0       11   55   83          FALSE              TRUE       11      46
#> 47     0       11   56   68          FALSE              TRUE       11      47
#> 48     0        6   56   73          FALSE             FALSE        6      48
#> 49     0       11   58   61          FALSE              TRUE       11      49
#> 50     0        5   61   73          FALSE             FALSE        5      50
#> 51     0        7   62   72          FALSE             FALSE        7      51
#> 52     0       11   62   84          FALSE              TRUE       11      52
#> 53     0       11   64   67          FALSE              TRUE       11      53
#> 54     0       11   64   93          FALSE              TRUE       11      54
#> 55     0       11   65   77          FALSE              TRUE       11      55
#> 56     0       11   67   72          FALSE              TRUE       11      56
#> 57     0       11   68   85          FALSE              TRUE       11      57
#> 58     0       11   68   96          FALSE              TRUE       11      58
#> 59     0       11   69   99          FALSE              TRUE       11      59
#> 60     0       11   70   79          FALSE              TRUE       11      60
#> 61     0       10   77   98          FALSE             FALSE       10      61
#> 62     0       11   79   84          FALSE              TRUE       11      62
#> 63     1       11   41   65          FALSE              TRUE       10      63
#> 64     1       11    8   46          FALSE              TRUE       10      64
#> 65     2       10   70   89          FALSE             FALSE        8      65
#> 66     2       11   68   92          FALSE              TRUE        9      66
#> 67     3       11   54   63          FALSE              TRUE        8      67
#> 68     5       11   24   40          FALSE              TRUE        6      68
#> 69     5       11   26   77          FALSE              TRUE        6      69
#> 70     5       11   11   30          FALSE              TRUE        6      70
#> 71     5       11   53   89          FALSE              TRUE        6      71
#> 72     6       11    8   27          FALSE              TRUE        5      72
#> 73     7       11    3   63          FALSE              TRUE        4      73
#> 74     8       11   44   45          FALSE              TRUE        3      74
#> 75     8       11   17   75          FALSE              TRUE        3      75
#> 76     8       11   21   47          FALSE              TRUE        3      76
#> 77     9       11   12   54          FALSE              TRUE        2      77
#> 78     9       11    5    7          FALSE              TRUE        2      78
#> 79     9       11   28   87          FALSE              TRUE        2      79
#> 80     9       10   48   98          FALSE             FALSE        1      80
#> 81     9       11   40   87          FALSE              TRUE        2      81
#> 82    10       11    6   78          FALSE              TRUE        1      82
#> 83    10       11   72   95          FALSE              TRUE        1      83

# Extract data from all simulations
as.data.frame(dx)
#>     sim onset terminus tail head onset.censored terminus.censored duration
#> 1     1     0        6    1   46          FALSE             FALSE        6
#> 2     1     0        7    1   86          FALSE             FALSE        7
#> 3     1     0       10    2   97          FALSE             FALSE       10
#> 4     1     0       11    3    8          FALSE              TRUE       11
#> 5     1     0        3    6   79          FALSE             FALSE        3
#> 6     1     0       11    7   86          FALSE              TRUE       11
#> 7     1     0        1    8   67          FALSE             FALSE        1
#> 8     1     0        7    9   14          FALSE             FALSE        7
#> 9     1     0       11    9   64          FALSE              TRUE       11
#> 10    1     0       11    9   82          FALSE              TRUE       11
#> 11    1     0        6   10   25          FALSE             FALSE        6
#> 12    1     0       11   11   15          FALSE              TRUE       11
#> 13    1     0       11   11   90          FALSE              TRUE       11
#> 14    1     0       11   12   81          FALSE              TRUE       11
#> 15    1     0       11   13   44          FALSE              TRUE       11
#> 16    1     0       11   14   16          FALSE              TRUE       11
#> 17    1     0       11   15   43          FALSE              TRUE       11
#> 18    1     0       11   15   76          FALSE              TRUE       11
#> 19    1     0       11   18   79          FALSE              TRUE       11
#> 20    1     0       11   19   98          FALSE              TRUE       11
#> 21    1     0       11   20   55          FALSE              TRUE       11
#> 22    1     0       11   22   73          FALSE              TRUE       11
#> 23    1     0        8   22   86          FALSE             FALSE        8
#> 24    1     0        9   23   54          FALSE             FALSE        9
#> 25    1     0        3   23   59          FALSE             FALSE        3
#> 26    1     0        2   24   57          FALSE             FALSE        2
#> 27    1     0        4   24   79          FALSE             FALSE        4
#> 28    1     0       11   25   40          FALSE              TRUE       11
#> 29    1     0       11   26   60          FALSE              TRUE       11
#> 30    1     0       11   30   50          FALSE              TRUE       11
#> 31    1     0        9   30   75          FALSE             FALSE        9
#> 32    1     0       11   31   67          FALSE              TRUE       11
#> 33    1     0        5   31   95          FALSE             FALSE        5
#> 34    1     0       11   34   54          FALSE              TRUE       11
#> 35    1     0       11   34   80          FALSE              TRUE       11
#> 36    1     0       11   38   39          FALSE              TRUE       11
#> 37    1     0        5   38   56          FALSE             FALSE        5
#> 38    1     0       11   39   85          FALSE              TRUE       11
#> 39    1     0       10   41   69          FALSE             FALSE       10
#> 40    1     0       11   42   76          FALSE              TRUE       11
#> 41    1     0       11   43   97          FALSE              TRUE       11
#> 42    1     0       11   45   50          FALSE              TRUE       11
#> 43    1     0        7   46   82          FALSE             FALSE        7
#> 44    1     0        1   51   67          FALSE             FALSE        1
#> 45    1     0        4   53   75          FALSE             FALSE        4
#> 46    1     0       11   55   83          FALSE              TRUE       11
#> 47    1     0       11   56   68          FALSE              TRUE       11
#> 48    1     0        6   56   73          FALSE             FALSE        6
#> 49    1     0       11   58   61          FALSE              TRUE       11
#> 50    1     0        5   61   73          FALSE             FALSE        5
#> 51    1     0        7   62   72          FALSE             FALSE        7
#> 52    1     0       11   62   84          FALSE              TRUE       11
#> 53    1     0       11   64   67          FALSE              TRUE       11
#> 54    1     0       11   64   93          FALSE              TRUE       11
#> 55    1     0       11   65   77          FALSE              TRUE       11
#> 56    1     0       11   67   72          FALSE              TRUE       11
#> 57    1     0       11   68   85          FALSE              TRUE       11
#> 58    1     0       11   68   96          FALSE              TRUE       11
#> 59    1     0       11   69   99          FALSE              TRUE       11
#> 60    1     0       11   70   79          FALSE              TRUE       11
#> 61    1     0       10   77   98          FALSE             FALSE       10
#> 62    1     0       11   79   84          FALSE              TRUE       11
#> 63    1     1       11   41   65          FALSE              TRUE       10
#> 64    1     1       11    8   46          FALSE              TRUE       10
#> 65    1     2       10   70   89          FALSE             FALSE        8
#> 66    1     2       11   68   92          FALSE              TRUE        9
#> 67    1     3       11   54   63          FALSE              TRUE        8
#> 68    1     5       11   24   40          FALSE              TRUE        6
#> 69    1     5       11   26   77          FALSE              TRUE        6
#> 70    1     5       11   11   30          FALSE              TRUE        6
#> 71    1     5       11   53   89          FALSE              TRUE        6
#> 72    1     6       11    8   27          FALSE              TRUE        5
#> 73    1     7       11    3   63          FALSE              TRUE        4
#> 74    1     8       11   44   45          FALSE              TRUE        3
#> 75    1     8       11   17   75          FALSE              TRUE        3
#> 76    1     8       11   21   47          FALSE              TRUE        3
#> 77    1     9       11   12   54          FALSE              TRUE        2
#> 78    1     9       11    5    7          FALSE              TRUE        2
#> 79    1     9       11   28   87          FALSE              TRUE        2
#> 80    1     9       10   48   98          FALSE             FALSE        1
#> 81    1     9       11   40   87          FALSE              TRUE        2
#> 82    1    10       11    6   78          FALSE              TRUE        1
#> 83    1    10       11   72   95          FALSE              TRUE        1
#> 84    2     0        6    2   62          FALSE             FALSE        6
#> 85    2     0       11    4   74          FALSE              TRUE       11
#> 86    2     0       11    5   37          FALSE              TRUE       11
#> 87    2     0        3    6   66          FALSE             FALSE        3
#> 88    2     0        6    7   62          FALSE             FALSE        6
#> 89    2     0       11    7   72          FALSE              TRUE       11
#> 90    2     0       11    9   71          FALSE              TRUE       11
#> 91    2     0       11   11   28          FALSE              TRUE       11
#> 92    2     0       11   12   29          FALSE              TRUE       11
#> 93    2     0       11   12   68          FALSE              TRUE       11
#> 94    2     0       11   12   91          FALSE              TRUE       11
#> 95    2     0        3   13   37          FALSE             FALSE        3
#> 96    2     0       11   13   78          FALSE              TRUE       11
#> 97    2     0        7   13   89          FALSE             FALSE        7
#> 98    2     0       11   14   19          FALSE              TRUE       11
#> 99    2     0        8   14   24          FALSE             FALSE        8
#> 100   2     0       11   16   82          FALSE              TRUE       11
#> 101   2     0       11   19   51          FALSE              TRUE       11
#> 102   2     0       11   22   64          FALSE              TRUE       11
#> 103   2     0       11   24   61          FALSE              TRUE       11
#> 104   2     0       11   25   42          FALSE              TRUE       11
#> 105   2     0       11   26   40          FALSE              TRUE       11
#> 106   2     0        2   26   79          FALSE             FALSE        2
#> 107   2     0       11   29   40          FALSE              TRUE       11
#> 108   2     0        8   29   58          FALSE             FALSE        8
#> 109   2     0       11   33   80          FALSE              TRUE       11
#> 110   2     0       11   34   51          FALSE              TRUE       11
#> 111   2     0       11   35   75          FALSE              TRUE       11
#> 112   2     0       11   36   81          FALSE              TRUE       11
#> 113   2     0        3   36   89          FALSE             FALSE        3
#> 114   2     0        5   37   80          FALSE             FALSE        5
#> 115   2     0       11   38   74          FALSE              TRUE       11
#> 116   2     0        1   38   80          FALSE             FALSE        1
#> 117   2     0       11   39   67          FALSE              TRUE       11
#> 118   2     0       11   41   90          FALSE              TRUE       11
#> 119   2     0       11   45   51          FALSE              TRUE       11
#> 120   2     0       11   49   70          FALSE              TRUE       11
#> 121   2     0        1   55   68          FALSE             FALSE        1
#> 122   2     0        5   62   63          FALSE             FALSE        5
#> 123   2     0        1   64   75          FALSE             FALSE        1
#> 124   2     0       11   64   97          FALSE              TRUE       11
#> 125   2     0        1   66   88          FALSE             FALSE        1
#> 126   2     0        5   74   89          FALSE             FALSE        5
#> 127   2     0       11   78   86          FALSE              TRUE       11
#> 128   2     0        8   82   94          FALSE             FALSE        8
#> 129   2     1       11   11   61          FALSE              TRUE       10
#> 130   2     1       11   72   91          FALSE              TRUE       10
#> 131   2     1        6   10   33          FALSE             FALSE        5
#> 132   2     1       11   57   78          FALSE              TRUE       10
#> 133   2     2        5    4   47          FALSE             FALSE        3
#> 134   2     2       11   48   91          FALSE              TRUE        9
#> 135   2     2       11   68   92          FALSE              TRUE        9
#> 136   2     2       11   63   91          FALSE              TRUE        9
#> 137   2     3       11   13   60          FALSE              TRUE        8
#> 138   2     3       11   25   36          FALSE              TRUE        8
#> 139   2     4       11   21   37          FALSE              TRUE        7
#> 140   2     5       11    5   73          FALSE              TRUE        6
#> 141   2     5       11    7  100          FALSE              TRUE        6
#> 142   2     5       11   12   39          FALSE              TRUE        6
#> 143   2     5       11   23   32          FALSE              TRUE        6
#> 144   2     6       10    1   27          FALSE             FALSE        4
#> 145   2     6       11   90   92          FALSE              TRUE        5
#> 146   2     6       11   75   83          FALSE              TRUE        5
#> 147   2     6       11   21   69          FALSE              TRUE        5
#> 148   2     6       11   41   70          FALSE              TRUE        5
#> 149   2     7       10   10   43          FALSE             FALSE        3
#> 150   2     7       11   27   64          FALSE              TRUE        4
#> 151   2     7       11   15   35          FALSE              TRUE        4
#> 152   2     8       11   61   87          FALSE              TRUE        3
#> 153   2     8       11   35   90          FALSE              TRUE        3
#> 154   2     8       11   79   92          FALSE              TRUE        3
#> 155   2     9       10   39   49          FALSE             FALSE        1
#> 156   2     9       11   69   94          FALSE              TRUE        2
#> 157   2     9       11   40   82          FALSE              TRUE        2
#> 158   2     9       11    9   66          FALSE              TRUE        2
#> 159   2     9       11   36   47          FALSE              TRUE        2
#> 160   2    10       11   26   33          FALSE              TRUE        1
#> 161   2    10       11   18   32          FALSE              TRUE        1
#> 162   3     0       11    3   92          FALSE              TRUE       11
#> 163   3     0        3    5   30          FALSE             FALSE        3
#> 164   3     0        5    5   65          FALSE             FALSE        5
#> 165   3     0       11    6   40          FALSE              TRUE       11
#> 166   3     0        1    7   18          FALSE             FALSE        1
#> 167   3     0       11    7   71          FALSE              TRUE       11
#> 168   3     0        8    7   96          FALSE             FALSE        8
#> 169   3     0        8    8   55          FALSE             FALSE        8
#> 170   3     0        1    9   33          FALSE             FALSE        1
#> 171   3     0        2   10   41          FALSE             FALSE        2
#> 172   3     0        7   11   48          FALSE             FALSE        7
#> 173   3     0       11   12   76          FALSE              TRUE       11
#> 174   3     0        8   12   77          FALSE             FALSE        8
#> 175   3     0       11   12   80          FALSE              TRUE       11
#> 176   3     0       11   14   32          FALSE              TRUE       11
#> 177   3     0       11   14   35          FALSE              TRUE       11
#> 178   3     0       11   15   53          FALSE              TRUE       11
#> 179   3     0       11   16   94          FALSE              TRUE       11
#> 180   3     0       11   17   27          FALSE              TRUE       11
#> 181   3     0       11   17   39          FALSE              TRUE       11
#> 182   3     0       11   17   74          FALSE              TRUE       11
#> 183   3     0        6   18   27          FALSE             FALSE        6
#> 184   3     0       11   18   83          FALSE              TRUE       11
#> 185   3     0       11   19   58          FALSE              TRUE       11
#> 186   3     0        8   20   48          FALSE             FALSE        8
#> 187   3     0        7   21   44          FALSE             FALSE        7
#> 188   3     0        1   21   79          FALSE             FALSE        1
#> 189   3     0       10   25   95          FALSE             FALSE       10
#> 190   3     0        6   28   67          FALSE             FALSE        6
#> 191   3     0       11   31   46          FALSE              TRUE       11
#> 192   3     0        7   32   53          FALSE             FALSE        7
#> 193   3     0       11   32   74          FALSE              TRUE       11
#> 194   3     0       11   34   37          FALSE              TRUE       11
#> 195   3     0       11   34   73          FALSE              TRUE       11
#> 196   3     0        4   36   56          FALSE             FALSE        4
#> 197   3     0       11   37   90          FALSE              TRUE       11
#> 198   3     0        2   40   69          FALSE             FALSE        2
#> 199   3     0        3   41   50          FALSE             FALSE        3
#> 200   3     0       11   42   98          FALSE              TRUE       11
#> 201   3     0        6   43   50          FALSE             FALSE        6
#> 202   3     0        4   47   80          FALSE             FALSE        4
#> 203   3     0       11   50   94          FALSE              TRUE       11
#> 204   3     0        1   53   99          FALSE             FALSE        1
#> 205   3     0       11   55   63          FALSE              TRUE       11
#> 206   3     0        2   58   96          FALSE             FALSE        2
#> 207   3     0        3   61   87          FALSE             FALSE        3
#> 208   3     0       11   61   99          FALSE              TRUE       11
#> 209   3     0        2   63   83          FALSE             FALSE        2
#> 210   3     0       11   63   90          FALSE              TRUE       11
#> 211   3     0       11   65   81          FALSE              TRUE       11
#> 212   3     0       11   65   82          FALSE              TRUE       11
#> 213   3     0       11   67   81          FALSE              TRUE       11
#> 214   3     0       11   70   85          FALSE              TRUE       11
#> 215   3     0        2   72   73          FALSE             FALSE        2
#> 216   3     0       11   89   94          FALSE              TRUE       11
#> 217   3     1       11   83   98          FALSE              TRUE       10
#> 218   3     1       11   37   71          FALSE              TRUE       10
#> 219   3     1       11   10   66          FALSE              TRUE       10
#> 220   3     1        7   53   60          FALSE             FALSE        6
#> 221   3     2       11   10   77          FALSE              TRUE        9
#> 222   3     2       11    7   67          FALSE              TRUE        9
#> 223   3     3       11   67   90          FALSE              TRUE        8
#> 224   3     3       11    5   50          FALSE              TRUE        8
#> 225   3     4        7   26   74          FALSE             FALSE        3
#> 226   3     4       11    8   59          FALSE              TRUE        7
#> 227   3     5       11    3   14          FALSE              TRUE        6
#> 228   3     6       10   15   57          FALSE             FALSE        4
#> 229   3     6       11   65   97          FALSE              TRUE        5
#> 230   3     6       11   28   64          FALSE              TRUE        5
#> 231   3     7       11   21   68          FALSE              TRUE        4
#> 232   3     7       11    2   83          FALSE              TRUE        4
#> 233   3     7       11   21   26          FALSE              TRUE        4
#> 234   3     7       11   89   91          FALSE              TRUE        4
#> 235   3     7       11   55   92          FALSE              TRUE        4
#> 236   3     7       11   10   99          FALSE              TRUE        4
#> 237   3     8       11   81   90          FALSE              TRUE        3
#> 238   3     8       10    5   89          FALSE             FALSE        2
#> 239   3     8       11   18   55          FALSE              TRUE        3
#> 240   3     8       11   17   64          FALSE              TRUE        3
#> 241   3     8       11   33   53          FALSE              TRUE        3
#> 242   3     9       11   27   62          FALSE              TRUE        2
#> 243   3     9       11   34   42          FALSE              TRUE        2
#> 244   3     9       11   13   60          FALSE              TRUE        2
#> 245   3     9       11    6   29          FALSE              TRUE        2
#> 246   3     9       11   31   77          FALSE              TRUE        2
#> 247   3     9       11   11   45          FALSE              TRUE        2
#> 248   3    10       11   27   45          FALSE              TRUE        1
#> 249   3    10       11   29   58          FALSE              TRUE        1
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
#> 84        1
#> 85        2
#> 86        3
#> 87        4
#> 88        5
#> 89        6
#> 90        7
#> 91        8
#> 92        9
#> 93       10
#> 94       11
#> 95       12
#> 96       13
#> 97       14
#> 98       15
#> 99       16
#> 100      17
#> 101      18
#> 102      19
#> 103      20
#> 104      21
#> 105      22
#> 106      23
#> 107      24
#> 108      25
#> 109      26
#> 110      27
#> 111      28
#> 112      29
#> 113      30
#> 114      31
#> 115      32
#> 116      33
#> 117      34
#> 118      35
#> 119      36
#> 120      37
#> 121      38
#> 122      39
#> 123      40
#> 124      41
#> 125      42
#> 126      43
#> 127      44
#> 128      45
#> 129      46
#> 130      47
#> 131      48
#> 132      49
#> 133      50
#> 134      51
#> 135      52
#> 136      53
#> 137      54
#> 138      55
#> 139      56
#> 140      57
#> 141      58
#> 142      59
#> 143      60
#> 144      61
#> 145      62
#> 146      63
#> 147      64
#> 148      65
#> 149      66
#> 150      67
#> 151      68
#> 152      69
#> 153      70
#> 154      71
#> 155      72
#> 156      73
#> 157      74
#> 158      75
#> 159      76
#> 160      77
#> 161      78
#> 162       1
#> 163       2
#> 164       3
#> 165       4
#> 166       5
#> 167       6
#> 168       7
#> 169       8
#> 170       9
#> 171      10
#> 172      11
#> 173      12
#> 174      13
#> 175      14
#> 176      15
#> 177      16
#> 178      17
#> 179      18
#> 180      19
#> 181      20
#> 182      21
#> 183      22
#> 184      23
#> 185      24
#> 186      25
#> 187      26
#> 188      27
#> 189      28
#> 190      29
#> 191      30
#> 192      31
#> 193      32
#> 194      33
#> 195      34
#> 196      35
#> 197      36
#> 198      37
#> 199      38
#> 200      39
#> 201      40
#> 202      41
#> 203      42
#> 204      43
#> 205      44
#> 206      45
#> 207      46
#> 208      47
#> 209      48
#> 210      49
#> 211      50
#> 212      51
#> 213      52
#> 214      53
#> 215      54
#> 216      55
#> 217      56
#> 218      57
#> 219      58
#> 220      59
#> 221      60
#> 222      61
#> 223      62
#> 224      63
#> 225      64
#> 226      65
#> 227      66
#> 228      67
#> 229      68
#> 230      69
#> 231      70
#> 232      71
#> 233      72
#> 234      73
#> 235      74
#> 236      75
#> 237      76
#> 238      77
#> 239      78
#> 240      79
#> 241      80
#> 242      81
#> 243      82
#> 244      83
#> 245      84
#> 246      85
#> 247      86
#> 248      87
#> 249      88
```
