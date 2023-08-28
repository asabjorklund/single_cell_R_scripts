Doublet rate prediction
================
Asa Bjorklund
28/08/23

- <a href="#2-samples-at-equal-distribution"
  id="toc-2-samples-at-equal-distribution">2 samples at equal
  distribution.</a>
- <a href="#test-different-scenarios"
  id="toc-test-different-scenarios">Test different scenarios</a>
  - <a href="#more-of-1-sample" id="toc-more-of-1-sample">More of 1
    sample</a>
  - <a href="#3-samples" id="toc-3-samples">3 samples</a>
  - <a href="#10-samples" id="toc-10-samples">10 samples</a>
  - <a href="#20-samples" id="toc-20-samples">20 samples</a>
- <a href="#real-data" id="toc-real-data">Real data</a>

``` r
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))

source("predict_hashed_doublet_rates.R")
```

Given the number of cells per hashtag in an experiment and the number of
double tags observed, the intra-sample doublet rate can be calculated
for each sample.

Assuming that the sample distribution in the hashed cells and in the
doublets is the same, e.g. that there is no sample bias in doublet
creation. Also, assumes that all multiplets consists of 2 cells.

The function `predict_hashed_doublet_rates` will create artificial
doublets from the observed distribution, randomly selected 100 times.
Simulations are done across different doublet rates and then identifies
the rate closest to the observed for across hashtag doublets.

## 2 samples at equal distribution.

Define cells as 500 of each sample and 75 multiplets.

``` r
nMulti = 75
nCells = c(500,500)
names(nCells) = c("A","B")

print(nCells)
```

    ##   A   B 
    ## 500 500

``` r
predict_hashed_doublet_rates(nCells, 75)
```

![](doublet_rate_predictor_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["nDoublets"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["PercentageTotal"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["PercentageSample"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"37.16","2":"3.716","3":"7.432","_rn_":"A"},{"1":"36.80","2":"3.680","3":"7.360","_rn_":"B"},{"1":"75.04","2":"7.504","3":"NA","_rn_":"HashDoublet"},{"1":"149.00","2":"14.900","3":"NA","_rn_":"totalDoublets"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

## Test different scenarios

### More of 1 sample

``` r
nCells[1] = 1000
nCells
```

    ##    A    B 
    ## 1000  500

``` r
predict_hashed_doublet_rates(nCells, 75)
```

![](doublet_rate_predictor_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["nDoublets"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["PercentageTotal"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["PercentageSample"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"64.29","2":"4.286000","3":"6.429","_rn_":"A"},{"1":"17.13","2":"1.142000","3":"3.426","_rn_":"B"},{"1":"67.58","2":"4.505333","3":"NA","_rn_":"HashDoublet"},{"1":"149.00","2":"9.933333","3":"NA","_rn_":"totalDoublets"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

### 3 samples

Different number of each sample, still 75 observed multiplets.

``` r
nCells = c(200,500,700)
names(nCells) = LETTERS[1:length(nCells)]
nCells
```

    ##   A   B   C 
    ## 200 500 700

``` r
predict_hashed_doublet_rates(nCells, 75)
```

![](doublet_rate_predictor_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["nDoublets"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["PercentageTotal"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["PercentageSample"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"2.67","2":"0.1907143","3":"1.335","_rn_":"A"},{"1":"17.12","2":"1.2228571","3":"3.424","_rn_":"B"},{"1":"32.13","2":"2.2950000","3":"4.590","_rn_":"C"},{"1":"75.08","2":"5.3628571","3":"NA","_rn_":"HashDoublet"},{"1":"127.00","2":"9.0714286","3":"NA","_rn_":"totalDoublets"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

### 10 samples

Different number of each sample

``` r
nCells = seq(100,1100,100)
names(nCells) = LETTERS[1:length(nCells)]
nCells
```

    ##    A    B    C    D    E    F    G    H    I    J    K 
    ##  100  200  300  400  500  600  700  800  900 1000 1100

``` r
predict_hashed_doublet_rates(nCells, 75)
```

![](doublet_rate_predictor_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["nDoublets"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["PercentageTotal"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["PercentageSample"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"0.02","2":"0.0003030303","3":"0.0200000","_rn_":"A"},{"1":"0.06","2":"0.0009090909","3":"0.0300000","_rn_":"B"},{"1":"0.18","2":"0.0027272727","3":"0.0600000","_rn_":"C"},{"1":"0.25","2":"0.0037878788","3":"0.0625000","_rn_":"D"},{"1":"0.41","2":"0.0062121212","3":"0.0820000","_rn_":"E"},{"1":"0.67","2":"0.0101515152","3":"0.1116667","_rn_":"F"},{"1":"0.95","2":"0.0143939394","3":"0.1357143","_rn_":"G"},{"1":"1.29","2":"0.0195454545","3":"0.1612500","_rn_":"H"},{"1":"1.39","2":"0.0210606061","3":"0.1544444","_rn_":"I"},{"1":"2.19","2":"0.0331818182","3":"0.2190000","_rn_":"J"},{"1":"2.37","2":"0.0359090909","3":"0.2154545","_rn_":"K"},{"1":"75.22","2":"1.1396969697","3":"NA","_rn_":"HashDoublet"},{"1":"85.00","2":"1.2878787879","3":"NA","_rn_":"totalDoublets"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

### 20 samples

Equal number of each sample

``` r
nCells = rep(500, 20)
names(nCells) = LETTERS[1:length(nCells)]
nCells
```

    ##   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T 
    ## 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500 500

``` r
predict_hashed_doublet_rates(nCells, 75)
```

![](doublet_rate_predictor_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["nDoublets"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["PercentageTotal"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["PercentageSample"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"0.22","2":"0.0022","3":"0.044","_rn_":"A"},{"1":"0.22","2":"0.0022","3":"0.044","_rn_":"B"},{"1":"0.20","2":"0.0020","3":"0.040","_rn_":"C"},{"1":"0.25","2":"0.0025","3":"0.050","_rn_":"D"},{"1":"0.20","2":"0.0020","3":"0.040","_rn_":"E"},{"1":"0.22","2":"0.0022","3":"0.044","_rn_":"F"},{"1":"0.20","2":"0.0020","3":"0.040","_rn_":"G"},{"1":"0.16","2":"0.0016","3":"0.032","_rn_":"H"},{"1":"0.13","2":"0.0013","3":"0.026","_rn_":"I"},{"1":"0.12","2":"0.0012","3":"0.024","_rn_":"J"},{"1":"0.22","2":"0.0022","3":"0.044","_rn_":"K"},{"1":"0.18","2":"0.0018","3":"0.036","_rn_":"L"},{"1":"0.23","2":"0.0023","3":"0.046","_rn_":"M"},{"1":"0.17","2":"0.0017","3":"0.034","_rn_":"N"},{"1":"0.22","2":"0.0022","3":"0.044","_rn_":"O"},{"1":"0.21","2":"0.0021","3":"0.042","_rn_":"P"},{"1":"0.17","2":"0.0017","3":"0.034","_rn_":"Q"},{"1":"0.19","2":"0.0019","3":"0.038","_rn_":"R"},{"1":"0.20","2":"0.0020","3":"0.040","_rn_":"S"},{"1":"0.23","2":"0.0023","3":"0.046","_rn_":"T"},{"1":"75.06","2":"0.7506","3":"NA","_rn_":"HashDoublet"},{"1":"79.00","2":"0.7900","3":"NA","_rn_":"totalDoublets"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

## Real data

Data from Nima with 4 samples, uneven distribution of cell numbers.

``` r
stats = read.csv("stats.txt")

nCells = stats$Cells_Called[-1]
names(nCells) = stats$Sample_Tag[-1]

barplot(nCells, las = 2)
```

![](doublet_rate_predictor_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
nCells
```

    ## SampleTag01_hs SampleTag02_hs SampleTag03_hs SampleTag04_hs SampleTag05_hs 
    ##              0              2              0              0              0 
    ## SampleTag06_hs SampleTag07_hs SampleTag08_hs SampleTag09_hs SampleTag10_hs 
    ##              0              0           3806            315           1291 
    ## SampleTag11_hs SampleTag12_hs      Multiplet   Undetermined 
    ##           1479              0            310            112

Have 4 samples with 315 to 3806 cells, 4% doublets and 1.5% undermined.

``` r
# Ignore all other sample tags
nMulti = nCells["Multiplet"]
nCells = nCells[c(8:11)]

predict_hashed_doublet_rates(nCells, nMulti)
```

![](doublet_rate_predictor_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["nDoublets"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["PercentageTotal"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["PercentageSample"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"155.90","2":"2.26237121","3":"4.0961640","_rn_":"SampleTag08_hs"},{"1":"1.15","2":"0.01668843","3":"0.3650794","_rn_":"SampleTag09_hs"},{"1":"17.39","2":"0.25235815","3":"1.3470178","_rn_":"SampleTag10_hs"},{"1":"24.19","2":"0.35103759","3":"1.6355646","_rn_":"SampleTag11_hs"},{"1":"310.37","2":"4.50399071","3":"NA","_rn_":"HashDoublet"},{"1":"509.00","2":"7.38644609","3":"NA","_rn_":"totalDoublets"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>
