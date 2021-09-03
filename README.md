# hccPIRS

<!-- badges: start -->
<!-- badges: end -->

This package calculates replication stress-related prognostic index (PIRS) for HBV-associated HCC patients, and estimates the enrichment of 21 replication stress signatures. If specified, a heatmap will be generated to show the landscape of the replication stress signatures in an ascending order of PIRS score.

## Citation

For now, you can cite the following bioRxiv preprint

## Installation

You may install this package with:

``` r
if (!require("devtools")) 
    install.packages("devtools")
devtools::install_github("xlucpu/hccPIRS")
```

## Example

``` r
library(hccPIRS)
## basic example code
library(hccPIRS)
load(system.file("extdata", "tpm.demo.RData", package = "function", mustWork = TRUE)) # load example data

res <- hccPIRS(expr = tpm.demo,
               scaleFlag  = FALSE,
               centerFlag = FALSE,
               doplot = TRUE,
               fig.path = getwd(),
               fig.name   = "heatmap of replication stress",
               enrich = "gsva",
               width = 6,
               height = 4)

pirs <- res$pirs # extract normalized PIRS score for each sample
head(pirs)

rsMat <- res$RS.sscore # extract enrichment score for replication stress signatures
head(rsMat)

res$hm # show the heatmap
```

