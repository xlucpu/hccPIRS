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
load(system.file("extdata", "tpm.demo.RData", package = "hccPIRS", mustWork = TRUE)) # load example data

res <- hccPIRS(expr       = tpm.demo,
               scaleFlag  = FALSE, # no scale for input data
               centerFlag = FALSE, # no center for input data
               doplot     = TRUE, # generate heatmap
               fig.path   = getwd(),
               fig.name   = "heatmap of replication stress",
               enrich     = "gsva", # use gsva to quantify enrichment
               width      = 6,
               height     = 4)

pirs <- res$pirs # extract normalized PIRS score for each sample
print(pirs)

rsMat <- res$RS.score # extract enrichment score for replication stress signatures
rsMat[1:21,1:3]

res$hm # show the heatmap
```
<img src="https://user-images.githubusercontent.com/57204704/131953321-cc06a99f-8f68-4505-ab5c-5bbbba886bbf.jpg height="230" align="right" />
