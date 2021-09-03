#' @name hccPIRS
#' @title Replication stress-related prognostic index in HCC
#' @description This function calculates replication stress-related prognostic index (PIRS) for HBV-associated HCC patients, and estimates the enrichment of 21 replication stress signatures.
#' @author Xiaofan Lu
#'
#' @param expr A numeric expression matrix with row features and sample columns.
#' @param scaleFlag A logic value to indicate if the expression data should be further scaled. FALSE by default.
#' @param centerFlag A logic value to indicate if the expression data should be further centered. TRUE by default.
#' @param doplot A logic value to indicate whether to generate heatmap of replication stress signatures and PIRS score; FALSE by default.
#' @param fig.path A string value to indicate the output path for storing the heatmap.
#' @param fig.name A string value to indicate the name of the heatmap.
#' @param enrich A string value to indicate the method for single-sample enrichment analysis. Allowed values contain c('gsva', 'ssgsea', 'zscore', 'plage'); 'gsva' by default.
#' @param width A numeric value to indicate the relative width for each cell in the heatmap; 1 by default.
#' @param height A numeric value to indicate the relative height for each cell in the heatmap; 4 by default.
#'
#' @return
#' @export
#' @import ComplexHeatmap
#' @import GSVA
#' @import circlize
#' @import gplots
#' @import grid
#' @references Dreyer, SB, Upstill-Goddard, R, Paulus-Hock, V, Paris, C, Lampraki, E-M, Dray, E, et al. (2021). Targeting DNA Damage Response and Replication Stress in Pancreatic Cancer. Gastroenterology 160: 362-377.e313.

#' @examples
#' library(hccPIRS)
#' load(system.file("extdata", "tpm.demo.RData", package = "hccPIRS", mustWork = TRUE)) # load example data
#' res <- hccPIRS(expr = tpm.demo,
#'                scaleFlag  = FALSE,
#'                centerFlag = FALSE,
#'                doplot = TRUE,
#'                fig.path = getwd(),
#'                fig.name   = "heatmap of replication stress",
#'                enrich = "gsva",
#'                width = 6,
#'                height = 4)
#' pirs <- res$pirs # extract normalized PIRS score for each sample
#' print(pirs)
#' rsMat <- res$RS.score # extract enrichment score for replication stress signatures
#' rsMat[1:21, 1:3]
#' res$hm # show the heatmap

hccPIRS <- function(expr = NULL,
                    scaleFlag  = FALSE,
                    centerFlag = FALSE,
                    doplot = TRUE,
                    fig.path = getwd(),
                    fig.name   = "heatmap of replication stress",
                    enrich = "gsva",
                    width = 1,
                    height = 4) {

  # customized function for min-max normalization
  range01 <- function(x){(x-min(x))/(max(x)-min(x)) * 10}

  # customized function for creating heatmap
  standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {
    outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
    if (!is.null(halfwidth)) {
      outdata[outdata>halfwidth]=halfwidth
      outdata[outdata<(-halfwidth)]= -halfwidth
    }
    return(outdata)
  }

  # initial checking
  if(max(expr) < 25 | (max(expr) >= 25 & min(expr) < 0)) {
    message("--expression profile seems to have veen standardised (z-score or log transformation), no more action will be performed.")
    gset <- expr
  }
  if(max(expr) >= 25 & min(expr) >= 0){
    message("--log2 transformation done for expression data.")
    gset <- log2(expr + 1)
  }
  if(!all(is.element(rownames(pirs.coeff),rownames(gset)))) {
    stop("--not all features matched in your expression profile, please check for the following missing:")
    print(setdiff(rownames(pirs.coeff), intersect(rownames(pirs.coeff),rownames(gset))))
  }

  # scale data
  emat <- t(scale(t(gset), scale = scaleFlag, center = centerFlag))

  # call PIRS
  pirs.raw <- apply(t(emat[rownames(pirs.coeff),]),1,function(x) {x %*% pirs.coeff$coeff})
  pirs <- range01(pirs.raw)

  # calculate enrichment for 21 replication stress signatures
  label <- c("Homologous recombination",
             "Activation of ART in response to replication stress",
             "E2F inhibition of pre replication complex",
             "Cell cycle",
             "Mitotic Spindle Checkpoint",
             "G2/M DNA damage checkpoint",
             "Homologous DNA Pairing/Strand Exchange",
             "HDR through Single Strand Annealing SSA",
             "Senescence Associated Heterochromatin Foci",
             "Cyclin A/B1 associated events during G2/M transition",
             "Chk1 Chk2 inactivation of Cyclin B/Cdk1 complex",
             "APC Cdc20 to APC Cdh1 conversion: late anaphase",
             "MMR",
             "G1/S DNA Damage Checkpoints",
             "Recognition of DNA damage: PCNA complex",
             "M/G1 Transition",
             "DNA Replication pre initiation",
             "Formation of TC NER Pre Incision Complex",
             "DNA Damage Recognition in GG NER",
             "Pyrimidine metabolism",
             "Fanconi Anemia Pathway")

  RS.score <- suppressWarnings(gsva(as.matrix(emat),
                    RS.signature,
                    method = enrich))
  rownames(RS.score) <- label

  # generate heatmap
  if(doplot) {
    indata <- standarize.fun(RS.score,
                             scaleFlag = TRUE,
                             centerFlag = TRUE,
                             halfwidth = )

    samOrder <- sort(pirs)
    col_fun = colorRamp2(c(-1,-0.5, 0, 0.5, 1), colorpanel(5,low="#44A0D5",mid = "white",high="#A94747"))
    top_anno <- anno_barplot(as.numeric(samOrder),
                             border = T,
                             bar_width = 1,
                             gp = gpar(fill = colorpanel(length(samOrder),low="#3B4E55",mid = "#81AF96",high="#D69044"),border =NA,lty="blank"),
                             height = unit(2, "cm"))
    hm <- Heatmap(as.matrix(indata),
                  border = TRUE,
                  col = col_fun,
                  cluster_columns = FALSE,
                  cluster_rows = FALSE,
                  show_row_names = TRUE,
                  row_names_side = "right",
                  show_column_names = FALSE,
                  width = ncol(indata)*unit(width, "mm"),
                  height = nrow(indata)*unit(height, "mm"),
                  row_names_gp = gpar(fontsize = 10, col = "black"),
                  column_names_gp = gpar(fontsize = 10),
                  column_names_side = "bottom",
                  name = "Enrichment",
                  top_annotation = HeatmapAnnotation(PIRS = top_anno))

    outFig <- paste0(fig.name,".pdf")
    pdf(file = file.path(fig.path, outFig))
    draw(hm, heatmap_legend_side = "left")
    invisible(dev.off())

    return(list(pirs.raw = pirs.raw,
                pirs = pirs,
                RS.score = RS.score,
                hm = hm))
  } else {
    return(list(pirs.raw = pirs.raw,
                pirs = pirs,
                RS.score = RS.score))
  }
}
