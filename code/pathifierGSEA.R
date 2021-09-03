#----------------------------------------------#
# Functional analysis using pathifier and GSEA #
#----------------------------------------------#

library(pathifier)
library(clusterProfiler)

## TCGA mRNA for pathfier
tpm.norm <- read.table(file.path(data.path,"TCGA_LIHCnorm_mRNA_TPM.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
msigdb.hm <- gmt2list(file.path(data.path,"h.all.v7.4.symbols.gmt")) 
indata <- cbind.data.frame(tpm[rownames(tpm.norm),names(pirs.tcga)],tpm.norm)
indata <- as.matrix(log2(indata + 1))
path.res <- quantify_pathways_deregulation(data = indata, 
                                           allgenes = rownames(indata), 
                                           syms = msigdb.hm, 
                                           pathwaynames = names(msigdb.hm), 
                                           normals = rep(c(FALSE,TRUE),c(103,50)), 
                                           attempts = 100, 
                                           min_exp = 0, 
                                           min_std = 0)

path.score <- as.data.frame(rbindlist(lapply(path.res$scores,as.data.frame))); dimnames(path.score) <- list(names(path.res$xs),colnames(indata))

path.score.cor <- NULL
for (i in rownames(path.score)) {
  tmp <- data.frame(pirs = as.numeric(pirs.tcga),
                    expr = as.numeric(path.score[i,names(pirs.tcga)]),
                    stringsAsFactors = F)
  path.score.cor <- rbind.data.frame(path.score.cor,
                                     data.frame(path = i,
                                                avg = mean(as.numeric(path.score[i,names(pirs.tcga)])),
                                                r = cor.test(tmp$pirs, tmp$expr, method = "pearson")$estimate,
                                                p = cor.test(tmp$pirs, tmp$expr, method = "pearson")$p.value,
                                                stringsAsFactors = F),
                                     stringsAsFactors = F)
}
rownames(path.score.cor) <- path.score.cor$path
rownames(path.score.cor) <- gsub("HALLMARK_","",rownames(path.score.cor))

## CN mRNA for pathfier
cn.expr2 <- read.delim(file.path(data.path,"LIHC-CN-TPM.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
cn.expr2 <- cn.expr2[rowSums(cn.expr2) > 0,]
cn.expr2 <- cn.expr2[apply(cn.expr2,1,var) > 0,]
path.res.cn <- quantify_pathways_deregulation(data = as.matrix(cn.expr2), 
                                              allgenes = rownames(cn.expr2), 
                                              syms = msigdb.hm, 
                                              pathwaynames = names(msigdb.hm), 
                                              normals = rep(c(TRUE,FALSE),c(159,159)), 
                                              attempts = 100, 
                                              min_exp = 0, 
                                              min_std = 0) 

path.score.cn <- as.data.frame(rbindlist(lapply(path.res.cn$scores,as.data.frame))); dimnames(path.score.cn) <- list(names(path.res.cn$xs),colnames(cn.expr2))
colnames(path.score.cn) <- c(cn.sinfo$`Adjacent liver tissue (N) sample ID`,rownames(cn.sinfo))
path.score.cor.cn <- NULL
for (i in rownames(path.score.cn)) {
  tmp <- data.frame(pirs = as.numeric(pirs.cn),
                    expr = as.numeric(path.score.cn[i,names(pirs.cn)]),
                    stringsAsFactors = F)
  path.score.cor.cn <- rbind.data.frame(path.score.cor.cn,
                                        data.frame(path = i,
                                                   avg = mean(as.numeric(path.score.cn[i,names(pirs.cn)])),
                                                   r = cor.test(tmp$pirs, tmp$expr, method = "pearson")$estimate,
                                                   p = cor.test(tmp$pirs, tmp$expr, method = "pearson")$p.value,
                                                   stringsAsFactors = F),
                                        stringsAsFactors = F)
}
rownames(path.score.cor.cn) <- path.score.cor.cn$path
rownames(path.score.cor.cn) <- gsub("HALLMARK_","",rownames(path.score.cor.cn))
compath <- intersect(rownames(path.score),rownames(path.score.cn))

## CN-LIHC protein for GSEA
MSIGDB.HM <- read.gmt(file.path(data.path,"h.all.v7.4.symbols.gmt"))
cor.res <- NULL
for (i in rownames(cn.prot)) {
  tmp <- data.frame(pirs = as.numeric(pirs.cn),
                    expr = as.numeric(cn.prot[i,names(pirs.cn)]),
                    stringsAsFactors = F)
  cor.res <- rbind.data.frame(cor.res,
                              data.frame(protein = i,
                                         r = cor.test(tmp$pirs, tmp$expr, method = "pearson")$estimate,
                                         p = cor.test(tmp$pirs, tmp$expr, method = "pearson")$p.value,
                                         stringsAsFactors = F),
                              stringsAsFactors = F)
}
geneList <- cor.res$r; names(geneList) <- cor.res$protein
geneList <- sort(geneList,decreasing = T)
gsea_pirs_cn <- GSEA(geneList = geneList,eps = 0,TERM2GENE = MSIGDB.HM,pvalueCutoff = 1,seed = T,minGSSize = 0,verbose = T)
gsea.res <- as.data.frame(gsea_pirs_cn)
rownames(gsea.res) <- gsub("HALLMARK_","",rownames(gsea.res))
gsea.res <- gsea.res[rownames(path.score.cor),]
gsea.res <- gsea.res[order(gsea.res$NES,decreasing = T),]
