#-------------------------------------------#
# Ridge regression to predict drug response #
#-------------------------------------------#

library(data.table)
library(pRRophetic)
library(impute)
library(SimDesign)

# ridge prediction
depmap.sinfo <- read.csv(file.path(data.path,"PRISM/sample_info.csv"),row.names = 1,check.names = F,stringsAsFactors = F,header = T)
depmap.sinfo <- depmap.sinfo[which(depmap.sinfo$primary_disease %in% c("Liver Cancer")),]

depmap <- read.csv(file.path(data.path,"PRISM/secondary-screen-dose-response-curve-parameters.csv"),row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
comdepmap <- intersect(rownames(depmap.sinfo),unique(depmap$depmap_id))
depmap <- depmap[which(depmap$depmap_id %in% comdepmap),]

depmap <- depmap[,c("depmap_id","broad_id","auc")]
depmap.auc <- reshape(depmap, idvar = "depmap_id", timevar = "broad_id", direction = "wide")
rownames(depmap.auc) <- depmap.auc$depmap_id
depmap.auc <- depmap.auc[,-1]
colnames(depmap.auc) <- gsub("auc.","",colnames(depmap.auc),fixed = T)

depmap.expr <- fread(file.path(data.path,"PRISM/PRISM_CCLE_expression_21Q2.csv"),check.names = F,stringsAsFactors = F,header = T,data.table = F)
rownames(depmap.expr) <- depmap.expr[,1]; depmap.expr <- depmap.expr[,-1]
colnames(depmap.expr) <- sapply(strsplit(colnames(depmap.expr)," ",fixed = T), "[",1)
comdepmap <- intersect(rownames(depmap.auc),rownames(depmap.expr))
depmap.expr <- depmap.expr[comdepmap,]
depmap.auc <- depmap.auc[comdepmap,]
comdepmap.drug <- unique (grep(paste0(substr(selected.depmap.drug$id,1,13),collapse="|"), 
                             colnames(depmap.auc), value=TRUE))
comdepmap.drug <- comdepmap.drug[1:5]

# applied to TCGA-LIHC
trainExpr <- t(depmap.expr)
trainExpr <- trainExpr[rowSums(trainExpr) > 0,]
trainExpr <- trainExpr[apply(trainExpr, 1, function(x) {sum(x > 0) > 0.9 * ncol(trainExpr)}),]
trainPtype <- as.data.frame(depmap.auc[,comdepmap.drug])
trainPtype <- as.data.frame(t(impute.knn(t(as.matrix(trainPtype)))$data))

testExpr <- log2(tpm[hbv.mt] + 1)
testExpr <- testExpr[rowSums(testExpr) > 0,]
testExpr <- testExpr[apply(testExpr, 1, function(x) {sum(x > 0) > 0.9 * ncol(testExpr)}),]
comgene <- intersect(rownames(trainExpr),rownames(testExpr))
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- as.matrix(testExpr[comgene,])

drugsen.tcga <- NULL
for (i in 1:ncol(trainPtype)) { 
  display.progress(index = i,totalN = ncol(trainPtype))
  d <- colnames(trainPtype)[i]
  tmp <- as.vector(trainPtype[,d])
  
  ptypeOut <- quiet(calcPhenotype(trainingExprData = as.matrix(trainExpr),
                                  trainingPtype = tmp,
                                  testExprData = as.matrix(testExpr),
                                  powerTransformPhenotype = T,
                                  selection = 1))
  drugsen.tcga <- rbind.data.frame(drugsen.tcga,ptypeOut)
  
}
dimnames(drugsen.tcga) <- list(colnames(trainPtype),colnames(testExpr))
drugsen.tcga <- as.data.frame(t(drugsen.tcga))
quantile3 <- quantileCut(drugsen.tcga$pirs,3)
drugsen.tcga$pirs <- ifelse(as.character(quantile3) == levels(quantile3)[3],"PIHRS","PILRS")
