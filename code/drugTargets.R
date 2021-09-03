#-----------------------------------------#
# identification of potential drug target #
#-----------------------------------------#
library(ggplot2)

# load drug target
drug.target <- read.delim(file = file.path(data.path,"targetable gene list.txt"),sep = "\t",row.names = NULL,header = T,check.names = F,stringsAsFactors = F)
drug.target <- unique(drug.target$`Target genes`)

# get candidates from TCGA-LIHC and CN-LIHC using transcriptome expression profiles using pearson correlation
cor.drug.tcga.res <- cor.drug.cn.res <- NULL
for (i in drug.target) {
  cor.res <- cor.test(as.numeric(tcgahbvmt.expr.rmbatch[i,names(pirs.tcga)]),pirs.tcga,method = "pearson")
  cor.drug.tcga.res <- rbind.data.frame(cor.drug.tcga.res,
                                        data.frame(target = i,
                                                   r = cor.res$estimate,
                                                   p = cor.res$p.value,
                                                   stringsAsFactors = F),
                                        stringsAsFactors = F)
  
  cor.res <- cor.test(as.numeric(cn.expr.rmbatch[i,names(pirs.cn)]),pirs.cn,method = "pearson")
  cor.drug.cn.res <- rbind.data.frame(cor.drug.cn.res,
                                      data.frame(target = i,
                                                 r = cor.res$estimate,
                                                 p = cor.res$p.value,
                                                 stringsAsFactors = F),
                                      stringsAsFactors = F)
  
}
potential.target.tcga <- cor.drug.tcga.res[which(cor.drug.tcga.res$r > 0.3 & cor.drug.tcga.res$p < 0.05),]
potential.target.cn <- cor.drug.cn.res[which(cor.drug.cn.res$r > 0.3 & cor.drug.cn.res$p < 0.05),]

rownames(potential.target.tcga) <- potential.target.tcga$target
rownames(potential.target.cn) <- potential.target.cn$target

# load CERES data
ceres <- read.csv(file.path(data.path,"CRISPR_gene_effect_21Q2.csv"),row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ceres <- as.data.frame(t(ceres))
rownames(ceres) <- sapply(strsplit(rownames(ceres)," ",fixed = T),"[",1)

dep.sam <- names(pirs.depmap) # get cell lines name
ceres.hcc <- ceres[,intersect(dep.sam,colnames(ceres))]
cor.ceres.hcc.res <- NULL
for(i in rownames(ceres.hcc)) {
  
  # spearman correlation between CERES and PIRS scores of HCCLs
  cor.res <- cor.test(as.numeric(ceres.hcc[i,]),pirs.depmap[colnames(ceres.hcc)],method = "spearman")
  cor.ceres.hcc.res <- rbind.data.frame(cor.ceres.hcc.res,
                                        data.frame(target = i,
                                                   r = cor.res$estimate,
                                                   p = cor.res$p.value,
                                                   stringsAsFactors = F),
                                        stringsAsFactors = F)
}
potential.target.ceres <- cor.ceres.hcc.res[which(cor.ceres.hcc.res$r < -0.3 & cor.ceres.hcc.res$p < 0.05),]
rownames(potential.target.ceres) <- potential.target.ceres$target
potential.target <- intersect(potential.target.tcga$target,potential.target.cn$target)
final.target <- intersect(potential.target,potential.target.ceres$target) # initial drug targets

# boxplot of CERES scores
pdf(file.path(fig.path,"boxplot for potential drug targets candidate.pdf"), width = 3,height = 3)
par(bty="o", mgp = c(2,0.5,0), mar = c(5.1,3.1,2.1,2.1),tcl=-.25, xpd = F, las = 1)
boxinput <- ceres.hcc[final.target,]
boxplot(t(boxinput[names(sort(rowMeans(boxinput))),]),
        ylab = "Mean of CERES scores",
        outline = F,
        xaxt = "n")
abline(h = 0, lwd = 1.5, lty = 2, col = "grey70")
axis(side = 1, at = 1:11, labels = F)
text(1:11, par("usr")[3]-0.05, labels = names(sort(rowMeans(boxinput))), pos = 2,srt = 90, xpd = TRUE, cex = 0.9)
invisible(dev.off())

# correlation plot
pdf(file = file.path(fig.path,"combined correlation between drug target and pirs in tcga and cn cohort.pdf"), width = 5, height = 5)
par(bty="o", mgp = c(2,0.5,0), mar = c(3.1,3.1,2.1,2.1),tcl=-.25, xpd = F, las = 1)
## TCGA-LIHC
plot(cor.drug.tcga.res$r,
     -log10(cor.drug.tcga.res$p),
     xlab = "Pearson's correlation coefficient",
     ylab = "-log10(P-value)",
     type = "p",
     pch = 19,
     cex = 0.8,
     xaxt = "n",
     col = "#BDBDBD")
axis(side = 1, at = c(-0.5,-0.3,0,0.3,0.5,1), labels = c(-0.5,-0.3,0,0.3,0.5,1))
abline(h = -log10(0.05), lwd = 2, lty = 2)
abline(v = 0.3, lwd = 2, lty = 2)
points(potential.target.tcga$r,
       -log10(potential.target.tcga$p),
       pch = 19,
       cex = 1.4,
       col = "white")
points(potential.target.tcga$r,
       -log10(potential.target.tcga$p),
       pch = 19,
       cex = 1.4,
       col = alpha(jama[5],0.2))
points(potential.target.tcga[which(potential.target.tcga$target %in% final.target),"r"],
       -log10(potential.target.tcga[which(potential.target.tcga$target %in% final.target),"p"]),
       pch = 19,
       cex = 1.7,
       col = jama[1])
text(potential.target.tcga[which(potential.target.tcga$target %in% final.target),"r"] - 0.05,
     -log10(potential.target.tcga[which(potential.target.tcga$target %in% final.target),"p"]) - 5,
     labels = potential.target.tcga[which(potential.target.tcga$target %in% final.target),"target"])

## CN-LIHC
points(cor.drug.cn.res$r,
       -log10(cor.drug.cn.res$p),
       pch = 19,
       cex = 0.8,
       col = "#BDBDBD")
points(potential.target.cn$r,
       -log10(potential.target.cn$p),
       pch = 19,
       cex = 1.4,
       col = "white")
points(potential.target.cn$r,
       -log10(potential.target.cn$p),
       pch = 19,
       cex = 1.4,
       col = alpha(jama[2],0.2))
points(cor.drug.cn.res[which(cor.drug.cn.res$target %in% final.target),"r"],
       -log10(cor.drug.cn.res[which(cor.drug.cn.res$target %in% final.target),"p"]),
       pch = 19,
       cex = 1.7,
       col = jama[4])
text(cor.drug.cn.res[which(cor.drug.cn.res$target %in% final.target),"r"] - 0.05,
     -log10(cor.drug.cn.res[which(cor.drug.cn.res$target %in% final.target),"p"]) + 5,
     labels = cor.drug.cn.res[which(cor.drug.cn.res$target %in% final.target),"target"])
legend("topleft",
       col = jama[1:2],
       legend = c("TCGA-LIHC","CN-LIHC"),
       pch = 19,
       bty = "n")
invisible(dev.off())

# remove targets with average CERES score greater than 0
tmp <- ceres.hcc[final.target,]
tmp <- tmp[rowMeans(tmp) < 0,]
final.target <- rownames(tmp) 

pdf(file = file.path(fig.path,"correlation between drug target and pirs in ccle cohort.pdf"), width = 5, height = 5)
par(bty="o", mgp = c(2,0.5,0), mar = c(3.1,3.1,2.1,2.1),tcl=-.25, xpd = F, las = 1)
plot(cor.ceres.hcc.res$r,
     -log10(cor.ceres.hcc.res$p),
     xlab = "Spearman's correlation coefficient",
     ylab = "-log10(P-value)",
     type = "p",
     pch = 19,
     xlim = c(-1,1),
     cex = 0.8,
     col = "#BDBDBD")
abline(h = -log10(0.05), lwd = 2, lty = 2)
abline(v = -0.3, lwd = 2, lty = 2)
points(potential.target.ceres$r,
       -log10(potential.target.ceres$p),
       pch = 19,
       cex = 1.4,
       col = "white")
points(potential.target.ceres$r,
       -log10(potential.target.ceres$p),
       pch = 19,
       cex = 1.4,
       col = alpha("#4AD0EB",0.2))
points(potential.target.ceres[which(potential.target.ceres$target %in% final.target),"r"],
       -log10(potential.target.ceres[which(potential.target.ceres$target %in% final.target),"p"]),
       pch = 19,
       cex = 2,
       col = alpha(c("red",jama[4],jama[5],jama[2],jama[6],jama[1],jama[7]),0.8))
invisible(dev.off())
