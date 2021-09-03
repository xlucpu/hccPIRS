#--------------------------------------------------------#
# REM-meta analysis for differential expression and GSEA #
#--------------------------------------------------------#

library(MetaVolcanoR)

# REM meta-analysis
meta_degs_rem <- rem_mv(diffexp = diffexplist, # a list of differential expression analysis results
                        pcriteria = "padj",
                        foldchangecol = 'log2fc',
                        genenamecol = 'id',
                        geneidcol = NULL,
                        collaps = FALSE,
                        llcol = 'CI.L',
                        rlcol = 'CI.R',
                        vcol = NULL,
                        cvar = TRUE,
                        metathr = 0.01,
                        jobname = "LIHC-REM",
                        outputfolder = res.path,
                        draw = 'HTML',
                        ncores = 1)
meta_degs_rem@MetaVolcano
meta.res <- meta_degs_rem@metaresult

# perform GSEA
MSIGDB.HM <- read.gmt(file.path(data.path,"h.all.v7.4.symbols.gmt"))
geneList <- meta.res$randomSummary; names(geneList) <- meta.res$id
geneList <- sort(geneList,decreasing = T)
gsea_pirs_meta <- GSEA(geneList = geneList,eps = 0,TERM2GENE = MSIGDB.HM,pvalueCutoff = 1,nPerm = 10000,seed = T,minGSSize = 0,verbose = T)

## DEG volcano plot
x <- meta.res[,c("id","randomSummary","randomCi.lb","randomCi.ub","randomP")]
x$padj <- p.adjust(x$randomP, "BH")
upx <- x[which(x$padj < 0.05 & x$randomSummary > log2(1.5)),]
dnx <- x[which(x$padj < 0.05 & x$randomSummary < -(log2(1.5))),]
upx2 <- x[which(x$padj < 0.05 & x$randomSummary > log2(2)),]
dnx2 <- x[which(x$padj < 0.05 & x$randomSummary < -(log2(2))),]

pdf(file.path(fig.path,"volcano plot for meta differential analysis.pdf"), width = 8, height = 8)
par(bty="o", mgp = c(1.9,.33,0), mar=c(3.1,3.1,2.1,2.1)+.1, las=1, tcl=-.25)
plot(NULL, NULL, ylim = c(0,80), xlim = c(-3,3),
     xlab = "Summary log2FoldChange", ylab = "-log10(Summary FDR)",col="white",
     main = "Meta-analysis (n=606)")
rect(par("usr")[1],
     par("usr")[3],
     par("usr")[2],
     par("usr")[4],
     col = "white",
     border = F)
grid(col = "grey90", lty = 2, lwd = 1.5)
points(x$randomSummary,-log10(x$randomP), col = alpha("grey80",0.6), pch = 19, cex = 0.8)
points(upx$randomSummary,-log10(upx$randomP), col = alpha(jama[2],0.6), pch = 19, cex = 1.2)
for (i in 1:nrow(upx)) {
  lines(x = c(upx$randomCi.lb[i], upx$randomCi.ub[i]), y = c(-log10(upx$randomP)[i],-log10(upx$randomP)[i]), col = alpha(jama[2],0.6), lwd = 1.2)
}
points(dnx$randomSummary,-log10(dnx$randomP), col = alpha(jama[5],0.6), pch = 19, cex = 1.2)
for (i in 1:nrow(dnx)) {
  lines(x = c(dnx$randomCi.lb[i], dnx$randomCi.ub[i]), y = c(-log10(dnx$randomP)[i],-log10(dnx$randomP)[i]), col = alpha(jama[5],0.6), lwd = 1.2)
}
points(upx2$randomSummary,-log10(upx2$randomP), col = alpha(jama[4],0.6), pch = 19, cex = 1.5)
for (i in 1:nrow(upx2)) {
  lines(x = c(upx2$randomCi.lb[i], upx2$randomCi.ub[i]), y = c(-log10(upx2$randomP)[i],-log10(upx2$randomP)[i]), col = alpha(jama[4],0.6), lwd = 1.2)
}
points(dnx2$randomSummary,-log10(dnx2$randomP), col = alpha(jama[1],0.6), pch = 19, cex = 1.5)
for (i in 1:nrow(dnx2)) {
  lines(x = c(dnx2$randomCi.lb[i], dnx2$randomCi.ub[i]), y = c(-log10(dnx2$randomP)[i],-log10(dnx2$randomP)[i]), col = alpha(jama[1],0.6), lwd = 1.2)
}
par(new = T, bty="o")
plot(-1, -1,
     col = "white",
     xlim = c(-3,3), ylim = c(0,80),
     xlab = "", ylab = "",
     xaxt = "n", yaxt = "n")
invisible(dev.off())

## GSEA meta barplot
df <- as.data.frame(gsea_pirs_meta); df <- df[,c("ID","NES","p.adjust")]
df$ID <- gsub("HALLMARK_","",df$ID)
df$group <- ifelse(df$NES > 0 & df$p.adjust < 0.05,"c", ifelse(df$NES < 1 & df$p.adjust < 0.05,"a","b"))
sortdf <- df[order(df$NES),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)

ggplot(sortdf, aes(ID, NES, fill = group)) + geom_bar(stat = 'identity') +
  coord_flip() +
  scale_fill_manual(values = c(jama[1], 'snow3', jama[2]), guide = FALSE) +

  geom_text(data = subset(df, NES > 0),
            aes(x=ID, y= -0.05, label= paste0(" ", ID), color = group),
            size = 4,
            hjust = "inward" ) +
  geom_text(data = subset(df, NES < 0),
            aes(x=ID, y= 0.05, label=ID, color = group),
            size = 4, hjust = "outward") +
  scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) +

  xlab("") +ylab("Normalized Enrichment Score\nMeta analysis of PIHRS vs. PILRS")+
  theme_bw() +
  theme(panel.grid =element_blank(),
        axis.text.x = element_text(size = 10, colour = "black"),
        panel.border = element_rect(size = 0.6),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
ggsave(file.path(fig.path,"gsea plot for meta analysis of PIHRS vs. PILRS.pdf"), width = 8,height = 8)
