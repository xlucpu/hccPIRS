#----------------------------------------------------------#
# Random survival forest with independent variable hunting #
#----------------------------------------------------------#

library(randomForestSRC)
res.rsf <- rfsrc(Surv(OS.time, OS) ~ ., surv, 
                 nodesize = 20, 
                 proximity = T, 
                 tree.err = T, 
                 forest = T, 
                 ntree = 1000,
                 splitrule = "logrankscore", 
                 importance = TRUE)

res.trc <- res.trcoob <- c()
topvars <- list()
for (j in 1:1000) { 
  print(paste("trying for", j, "times"))
  vars<-var.select(object = res.rsf,
                   cause = 1,
                   method = "vh.vimp", 
                   conservative = "high", 
                   ntree = 1000,
                   nodesize = 20,
                   splitrule = "logrankscore", 
                   nsplit = 10, 
                   xvar.wt = NULL,
                   refit = T, 
                   fast = T,
                   na.action = "na.impute", 
                   always.use = NULL, 
                   nrep = 10,
                   prefit =  list(action = T, 
                                  ntree = 1000,
                                  nodesize = 20, 
                                  nsplit = 10),
                   verbose = TRUE)
  
  trc <- rcorr.cens(-vars$rfsrc.refit.obj$predicted, 
                  Surv(surv$OS.time, surv$OS))["C Index"]
  trcoob <- rcorr.cens(-vars$rfsrc.refit.obj$predicted.oob, 
                     Surv(surv$OS.time, surv$OS))["C Index"]
  
  res.trc <- rbind(res.trc, trc)
  res.trcoob <- rbind(res.trcoob, trcoob)
  if(length(vars$topvars) == 0) {
    topvars[[j]] <- "[Not Available]"
  } else {
    topvars[[j]] <- vars$topvars
  }
}
result <- data.frame(res.trc,res.trcoob,row.names = 1:nrow(res.trc))
colnames(result) <- c("res.trc.cindex","res.trcoob.cindex")
bestresult <- result[result$res.trcoob.cindex == max(result$res.trcoob.cindex),] 
bestvars <- unique(topvars[[as.numeric(rownames(bestresult))]]) 
rsf.res <- cbind(surv[,1:2],surv[,as.character(bestvars)])
