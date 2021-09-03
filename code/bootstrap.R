#------------------------------------------------------------------#
# bootstrap approach to enhance the robustness of prognostic value #
#------------------------------------------------------------------#

library(survival)
outTab <- NULL
for(i in 3:ncol(surv)){ # survival information (OS in this case)
  
  # customized function
  display.progress = function (index, totalN, breakN=20) {
    if ( index %% ceiling(totalN/breakN)  ==0  ) {
      cat(paste(round(index*100/totalN), "% ", sep=""))
    }
  }    
  
  display.progress(index = i, totalN = ncol(surv)) # show running progression
  gene <- colnames(surv)[i]
  Mboot <- replicate(1000, expr = { # bootstrap for 1,000 times
    indices <- sample(rownames(surv), size = nrow(surv) * 0.8, replace = F) # extract 80% samples at each bootsratp
    data <- surv[indices,]
    fmla1 <- as.formula(Surv(data[,"OS.time"],data[,"OS"]) ~ data[,gene])
    mycox <- coxph(fmla1,data = data)
    coxResult <- summary(mycox)
    P <- coxResult$coefficients[,"Pr(>|z|)"]
  }
  )
  times <- length(Mboot[which(Mboot < 0.01)])
  outTab <- rbind(outTab,
                  cbind(gene = gene,
                        times = times))
}
outTab <- as.data.frame(outTab)
bootGene <- outTab[as.numeric(as.character(outTab$times)) > 800,] # get final genes which present in more than 800 times