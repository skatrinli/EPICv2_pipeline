################################################################################
## Purpose: QQ and Manhattan Plots
################################################################################
install.packages("qqman")
library(qqman)

PROJECTNAME <- "PREG_TOLD_Dystocia_Binary"

## Load results file
df <- read.csv("/Volumes/som-ts/GynOb/SmithLab/NewData/People/Abby/PREG_TOLD_EWAS/PREG_TOLD_QC_FINAL_Dystocia_Binary_EWAS.csv")

################################################################################
## QQ plot
################################################################################
lambda <- qchisq(median(as.numeric(df$P.value),na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)

jpeg(filename = paste(PROJECTNAME,"QQplot.jpeg", sep = "_"), width = 2, height = 2, 
     units = "in", res = 300, pointsize = 4)
qq(df$P.value)
leg<-parse(text=paste("lambda==",round(lambda,2),sep=""))
legend("topleft", legend=leg)
dev.off()

################################################################################
## Manhattan plot
################################################################################
df$CpG_chrm <- gsub("chr","",df$CpG_chrm)
df$CpG_chrm[df$CpG_chrm == "X"] = 23
df$CpG_chrm[df$CpG_chrm == "Y"] = 24
df$CpG_chrm<-as.numeric(df$CpG_chrm)

jpeg(filename = paste(PROJECTNAME,"_Manhattan_plot.jpeg", sep = "_"), width = 4, height = 2, 
     units = "in", res = 300, pointsize = 4)
manhattan(df, chr = "CpG_chrm", bp = "CpG_beg", p = "P.value", snp = "CPG.Labels", col = c("darkred", "darkblue"), 
          chrlabs = c(1:22, "X", "Y"), genomewideline = -log10(9e-08), suggestiveline = F, logp = TRUE)
dev.off()
