################################################################################
## Purpose: Plotting CpGs for Categorical Variables
################################################################################
## Install and Load packages
install.packages("ggpubr")
library(ggpubr)

## Define variables
VarName <- "Dystocia_Binary" # the name of the variable in the phenotype file
VarLabel <- "Dystocia" # how you want the variable name to be shown in the x-axis of the plot
cg <- c("cg18718657",
        "cg01482375",
        "cg17468553",
        "cg18603700",
        "cg06175578") # The top CpGs that you want to plot (5 CpGs as an example here)

## Load Phenotype
pheno <- read.csv("/Volumes/som-ts/GynOb/SmithLab/NewData/People/Abby/PREG_TOLD_EWAS/QC_Results/PREG_TOLD_QC_FINAL_pheno_ENmix_passQC_wEstimates_AJB.csv")
  
## Recode phenotype as case/control
pheno[,VarName] <- ifelse(pheno[,VarName] == 0, "Control","Case")
pheno[,VarName] <- as.factor(pheno[,VarName])
table(pheno[,VarName])

## Load DNAm beta values - write the correct path and file name 
beta <- get(load("/Volumes/som-ts/GynOb/SmithLab/NewData/People/Abby/qc_directory/PREG_TOLD_QC/PREG_TOLD_QC_ENmix_qcd_combat_CP_wcovar_pheno_age_sex.RData"))
rm(combat_beta)
gc()

## Select CpGs and align
beta <- beta[cg, ]
beta <- t(beta)

## Merge
df <- merge(pheno, beta, by.x = "SampleID", by.y = "row.names")

## Plot
for (i in cg) {
  ggboxplot(df, x = VarName, y = i, notch = TRUE, fill = "gray", xlab = VarLabel) %>%
  ggexport(filename = paste0(VarName,"_",i,".png"), width = 600, height = 600, res = 150)
}


