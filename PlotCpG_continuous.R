################################################################################
## Purpose: Plotting CpGs for Continuous Variables
################################################################################
## Install and Load packages
install.packages("ggpubr")
library(ggpubr)

## Define variables
VarName <- "CTQ_total_score" # the name of the variable in the phenotype file
VarLabel <- "CTQ total score" # how you want the variable name to be shown in the x-axis of the plot
cg <- c("cg26188218",
        "cg20399509",
        "cg12571007",
        "cg22504528",
        "cg03648500") # The top CpGs that you want to plot (5 CpGs as an example here)

## Load Phenotype
pheno <- read.csv("/Volumes/som-ts/GynOb/SmithLab/NewData/People/Abby/PREG_TOLD_EWAS/QC_Results/PREG_TOLD_QC_FINAL_pheno_ENmix_passQC_wEstimates_AJB.csv")

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
  ggscatter(df, x = VarName, y = i, xlab = VarLabel, 
            shape = 20, add = "reg.line", conf.int = TRUE, add.params = list(color = "blue", fill = "lightgray")) %>%
    ggexport(filename = paste0(VarName,"_",i,".png"), width = 600, height = 600, res = 150)
}


