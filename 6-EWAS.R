#########################################################################
## Seyma Katrinli Wise
## Purpose: EWAS
#########################################################################
#' Clean
rm(list=ls())
gc()

#' Load packages
install.packages("CpGassoc")
library(CpGassoc)
library(data.table)
#########################################################################
## Define the paths here.
OUTPUTPATH <- "/my/output/path/"
PROJECTNAME <- "my_project_name"
OTHERPATH <- "/my/otherfile/path/"

## Load Methylation Data
beta <- get(load(paste0(OUTPUTPATH,PROJECTNAME,"_ENmix_qcd_combat_CP_wcovar_pheno_age_sex.Rdata")))

## Load phenotype file with MethylationID (combination of Sentrix ID and Sentrix Position) as the 1st column
pheno <- read.csv(paste0(OUTPUTPATH, PROJECTNAME, "_pheno_ENmix_passQC_wEstimates.csv"), row.names = 1)

## Define variables
mainVar <- "PTSDpm" # name of the ptsd variable, coded as: cases = 1 and controls = 0

## Define covariates to be adjusted for EWAS
## Covariates to be included:
##  - cell types from step 3 ("CD8T","CD4T","NK","Bcell","Mono")
##  - GWAS PCs PC1 and PC2 (if available), if not mPC2 (Comp.2) and mPC3 (Comp.3) from step 3.1
##  - age
##  - sex (if applicable)
covar <- data.frame(pheno[,c("Comp.2","Comp.3","CD8T","CD4T","NK","Bcell","Mono","female","Age")])

#########################################################################

# A function here to get and order the required data
clean_order <- function(beta, pheno){
  cpg <- beta[, colnames(beta) %in% row.names(pheno)]
  cpg <- cpg[, order(colnames(cpg))]
  pheno <- pheno[rownames(pheno) %in% colnames(cpg), ]
  pheno <- pheno[order(rownames(pheno)), ]
  print(table(rownames(pheno) == colnames(cpg))) # should be TRUE
  stopifnot(all(rownames(pheno) == colnames(cpg)))
  return(list(pheno = pheno, cpg = cpg))
}

cleaned_df <- clean_order(beta = beta, pheno = pheno)
pheno <- cleaned_df$pheno


# Function to run EWAS with CpGAssoc
cpg_assoc_test <- function(cpg, pheno, covar){
  message("Running test, patience ...")
  test <- cpg.assoc(cpg, pheno, covar, logit.transform = T, large.data=TRUE)
  assoc <- test$results
  eff <- test$coefficients
  results <- cbind(assoc, eff)
}

# Run EWAS with CpGAssoc
results <- cpg_assoc_test(cpg = cleaned_df$cpg,
                          pheno = cleaned_df$pheno[, mainVar],
                          covar = covar)
#########################################################################
# Annotate
epicv2manifestfilename = paste0(OTHERPATH, "pidsley2024.csv")
annot = data.table::fread(epicv2manifestfilename, header=T,
                               stringsAsFactors = F, sep=",", fill=T)

annot <- annot[,c("Name","CpG_chrm","CpG_beg","UCSC_RefGene_Name")]

# Merge ewas with annot, keeping all rows from ewas
ewas_merged <- merge(results, annot, by = 1)
write.csv(ewas_merged, file = paste0(OUTPUTPATH, PROJECTNAME,"_",mainVar,"_EWAS.csv"), row.names = FALSE) # Save the phenotype file with cell types
