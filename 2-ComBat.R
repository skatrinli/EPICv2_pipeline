################################################################################
## Seyma Katrinli Wise
## Purpose: ComBAT Normalization to remove batch effects
################################################################################
## - Clean environ start.
rm(list = ls())
gc()
################################################################################
## Define the paths here.
OUTPUTPATH <- "/my/output/path/"
PROJECTNAME <- "my_project_name"
OTHERPATH <- "/my/otherfile/path/"

#' See the example sample sheet. Must have columns: 
#' SampleID - Methylation ID as SentrixID_SentrixPosition
#' Sample_Name - Study specific sample name
#' CONTROL - Control Samples, TRUE for control samples, FALSE for real samples
#' Sex - MALE or FEMALE
#' Age
#' Pheno - Phenotype of Interest, 1 for cases and 0 for controls

################################################################################
# path to source file required to install the packages
source(paste0(OTHERPATH, "AncillaryFiles/install_needed_packages.R")) # PATH TO SOURCE FILE, this is an example path, change accordingly

#' Lets install the required packages
bioc_packages <- c("BiocManager", "sva", "impute", "BiocParallel")

#' Call functions
install_bioconductor_pkgs(pkgs = bioc_packages)

#' Run test
check_installed(pkgs = c(bioc_packages))

#' Load all packages, if needed
lapply(c(bioc_packages), require, character.only = TRUE)

################################################################################
## Define Variables for Combat step
## Define model.matrix, which includes the variables that you want to protect when adjusting for chip and position
## Generally variables that we use as covariates in the EWAS (sex, age, main phenotype -PTSD-, smoking) are included in the model.matrix
# phen$Sentrix_ID <- substr(phen$SampleID,1,12)
# phen$Sentrix_Position <- substr(phen$SampleID,14,20)
sex <- "Sex"         # Name of Sex Variable
age <- "Age_Draw"         # Name of Age Variable
mainvar <- "Response"       # Variable of interest: Cases coded as 1, controls coded as 0
chip <- "Sentrix_ID"
position <- "Sentrix_Position"
################################################################################
# Saving the output
sink(paste0(OUTPUTPATH, PROJECTNAME, "_",Sys.Date(), "_ComBAT_Output.txt"), split = TRUE)
################################################################################
# Load data
beta <- get(load(paste0(OUTPUTPATH, PROJECTNAME, "_betas_without_suffix_passQC.RData")))
rm(BETAFILE_PASS)
gc()

pheno <- read.csv(paste0(OUTPUTPATH, PROJECTNAME, "_pheno_ENmix_passQC.csv"), row.names = 1)

# Methylation IDs (Row names of phenotype file and Column names of beta matrix) should match
phen <- pheno[row.names(pheno) %in% colnames(beta), ]
beta <- beta[, which(colnames(beta) %in% row.names(pheno))]
beta <- beta[, order(colnames(beta))]
phen <- phen[order(row.names(phen)), ]
stopifnot(all(colnames(beta) == rownames(phen)))

########################################################################################################
############# ComBAT normalization to adjust for batch effects (chip and array position) ############
########################################################################################################
print("Running ComBat")

## You should not have NAs in model matrix, so we remove subjects with no phenotype info
print(paste0("Samples with no Main Variable information = ", sum(is.na(phen[,mainvar])))) # 0
print(paste0("Samples with no Sex information = ", sum(is.na(phen[,sex])))) # 0
print(paste0("Samples with no Age information = ", sum(is.na(phen[,age])))) # 0

phen <- phen[,c(chip,position,mainvar,age,sex)]
phen <- phen[complete.cases(phen), ]
beta <- beta[,colnames(beta)%in%rownames(phen)]
stopifnot(all(colnames(beta) == rownames(phen)))

chip <- as.factor(phen[,chip])
position <- as.factor(phen[,position])
mainvar <- as.factor(phen[,mainvar])
sex <- as.factor(phen[,sex])
age <- as.numeric(phen[,age])

## Define model matrix
if (length(levels(sex)) == 1) {
  moddata <- model.matrix(~mainvar+age)
} else {
  moddata <- model.matrix(~mainvar+age+sex)
}

# Remaining samples
print(paste0("Remaining Samples = ", nrow(phen))) 

## ComBAT does not handle NAs in the methylation file,
## So we have to impute NAs in the methylation beta matrix
beta.log <- log((beta/(1-beta)))
beta.imputed <- impute.knn(as.matrix(beta.log))
beta.imputed <- beta.imputed$data

rm(beta.log)
gc()

# Run ComBat
SerialParam()
combat_beta <- ComBat(dat = beta.imputed, mod = moddata, batch = chip, BPPARAM = SerialParam())
combat_beta <- ComBat(dat = combat_beta, mod = moddata, batch = position, BPPARAM = SerialParam())

# Reverse Beta Values
combat_beta <- 1/(1+(1/exp(combat_beta)))

## We need to put NAs back (from the original matrix) to ComBAT adjusted beta matrix
## We don't want to use imputed beta values for missing data
beta <- beta[, colnames(beta) %in% colnames(combat_beta)]
beta <- beta[, order(colnames(beta))]
combat_beta <- combat_beta[, order(colnames(combat_beta))]
table(colnames(beta) == colnames(combat_beta))

combat_beta[is.na(beta)] <- NA
range(combat_beta, na.rm = T) # 1.860516e-05 0.9991135
combat_beta[combat_beta >= 0.9999998] <- 0.9999999
save(combat_beta,file=paste0(OUTPUTPATH,PROJECTNAME,"_ENmix_qcd_combat_CP_wcovar_pheno_age_sex.Rdata"))

################################################################################
sessionInfo()

print("Script complete. Exiting.")

################################################################################
# Close sink
sink()

