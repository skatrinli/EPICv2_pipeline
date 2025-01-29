################################################################################
## Seyma Katrinli
## Purpose: Smoking Score Estimation
################################################################################
#' remove all objects from your workspace
rm(list=ls())
gc()

# Load packages
library(impute)

################################################################################
# Define Variables and Load Data
PROJECTNAME <- "my_project_name"
OUTPUTPATH <- "/my/output/path/"
OTHERPATH <- "/my/otherfile/path/"

## Load beta values
BETAFILE_PASS <- get(load(paste0(OUTPUTPATH,PROJECTNAME,"_ENmix_qcd_combat_CP_wcovar_pheno_age_sex.Rdata")))

## Load the phenotype file
pheno <- read.csv(paste0(OUTPUTPATH, PROJECTNAME, "_pheno_ENmix_passQC_wEstimates.csv"))

## Load smoking score coefficients shared in the pipeline package
Scoefs <- read.csv(paste0(OTHERPATH, "AncillaryFiles/Smoking_probes_betas.csv"))
names(Scoefs)<-c("Marker","Coefficient")

################################################################################
#' This function will calculate the smoking score and combine it with the phenotype information
#' input: Combat adjusted beta values, phenotype file, smoking coefficients,
#' Output: Phenotye information with smoking scores
#'
smoking_score <- function(beta, pheno, Scoefs){
  # Changing beta values of 0 to 0.0001
  beta[beta < 0.0001] <- 0.0001
  beta[beta > 0.9999] <- 0.9999

  # Convert to M-vals
  # --------------------------------
  # Do we need to convert beta to m values here
  # because they are not used after conversion
  # SK: They are used for matrix multiplication (only the selected probes)
  # SK: Technically we can convert to M-values after selecting the 39 smoking probes (Sprobes2), but I kept Mark's order here.
  beta <- log2(beta/(1-beta))

  # How many probes are we using?
  Sprobes <- beta[which(row.names(beta) %in% as.character(Scoefs$Marker)),]
  message("Number of probes used (out of 39): ", nrow(Sprobes))
  message(paste0("Missing probes: "),
          paste0(Scoefs$Marker[!Scoefs$Marker%in%row.names(Sprobes)], collapse = ","))
  Scoefs2 <- Scoefs[as.character(Scoefs$Marker) %in% row.names(Sprobes),]
  Sprobes2 <- Sprobes[as.character(Scoefs2$Marker), ]
  stopifnot(all(rownames(Sprobes2)==Scoefs2$Marker))

  ## Since EWAS was computed for M values, we need to transform.
  ## SK: This line is from Mark, I just copied and pasted, but I think he refers to previous M-value transformation
  ## SK: We should delete this line, it's confusing :) 
  # ----------------------------------------- Not sure what you mean by the comment here
  Smo <- t(Scoefs2$Coefficient) %*% as.matrix(Sprobes2)
  Smo2 <- t(Smo)
  SmoScore <- data.frame(row.names(Smo2),Smo2) 
  names(SmoScore) <- c("ID", "SmoS") 

  # also removed all.x = True because by default merge will combine the ids that are common
  ## SK: The reason I keep all.x = T is to not lose phenotype information, and get NA for samples don't have methylation data
  ## SK: Because, I don't remove QC'ed out samples from my main phenotype document, just flag them; but remove them from methylation file.
  
  dat <- merge(pheno, SmoScore, by = 1)
  return(dat)
}
################################################################################
## Impute
beta.imputed <- impute.knn(as.matrix(BETAFILE_PASS))
beta.imputed <- beta.imputed$data
################################################################################
## Run the function to compute smoking score
dat <- smoking_score(beta = beta.imputed, pheno = pheno, Scoefs = Scoefs) # 2 probes missing out of 39
write.csv(dat, file = paste0(OUTPUTPATH, PROJECTNAME, "_pheno_ENmix_passQC_wEstimates.csv"), row.names = FALSE) # Save the phenotype file with smoking scores




