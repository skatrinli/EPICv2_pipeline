################################################################################
## Seyma Katrinli Wise
## Purpose: Cell Proportion Estimation
################################################################################
#' remove all objects from your workspace
rm(list=ls())
gc()

################################################################################
# Define Variables and Load Data
PROJECTNAME <- "my_project_name"
OUTPUTPATH <- "/my/output/path/"
OTHERPATH <- "/my/otherfile/path/"

## Load beta values (not ComBAT adjusted beta values, just normalized beta values)
load(paste0(OUTPUTPATH, PROJECTNAME, "_betas_without_suffix_passQC.RData"))

## Load ref dataset shared in the pipeline package
load(paste0(OTHERPATH, "AncillaryFiles/IDOLOptimizedCpGs_REF.Rdata"))

## Load the phenotype file
pheno <- read.csv(paste0(OUTPUTPATH, PROJECTNAME, "_pheno_ENmix_passQC.csv"), row.names = 1)

################################################################################
# path to source file required to install the packages
source(paste0(OTHERPATH, "AncillaryFiles/install_needed_packages.R")) # PATH TO SOURCE FILE, this is an example path, change accordingly

#' Install packages if not already installed
need_pkgs <- c("FlowSorted.Blood.EPIC", "EpiDISH")
install_bioconductor_pkgs(pkgs = need_pkgs)

#' Run test
check_installed(pkgs = need_pkgs)

#' Load all packages, if needed
load_pkgs <- c("FlowSorted.Blood.EPIC", "EpiDISH")
lapply(load_pkgs, require, character.only = TRUE)

#' This function will calculate the cell proportions and combine it with the phenotype information
#' input: beta values before Combat, phenotype file, reference cell proportions,
#' xcol- column name on which the data should be merged
#' ycol- column name
#' Output: Phenotye information with cell types
estimate_cell_proportion <- function(betavals, phenotyps, cell_prop_ref, xcol, ycol ){
  
  ## Calculate cell types using RPC method
  RPC <- epidish(betavals, as.matrix(cell_prop_ref), method = "RPC")
  
  cellTypes <- as.data.frame(RPC$estF) #RPC count estimates
  
  phen <- merge(phenotyps, cellTypes, by.x = xcol, by.y = ycol, all.x = T)
  
  return(phen)
}

## Run the function
output <- estimate_cell_proportion(betavals = BETAFILE_PASS, phenotyps = pheno, cell_prop_ref = ref, xcol = "SampleID", ycol = 'row.names' )
write.csv(output, file = paste0(OUTPUTPATH, PROJECTNAME, "_pheno_ENmix_passQC_wEstimates.csv"), row.names = FALSE) # Save the phenotype file with cell types


