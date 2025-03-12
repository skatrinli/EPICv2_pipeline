################################################################################
## Seyma Katrinli Wise
## Purpose: Calculate methylation PCs
################################################################################

# Define Variables and Load Data
PROJECTNAME <- "my_project_name"
OUTPUTPATH <- "/my/output/path/"
OTHERPATH <- "/my/otherfile/path/"

### Modify this line to indicate correct pathname, if the annotation file that we've shared is in the working dir, the code will be:
file.pathname = paste0(OTHERPATH, "AncillaryFiles/cpgs_within_0bp_of_TGP_SNP.rdata")

## Load beta values (not ComBAT adjusted beta values, just normalized beta values)
BETAFILE_PASS <- get(load(paste0(OUTPUTPATH, PROJECTNAME, "_betas_without_suffix_passQC.RData")))

### Load and call location.based.pc() function that we've shared to compute principal components
### This function will select only CpG sites in the location-based annotation file,
### set missing values to the CpG-site-average, and compute principal components

source(paste0(OTHERPATH, "AncillaryFiles/location.based.pc.R")) # PATH TO SOURCE FILE, this is an example path, change accordingly
pc <- location.based.pc(BETAFILE_PASS,file.pathname)

### The returned result will be a princomp object.  The principal components will be available as pc$loadings

### Principal components can be incorporated as covariates in regression analysis
### Samples will be sorted in the same order as beta.obj, but if using with other data, make sure IDs match up!

top10pc <- pc$loadings[,1:10]

## You can merge the mPCs with phenotype file
# Load phenotype file: Samples as rows, the first column should be methylation IDs (same as column names of beta matrix)
pheno <- read.csv(paste0(OUTPUTPATH, PROJECTNAME, "_pheno_ENmix_passQC_wEstimates.csv"))
pheno <- merge(pheno,top10pc,by.x = 1,by.y = "row.names",all.x = T)
write.csv(pheno, file = paste0(OUTPUTPATH, PROJECTNAME, "_pheno_ENmix_passQC_wEstimates.csv"), row.names = FALSE) # Save the phenotype file with cell types


