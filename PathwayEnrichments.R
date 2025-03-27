# BiocManager::install("missMethyl") # install the package
library(missMethyl)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)


# Load summary statistics from EWAS and order based on p-values
beta<-read.csv("PREG_TOLD_QC_FINAL_Dystocia_Binary_EWAS_final_annotated.csv")
beta<-beta[order(beta$P.value),]

## Select top 500 CpGs
topCpGs<-beta[1:500,] 
range(topCpGs$P.value) # make sure the highest p < 0.05, if not try first 100 CpGs

## Get the names of selected CpGs
sigCpGs<-topCpGs$CPG.Labels

## GO pathway
gst<-gometh(sig.cpg = sigCpGs, collection = "GO", array.type = "EPIC", prior.prob = T, sig.genes = T)

# Subset to BP
gsa <- subset(gst, ONTOLOGY == "BP")
gsa$FDR <- p.adjust(gsa$P.DE, method = "fdr")
gsa <- gsa[order(gsa$P.DE),]
write.csv(gsa, file = "PREG_TOLD_QC_FINAL_Dystocia_Binary_top500_GO_BP.csv")

## KEGG
gst.kegg<-gometh_patch(sig.cpg = sigCpGs, collection = "KEGG", array.type = "EPIC", prior.prob = T, sig.genes = T)
gst.kegg <- gst.kegg[order(gst.kegg$P.DE),]
write.csv(gst.kegg,file = "PREG_TOLD_QC_FINAL_Dystocia_Binary_top500_KEGG.csv") 
