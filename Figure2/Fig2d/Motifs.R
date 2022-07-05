##loading libraries
library(RcisTarget)
library(DT)

###loading genes in list form
geneList1= read.table("Module1.txt")[,1]
geneLists <- list(geneListName=geneList1)


# Select motif database to use which can be downloaded from: https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/
data(motifAnnotations_hgnc)
motifRankings <- importRankings("hg19-tss-centered-10kb-7species.mc9nr.feather")

# Motif enrichment analysis
motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                         motifAnnot=motifAnnotations_hgnc)
motifs_AUC <- calcAUC(geneLists, motifRankings)

#Select significant motifs
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, 
                                           motifAnnot=motifAnnotations_hgnc)


# Identify significant genes for each motif
motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable, 
                                                   geneSets=geneLists,
                                                   rankings=motifRankings, 
                                                   nCores=1,
                                                   method="aprox")




motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)
resultsSubset <- motifEnrichmentTable_wGenes_wLogo

datatable(resultsSubset[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], 
          escape = FALSE, # To show the logo
          filter="top", options=list(pageLength=5))
signifMotifNames <- motifEnrichmentTable$motif[1:20]
