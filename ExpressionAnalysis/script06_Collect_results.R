#==================================
# make tables of results
#==================================

# 2022-05-26
# corrected one line DER code bug for calculating total DE genes

# Clean workspace
rm(list=ls())

library(GO.db)
library(AnnotationDbi)
library("data.table")
library(data.tree) #plot list structure

DIR_Output = "../output/"
method = 'mash'

# -------------------------------------------------------------------------
# SECTION 01: Collect DE/DASE individual genes + gene set test results
# -------------------------------------------------------------------------

# -----------------------------------------------------------
# DE 
# - collect single gene test results for total expression
# - collect GSVA results for total expression
expression = 'DE'
if(method == 'mash') {
  results = readRDS(file.path(DIR_Output,sprintf('%s_mash_results.rds',expression)))
} else {
  results = readRDS(file.path(DIR_Output,sprintf('%s_results_tissues.rds',expression)))
}
FromListSimple(results)

effectNames = names(results[[1]])
traits = names(results) #"DE"   "GSVA"
tissues = colnames(results[[1]][[1]]$PosteriorMeans)

results_table = list()
for(trait in traits) {
  IDs = data.frame(ID = rownames(results[[trait]][[1]]$PosteriorMeans),stringsAsFactors = F)
  if(trait %in% c('GO','GSVA')) IDs = select(GO.db,IDs$ID,c('GOID','ONTOLOGY','TERM'))   
  results_table[[trait]] = data.frame(IDs,
                              do.call(cbind,lapply(effectNames,function(effect) {
                                do.call(cbind,lapply(tissues,function(tissue) {
                                  table = data.frame(results[[trait]][[effect]]$PosteriorMeans[,tissue],
                                                     results[[trait]][[effect]]$lfsrs[,tissue])
                                  colnames(table) = sprintf('%s_%s_%s',effect,tissue,c('effect','lfsr'))
                                  table
                                }))
                              })),check.names=F)
  rownames(results_table[[trait]]) = NULL
  write.csv(results_table[[trait]],file = file.path(DIR_Output,sprintf('%s_results_table_%s.csv',trait,method)))
}
# Filenams in DIR_Output:
# - DE_results_table_mash.csv
# - GSVA_results_table_mash.csv

# --------------------------------------------------------------
# DASE 
# - collect single gene test results for ASE
# - DER: no significant GSVA results for ASE

expression = 'DASE'
if(method == 'mash') {
  results = readRDS(file.path(DIR_Output,sprintf('%s_mash_results.rds',expression)))
} else {
  results = readRDS(file.path(DIR_Output,sprintf('%s_results_tissues.rds',expression)))
}

effectNames = names(results[[1]])
traits = names(results)
averageExpression = readRDS('../processed/averageASECounts.rds')
IDs = unique(unlist(lapply(averageExpression,names)))
averageExpression$Joint = rowMeans(sapply(averageExpression,function(x) x[IDs]),na.rm=T)
names(averageExpression$Joint) = IDs
tissues = c("Joint")
for(effect in effectNames) {
  results$DASE[[effect]] = list(
    PosteriorMeans = as.matrix(data.frame(Joint = rowMeans(results$DASE[[effect]]$PosteriorMeans,na.rm=T))),
    lfsrs = as.matrix(data.frame(Joint = apply(results$DASE[[effect]]$lfsrs,1,min,na.rm=T)))
  )
}

trait = 'DASE'
IDs = data.frame(ID = rownames(results[[trait]][[1]]$PosteriorMeans))
results_table[[trait]] = data.frame(IDs,
                                    do.call(cbind,lapply(effectNames,function(effect) {
                                      do.call(cbind,lapply(tissues,function(tissue) {
                                        table = data.frame(results[[trait]][[effect]]$PosteriorMeans[,tissue],
                                                           results[[trait]][[effect]]$lfsrs[,tissue])
                                        colnames(table) = sprintf('%s_%s_%s',effect,tissue,c('effect','lfsr'))
                                        table
                                      }))
                                    })),check.names=F)
rownames(results_table[[trait]]) = NULL
write.csv(results_table[[trait]],file = file.path(DIR_Output,sprintf('%s_results_table_%s.csv',trait,method)))
# Filenams in DIR_Output:
# - DASE_results_table_mash.csv
# - ASE has no GSVA results

#--------------------------------------------------------------------
# traits significant by mash
# print out how many traits significant for above DE and DASE results
get_significant_results_NA = function(m,...){
  m$result$lfsr[is.na(m$result$lfsr)] = 1
  get_significant_results(m,...)
}

sig_table = function(table,trait) {
  tissues = c('MetLeaftip','MetLeafbase','PvLeaftip','Joint')
  tissues = tissues[sapply(tissues,function(tissue) length(grep(tissue,colnames(table)))>0)]
  sig_table = sapply(tissues,function(tissue){
    sig_table = table[,sprintf('%s_%s_%s',effectNames,tissue,'lfsr')] < thresh
    colnames(sig_table) = effectNames
    sum_sig = colSums(sig_table,na.rm=T)
    sum_sig['Total'] = sum(!is.na(sig_table[,1]))
    sum_sig
  })
  Total = sapply(effectNames,function(effect) {
    sig_table = table[,sprintf('%s_%s_%s',effect,tissues,'lfsr'),drop=F] < thresh
    colnames(sig_table) = tissues
    # sum(rowSums(sig_table)>0,na.rm=T) #DER old code 04/24/22
    sum(rowSums(sig_table,na.rm = T)>0) #HH corrected code 05/26/22
  })
  Total['Total'] = nrow(table)
  rbind(t(sig_table),Total)
}

thresh = 0.05
sig_tables = list()
for(trait in names(results_table)) { #"DE"   "GSVA" "DASE"
  print(trait)
  print(sig_table(results_table[[trait]],trait))
}
# [1] "DE"
# SA Mex:High SA:High Mex:Lat_std SA:Lat_std Total
# MetLeaftip  114     1278     715          50         60 18369
# MetLeafbase 106     3716    1626          22         50 20401
# PvLeaftip   107      319     368          10        122 18079
# Total       124     4432    1816          60        131 21599
# [1] "GSVA"
# SA Mex:High SA:High Mex:Lat_std SA:Lat_std Total
# MetLeaftip   0       14     167           1          0  1455
# MetLeafbase  0      303     185          44          0  1403
# PvLeaftip   19       11       1           0          0  1455
# Total       19      303     216          45          0  1455
# [1] "DASE"
#        SA Mex:High SA:High Mex:Lat_std SA:Lat_std Total
# Joint 249      341     260          17         23 13632
# Total 249      341     260          17         23 13632

dim(results_table$DE)
head(results_table$DE)
tail(results_table$DE)
sum(results_table$DE$`Mex:High_MetLeaftip_lfsr`<thresh,na.rm = T)
length(results_table$DE$ID[results_table$DE$`Mex:High_MetLeaftip_lfsr`<thresh])

#--------------------------------------------------------------------------
# GO by Population (GO Enrichment test by goseq for each Tissue/DB separately)
# goseq analysis was applied to total expression data 
#  - for each Loc:Tiss separately
#  - for each geneset DB separately
#  - for up/down regulated genes separately
#  - see the list structure visualized by the FromListSimple() function below

expression = 'DE'
thresh = 0.1
Enrichment_results = readRDS(file = file.path(DIR_Output,sprintf('%s_%s_GO_enrichment_sign.rds',expression,method)))
# view list structure
FromListSimple(Enrichment_results)
FromListSimple(Enrichment_results$`0.1`)

types = names(Enrichment_results[[as.character(thresh)]])
tissues = names(Enrichment_results[[1]][[1]])
effectNames = names(Enrichment_results[[1]][[1]][[1]])
dim(Enrichment_results$`0.05`$GO$MetLeaftip$`Mex:High`$up)

for(type in types) {
  table = c()
  for(effect in effectNames) {
    for(tissue in tissues) {
      for(sign in c('up','down')) {
        table = rbind(table,Enrichment_results[[as.character(thresh)]][[type]][[tissue]][[effect]][[sign]])
      }
    }
  }
  results_table[[type]] = table
  write.csv(results_table[[type]],file = file.path(DIR_Output,sprintf('Enrichment_%s_results_table_%s.csv',type,method,thresh)))
}
sapply(results_table[types],function(table) table(subset(table,BH<thresh)$Tissue))
# $GO
# MetLeafbase  MetLeaftip 
# 1331          54 

# $CornCyc
# MetLeafbase  MetLeaftip   PvLeaftip 
# 12           3           1 

# $KEGG
# MetLeafbase  MetLeaftip   PvLeaftip 
# 40          20           1 

# $MaizeGRN
# MetLeafbase  MetLeaftip   PvLeaftip 
# 239         153          43 

# Filenams in DIR_Output:
# - Enrichment_CornCyc_results_table_mash.csv  
# - Enrichment_GO_results_table_mash.csv  
# - Enrichment_KEGG_results_table_mash.csv  
# - Enrichment_MaizeGRN_results_table_mash.csv

# summary of Results in SECTION 01:
A <- list()
A[["totalExpr"]]=list(); A[["ASE"]]=list()
A[["totalExpr"]][["SingleGeneAnalysis"]]="DE_results_table_mash.csv"
A[["totalExpr"]][["GeneSetAnalysis"]]=c("GSVA_results_table_mash.csv",
                                           "Enrichment_CornCyc_results_table_mash.csv",
                                           "Enrichment_GO_results_table_mash.csv",
                                           "Enrichment_KEGG_results_table_mash.csv",
                                           "Enrichment_MaizeGRN_results_table_mash.csv")
A[["ASE"]][["SingleGeneAnalysis"]]="DASE_results_table_mash.csv"
A[["ASE"]][["GeneSetAnalysi"]]=NA
A
FromListSimple(A)
str(A)

# --------------------------------------------------------
# SECTION 02: Collect DE/DASE Convergent results
# --------------------------------------------------------

#------------------------------------------------------------
# Convergent lists
all_categories_info = fread('../input/all_category_info.csv')
expression = 'DE'
convergent_lists = readRDS(file = file.path(DIR_Output,sprintf('%s_%s_convergent_lists.rds',expression,method)))
names(convergent_lists) #"0.05" "0.1"  "0.2" 
names(convergent_lists$`0.05`) #"DE"         "GSVA"       "Enrichment"
FromListSimple(convergent_lists)
convergent_lists$`0.05`$GSVA

# DE - convergent genes for each Loc:Tiss
names(convergent_lists$`0.05`$DE) #"MetLeaftip"  "MetLeafbase" "PvLeaftip" 
thresh = 0.05
tissues = names(convergent_lists[[as.character(thresh)]]$DE)
convergent_table = c()
for(tissue in tissues) {
  try({convergent_table = rbind(convergent_table,data.frame(Tissue = tissue,Gene = convergent_lists[[as.character(thresh)]]$DE[[tissue]]$sig))})
}
write.csv(convergent_table,file = file.path(DIR_Output,sprintf('%s_Convergent_genes_%s_%0.2f.csv',expression,method,thresh)))
dim(convergent_table) #567   2
head(convergent_table)
table(convergent_table$Tissue)
# MetLeaftip MetLeafbase   PvLeaftip 
# 126         411          30 

# -----------------------------------------------
# number of common gene IDs between Mex and SA
numComID=function(table,tissue1,tissue2){
  t1=results_table$DE$ID[results_table$DE[,tissue1]<thresh & !is.na(results_table$DE[,tissue1])]
  t2=results_table$DE$ID[results_table$DE[,tissue2]<thresh & !is.na(results_table$DE[,tissue2])]
  return(length(intersect(t1,t2)))
}
numComID(results_table$DE,"Mex:High_MetLeaftip_lfsr","SA:High_MetLeaftip_lfsr") #131
numComID(results_table$DE,"Mex:High_MetLeafbase_lfsr","SA:High_MetLeafbase_lfsr") #429
numComID(results_table$DE,"Mex:High_PvLeaftip_lfsr","SA:High_PvLeaftip_lfsr" )    #30

# per-tissue overlap Testing (adapted from Detect_convergence.R)
test_overlap = function(lfsr_1,lfsr_2,thresh) {
  lfsr_1=lfsr_1[!is.na(lfsr_1)]
  lfsr_2=lfsr_2[!is.na(lfsr_2)]
  stat = sprintf('%d of %d',sum(lfsr_1<thresh & lfsr_2 < thresh),length(lfsr_1))
  p = phyper(sum(lfsr_1<thresh & lfsr_2 < thresh)-1,
             sum(lfsr_1<thresh),
             length(lfsr_1)-sum(lfsr_1<thresh),
             sum(lfsr_2<thresh),lower.tail=F)
  c(stat,signif(p,3))
}
tissue1="Mex:High_MetLeaftip_lfsr"; tissue2="SA:High_MetLeaftip_lfsr"
test_overlap(results_table$DE[,tissue1],results_table$DE[,tissue2],
             thresh = 0.05) #"131 of 18369" "2.74e-25"
tissue1="Mex:High_MetLeafbase_lfsr"; tissue2="SA:High_MetLeafbase_lfsr"
test_overlap(results_table$DE[,tissue1],results_table$DE[,tissue2],
             thresh = 0.05) #"429 of 20401" "1.12e-17" 
tissue1="Mex:High_PvLeaftip_lfsr"; tissue2="SA:High_PvLeaftip_lfsr"
test_overlap(results_table$DE[,tissue1],results_table$DE[,tissue2],
             thresh = 0.05) #"30 of 18079" "3.28e-12"

# GSVA - total expression
names(convergent_lists$`0.05`$GSVA) #"MetLeaftip"  "MetLeafbase" "PvLeaftip" 
thresh = 0.1
tissues = names(convergent_lists[[as.character(thresh)]]$GSVA) ##"MetLeaftip"  "MetLeafbase" "PvLeaftip" 
convergent_table = c()
(tissue = tissues[2])
for(tissue in tissues) {
  try({convergent_table = rbind(convergent_table,data.frame(Tissue = tissue,ID = convergent_lists[[as.character(thresh)]]$GSVA[[tissue]]$sig))})
}
head(convergent_table)
dim(convergent_table)

convergent_table = data.frame(convergent_table,all_categories_info[match(convergent_table$ID,all_categories_info$ID),])
types = unique(convergent_table$Type)
convergent_tables = lapply(types,function(type) {
  table = subset(convergent_table,Type == type)
  table = table[,colSums(!is.na(table))>0]
  write.csv(convergent_table,file = file.path(DIR_Output,sprintf('%s_Convergent_GSVA_%s_%s_%0.2f.csv',expression,type,method,thresh)))
  table
})
sapply(convergent_tables,function(table) table(table$Tissue))
unique(convergent_table$Tissue)

# Enrichment (i.e. geseq results)
names(convergent_lists$`0.05`$Enrichment) #"GO"       "CornCyc"  "KEGG"     "MaizeGRN"
sapply(convergent_lists$`0.05`$Enrichment,function(table) identical(colnames(table), colnames(results_table$GO)))
convergent_lists$`0.05`$Enrichment$GO
convergent_lists$`0.1`$Enrichment$CornCyc
convergent_lists$`0.05`$Enrichment$KEGG

thresh = 0.1
types = names(convergent_lists[[as.character(thresh)]]$Enrichment) #"GO"       "CornCyc"  "KEGG"     "MaizeGRN"
for(type in types) {
  table = do.call(rbind,lapply(tissues,function(tissue) convergent_lists[[as.character(thresh)]]$Enrichment[[type]][[tissue]]))
  if(!is.null(table)) {
    # break
    print(sprintf('%s %d',type,nrow(table)))
    table = table[,colSums(!is.na(table))>0]
    write.csv(convergent_table,file = file.path(DIR_Output,sprintf('%s_Convergent_Enrichment_%s_%s_%0.2f.csv',expression,type,method,thresh)))
  }
}
dim(convergent_table) #567   2
head(convergent_table)

# DASE - convergent genes
thresh = 0.05
expression = 'DASE'
convergent_lists = readRDS(file = file.path(DIR_Output,sprintf('%s_%s_convergent_lists.rds',expression,method)))
tissues = names(convergent_lists[[as.character(thresh)]]$DASE)
convergent_table = c()
for(tissue in tissues) {
  try({convergent_table = rbind(convergent_table,data.frame(Tissue = tissue,Gene = convergent_lists[[as.character(thresh)]]$DASE[[tissue]]$sig))})
}
table(convergent_table$Tissue)
write.csv(convergent_table,file = file.path(DIR_Output,sprintf('%s_Convergent_genes_%s_%0.2f.csv',expression,method,thresh)))

# Enrichment of convergent genes
expression = 'DE'
GO_convergent = readRDS(file = file.path(DIR_Output,sprintf('%s_%s_GO_of_convergentGenes.rds',expression,method)))
FromListSimple(GO_convergent)

tissues = names(GO_convergent[[as.character(thresh)]]) # "MetLeaftip"  "MetLeafbase" "PvLeaftip"  
convergent_table = c()
thresh = 0.05
for(tissue in tissues) {
  try({convergent_table = rbind(convergent_table,data.frame(Tissue = tissue,GO_convergent[[as.character(thresh)]][[tissue]]))})
}
head(convergent_table)
FDR=0.1
convergent_table = subset(convergent_table,BH < FDR)
table(convergent_table$Tissue)

write.csv(convergent_table,file = file.path(DIR_Output,sprintf('%s_Enrichment_of_Convergent_genes_%s_%0.2f.csv',expression,method,thresh)))

names(convergent_lists$`0.05`$Enrichment$GO)
