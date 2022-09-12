library(data.table)
library(goseq)

effectNames = c('SA','Mex:High','SA:High','Mex:Lat_std','SA:Lat_std')

# set up dirs
DIR_Input = "../input"
DIR_Output = "../output"

all_categories = fread(file = '../input/all_categories.csv',data.table=F)
all_category_info = fread(file = '../input/all_category_info.csv',data.table=F)
types = unique(all_category_info$Type)

# results = list()

v_Totals = readRDS('../processed/v_Totals.rds')
averageExpression = readRDS('../processed/averageExpression.rds')

method = 'mash'
expression = 'DE'

for(expression in c('DE','DASE')) {  
for(method in c('mash','raw')) {
  
  if(method == 'mash') {
    results = readRDS(file.path(DIR_Output,sprintf('%s_mash_results.rds',expression)))
  } else {
    results = readRDS(file.path(DIR_Output,sprintf('%s_results_tissues.rds',expression)))
  }
  if(expression == 'DE') {
    averageExpression = readRDS('../processed/averageExpression.rds')
    tissues = c("MetLeaftip","MetLeafbase","PvLeaftip")
  } else {
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
  }
  # IDs = unique(unlist(lapply(voom_Elists,function(x) rownames(x$E))))
  # for(effect in effectNames) {
  #   for(i in names(results[[expression]][[effect]])) {
  #     stopifnot(all(rownames(results[[expression]][[effect]][[i]]) == IDs,na.rm=T))
  #     rownames(results[[expression]][[effect]][[i]]) = IDs
  #   }
  # }
  # if(method == 'mash') {
  #   saveRDS(results,file.path(DIR_Output,sprintf('%s_mash_results.rds',expression)))
  # } else {
  #   saveRDS(results,file.path(DIR_Output,sprintf('%s_results_tissues.rds',expression)))
  # }
  
    library(goseq)
    Enrichment_sign = list()
    Loc_Tiss = tissue = "MetLeaftip"
    # Loc_Tiss = "MetLeafbase"
    # Loc_Tiss = "PvLeaftip"
    Loc_Tiss = tissue = tissues[2]
    for(thresh in  c(0.05,0.1,0.2)) {
      for(type in types) {
        all_categories_type = subset(all_categories,ID %in% subset(all_category_info,Type == type)$ID)
        Enrichment_sign[[as.character(thresh)]][[type]] = list()
        for(Loc_Tiss in tissues){
          Enrichment_sign[[as.character(thresh)]][[type]][[Loc_Tiss]] = list()
          for(effect in c('SA','Mex:High','SA:High')[2:3]) {
            Enrichment_sign[[as.character(thresh)]][[type]][[Loc_Tiss]][[effect]] = list()
            for(sign in c('up','down')) {
              wrong_sign = ifelse(sign == 'up',1,-1) != sign(results[[expression]][[effect]]$PosteriorMeans[,Loc_Tiss])
              lfsr = results[[expression]][[effect]]$lfsrs[,Loc_Tiss]
              lfsr[wrong_sign]=1
              lfsr = na.omit(lfsr)
              categories_tissue = subset(all_categories_type,Gene %in% names(lfsr))
              map_cat2gene <- split(categories_tissue$Gene, f = categories_tissue$ID)
              map_cat2gene = map_cat2gene[lengths(map_cat2gene)>=10 & lengths(map_cat2gene) <= 1000]
              DE_genes = ifelse(lfsr<thresh,1,0)
              names(DE_genes) = names(lfsr)
              np = nullp(DE_genes,bias.data = log(averageExpression[[Loc_Tiss]][names(DE_genes)]))
              res = goseq(np,gene2cat = map_cat2gene,use_genes_without_cat = F)
              res$BH = p.adjust(res$over_represented_pvalue,method = 'BH')
              
              print(sprintf('%s %s %s %s %0.2f %d',type,Loc_Tiss,effect,sign, thresh,sum(res$BH < 0.05)))
              Enrichment_sign[[as.character(thresh)]][[type]][[Loc_Tiss]][[effect]][[sign]] = data.frame(Tissue = Loc_Tiss, Effect= effect, Sign = sign, res)
            }
          }
        }
      }
    }
    saveRDS(Enrichment_sign,file = file.path(DIR_Output,sprintf('%s_%s_GO_enrichment_sign.rds',expression,method)))
  }
}

