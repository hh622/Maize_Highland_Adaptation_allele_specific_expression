library(goseq)
test_overlap = function(lfsr_1,lfsr_2,thresh) {
  stat = sprintf('%d of %d',sum(lfsr_1<thresh & lfsr_2 < thresh),length(lfsr_1))
  p = phyper(sum(lfsr_1<thresh & lfsr_2 < thresh)-1,
             sum(lfsr_1<thresh),
             length(lfsr_1)-sum(lfsr_1<thresh),
             sum(lfsr_2<thresh),lower.tail=F)
  c(stat,signif(p,3))
}



expression='DE'
method = 'mash'
DIR_Output="../output/"
all_category_info = fread('../input/all_category_info.csv') #HH added on 4/26
all_categories = fread(file = '../input/all_categories.csv',data.table=F) #HH added on 4/26

expression ='DE'; method ='mash' #for testing
for(expression in c('DE','DASE')) {  
  for(method in c('mash','raw')) {
    print(paste(expression,method))
    
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
    traits = names(results)
    # comparison ='High'
    for(comparison in c('High')) {
      print(comparison)
      comparison_effects = paste(c('Mex','SA'),comparison,sep=':')
      
      # Hypergeometric test for # overlaps
      # HH 2022-04-30: this performs across-tissue/joint overlapping test
      for(thresh in c(0.05,0.1,0.2)) {
        print(thresh)
        table = sapply(traits,function(trait) {
          test_overlap(apply(results[[trait]][[comparison_effects[1]]]$lfsrs,1,min,na.rm=T),
                       apply(results[[trait]][[comparison_effects[2]]]$lfsrs,1,min,na.rm=T),
                       thresh)
        })
        print(table)
      }
      # length(apply(results[[trait]][[comparison_effects[1]]]$lfsrs,1,min,na.rm=T))
      # length(apply(results[[trait]][[comparison_effects[2]]]$lfsrs,1,min,na.rm=T))
      # sum(is.na(apply(results[[trait]][[comparison_effects[2]]]$lfsrs,1,min,na.rm=T)))
      # sum(is.na(apply(results[[trait]][[comparison_effects[1]]]$lfsrs,1,min,na.rm=T)))
      
      # per-tissue correlation of "Convergent" genes
      for(trait in traits) {
        print(trait)
        for(thresh in c(0.05,0.1,0.2)) {
          print(thresh)
          stats=sapply(tissues,function(Loc_Tiss) {
            lfsrs = cbind(results[[trait]][[comparison_effects[1]]]$lfsrs[,Loc_Tiss],results[[trait]][[comparison_effects[2]]]$lfsrs[,Loc_Tiss])
            PosteriorMeans = cbind(results[[trait]][[comparison_effects[1]]]$PosteriorMeans[,Loc_Tiss],results[[trait]][[comparison_effects[2]]]$PosteriorMeans[,Loc_Tiss])
            sig = apply(lfsrs,1,max) < thresh & rowSums(is.na(lfsrs))==0
            # sig = na.omit(sig)
            if(sum(sig)<2) return(c(Correlation=NA,Signs=NA))
            c(Correlation=cor(PosteriorMeans[sig,],use='p')[2],Signs=mean(sign(PosteriorMeans[sig,1])==sign(PosteriorMeans[sig,2])))
          })
          print(stats)
        }
      }
      
      # get lists of convergent genes
      
      # find overlapping traits between Mex:High and SA:High
      convergent_lists = list()
      for(thresh in c(.05,.1,.2)) {
        print(thresh)
        convergent_lists[[as.character(thresh)]] = list()
        for(trait in traits) {
          convergent_lists[[as.character(thresh)]][[trait]] = list()
          for(Loc_Tiss in tissues) {
            lfsrs = cbind(results[[trait]][[comparison_effects[1]]]$lfsrs[,Loc_Tiss],results[[trait]][[comparison_effects[2]]]$lfsrs[,Loc_Tiss])
            PosteriorMeans = cbind(results[[trait]][[comparison_effects[1]]]$PosteriorMeans[,Loc_Tiss],results[[trait]][[comparison_effects[2]]]$PosteriorMeans[,Loc_Tiss])
            max_lfsrs = apply(lfsrs,1,max)
            lfsrs = lfsrs[order(max_lfsrs),]
            PosteriorMeans = PosteriorMeans[order(max_lfsrs),]
            conv = apply(lfsrs,1,max)<thresh & sign(PosteriorMeans[,1])==sign(PosteriorMeans[,2]) 
            # HH notes on 2022-04-30: sign checking between Mex and SA grantees regulation direction the same in highland
            convergent_lists[[as.character(thresh)]][[trait]][[Loc_Tiss]]$sig = names(conv[!is.na(conv) & conv])
            convergent_lists[[as.character(thresh)]][[trait]][[Loc_Tiss]]$background = names(conv[!is.na(conv)])
            print(sprintf('%s %s %d',trait,Loc_Tiss,sum(conv,na.rm=T)))
          }
        }
      }
      
      # convergent Enrichment
      Enrichment_sign = readRDS(file.path(DIR_Output,sprintf('%s_%s_GO_enrichment_sign.rds',expression,method)))
      # for testing
      # Enrichment_sign = readRDS(file.path("../output/",sprintf('%s_%s_GO_enrichment_sign.rds',expression,method)))
      # FromListSimple(Enrichment_sign); Loc_Tiss = tissues[1]
      types = names(Enrichment_sign[[1]])
      for(thresh in c(.05,.1,.2)) {
        print(thresh)
        convergent_lists[[as.character(thresh)]]$Enrichment = list()
        for(type in types) {
          convergent_lists[[as.character(thresh)]]$Enrichment[[type]] = list()
            for(Loc_Tiss in tissues) {
              convergent = c()
              for(sign in c('up','down')) {
                res_Mex = Enrichment_sign[[thresh]][[type]][[Loc_Tiss]][[comparison_effects[1]]][[sign]]
                res_SA = Enrichment_sign[[thresh]][[type]][[Loc_Tiss]][[comparison_effects[2]]][[sign]]
                #key code for Enrichment convergent
                categories = intersect(subset(res_Mex,BH<thresh)$category,subset(res_SA,BH<thresh)$category) 
                # because goseq enrichment was tested for up/down regulated genes separately, so can simply do intersection set 
                if(length(categories) == 0) {
                  print(sprintf('%s %s %d',type, Loc_Tiss,0))
                  next
                }
                convergent = rbind(convergent,data.frame(Tissue = Loc_Tiss, Sign = sign,
                                                         all_category_info[match(categories,all_category_info$ID),]))
              }
              convergent_lists[[as.character(thresh)]]$Enrichment[[type]][[Loc_Tiss]] = convergent
              print(sprintf('%s %s %d',type, Loc_Tiss,nrow(convergent)))
          }
        }
      }
      
      saveRDS(convergent_lists,file = file.path(DIR_Output,sprintf('%s_%s_convergent_lists.rds',expression,method)))
      
      # Enrichment of convergent genes
      GO_convergent = list()
      for(thresh in c(.05,.1,.2)){
        GO_convergent[[as.character(thresh)]] = list()
        for(tissue in tissues) {
          all_genes = unique(unlist(convergent_lists[[as.character(thresh)]][[expression]][[tissue]]))
          
          categories_tissue = subset(all_categories,Gene %in% all_genes)
          map_cat2gene <- split(categories_tissue$Gene, f = categories_tissue$ID)
          map_cat2gene = map_cat2gene[lengths(map_cat2gene)>=10 & lengths(map_cat2gene) <= 1000]
          
          DE_genes = ifelse(all_genes %in% convergent_lists[[as.character(thresh)]][[expression]][[tissue]]$sig,1,0)
          names(DE_genes) = all_genes
          np = nullp(DE_genes,bias.data = log(averageExpression[[tissue]][names(DE_genes)]))
          res = goseq(np,gene2cat = map_cat2gene,use_genes_without_cat=T)
          res$BH = p.adjust(res$over_represented_pvalue,method='BH')
          # res$Name = kegg_p$PathwayName[match(res$category,kegg_p$PathwayID)]
          head(res)
          GO_convergent[[as.character(thresh)]][[tissue]] = data.frame(Thresh = thresh, Tissue = tissue, res)
          print(sprintf('%0.2f %s %d',thresh, tissue, sum(res$BH<0.05)))
        }
      }
      saveRDS(GO_convergent,file = file.path(DIR_Output,sprintf('%s_%s_GO_of_convergentGenes.rds',expression,method)))
    }
  }
}


