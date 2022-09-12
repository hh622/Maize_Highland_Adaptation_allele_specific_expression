library(limma)
library(data.table)
library(ashr)
library(mashr)
library(parallel)
source('gsva_voom.R')
cores=9

tissues = c("MetLeaftip","MetLeafbase","PvLeaftip")
Populations = c('Mex','SA')

# set up dirs
DIR_Input = "../input"
DIR_Output = "../output"

all_categories = fread(file = '../input/all_categories.csv',data.table=F)
all_category_info = fread(file = '../input/all_category_info.csv',data.table=F)

# results = list()

# for(expression in c('DE','DASE')) {
  
  if(expression == 'DE') {
    voom_Elists = readRDS('../processed/v_Totals.rds')
    averageExpression = readRDS('../processed/averageExpression.rds')
  } else {
    voom_Elists = readRDS('../processed/v_ASEs_10_10.rds')
    averageExpression = readRDS('../processed/averageASECounts.rds')
  }
  
  results = list()
  
  # expression
  results[[expression]] = list()
  
  
  # GSVA
  # Loc_Tiss = tissue = tissues[2]
  trait = 'GSVA'
  gsvaScores = list()
  for(Loc_Tiss in tissues){
    print(Loc_Tiss)
    gsvaScores[[Loc_Tiss]] = list()
    v_Total = voom_Elists[[Loc_Tiss]]
    design = v_Total$design
    cols = v_Total$cols

    # v_Total$E = v_Total$E[1:5000,]
    # v_Total$weights = v_Total$weights[1:5000,]
    # averageExpression[[Loc_Tiss]] = averageExpression[[Loc_Tiss]][1:5000]
    
    # fit covariates
    design0 = design[,-cols[-1]]
    fit0 = lmFit(v_Total,design0)
    if(!weights) v_Total$weights[] = 1
    E_resid = v_Total$E - fit0$coefficients %*% t(design0)
    expr = do_kdcf(E_resid,v_Total$weights,cores=cores)
    
    if(reduce) voom_Elists[[Loc_Tiss]]$design = design[,cols]
    
    for(type in unique(all_category_info$Type)) {
      print(type)
  
      categories_tissue = subset(all_categories,Gene %in% rownames(v_Total$E) & ID %in% subset(all_category_info,Type == type)$ID)
      map_cat2gene <- split(categories_tissue$Gene, f = categories_tissue$ID)
      map_cat2gene = map_cat2gene[lengths(map_cat2gene)>=10 & lengths(map_cat2gene) <= 2000]
  
      sets = get_sets(v_Total$E,map_cat2gene)
      genes = rowSums(sets,na.rm=T)>0
      gsvaScore = do_ES(expr[genes,],v_Total$weights[genes,],sets[genes,],cores=cores)
  
      # null distributions
      sets_null = scramble_sets(sets,averageExpression[[Loc_Tiss]][rownames(v_Total$E)],seed=1)
      genes_null = rowSums(sets_null,na.rm=T)>0
      gsvaScore_null = do_ES(expr[genes_null,],v_Total$weights[genes_null,],sets_null[genes_null,],cores=cores)
  
      gsvaScores[[Loc_Tiss]][[type]] = list(
        gsvaScore = gsvaScore,
        gsvaScore_null = gsvaScore_null
      )
    }
  }
  saveRDS(gsvaScores,file = file.path(DIR_Output,sprintf('%s_gsvaScores_residuals_%s_%s.rds',expression,weights,reduce)))
  gsvaScores_all = gsvaScores
  # 
  # gsvaScores_all = readRDS(file.path(DIR_Output,sprintf('%s_gsvaScores_residuals.rds',expression)))
  results = list()
  for(Loc_Tiss in tissues){
    results[[Loc_Tiss]] = list()
    for(type in unique(all_category_info$Type)) {
      gsvaScores = gsvaScores_all[[Loc_Tiss]][[type]]
      print(sprintf('%s %s',Loc_Tiss,type))
      
      null_scores = gsvaScores$gsvaScore_null$ES
      scores = gsvaScores$gsvaScore$ES
      null_threshold = quantile(na.omit(abs(null_scores)),0.95)
      print(null_threshold)
      
      sig_categories = rowSums(abs(scores)>null_threshold,na.rm=T) >= .1*ncol(scores)
      sig_categories_null = rowSums(abs(null_scores)>null_threshold,na.rm=T) >= .1*ncol(scores)
      print(c(sum(sig_categories),sum(sig_categories_null)))
      print(sum(sig_categories_null)/(sum(sig_categories)))
  #   }
  # }
      
      try({
        design = voom_Elists[[Loc_Tiss]]$design
        cols = voom_Elists[[Loc_Tiss]]$cols
        if(reduce) cols = 1:ncol(design)
        gsvaScore = gsvaScores$gsvaScore
        ES = gsvaScore$ES
        norm_factors = gsvaScore$norm_factors
  
        # do voom
        span = 0.5
        fit1 = lmFit(ES,design)
        sy <- sqrt(fit1$sigma)
        sx = rowMeans(log(norm_factors),na.rm=T)
        plot(sy~sx,main = Loc_Tiss)
        l <- lowess(sx, sy, f = span)
        lines(l,col=2)
        f <- approxfun(l, rule = 2, ties = list("ordered", mean))
        w = matrix(1/f(log(norm_factors)),nrow = nrow(norm_factors))^4
  
        sig_categories_tissue = sig_categories[rownames(ES)]
        gsva_Elist = new('EList',list(E = ES[sig_categories_tissue,],weights=w[sig_categories_tissue,]))
  
  
        categories_info_tissue = all_category_info[match(rownames(gsva_Elist$E),all_category_info$ID),]
  
        # do final fit
        fit2 = lmFit(gsva_Elist,design)
        # fit2 = lmFit(ES,design)
        efit2 <- eBayes(fit2)
  
        effects = efit2$coefficients[,cols]
        ses = effects/efit2$t[,cols]
        ps = 2*pt(abs(effects)/ses,efit2$df.residual,lower.tail=F)
        BH = apply(ps,2,p.adjust,method = 'BH')
        ash = lapply(1:length(cols),function(x) ash(effects[,x],ses[,x],df=median(efit2$df.residual))$result)
        ash_PosteriorMean = sapply(ash,function(x) x$PosteriorMean)
        ash_lfsr = sapply(ash,function(x) x$lfsr)
        colnames(ash_PosteriorMean) = colnames(ash_lfsr) = colnames(effects)
        rownames(ash_PosteriorMean) = rownames(ash_lfsr) = rownames(effects)
        print(colSums(ash_lfsr<0.05))
  
        results[[Loc_Tiss]][[type]] = list(
                                            effects = effects,
                                            SEs = ses,
                                            ps = ps,
                                            BH = BH,
                                            PosteriorMean = ash_PosteriorMean,
                                            lfsr = ash_lfsr
                                          )
      })
    }
  }
  saveRDS(results,file = file.path(DIR_Output,sprintf('%s_results_gsvaScores_residuals_%s_%s.rds',expression,weights,reduce)))
  
  print(c(weights,reduce))
  counts = c()
  for(Loc_Tiss in tissues){
  for(type in unique(all_category_info$Type)) {
      # print(sprintf('%s %s %d',Loc_Tiss,type,nrow(results[[Loc_Tiss]][[type]]$effects)))
      # print(colSums(results[[Loc_Tiss]][[type]]$lfsr<0.05 ))
      counts = rbind(counts,data.frame(Tissue = Loc_Tiss, Type = type, N = nrow(results[[Loc_Tiss]][[type]]$effects), t(colSums(results[[Loc_Tiss]][[type]]$lfsr<0.05 ))))
    }
  }
  print(c(weights,reduce))
  counts
  
  
  



  # effectNames = colnames(results[[1]]$MetLeaftip$effects)
  # fill_matrix = function(ls,names=NULL) {
  #   if(is.null(names)) names = names(ls)
  #   IDs = unique(unlist(lapply(ls,function(x) names(x))))
  #   mat = sapply(ls,function(x) x[IDs])
  #   colnames(mat) = names
  #   rownames(mat) = IDs
  #   mat
  # }
  # results_tissues = list()
  # for(trait in names(results)) {
  #   results_tissues[[trait]] = list()
  #   for(effect in effectNames) {
  #     results_tissues[[trait]][[effect]] = list(
  #       effects = fill_matrix(lapply(tissues,function(tissue) results[[trait]][[tissue]]$effects[,effect]),tissues),
  #       SEs = fill_matrix(lapply(tissues,function(tissue) results[[trait]][[tissue]]$SE[,effect]),tissues),
  #       PosteriorMeans = fill_matrix(lapply(tissues,function(tissue) results[[trait]][[tissue]]$PosteriorMean[,effect]),tissues),
  #       lfsrs = fill_matrix(lapply(tissues,function(tissue) results[[trait]][[tissue]]$lfsr[,effect]),tissues)
  #     )
  #   }
  # }
#   
#   saveRDS(results_tissues,file = file.path(DIR_Output,sprintf('%s_results_tissues.rds',expression)))
#     
#   results_tissues = readRDS(file.path(DIR_Output,sprintf('%s_results_tissues.rds',expression)))
#   effectNames = names(results_tissues[[expression]])
#   traits = names(results_tissues)
#   
#   ##mash
#   mash_results = list()
#   library(mashr)
#   for(trait in traits) {
#     mash_results[[trait]] = list()
#     for(effect in effectNames[-1]) {
#       print(c(effect,trait))
#       effects = results_tissues[[trait]][[effect]]$effects
#       SEs = results_tissues[[trait]][[effect]]$SEs
#       lfsrs = results_tissues[[trait]][[effect]]$lfsrs
#       
#       data = mash_set_data(as.matrix(effects), as.matrix(SEs))
#       U.c = cov_canonical(data)
#       
#       # Step 2: Obtain initial data-driven covariance matrices
#       
#       # select strong signals
#       m.1by1 = mash_1by1(data)
#       strong = get_significant_results(m.1by1,0.05)
#       if(length(strong)<20) {
#         strong = order(apply(m.1by1$result$lfsr,1,min))[1:min(20,nrow(effects))]
#       }
#       
#       # Perform PCA on data and return list of candidate covariance matrices
#       U.pca = cov_pca(data,npc=ncol(effects),subset=strong)
#       # npc:	the number of PCs to use, should be less than or equal to n_conditions(data)
#       # subset: indices of the subset of data to use (set to NULL for all data)
#       # print(names(U.pca))
#       
#       ## ----------------------------------------------------------------
#       # Step 3: prepare canonical/data-driven covariance matrices
#       # Perform "extreme deconvolution" (Bovy et al) on a subset of the data
#       U.ed = cov_ed(data, U.pca, subset=strong)
#       # subset: a subset of data to be used when ED is run (set to NULL for all the data)
#       # The function cov_ed is used to apply the ED algorithm from a specified initialization
#       # (here U.pca) and to a specified subset of signals.
#       
#       
#       ## ------------------------------------------------------------------
#       # Step 4: Run mash for different models
#       # using just the simple canonical covariances as in the initial introductory vignette.
#     
#       V.em_c = mash_estimate_corr_em(data, c(U.c,U.ed), details = TRUE)
#       m.Vem_c = V.em_c$mash.model
#       m.Vem_c$result$NAs = is.na(effects)
#       m.Vem_c$V = V.em_c$V
#       PosteriorMeans = m.Vem_c$result$PosteriorMean
#       lfsrs = m.Vem_c$result$lfsr
#       PosteriorMeans[is.na(effects)] = NA
#       lfsrs[is.na(effects)] = NA
#       mash_results[[trait]][[effect]] = list(
#         PosteriorMeans = PosteriorMeans,
#         lfsrs = lfsrs
#       )
#       saveRDS(mash_results,file = file.path(DIR_Output,sprintf('%s_mash_results.rds',expression)))
#     }
#   }
#   
#   saveRDS(mash_results,file = file.path(DIR_Output,sprintf('%s_mash_results.rds',expression)))
# 
# # }
# 
# 
# 
















