wtd.var = function (x, weights = NULL, normwt = FALSE, na.rm = TRUE, method = c("unbiased","ML")){
  # from Hmisc package, can't install on Farm
  method <- match.arg(method)
  if (!length(weights)) {
    if (na.rm)
      x <- x[!is.na(x)]
    return(var(x))
  }
  if (na.rm) {
    s <- !is.na(x + weights)
    x <- x[s]
    weights <- weights[s]
  }
  if (normwt)
    weights <- weights * length(x)/sum(weights)
  if (normwt || method == "ML")
    return(as.numeric(stats::cov.wt(cbind(x), weights, method = method)$cov))
  sw <- sum(weights)
  if (sw <= 1)
    warning("only one effective observation; variance estimate undefined")
  xbar <- sum(weights * x)/sw
  sum(weights * ((x - xbar)^2))/(sw - 1)
}

scramble_sets = function(sets,gene_stat,seed=1) {
  ordered_genes = names(gene_stat[order(gene_stat)])
  names(ordered_genes) = ordered_genes[(seed-1+1:length(ordered_genes)) %% length(ordered_genes)+1]
  sets0 = sets[ordered_genes[rownames(sets)],]
  sets0
}


get_sets = function(expr,gs) {
  sets = do.call(cbind,mclapply(gs,function(set) rownames(expr) %in% set,mc.cores=cores))
  rownames(sets) = rownames(expr)
  colnames(sets) = names(gs)
  sets
}

do_kdcf = function(Y,weights,cores=1) {
  # require(Hmisc)
  expr_kernel = do.call(rbind,mclapply(1:nrow(Y),function(i) {
    w = weights[i,]
    w = w/mean(w,na.rm=T)
    sd = sqrt(wtd.var(Y[i,],weights = w,na.rm = T))
    e = rowMeans(sweep(pnorm(outer(Y[i,],Y[i,],'-'),0,sd/4),2,w,'*'),na.rm=T)
    e
  },mc.cores = cores))
  expr_kernel
}

do_ES = function(expr_kernel,weights,sets,tau = 1,standardize_scores=T,cores=1) {
  # require(Hmisc)
  require(Rcpp)
  require(parallel)
  sourceCpp('ES_score_matrix.cpp')
  # recover()
  ES = do.call(cbind,mclapply(1:ncol(expr_kernel),function(j) {
    # if(verbose) print(j)
    d = !is.na(expr_kernel[,j])
    x = expr_kernel[d,j]
    w = weights[d,j]#/sum(weights[d,j])
    #
    o = order(x)
    x = x[o]
    w = w[o]
    z = 1:length(x)
    r = abs(length(z)/2 - z)^tau
    # r[] = 1
    sets_j = sets[which(d)[o],,drop=F]
    set_stats = get_set_stats(w,r,sets_j)
    es = ES_score_matrix(r*w,set_stats$sum_rw_sets,sets_j,w,set_stats$dec)
    if(standardize_scores) es = c(es/set_stats$norm_factor,set_stats$norm_factor)
    # n_sets = colSums(w*sets_j)^2/colSums(w^2*sets_j)
    # n_not_sets = colSums(w*(1-sets_j))^2/colSums(w^2*(1-sets_j))
    # nf2 = 1/sqrt(n_sets*n_not_sets/(n_sets+n_not_sets))
    # nr_sets = colSums(r*w*sets_j)^2/colSums((r*w)^2*sets_j)
    # n_not_sets = colSums(w*(1-sets_j))^2/colSums(w^2*(1-sets_j))
    # nf3 = 1/sqrt(nr_sets*n_not_sets/(nr_sets+n_not_sets))
    # n_sets = colSums(w*sets_j)^2/colSums(w^2*sets_j)
    # r_sets = colSums(r*sets_j)^2/colSums(r^2*sets_j)
    # n_sets2 = sqrt(n_sets*r_sets)
    # n_not_sets = colSums(w*(1-sets_j))^2/colSums(w^2*(1-sets_j))
    # nf4 = 1/sqrt(n_sets2*n_not_sets/(n_sets2+n_not_sets))
    # es = c(es,set_stats$norm_factor,nf2,nf3,nf4)
    # },{
    # # 
    # # #
    # rw = r*w
    # # rset_w = abs(rw)*sets_j
    # # sum_rset_w = colSums(rset_w)
    # sum_rset_w = get_sum_rset_w(rw,sets_j)
    # #
    # # sum_w = sum(w)
    # #
    # # dec = 1/(sum_w - colSums(w*sets_j))
    # dec = get_dec(w,sets_j)
    # # cumsum = mx_pos = mx_neg = rep(0,ncol(sets_j))
    # # for(i in 1:nrow(sets_j)) {
    # #   ind = sets_j[i,]
    # #   cumsum[ind] = cumsum[ind] + abs(rw[i])/sum_rset_w[ind]
    # #   cumsum[!ind] = cumsum[!ind] - w[i]*dec[!ind]
    # #   mx_pos[cumsum>mx_pos] = cumsum[cumsum>mx_pos]
    # #   mx_neg[cumsum<mx_neg] = cumsum[cumsum<mx_neg]
    # # }
    # # es = -1*(mx_pos+mx_neg)
    # es = ES_score_matrix(rw,sum_rset_w,sets_j,w,dec)
    # if(standardize_scores) es = es*colSums(sets_j)/sqrt(get_sum_rset_w(w,sets_j))
    # },times=10)
    es
  },mc.cores = cores))
  # return(list(ES = ES[1:ncol(sets),],nf1 = ES[ncol(sets)+1:ncol(sets),],nf2 = ES[2*ncol(sets)+1:ncol(sets),]
  #             ,nf3 = ES[3*ncol(sets)+1:ncol(sets),],nf4 = ES[4*ncol(sets)+1:ncol(sets),]))
  if(standardize_scores) {
    norm_factors = ES[-c(1:ncol(sets)),]
    ES = ES[1:ncol(sets),]
  } else{
    norm_factors = matrix(1,nrow = ncol(sets),ncol = ncol(expr))
  }
  rownames(ES) = rownames(norm_factors) = colnames(sets)
  colnames(ES) = colnames(norm_factors) = colnames(expr)
  return(list(
    ES = ES,
    norm_factors = norm_factors
  ))
}

gsva_voom = function(expr,gs,weights = NULL,cores = 1,verbose=F) {
  # require(Hmisc)
  require(Rcpp)
  sourceCpp('ES_score_matrix.cpp')
  if(is.null(weights)) weights = matrix(1,nrow(expr),ncol(expr))
  
  # make gs into a binary matrix
  if(verbose) print('making sets')
  sets = get_sets(expr,gs)
  
  # weights = 1+0*weights
  
  # Weighted Gaussian kernel
  if(verbose) print('Kernel')
  expr_kernel = do_kdcf(expr,weights,cores)
  
  # weights = 1+0*weights
  
  # modified KS statistic max.diff = T
  if(verbose) print('ES')
  # recover()
  ES = do_ES(expr_kernel,weights,sets,cores)
  
  return(ES)
}


