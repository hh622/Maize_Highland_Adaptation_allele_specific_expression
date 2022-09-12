
# # set up dirs
# DIR_GSVA_Input = "HiLoASEGSVA/input/"
# DIR_GSVA_Output = "HiLoASEGSVA/output/"
# if(!dir.exists(DIR_GSVA_Output)) dir.create(DIR_GSVA_Output)
# list.files(DIR_GSVA_Input)
# 
# 
# # Collect the ASE results together in a list to make processing easier.
# ASE_objects = list()
# Loc_Tiss = "MetLeaftip"
# # Loc_Tiss = "MetLeafbase"
# 
# ls=load(file = paste0(DIR_GSVA_Input,"HiLo_DASE_SingleTissAnalysis_DGEList_",Loc_Tiss,".RData"))
# 
# head(xRef$samples)
# str(xRef$samples)
# #this is the design matrix confirmed by Dan in meeting 01/26/21
# if(Loc_Tiss == "PvLeaftip"){
#   design = model.matrix(~Population+Population:Elevation_class+Population:Lat,xRef$samples)
# }
# if(Loc_Tiss == "MetLeaftip" | Loc_Tiss == "MetLeafbase"){
#   design = model.matrix(~Block+Population+Population:Elevation_class+Population:Lat,xRef$samples)
# }
# 
# #modify colnames of design matrix
# colnames(design) <- gsub("Block|Population|Elevation_class", "", colnames(design))
# head(design)
# 

voom_ASE = function(counts_REF,counts_ALT, design = NULL, lib.size = NULL, weights = NULL, span = 0.5, 
                    plot = FALSE, save.plot = FALSE) {
  out <- list()
  counts_REF = as.matrix(counts_REF)
  counts_ALT = as.matrix(counts_ALT)
  n <- nrow(counts_REF)
  if (n < 2L) 
    stop("Need at least two genes to fit a mean-variance trend")
  m <- min(counts_REF+counts_ALT)
  if (is.na(m)) 
    stop("NA counts not allowed")
  if (m < 0) 
    stop("Negative counts now allowed")
  if (is.null(design)) {
    design <- matrix(1, ncol(counts), 1)
    rownames(design) <- colnames(counts)
    colnames(design) <- "GrandMean"
  }
  if (is.null(lib.size)) 
    lib.size <- colSums(counts_REF+counts_ALT)
  ASE = t(log2(t(counts_ALT + 0.5)/(lib.size/2 + 1) * 1e+06)) - t(log2(t(counts_REF + 0.5)/(lib.size/2 + 1) * 1e+06))
  weights = matrix(1,nrow(ASE),ncol(ASE))
  weights[is.infinite(ASE)] = 0
  weights[is.na(ASE)] = 0
  weights[counts_ALT+counts_REF == 0] = 0
  ASE[weights==0] = NA
  fit <- lmFit(ASE, design, weights = weights)
  total_counts = counts_REF+counts_ALT
  normalized_totalcounts = sweep(log2(total_counts + 0.5),2,log2(lib.size+1)-log2(1e6),'-')
  sx <- rowMeans(normalized_totalcounts)
  # sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
  sy <- sqrt(fit$sigma)
  allzero <- rowSums(weights) == 0
  if (any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- lowess(sx, sy, f = span)
  if (plot) {
    plot(sx, sy, xlab = "log2( total counts + 0.5 )", ylab = "Sqrt( standard deviation )", 
         pch = 16, cex = 0.25)
    title("voom: Mean-variance trend for ASE")
    lines(l, col = "red")
  }
  f <- approxfun(l, rule = 2, ties = list("ordered", mean))
  
  # input actual counts into f, on the log2 scale, normalized by log2(mean(lib.size)/1e6)
  w <- 1/f(log2(total_counts) - log2(mean(lib.size)) + log2(1e6))^4
  dim(w) <- dim(ASE)
  w[weights == 0] = 0 # set weights for 0-counts to 0
  out$E <- ASE
  out$weights <- w
  rownames(out$weights) = rownames(out$E)
  colnames(out$weights) = colnames(out$E)
  out$design <- design
  
  if (save.plot) {
    out$voom.xy <- list(x = sx, y = sy, xlab = "log2( count size + 0.5 )", 
                        ylab = "Sqrt( standard deviation )")
    out$voom.line <- l
  }
  new("EList", out)
}

# x_ASE = voom_ASE(xRef,xAlt,design,plot=T)
