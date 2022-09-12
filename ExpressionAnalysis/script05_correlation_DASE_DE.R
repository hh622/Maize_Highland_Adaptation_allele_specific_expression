library(ggplot2)

tissues=c("MetLeaftip","MetLeafbase","PvLeaftip" )
Loc_Tiss = tissues[2]

v_Total = readRDS('../processed/v_Totals.rds')[[Loc_Tiss]]
v_ASE = readRDS('../processed/v_ASEs_10_10.rds')[[Loc_Tiss]]

ID = intersect(rownames(v_Total$E),rownames(v_ASE$E))
v_Total = v_Total[ID,]
v_ASE = v_ASE[ID,]
samples = v_Total$samples

m = rowMeans(v_Total$E)
i = order(-m)[1:300]
cors = sapply(1:nrow(v_Total$E),function(i) cor(v_Total$E[i,],v_ASE$E[i,],use='p'))
slope = sapply(1:nrow(v_Total$E),function(i) coef(lm(v_Total$E[i,]~v_ASE$E[i,],weights = sqrt(v_Total$weights[i,]*v_ASE$weights[i,])))[2])
slope_w = sapply(1:nrow(v_Total$E),function(i) coef(lm(v_Total$E[i,]~v_ASE$E[i,],weights = sqrt(v_Total$weights[i,]*v_ASE$weights[i,])))[2])


samples$Lat_std = resid(lm(abs(Lat)~Population:Elevation_class,samples))
samples$Row_std = resid(lm(Row~SamplingGroup,samples))
samples$Row_group = samples$Row-tapply(samples$Row,samples$SamplingGroup,min)[samples$SamplingGroup]
design = model.matrix(~Plate + SamplingGroup+poly(Row_std,3)+SamplingGroup:poly(Row_std,3) + Block+Population+Population:Elevation_class+Population:Lat_std,samples)
# design = model.matrix(~Plate + SamplingGroup+SamplingGroup:bs(Row_std,knots=1) + Block+Population+Population:Elevation_class+Population:Lat_std,samples)
colnames(design) <- gsub("Block|Population|Elevation_class", "", colnames(design))
qr_design = qr(design)
design = design[,qr_design$pivot[1:qr_design$rank]]
cols = match(c('(Intercept)','SA','Mex:High','SA:High','Mex:Lat_std','SA:Lat_std'),colnames(design))[3:4]

design_technical = design[,c(1:10,12:17)]
fitted = lmFit(v_Total,design_technical)$coeff %*% t(design_technical)
cors_fitted = sapply(1:nrow(v_Total$E),function(i) cor(fitted[i,],v_ASE$E[i,],use='p'))
ase_scrambled = v_ASE$E[,sample(1:ncol(v_ASE))]
cors_scrambled = sapply(1:nrow(v_Total$E),function(i) cor(v_Total$E[i,],ase_scrambled[i,],use='p'))
cors_fitted_scrambled = sapply(1:nrow(v_Total$E),function(i) cor(fitted[i,],ase_scrambled[i,],use='p'))

plot(density(cors))
lines(density(cors_fitted),col=2)
lines(density(cors_scrambled),col=3)
lines(density(cors_fitted_scrambled),col=4)


vfit = lmFit(v_Total,design)
efit <- eBayes(vfit)

r = topTable(efit,coef = c(7:9,12:17),n=100)
j=3542;
d = data.frame(samples,y=v_Total$E[rownames(r)[j],],y2=v_ASE$E[rownames(r)[j],]);ggplot(d,aes(x=Row_group,y=y)) + geom_point(aes(color=SamplingGroup)) + 
  geom_smooth(aes(color=SamplingGroup),se=F,formula = y~poly(x,3),method='lm') + facet_wrap(~SamplingGroup,ncol=2)
j=order(-abs(cors))[3];
d = data.frame(samples,y=v_Total$E[rownames(r)[j],],y2=v_ASE$E[rownames(r)[j],]);
ggplot(d,aes(x=y2,y=y)) + geom_point(aes(color = Row_group))

r = topTable(efit,coef = c(7:9,12:17),n=Inf,sort='n')
r = topTable(efit,coef = cols[-1],n=Inf,sort='n')
