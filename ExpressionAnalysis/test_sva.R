library(sva)
library(limma)
library(ashr)
library(goseq)
all_categories = fread(file = '../input/all_categories.csv',data.table=F)
all_category_info = fread(file = '../input/all_category_info.csv',data.table=F)

tissue = Loc_Tiss = tissues[2]

v_Total = v_Totals[[tissue]]
E = v_Total$E
samples = v_Total$samples
design = v_Total$design
samples$Row_std = resid(lm(Row~Block,samples))
samples$Row_std2 = resid(lm(Row~SamplingGroup,samples))

ind = samples$Tissue == tissue
samples$X1[ind] = Met_field_info$X1[match(paste(samples$MaleParent[ind],samples$Block[ind],sep='_'),Met_field_info$ID)]
samples$Y1[ind] = Met_field_info$Y1[match(paste(samples$MaleParent[ind],samples$Block[ind],sep='_'),Met_field_info$ID)]
samples$X1 = factor(samples$X1)
samples$Row_std3 = resid(lm(Row~factor(X1),samples))
samples$Row_std4 = resid(lm(Row~factor(X1):Block+Plate,samples))
# samples$SamplingGroup2 = factor(samples$SamplingGroup)
samples$SamplingGroup2 = 1
samples$SamplingGroup2[samples$Row > 4080] = 2
samples$SamplingGroup2[samples$Row >= 4163] = 3
samples$SamplingGroup2 = factor(samples$SamplingGroup2)
samples$Row_std5 = resid(lm(Row~SamplingGroup2,samples))

# 4162
# sample

m = rowMeans(E)
v = rowMeans((E-m)^2)
i = m>5 & m < 8
plot(v~m)

pca = prcomp(t(E[i,]))
plot(pca$x[,1:2])


d = data.frame(x=pca$x[,1],y=pca$x[,2],samples)
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = Row_std)) + facet_wrap(~Block)
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = Row)) + facet_wrap(~SamplingGroup)
ggplot(d,aes(x=x,y=Row_std2)) + geom_point(aes(color = y)) + facet_wrap(~SamplingGroup)
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = Row_std3)) + facet_grid(X1~.)
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = Plate)) + facet_grid(X1~.)
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color = Row_std3)) + facet_grid(X1~group)
ggplot(d,aes(x=group,y=x)) + geom_boxplot(aes(color = Row_std3)) + facet_grid(X1~.)
ggplot(d,aes(x=x,y=Row_std3)) + geom_point(aes(color = y)) + facet_grid(X1~.)
ggplot(d,aes(x=x,y=Y1)) + geom_point(aes(color = y)) + facet_grid(X1~.)
ggplot(d,aes(x=y,y=Y1)) + geom_point(aes(color = x)) + facet_grid(X1~.)
ggplot(d,aes(x=group,y=y)) + geom_boxplot(aes(color = Row_std3)) + facet_grid(X1~.)

ggplot(d,aes(x=))

boxplot(d$x~d$X1)
boxplot(d$x~d$Plate)
plot(d$Row~jitter(as.numeric(as.factor(d$Plate))),col = factor(d$Block))
plot(d$x~jitter(as.numeric(as.factor(d$Plate))),col = factor(d$Block))
plot(d$x~jitter(as.numeric(as.factor(d$X1))),col = factor(d$Block))
ggplot(d,aes(x=X1,y=x)) + geom_boxplot(aes(color = group),position = position_dodge())
ggplot(d,aes(x=X1,y=x)) + geom_jitter(aes(color = group),width=.2) + facet_wrap(~Block+Plate)
ggplot(d,aes(x=X1,y=x)) + geom_jitter(aes(color = Row_std4),width=.2) + facet_wrap(~Block+Plate)
ggplot(d,aes(x=interaction(Block,Plate),y=x)) + geom_jitter(aes(color = Row_std3),width=.2) + facet_wrap(~X1)
ggplot(d,aes(x=Row,y=x)) + geom_point(aes(color=SamplingGroup2)) + facet_wrap(~X1+Block+Plate,scales='free')

ggplot(d,aes(x=Row,y=y)) + geom_point(aes(color=Plate)) + facet_wrap(SamplingGroup2~.,scales='free',ncol = 1)
ggplot(d,aes(x=Row,y=y)) + geom_point(aes(color=factor(SamplingGroup2))) + facet_wrap(Plate+Block~.,scales='free',ncol = 1)
ggplot(d,aes(x=x,y=lib.size)) + geom_point(aes(color=group)) + facet_grid(SamplingGroup2~.,scales='free')


ggplot(d,aes(x=Row,y=y)) + geom_point(aes(color=X1)) + facet_wrap(~SamplingGroup2,scales='free')
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color=Row_std5)) + facet_wrap(~SamplingGroup2,scales='free')
ggplot(d,aes(x=x,y=Row)) + geom_point(aes(color=Row_std5)) + facet_grid(SamplingGroup2~.,scales='free')
ggplot(d,aes(x=y,y=Row)) + geom_point(aes(color=Row_std5)) + facet_grid(SamplingGroup2~.,scales='free')
ggplot(d,aes(x=x,y=y)) + geom_point(aes(color=SamplingGroup2))# + facet_wrap(~SamplingGroup2,scales='free')
ggplot(d,aes(x=Row,y=(x^2+y^2))) + geom_point(aes(color=Row)) + facet_wrap(~SamplingGroup2,scales='free')

anova(lm(pca$x[,2]~Plate+Row:SamplingGroup2 + I(Row^2):SamplingGroup2 + Block + group,d))
anova(lm(pca$x[,2]~Plate+Row:SamplingGroup2 + I(Row^2):SamplingGroup2 + Block + group,d))
anova(lm(pca$x[,2]~Plate+Row:SamplingGroup2 + Block + group,d))

# geom_vline(xintercept = c(4163,4080))
anova(lm(y~X1*Plate,d))
anova(lm(pca$x[,1]~Row+SamplingGroup2:Row+X1+Block*Plate*group,d))
a=sapply(1:10,function(i) anova(lm(pca$x[,i]~Row+SamplingGroup2:Row+X1+Block*Plate*group,d))$P)[c(5,7,8),];a
anova(lm(pca$x[,1]~Row+SamplingGroup2:Row+Block*Plate*group,d))
b=sapply(1:10,function(i) anova(lm(pca$x[,i]~Row:SamplingGroup2+X1+Block*Plate*group,d))$P)[c(4,6,7),];b
plot(a[1,],b[1,],log='xy');abline(0,1)
anova(lm(pca$x[,1]~Row+SamplingGroup2:Row+Plate+Block*group,d))
c=sapply(1:10,function(i) anova(lm(pca$x[,i]~Row+SamplingGroup2:Row+Plate+Block*group,d))$P)[c(4,6),];c
plot(a[1,],c[1,],log='xy');abline(0,1)
plot(a[2,],c[2,],log='xy');abline(0,1)
anova(lm(pca$x[,1]~Row+SamplingGroup2:Row+Block*group,d))
c2=sapply(1:10,function(i) anova(lm(pca$x[,i]~Row+SamplingGroup2:Row+Block*group,d))$P)[c(3,5),];c2
plot(a[1,],c2[1,],log='xy');abline(0,1)
plot(a[2,],c2[2,],log='xy');abline(0,1)
anova(lm(pca$x[,1]~Row+SamplingGroup2:Row+Block+group,d))
c3=sapply(1:10,function(i) anova(lm(pca$x[,i]~Row+SamplingGroup2:Row+Block+group,d))$P)[c(3),];c3
plot(a[1,],c3,log='xy');abline(0,1)
anova(lm(pca$x[,1]~Row+SamplingGroup2:Row+group,d))
c4=sapply(1:10,function(i) anova(lm(pca$x[,i]~Row+SamplingGroup2:Row+group,d))$P)[c(2),];c4
plot(a[1,],c4,log='xy');abline(0,1)


design = model.matrix(~Block+Population+Population:Elevation_class+Population:Lat_std,samples)
colnames(design)
cols = 4:7
design = model.matrix(~Plate+Block+Population+Population:Elevation_class+Population:Lat_std,samples)
colnames(design)
cols = 7:10

design = model.matrix(~Plate+X1+Block+Population+Population:Elevation_class+Population:Lat_std,samples)
colnames(design)
cols = 11:14

design = model.matrix(~SamplingGroup2:Row+Plate+Block+Population+Population:Elevation_class+Population:Lat_std,samples)
cols = c(6,10:13)
design = model.matrix(~SamplingGroup2:poly(Row,2)+Plate+Block+Population+Population:Elevation_class+Population:Lat_std,samples)[,-c(12)]
cols = c(4,12:15)

vfit = lmFit(v_Total,design)
efit <- eBayes(vfit)

effects = efit$coefficients[,cols]
ts = efit$t[,cols]
ses = effects/ts
ps = 2*pt(abs(ts),efit$df.residual,lower.tail=F)
BH = apply(ps,2,p.adjust,method = 'BH')
ash = sapply(1:length(cols),function(x) ash(effects[,x],ses[,x],df=median(efit$df.residual))$result$lfsr)
colnames(ash) = colnames(effects)
rownames(ash) = rownames(effects)
print(colSums(ash<0.05))

thresh = 0.1
genes = rownames(ash)
categories_tissue = subset(all_categories,Gene %in% genes)
map_cat2gene <- split(categories_tissue$Gene, f = categories_tissue$ID)
map_cat2gene = map_cat2gene[lengths(map_cat2gene)>=10 & lengths(map_cat2gene) <= 1000]
DE_genes = ifelse(ash[,3]<thresh,1,0)
names(DE_genes) = genes
np = nullp(DE_genes,bias.data = log(averageExpression[[Loc_Tiss]][names(DE_genes)]))
res = goseq(np,gene2cat = map_cat2gene)
res$BH = p.adjust(res$over_represented_pvalue,method = 'BH')
head(res,n=20)
sum(res$BH < 0.05)

mm = model.matrix(~poly(Row,3):SamplingGroup2,samples)
fit = lmFit(v_Total,mm)
E_resid = E - fit$coeff %*% t(mm)
pca = prcomp(t(E_resid[i,]))
plot(pca$x[,1:2],col = factor(samples$Plate))
plot(pca$x[,1:2],col = factor(samples$group))

design = model.matrix(~SamplingGroup2:Row+Plate+Block+Population+Population:Elevation_class+Population:Lat_std,samples)
res_sva = sva(E,design,design[,-c(11:14)])

  