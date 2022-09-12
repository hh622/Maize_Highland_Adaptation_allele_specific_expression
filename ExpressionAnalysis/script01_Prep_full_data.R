library(edgeR)
library(data.table)
library(tidyverse)
library(GO.db)
source('voom_ASE.R')
# accession_info = fread('../input/Accession_selection_for_F1s_3_2_16.csv',data.table=F)
# randomization = fread('../../../F1_selection/F1_seed_randomization_v2.csv',data.table=F)
# randomization_substitute = unique(subset(randomization,!is.na(old_Pedigree) & old_Pedigree != '')[,c(4,5,8)])
# accession_info$MaleParent = paste(accession_info$Pair,accession_info$Elevation_class,sep='-')
# accession_info$MaleParent[match(randomization_substitute$old_Pedigree,accession_info$MaleParent)] = randomization_substitute$Pedigree
# accession_info_SA_16 = subset(accession_info,MaleParent == 'SA_32-High')
# accession_info_SA_16$MaleParent = 'SA_16-High'
# accession_info = rbind(accession_info,accession_info_SA_16)
# accession_info$MaleParent = sub('_','-',accession_info$MaleParent)
# write.csv(accession_info,file = '../input/Accession_info_final.csv',row.names=F)
accession_info = fread('../input/Accession_info_final.csv',data.table=F)

total_counts = fread('../Raw_data/HiLo_featurecounts_B73v4wholegenome_totalCounts_624samples.table',data.table=F)
rownames(total_counts) = sub('gene:','',total_counts$Geneid)
total_counts = as.matrix(total_counts[,-1])

samples = data.frame(file_name = colnames(total_counts), stringsAsFactors = F)

# correct some miss-named files
samples$name = samples$file_name
samples$name[samples$file_name == "Mex-18-High_met_blk2_V4leafbase_plate7_wellC1"] = "SA-28-High_met_blk2_V4leafbase_plate7_wellC1"
samples$name[samples$file_name == "Mex-18-High_met_blk2_V4leaftip_plate7_wellD1"] = "SA-28-High_met_blk2_V4leaftip_plate7_wellD1"
samples$name[samples$file_name == "SA-28-High_met_blk2_V4leafbase_plate6_wellA7"] = "Mex-18-High_met_blk2_V4leafbase_plate6_wellA7"
samples$name[samples$file_name == "SA-28-High_met_blk2_V4leaftip_plate6_wellB7"] = "Mex-18-High_met_blk2_V4leaftip_plate6_wellB7"

# drop some bad samples

# drop likely bad samples and checks
samples = samples[-grep('Check',samples$name),]
samples = samples[-grep('Mex-9-High_met_blk2',samples$name),]
samples = samples[-grep('SA-22-Low_met_blk1',samples$name),]
samples = samples[-grep('Mex-14-High',samples$name),]

total_counts = total_counts[,samples$file_name]

samples = separate(samples,'name',into = c('MaleParent','Loc','Block','Tissue','Plate','Well'),sep='_',remove = F)
samples$Population = ifelse(grepl('Mex',samples$MaleParent),'Mex','SA')
samples$Population[grepl('Check',samples$MaleParent)] = NA
samples$Elevation_class = factor(ifelse(grepl('High',samples$MaleParent),'High','Low'),levels = c('Low','High'))
samples$Loc = as.character(factor(samples$Loc,levels = c('met','pv'),labels = c('Met','Pv')))
samples$Tissue = as.character(factor(sub('V4','',samples$Tissue),levels = c('leaftip','leafbase'),labels = c('Leaftip','Leafbase')))
samples$Tissue = paste0(samples$Loc,samples$Tissue)
samples$Lat = accession_info$Lat[match(samples$MaleParent,accession_info$MaleParent)]
samples$Elevation = accession_info$Elevation[match(samples$MaleParent,accession_info$MaleParent)]
samples$lib.size = colSums(total_counts)
samples$group = paste(samples$Population,samples$Elevation_class,sep=':')

Met_field_info = readxl::read_excel('../Field_data/met_field_layout.xlsx')
Pv_field_info = readxl::read_excel('../Field_data/PV_field_layout.xlsx')
Met_field_info$ID = paste(sub('_','-',Met_field_info$Parent),Met_field_info$Block,sep='_blk')
Pv_field_info$ID = paste(sub('_','-',Pv_field_info$Parent),Pv_field_info$Block,sep='_blk')
ind = samples$Tissue == 'MetLeaftip'
samples$Row[ind] = Met_field_info$Row[match(paste(samples$MaleParent[ind],samples$Block[ind],sep='_'),Met_field_info$ID)]
samples$SamplingGroup[ind] = Met_field_info$SamplingGroup[match(paste(samples$MaleParent[ind],samples$Block[ind],sep='_'),Met_field_info$ID)]
ind = samples$Tissue == 'MetLeafbase'
samples$Row[ind] = Met_field_info$Row[match(paste(samples$MaleParent[ind],samples$Block[ind],sep='_'),Met_field_info$ID)]
samples$SamplingGroup[ind] = Met_field_info$SamplingGroup[match(paste(samples$MaleParent[ind],samples$Block[ind],sep='_'),Met_field_info$ID)]
ind = samples$Tissue == 'PvLeaftip'
samples$Row[ind] = Pv_field_info$Row[match(paste(samples$MaleParent[ind],samples$Block[ind],sep='_'),Pv_field_info$ID)]
samples$SamplingGroup[ind] = factor(Pv_field_info$SamplingGroup[match(paste(samples$MaleParent[ind],samples$Block[ind],sep='_'),Pv_field_info$ID)])

samples = subset(samples,lib.size >= 2e6)
total_counts = total_counts[,samples$file_name]
dim(total_counts)
total_counts[1:3,1:4]
write.csv(samples,file = '../processed/sample_info.csv',row.names=F)
samples = fread(file = '../processed/sample_info.csv',data.table = F)
samples$Elevation_class = factor(ifelse(grepl('High',samples$MaleParent),'High','Low'),levels = c('Low','High'))
samples$SamplingGroup = factor(samples$SamplingGroup)

# tissue = tissues[1]
# p=20
# i = c(order(results$DE[[tissue]]$DE_ash[,4])[1:p],order(results$DE[[tissue]]$DE_ash[,3])[1:p],order(results$DE[[tissue]]$DE_ash[,2])[1:p])
# j = grep("SA-3-Low",colnames(v_Totals[[tissue]]));j
# X = v_Totals[[tissue]]$E[i,];X = sweep(X,1,rowMeans(X),'-')
# c = cor(X,X[,j[2]])
# boxplot(c~v_Totals[[tissue]]$samples$Elevation_class+v_Totals[[tissue]]$samples$Population)
# c = cor(X)
# rownames(c) = paste(v_Totals[[tissue]]$samples$MaleParent,v_Totals[[tissue]]$samples$Block)
# h = hclust(as.dist(1-c),method = 'ward.D2')
# plot(h,cex=.5)
# Image(c[h$order,h$order])
# plot(h$order~v_Totals[[tissue]]$samples$Elevation_class+v_Totals[[tissue]]$samples$Population)


full_matrix = DGEList(counts = total_counts,samples = samples)
dim(full_matrix$counts); dim(full_matrix$samples)

tissues = c("MetLeaftip","MetLeafbase","PvLeaftip")

tissue_matrices = lapply(tissues,function(tissue) {
  samples = full_matrix$samples
  sub_matrix = full_matrix[,samples$Tissue == tissue]
  # filter using defaults of filterByExpr()
  filtered_matrix = sub_matrix[filterByExpr(sub_matrix,model.matrix(~group,sub_matrix$samples),min.count = 32),]
  # filtered_matrix2 = sub_matrix[apply(sub_matrix$counts,1,median)>=32 & rowMeans(sub_matrix$counts==0)<0.2,]
  filtered_matrix$samples$Lat_std = resid(lm(abs(Lat)~Population:Elevation_class,filtered_matrix$samples))
  filtered_matrix$samples$Row_std = resid(lm(Row~SamplingGroup,filtered_matrix$samples))
  
  filtered_matrix = calcNormFactors(filtered_matrix)
  print(dim(filtered_matrix))
  filtered_matrix
})
names(tissue_matrices) = tissues
tissue_matrices$MetLeaftip$counts[1:3,1:4]

v_Totals = lapply(tissues,function(tissue) {
  # design = model.matrix(~Plate + SamplingGroup + Row_std:SamplingGroup + Block+Population+Population:Elevation_class+Population:Lat_std,tissue_matrices[[tissue]]$samples)
  design = model.matrix(~Plate + SamplingGroup + poly(Row_std,3)+poly(Row_std,3):SamplingGroup + Block+Population+Population:Elevation_class+Population:Lat_std,tissue_matrices[[tissue]]$samples)
  #modify colnames of design matrix and drop non-estimable columns
  colnames(design) <- gsub("Block|Population|Elevation_class", "", colnames(design))
  qr_design = qr(design)
  design = design[,qr_design$pivot[1:qr_design$rank]]
  cols = match(c('(Intercept)','SA','Mex:High','SA:High','Mex:Lat_std','SA:Lat_std'),colnames(design))
 
  # do voom
  v_Total = voom(tissue_matrices[[tissue]],design = design,plot=T)
  v_Total$cols = cols
  v_Total$samples = tissue_matrices[[tissue]]$samples
  
  v_Total
})
names(v_Totals) = tissues
v_Totals$MetLeaftip$E[1:3,1:4]

averageExpression = lapply(tissues,function(tissue) {
  rowMeans(tissue_matrices[[tissue]]$counts)
})
names(averageExpression) = tissues

saveRDS(v_Totals,file = '../processed/v_Totals.rds')
saveRDS(averageExpression,file = '../processed/averageExpression.rds')


ASE_counts = fread('../Raw_data/HiLo_featurecounts_ALL_B73v4wholegenome_Thr10_ASE2.table',data.table=F)
rownames(ASE_counts) = sub('gene:','',ASE_counts$Geneid)
ASE_counts = as.matrix(ASE_counts[,-1])
ASE_counts = ASE_counts[,-grep('Check',colnames(ASE_counts))]
Ref_counts = ASE_counts[,grep('vA1',colnames(ASE_counts))]
Alt_counts = ASE_counts[,grep('vA2',colnames(ASE_counts))]
colnames(Ref_counts) = sub('.vA1','',colnames(Ref_counts))
colnames(Alt_counts) = sub('.vA2','',colnames(Alt_counts))
all(samples$file_name %in% colnames(Ref_counts))
Ref_counts = Ref_counts[,samples$file_name]
Alt_counts = Alt_counts[,samples$file_name]

tissue_ASE_matrices = lapply(tissues,function(tissue) {
  sub_samples = subset(samples,Tissue == tissue)
  sub_samples$Lat_std = resid(lm(abs(Lat)~Population:Elevation_class,sub_samples))
  sub_samples$Row_std = resid(lm(Row~SamplingGroup,sub_samples))
  groups = sub_samples$group
  # count per gene number of samples per group with at least 10 counts. Should be the same as number >0.
  non_zero_counts = sapply(unique(sub_samples$group),function(group) rowSums((Ref_counts+Alt_counts)[,sub_samples$file_name][,sub_samples$group == group] >= 32))
  genes = rownames(non_zero_counts[apply(non_zero_counts,1,min)>=10,])
  print(length(genes))
  list(xRef = Ref_counts[genes,sub_samples$file_name],
       xAlt = Alt_counts[genes,sub_samples$file_name],
       samples = sub_samples)
})
names(tissue_ASE_matrices) = tissues

v_ASEs = lapply(tissues,function(tissue) {
  ASE_matrices = tissue_ASE_matrices[[tissue]]
  # design = model.matrix(~Plate + SamplingGroup + Row_std:SamplingGroup + Block+Population+Population:Elevation_class+Population:Lat_std,ASE_matrices$samples)
  design = model.matrix(~Plate + SamplingGroup + poly(Row_std,3)+poly(Row_std,3):SamplingGroup + Block+Population+Population:Elevation_class+Population:Lat_std,tissue_matrices[[tissue]]$samples)
  #modify colnames of design matrix and drop non-estimable columns
  colnames(design) <- gsub("Block|Population|Elevation_class", "", colnames(design))
  qr_design = qr(design)
  design = design[,qr_design$pivot[1:qr_design$rank]]
  cols = match(c('(Intercept)','SA','Mex:High','SA:High','Mex:Lat_std','SA:Lat_std'),colnames(design))
   
  # do voom
  v_ASE = voom_ASE(ASE_matrices$xRef,ASE_matrices$xAlt,design = design,plot = T)
  
  v_ASE$samples = ASE_matrices$samples
  v_ASE$cols = cols
  v_ASE
})
names(v_ASEs) = tissues

averageASECounts = lapply(tissues,function(tissue) {
  rowMeans(tissue_ASE_matrices[[tissue]]$xRef+tissue_ASE_matrices[[tissue]]$xAlt)
})
names(averageASECounts) = tissues
saveRDS(v_ASEs,file = '../processed/v_ASEs_10_10.rds')
saveRDS(averageASECounts,file = '../processed/averageASECounts.rds')

# load categories
library(GO.db)
load(file = '../Raw_data/Modules/maizeGAMER_B73AGPv4_aggregate_gomap.RData')
corncyc = fread('../Raw_data/Modules//corncyc_pathways.20180702',data.table=F)
corncyc_info = fread('../Raw_data/Modules/2022-04-11MaizeGeneSet_Sharing/Input/CORNCYC_PathwayList_from_All_pathways_of_Z._mays_mays.csv',data.table=F)
kegg = fread('../Raw_data/Modules/2022-04-11MaizeGeneSet_Sharing/Input/KEGG_maize_pathways.csv',data.table=F)
kegg_info = fread('../Raw_data/Modules/2022-04-11MaizeGeneSet_Sharing/Input/KEGG_maize_pathway_list.csv',data.table=F)
MaizeGRN_info = readxl::read_excel('../Raw_data/Modules/MaizeGRN/studies.xlsx')
MaizeGRN_info = subset(MaizeGRN_info,Leaves == 'T')
MaizeGRN = do.call(rbind,lapply(MaizeGRN_info$nid,
                                function(file) data.frame(File = file,fread(file.path('../Raw_data/Modules/MaizeGRN/rf_100k',paste(file,'.tsv',sep='')),data.table=F))))


all_categories = c()
all_category_info = c()
for(category in c('GO','CornCyc','KEGG','MaizeGRN')) {
  if(category == 'GO') {
    categories = data.frame(Gene = gomap$Gene,ID = gomap$GO)
    category_names = AnnotationDbi::select(GO.db,unique(gomap$GO),c('TERM','ONTOLOGY'))
    category_names = subset(category_names,!is.na(category_names$TERM))
    category_names$Type = 'GO'
    categories = subset(categories,ID %in% category_names$GOID)
    colnames(category_names)[1] = 'ID'
  } else if(category == 'CornCyc') {
    category_names = data.frame(ID = corncyc_info$`Object ID`,Pathway = corncyc_info$Pathways,Type = 'CornCyc')
    categories = data.frame(Gene = corncyc$`Gene-name`,ID = corncyc$`Pathway-id`)
    categories = subset(categories,ID %in% category_names$ID)
  } else if(category == 'KEGG') {
    category_names = data.frame(ID = kegg_info$PathwayID,kegg_info,Type='KEGG')
    categories = data.frame(Gene = kegg$geneID,ID = kegg$PathwayID)
    categories = subset(categories,ID %in% category_names$ID)
  } else if(category == 'MaizeGRN') {
    category_names = data.frame(ID = unique(MaizeGRN$regulator),Type='MaizeGRN')
    categories = data.frame(Gene = MaizeGRN$target,ID = MaizeGRN$regulator)
  }
  all_categories = dplyr::bind_rows(all_categories,categories)
  all_category_info = dplyr::bind_rows(all_category_info,category_names)
}
write.csv(all_categories,file = '../input/all_categories.csv',row.names=F)
write.csv(all_category_info,file = '../input/all_category_info.csv',row.names=F)


# 
# 
# 
# tissue = tissues[1]
# samples = tissue_matrices[[tissue]]$samples
genes = c(
'Zm00001d043461', #ZCN12
'Zm00001d010752', #ZCN8
'Zm00001d038725', #ZCN7
'Zm00001d010987', #rap2.7
'Zm00001d048474', #ZmMADS1
'Zm00001d042315', #MADS69
'Zm00001d024909', #cct1 == zmcct10
'Zm00001d003162', #col11 == zmcct9
'Zm00001d000176', #cct2
'Zm00001d053858' #ehd1
)
# a=sapply(genes,function(gene) {
#   try(boxplot(tissues_voom[[tissue]]$E[gene,]~samples$Elevation_class+samples$Population,main = gene))
# })
# 
# 
samples = v_Totals[[tissue]]$samples
a=sapply(genes,function(gene) {
  try(boxplot(v_Totals[[tissue]]$E[gene,]~samples$Elevation_class+samples$Population,main = gene))
})
