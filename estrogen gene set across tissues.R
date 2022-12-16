# Set working directory (setwd) where preprocessed files are saved
setwd('')
estrogen_set = read.csv('datasets used/GSEA estrogen gene set.csv', check.names = F)
colnames(estrogen_set)= c('original', 'Entrez_id', 'Gene', 'Gene1')
#Load preprocessed GTEx data
load('datasets used/GTEx NA included env.RData')
GTEx_full=NULL
library(dplyr)
library(reshape2)
library(pheatmap)
working_dataset=GTEx_subfiltered
GTEx_subfiltered = NULL
row.names(working_dataset) = working_dataset$gene_tissue
working_dataset$gene_tissue=NULL
working_dataset = as.data.frame(t(working_dataset))
# load male-female GTEx metadata
sex_table = read.delim('datasets used/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
sex_table$GTEx_ID = gsub('GTEX-', '', sex_table$SUBJID)
sex_table$sexMF = ifelse(sex_table$SEX==1, 'M', 'F')
table(sex_table$sexMF)
# Keep only metadata samples matching the ones in working dataset
new_trts = sex_table[sex_table$GTEx_ID %in% row.names(working_dataset),]
table(new_trts$sexMF)
#make a list of estrogen gene set
est_genes = estrogen_set$Gene
# reshape working dataset to a long format
gg1 = reshape2::melt(as.matrix(working_dataset))
head(gg1)
colnames(gg1) = c('ID', 'gene_tissue', 'value')
gg1$gene_symbol = gsub("\\_.*","",gg1$gene_tissue)
gg1$gene_symbol[1:5]
gg1 = gg1[gg1$gene_symbol %in% est_genes,]
gg1$tissue = gsub(".*_","",gg1$gene_tissue)

#subsetting male and female working datasets
mm1 = gg1[gg1$ID %in% new_trts$GTEx_ID[new_trts$sexMF=='M'],]
ff1 = gg1[gg1$ID %in% new_trts$GTEx_ID[new_trts$sexMF=='F'],]
gg1=NULL
library(RColorBrewer)

test = mm1 %>% dplyr::select(gene_symbol, tissue, value) %>% 
  dplyr::group_by(gene_symbol, tissue) %>% 
  dplyr::mutate(z_score = scale(value))

test = as.data.frame(test)
head(test)
mheat = reshape2::dcast(test, gene_symbol ~ tissue, value.var = 'z_score', fun.aggregate = mean, na.rm=T)
row.names(mheat) = mheat$gene_symbol
mheat$gene_symbol=NULL
m2 = as.matrix(mheat)
m2[is.na(m2)] = 0
range(m2)
breaksList = seq(min(m2), max(m2), by = (max(m2)-min(m2))/1000)
m2 = na.omit(m2)
pdf(file = 'male counts of estrogen genes across tissues.pdf')
pheatmap(m2, color = colorRampPalette(brewer.pal(n = 7, name = "GnBu"))(length(breaksList)), breaks = breaksList, main='male counts of estrogen genes across tissues - row Z-score', fontsize_row = 1, fontsize_col = 2)
dev.off()


test = ff1 %>% dplyr::select(gene_symbol, tissue, value) %>% 
  dplyr::group_by(gene_symbol, tissue) %>% 
  dplyr::mutate(z_score = scale(value))
head(test)
test$z_score
mheat = dcast(test, gene_symbol ~ tissue, value.var = 'z_score', fun.aggregate = mean, na.rm=T)
row.names(mheat) = mheat$gene_symbol
mheat$gene_symbol=NULL
m2 = as.matrix(mheat)
m2[is.na(m2)] = 0
range(m2)
breaksList = seq(min(m2), max(m2), by = (max(m2)-min(m2))/1000)
m2 = na.omit(m2)
pdf(file = 'Female counts of estrogen genes across tissues.pdf')
pheatmap(m2, color = colorRampPalette(brewer.pal(n = 7, name = "GnBu"))(length(breaksList)), breaks = breaksList, main='Female counts of estrogen genes across tissues - row Z-score', fontsize_row = 1, fontsize_col = 2)
dev.off()



#################################
#specific tissues
brain_tissues = unique(mm1$tissue[grepl('Brain', mm1$tissue)])
tissue_list = c('Adipose - Subcutaneous', 'Adipose - Visceral (Omentum)', brain_tissues, 'Breast - Mammary Tissue', 'Muscle - Skeletal', 'Ovary', 'Pituitary', 'Testis', 'Thyroid', 'Uterus', 'Vagina')

mm2 = mm1[mm1$tissue %in% tissue_list,]
test = mm2 %>% dplyr::select(gene_symbol, tissue, value) %>% 
  dplyr::group_by(gene_symbol, tissue) %>% 
  dplyr::mutate(z_score = scale(value))
mheat = dcast(test, gene_symbol ~ tissue, value.var = 'z_score', fun.aggregate = mean, na.rm=T)
row.names(mheat) = mheat$gene_symbol
mheat$gene_symbol=NULL
m2 = as.matrix(mheat)
m2[is.na(m2)] = 0
range(m2)
breaksList = seq(min(m2), max(m2), by = (max(m2)-min(m2))/1000)
m2 = na.omit(m2)
pdf(file = 'male counts of estrogen genes across tissues - filtered list.pdf')
pheatmap(m2, color = colorRampPalette(brewer.pal(n = 7, name = "GnBu"))(length(breaksList)), breaks = breaksList, main='male counts of estrogen genes across tissues - row Z-score filtered list', fontsize_row = 1, fontsize_col = 2)
dev.off()


mm2 = ff1[ff1$tissue %in% tissue_list,]
test = mm2 %>% dplyr::select(gene_symbol, tissue, value) %>% 
  dplyr::group_by(gene_symbol, tissue) %>% 
  dplyr::mutate(z_score = scale(value))
mheat = dcast(test, gene_symbol ~ tissue, value.var = 'z_score', fun.aggregate = mean, na.rm=T)
row.names(mheat) = mheat$gene_symbol
mheat$gene_symbol=NULL
m2 = as.matrix(mheat)
m2[is.na(m2)] = 0
range(m2)
breaksList = seq(min(m2), max(m2), by = (max(m2)-min(m2))/1000)
m2 = na.omit(m2)
pdf(file = 'Female counts of estrogen genes across tissues - filtered list.pdf')
pheatmap(m2, color = colorRampPalette(brewer.pal(n = 7, name = "GnBu"))(length(breaksList)), breaks = breaksList, main='Female counts of estrogen genes across tissues - row Z-score filtered list', fontsize_row = 1, fontsize_col = 2)
dev.off()

###################################################################
#come up with thresholds
threshold  = 0
mm2 = mm1[mm1$tissue %in% tissue_list,]
test = mm2 %>% dplyr::select(gene_symbol, tissue, value) %>% 
  dplyr::group_by(gene_symbol, tissue) %>% 
  dplyr::mutate(z_score = scale(value))
males1 = test
cc1 = males1 %>% dplyr::select(gene_symbol, tissue, z_score) %>% dplyr::group_by(gene_symbol, tissue) %>% dplyr::summarise(avg=mean(z_score, na.rm=T))
head(cc1)
pdf(file = 'distribution of z-scores of estrogen genes (males).pdf')
hist(cc1$avg, main= 'distribution of row z-score - estrogen genes (males)', col = 'darkorange3')
abline(v=threshold, col= 'dodgerblue', lty=c(4), lwd=c(6))
dev.off()

test = mm2 %>% dplyr::select(ID, gene_symbol, tissue, value) %>% 
  dplyr::group_by(gene_symbol, tissue) %>% 
  dplyr::mutate(z_score = scale(value))
test$gene_tissue = paste0(test$gene_symbol, '_', test$tissue)
test = as.data.frame(test)
new_gene_tissue = test[test$z_score>threshold,]
new_gene_tissue = na.omit(new_gene_tissue)

test1 = test[test$gene_tissue %in% new_gene_tissue$gene_tissue,]
test1 = na.omit(test1)
test1$expr_category = ifelse(test1$z_score > threshold, 'High_expr', 'Low_expr')
wght_scheme = as.data.frame(table(test1$expr_category))
ratio1 = wght_scheme$Freq[wght_scheme$Var1=='High_expr'] / wght_scheme$Freq[wght_scheme$Var1=='Low_expr']

ind_list = test1 %>% dplyr::select(ID, expr_category) %>% dplyr::group_by(ID) %>% dplyr::count(expr_category)
ind_list = reshape2::dcast(ind_list, ID ~expr_category, fun.aggregate = sum, value.var = 'n')
ind_list$Low_expr = ind_list$Low_expr*ratio1
ind_list$final_expr_cat = ifelse(ind_list$High_expr > ind_list$Low_expr, 'High_expr', 'Low_expr')
male_estro_annots = ind_list

head(ind_list)

threshold  = 0
mm2 = ff1[ff1$tissue %in% tissue_list,]
test = mm2 %>% dplyr::select(gene_symbol, tissue, value) %>% 
  dplyr::group_by(gene_symbol, tissue) %>% 
  dplyr::mutate(z_score = scale(value))
males1 = test
cc1 = males1 %>% dplyr::select(gene_symbol, tissue, z_score) %>% dplyr::group_by(gene_symbol, tissue) %>% dplyr::summarise(avg=mean(z_score, na.rm=T))
head(cc1)
pdf(file = 'distribution of z-scores of estrogen genes (females).pdf')
hist(cc1$avg, main= 'distribution of row z-score - estrogen genes (females)', col = 'darkorange3')
abline(v=threshold, col= 'dodgerblue', lty=c(4), lwd=c(6))
dev.off()

test = mm2 %>% dplyr::select(ID, gene_symbol, tissue, value) %>% 
  dplyr::group_by(gene_symbol, tissue) %>% 
  dplyr::mutate(z_score = scale(value))
test$gene_tissue = paste0(test$gene_symbol, '_', test$tissue)
test = as.data.frame(test)
new_gene_tissue = test[test$z_score>threshold,]
new_gene_tissue = na.omit(new_gene_tissue)

test1 = test[test$gene_tissue %in% new_gene_tissue$gene_tissue,]
test1 = na.omit(test1)
test1$expr_category = ifelse(test1$z_score > threshold, 'High_expr', 'Low_expr')
wght_scheme = as.data.frame(table(test1$expr_category))
ratio1 = wght_scheme$Freq[wght_scheme$Var1=='High_expr'] / wght_scheme$Freq[wght_scheme$Var1=='Low_expr']

ind_list = test1 %>% dplyr::select(ID, expr_category) %>% dplyr::group_by(ID) %>% dplyr::count(expr_category)
ind_list = reshape2::dcast(ind_list, ID ~expr_category, fun.aggregate = sum, value.var = 'n')
ind_list$Low_expr = ind_list$Low_expr*ratio1
ind_list$final_expr_cat = ifelse(ind_list$High_expr > ind_list$Low_expr, 'High_expr', 'Low_expr')
female_estro_annots = ind_list



table(sex_table$sexMF)
new_trts = sex_table[sex_table$GTEx_ID %in% row.names(working_dataset),]
table(new_trts$sexMF)
males = new_trts[new_trts$sexMF=='M',]
head(males)
males$estrogen_cutoff = male_estro_annots$final_expr_cat[match(males$GTEx_ID, male_estro_annots$ID)]
table(males$estrogen_cutoff)

females = new_trts[new_trts$sexMF=='F',]
females$estrogen_cutoff = female_estro_annots$final_expr_cat[match(females$GTEx_ID, female_estro_annots$ID)]
head(females)
new_sex_table = as.data.frame(rbind(males, females))
table(females$estrogen_cutoff)
write.csv(new_sex_table, file = 'GTEx Subject IDs with estrogen annots.csv', row.names = F)
