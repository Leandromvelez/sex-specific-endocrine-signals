setwd('') # Set working directory where all files are stored
load('datasets used/GTEx NA included env.RData')

deg_table = read.csv('datasets used/hypothalamus DEGs FoverM.csv')
orths = read.delim('datasets used/ms_human_orthology.txt')
sex_annots = read.csv('datasets used/GTEx Subject IDs with estrogen annots.csv')

# Load filtered list of human secreted proteins
Secreted_proteins <- read.delim("datasets used/uniprot-secreted-filtered-organism__Homo+sapiens+(Human)+[9606]_.tab", header = T, check.names = F)


orth_table =  orths[orths$Common.Organism.Name =='mouse, laboratory',]
sse2 = orths[orths$Common.Organism.Name =='human',]
orth_table$human_gene_symbol = sse2$Symbol[match(orth_table$DB.Class.Key, sse2$DB.Class.Key)]
deg_table$human_orth = orth_table$human_gene_symbol[match(deg_table$Gene, orth_table$Symbol)]

library(reshape2)
library(dplyr)
library(ggplot2)
library(WGCNA)
library(qvalue)
library(ggpubr)

#Extract the dataset which has been filtered
working_dataset=GTEx_subfiltered
row.names(working_dataset) = working_dataset$gene_tissue
working_dataset$gene_tissue=NULL
working_dataset = as.data.frame(t(working_dataset))

mm1 = as.data.frame(colnames(working_dataset))
colnames(mm1) = 'gene_tissue'
mm1$gene_symbol = gsub("\\_.*","",mm1$gene_tissue)
mm1$tissue = gsub(".*_","",mm1$gene_tissue)
head(mm1)
brain_tissues = unique(mm1$tissue[grepl('Brain', mm1$tissue)])
tissue_list = c('Adipose - Subcutaneous', 'Adipose - Visceral (Omentum)', brain_tissues, 'Breast - Mammary Tissue', 'Muscle - Skeletal', 'Ovary', 'Stomach', 'Pituitary', 'Small Intestine - Terminal Ileum', 'Testis', 'Thyroid', 'Uterus', 'Vagina')
mm1 = mm1[mm1$tissue %in% tissue_list,]
new_working = working_dataset[,colnames(working_dataset) %in% mm1$gene_tissue]




##########################################################################################
#start here for pathhways

#Now crosstissue - female high estr

male_biased = deg_table[deg_table$adjP<0.05,]

tissue_list = c('Adipose - Subcutaneous', 'Adipose - Visceral (Omentum)', 'Muscle - Skeletal', 'Stomach')

ind_list = sex_annots[sex_annots$estrogen_cutoff=='High_expr',]
mm2 = mm1[mm1$gene_symbol %in% male_biased$human_orth,]
mm2 = mm2[mm2$tissue=='Brain - Hypothalamus',]
mm3 = mm1[mm1$tissue %in% tissue_list,]
working_1 = new_working[row.names(new_working) %in% ind_list$GTEx_ID,]
origin = working_1[,colnames(working_1) %in% mm3$gene_tissue]
targets = working_1[,colnames(working_1) %in% mm2$gene_tissue]

### Obtain bicor and p values
tissue.tissue.p = bicorAndPvalue(origin, targets, use='pairwise.complete.obs')

tt1 = tissue.tissue.p$p
tt1[is.na(tt1)] = 0.5
tt1[tt1==0] = 0.5
cc3 = as.data.frame(rowMeans(-log10(tt1)))


colnames(cc3) = 'Ssec_score'
cc3$gene_tissue = row.names(cc3)
cc3 = cc3[order(cc3$Ssec_score, decreasing = T),]
summary(cc3$Ssec_score)
scores=cc3
scores$gene_symbol =  gsub("\\_.*","",scores$gene_tissue)
scores$tissue = gsub(".*_","",scores$gene_tissue)
scores$secreted = ifelse(scores$gene_symbol %in% Secreted_proteins$`Gene names  (primary )`, 'Secreted', 'Non-secreted')
scores$Sig_P1e3 = ifelse(scores$Ssec_score > mean(scores$Ssec_score + (sd(scores$Ssec_score)*3)), 'Significant', 'NS')
scores$Sig_P1e2 = ifelse(scores$Ssec_score > mean(scores$Ssec_score + (sd(scores$Ssec_score)*2)), 'Significant', 'NS')
scores$Sig_P1e4 = ifelse(scores$Ssec_score > mean(scores$Ssec_score + (sd(scores$Ssec_score)*4)), 'Significant', 'NS')
scores$est_cat = paste0('high_estrogen')
ful_scores = scores[order(scores$Ssec_score, decreasing=T),]

ind_list = sex_annots[sex_annots$estrogen_cutoff=='Low_expr',]
mm2 = mm1[mm1$gene_symbol %in% male_biased$human_orth,]
mm2 = mm2[mm2$tissue=='Brain - Hypothalamus',]
mm3 = mm1[mm1$tissue %in% tissue_list,]
working_1 = new_working[row.names(new_working) %in% ind_list$GTEx_ID,]
origin = working_1[,colnames(working_1) %in% mm3$gene_tissue]
targets = working_1[,colnames(working_1) %in% mm2$gene_tissue]

tissue.tissue.p = bicorAndPvalue(origin, targets, use='pairwise.complete.obs')

tt1 = tissue.tissue.p$p
tt1[is.na(tt1)] = 0.5
tt1[tt1==0] = 0.5
cc3 = as.data.frame(rowMeans(-log10(tt1)))

table(tt1[tt1==0.5])

colnames(cc3) = 'Ssec_score'
cc3$gene_tissue = row.names(cc3)
cc3 = cc3[order(cc3$Ssec_score, decreasing = T),]
summary(cc3$Ssec_score)
scores=cc3
scores$gene_symbol =  gsub("\\_.*","",scores$gene_tissue)
scores$tissue = gsub(".*_","",scores$gene_tissue)
scores$secreted = ifelse(scores$gene_symbol %in% Secreted_proteins$`Gene names  (primary )`, 'Secreted', 'Non-secreted')
scores$Sig_P1e3 = ifelse(scores$Ssec_score > mean(scores$Ssec_score + (sd(scores$Ssec_score)*3)), 'Significant', 'NS')
scores$Sig_P1e2 = ifelse(scores$Ssec_score > mean(scores$Ssec_score + (sd(scores$Ssec_score)*2)), 'Significant', 'NS')
scores$Sig_P1e4 = ifelse(scores$Ssec_score > mean(scores$Ssec_score + (sd(scores$Ssec_score)*4)), 'Significant', 'NS')
scores$est_cat = paste0('low_estrogen')
ful_scores1 = as.data.frame(rbind(ful_scores, scores))

zz1 = ful_scores1[ful_scores1$Sig_P1e2=='Significant',]
zz1 = na.omit(zz1)
table(zz1$est_cat)

pdf(file = paste0('Cumulative Enrichments of all significant cors.pdf'))
ggplot(zz1, aes(x=tissue, y=Ssec_score, fill=est_cat)) + theme_classic() + geom_violin(width=0.6) + geom_boxplot(width=0.2, position = position_dodge(width=0.6), alpha=0.3, color='grey') + scale_fill_manual(values=c('darkorange3', 'darkorchid4')) + xlab('') + ylab('Cumulative Enrichments of all significant cors -log(pvalue)')
dev.off()

write.csv(ful_scores1, file ='crosstissue enrichments all genes.csv')
colnames(targets)[1:30]
head(ful_scores)
mean_table = ful_scores1 %>%  dplyr::select(tissue, est_cat, Ssec_score) %>% dplyr::group_by(est_cat, tissue) %>% dplyr::summarise(mean=mean(Ssec_score))
mean_table$est_tissue = paste0(mean_table$est_cat, '_', mean_table$tissue)
ful_scores1$est_tissue = paste0(ful_scores1$est_cat, '_', ful_scores1$tissue)
ful_scores1$ssec_mean = mean_table$mean[match(ful_scores1$est_tissue, mean_table$est_tissue)]
ful_scores1$normalized_ssec = ful_scores1$Ssec_score/ful_scores1$ssec_mean

low_co = mean(ful_scores1$normalized_ssec[ful_scores1$est_cat=='low_estrogen'] + (sd(ful_scores1$normalized_ssec[ful_scores1$est_cat=='low_estrogen']*2)))

high_co = mean(ful_scores1$normalized_ssec[!ful_scores1$est_cat=='low_estrogen'] + (sd(ful_scores1$normalized_ssec[!ful_scores1$est_cat=='low_estrogen']*2)))

ful_scores1$norm_CO = ifelse(ful_scores1$est_cat=='low_estrogen', paste0(low_co), paste0(high_co))
ful_scores1$norm_sig = ifelse(ful_scores1$normalized_ssec> ful_scores1$norm_CO, 'Significant', 'NS')

zz1 = ful_scores1[ful_scores1$norm_sig=='Significant',]

zz1 = na.omit(zz1)
table(zz1$est_cat)

pdf(file = paste0('Cumulative Enrichments of all significant cors - tissue_cond normalized.pdf'))
ggplot(zz1, aes(x=tissue, y=Ssec_score, fill=est_cat)) + theme_classic() + geom_violin(width=0.6) + geom_boxplot(width=0.2, position = position_dodge(width=0.6), alpha=0.3, color='grey') + scale_fill_manual(values=c('darkorange3', 'darkorchid4')) + xlab('') + ylab('')
dev.off()


write.csv(zz1, file ='crosstissue enrichments only sig genes.csv')

zz1 = ful_scores1[ful_scores1$Sig_P1e2=='Significant',]
zz1 = na.omit(zz1)

############################################################
#run for specific pathways

pathway_annots = read.delim('datasets used/uniprot-human-genes and goterms mapping.tab')
path_plots = function(pathway_term_set) {
  #path_term = 'ligand'
  tt2 = pathway_annots[grepl(pathway_term_set, pathway_annots$Gene.ontology..biological.process.),]
  tt2 = na.omit(tt2)
  zz1 = ful_scores1[ful_scores1$gene_symbol %in% tt2$Gene.names,]
  zz1 = na.omit(zz1)
  table(zz1$est_cat)
  zz1$est_cat = factor(zz1$est_cat, levels=c('low_estrogen','high_estrogen'))
  table(zz1$tissue)
  
  zz1$tissue= gsub('Adipose - Visceral (Omentum)', 'Adipose', zz1$tissue, fixed = T)
  zz1$tissue= gsub('Adipose - Subcutaneous', 'Adipose', zz1$tissue, fixed = T)
  
  
  head(zz1)
  zz1$tissue_cat = paste0(zz1$tissue, '_', zz1$est_cat)
  my_comparisons <- list( unique(zz1$tissue_cat))
  pdf(file = paste0('Cumulative Enrichments for ',pathway_term_set,  ' - all significant cors.pdf'))
  
  g1 = ggplot(zz1, aes(x=tissue, y=Ssec_score, fill=est_cat)) +  geom_violin(width=0.6) + geom_boxplot(width=0.1, position = position_dodge(width=0.6), alpha=0.3, color='grey') + scale_fill_manual(values=c('darkorange3', 'darkorchid4')) + xlab('') + ylab('') + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + stat_compare_means( method = "wilcox.test") + ggtitle(paste0(pathway_term_set, ' comparison by estrogen signaling category')) +  theme(axis.line.x.bottom = element_line(color = "black"), axis.line.y.left = element_line(color = "black"), panel.background = element_rect(fill = "white", color = "white"))
  print(g1)
  dev.off()
  
} ## Following plots are going to be created as pdf's in the working directory
path_plots('peptide hormone')
path_plots('feeding behavior')
path_plots('ligand')


pathway_term_set = 'All secreted proteins'
zz1 = ful_scores1[ful_scores1$gene_symbol %in% Secreted_proteins$`Gene names  (primary )`,]
zz1 = na.omit(zz1)
table(zz1$est_cat)
zz1$est_cat = factor(zz1$est_cat, levels=c('low_estrogen','high_estrogen'))
table(zz1$tissue)

zz1$tissue= gsub('Adipose - Visceral (Omentum)', 'Adipose', zz1$tissue, fixed = T)
zz1$tissue= gsub('Adipose - Subcutaneous', 'Adipose', zz1$tissue, fixed = T)


head(zz1)
zz1$tissue_cat = paste0(zz1$tissue, '_', zz1$est_cat)
my_comparisons <- list( unique(zz1$tissue_cat))
pdf(file = paste0('Cumulative Enrichments for ',pathway_term_set,  ' - all significant cors.pdf'))

g1 = ggplot(zz1, aes(x=tissue, y=Ssec_score, fill=est_cat)) + theme_classic() + geom_violin(width=0.6) + geom_boxplot(width=0.1, position = position_dodge(width=0.6), alpha=0.3, color='grey') + scale_fill_manual(values=c('darkorange3', 'darkorchid4')) + xlab('') + ylab('') + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + stat_compare_means( method = "wilcox.test") + ggtitle(paste0(pathway_term_set, ' comparison by estrogen signaling category') )
print(g1)
dev.off()


###################################################################################################################
#calculate interaction term for given pathways
pathway_term_set = 'feeding behavior'
tt2 = pathway_annots[grepl(pathway_term_set, pathway_annots$Gene.ontology..biological.process.),]
tt2 = na.omit(tt2)

#toggle here for all secreted proteins
#zz1 = ful_scores1[ful_scores1$gene_symbol %in% Secreted_proteins$`Gene names  (primary )`,]
zz1 = ful_scores1[ful_scores1$gene_symbol %in% tt2$Gene.names...primary..,]

zz1 = na.omit(zz1)
table(zz1$est_cat)
zz1$est_cat = factor(zz1$est_cat, levels=c('low_estrogen','high_estrogen'))
table(zz1$tissue)

zz1$tissue= gsub('Adipose - Visceral (Omentum)', 'Adipose', zz1$tissue, fixed = T)
zz1$tissue= gsub('Adipose - Subcutaneous', 'Adipose', zz1$tissue, fixed = T)


zz1$tissue_cat = paste0(zz1$tissue, '_', zz1$est_cat)
my_comparisons <- list( unique(zz1$tissue_cat))
interAB<-interaction(zz1$est_cat, zz1$tissue)
kw_test = kruskal.test(Ssec_score ~ interAB, data = zz1)
