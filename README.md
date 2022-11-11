# sex-specific-endocrine-signals
## This repository provides a detailed walk-through and scripts for GTEx integration work contributing to Massa et al., 2022
### All datasets for analyses are available here: https://drive.google.com/drive/folders/1FldQDp-I9NGqgBo0cW9WS3KuYCSeHDM7?usp=sharing
### Briefly, datasets include:
#### 1. Gene sets annotated to be involved in estrogen signaling: GSEA estrogen gene set.csv
###### Source: https://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_ESTROGEN_RESPONSE_EARLY.html

#### 2. Filtered GTEx gene expression data (V8) to compare genetic correlation structure of genes across tissues: GTEx NA included env
###### Source: https://elifesciences.org/articles/76887 and https://gtexportal.org/home/

#### 3. GTEx subject annotations including reported sex: GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
###### Source: https://elifesciences.org/articles/76887 and https://gtexportal.org/home/

#### 4. List of all known secreted proteins supplied by UniProt: human secreted proteins.tab
###### Source: https://www.uniprot.org/ and https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5935137/

#### 5. List of differential expression results by sex from mouse hypothalamus: hypothalamus DEGs FoverM.csv
###### Source: This study

#### 6. List of annotated mouse and human orthologue genes: ms_human_orthology.txt
###### Source: https://www.alliancegenome.org/ and http://www.informatics.jax.org/homology.shtml


#### 7. List of gene ontology terms for human genes: uniprot-human-genes and goterms mapping.tab
###### Source: http://geneontology.org/

## The scripts to analyze are broken up into two parts:
#### 1: Computation of estrogen signaling proxies per individual in GTEx to infer "low" vs "high" estrogen action: estrogen gene set across tissues.R
#### 2: Crosstissue correlation scores from Adipose, Skeletal Muscle and Stomach to Hypothalamus female-specific DEGs derived from mice: pathways and tissue-specific binning for low vs high estrogen.R
