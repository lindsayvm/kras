library(stringr)
library(data.table)

mut = "KRAS"

##### Load WGS clinical data: 
clinical.df = fread("//DATA/share/Voesties/data/HMF/update_5/metadata/metadata.tsv")
clinical.df$treatment[str_detect(string = tolower(clinical.df$treatment), pattern = "trifluridine")]
clinical.df$treatment[clinical.df$treatment %in% c("Trifluridine", "Trifluridine/Trifluridine")]
clinical.df = clinical.df[clinical.df$treatment %in% c("Trifluridine", "Trifluridine/Trifluridine")]

##### Load WGS
#list all file names of mutational drivers
fn = list.files(path = "/DATA/share/Voesties/data/HMF/update_5/somatics/",
                pattern = ".driver.catalog.tsv$",
                full.names = TRUE,  
                recursive = TRUE)
tmp = gsub(".*\\/", "", fn)
tmp = gsub(".driver.catalog.tsv","", tmp)
tmp = gsub("[a-zA-Z.]","", tmp)
fn = fn[tmp %in% gsub("[a-zA-Z.]","", clinical.df$`#patientId`)]
#keep only first biopsy
gsub(".*\\/|.driver.catalog.tsv$", "", fn)
gsub(".*\\/|.driver.catalog.tsv$", "", fn)[str_detect(fn, "TII", negate = T)]
fn = fn[str_detect(fn, "TII", negate = T)]

# read
ls = lapply(fn, function(x) {
  df = read.table(file = x, sep = '\t', header = TRUE, stringsAsFactors = TRUE)
  return(df)})
#Retrieve IDs
names(ls) = gsub(".*\\/|.driver.catalog.tsv$", "", fn)

#concatenate all files
nrow = 0
ncol = ncol(ls[[1]]) +1 # +1 because want to add patientID
mut.df = as.data.frame(matrix(rep(NA, ncol * nrow), ncol = ncol, nrow = nrow))
colnames(mut.df) = colnames(ls[[1]])
for(i in 1:length(fn)){
  tmp.df = ls[[i]]
  tmp.df = tmp.df[tmp.df$gene == "KRAS", ]
  if (nrow(tmp.df) != 0){
    tmp.df = cbind(sampleId = names(ls[i]), tmp.df)
  }
  mut.df = rbind(mut.df, tmp.df)
}
colnames(mut.df)
summary(mut.df)


clinical.df$mut = rep("WT")
clinical.df$mut[clinical.df$sampleId %in% mut.df$sampleId] = "MUT"
clinical.df$mut

#RNAseq: list all file names of fusion genes (n = 328)
#rna_fn = list.files(path = "//DATA/share/Voesties/data/HMF/update_5/isofox/",
#                    pattern = ".isf.gene_data.csv$",
#                    full.names = TRUE, recursive = TRUE)
#clinical files of TANGO RNAseq patients
#rna.names = gsub(pattern = ".*\\/|.isf.gene_data.csv", "", rna_fn)

#RNA_ID = rna.names[rna.names %in% clinical.df$sampleId]
#RNAseq available?
#clinical.df$RNAseq = rep("No", nrow(clinical.df))
#clinical.df$RNAseq[clinical.df$sampleId %in% RNA_ID] = "Yes"
#clinical.df = clinical.df[ ,c("#patientId","sampleId","RNAseq","tumorPurity","primaryTumorLocation","treatment")]
#write.table(clinical.df, file = "trifluridine_WGS_RNAseq.tsv", sep = "\t", row.names = F)
#tmp = fread("/home/l.leek/data/trifluridine_WGS_RNAseq.tsv")
#rsync -avhzP -e 'ssh -o "ProxyCommand ssh -A l.leek@rhpc.nki.nl -W %h:%p"' l.leek@darwin:/home/l.leek/data/trifluridine_WGS_RNAseq.tsv ~/Desktop/


#######RNA DATA
rna.df =  fread("/DATA/share/Voesties/data/rnaseq_hmf/HMF_RNAseq_after_QC/counts_full_QC_pass_vst_normalized_with_clinical_data_batch_corrected_for_biopsyLocation.txt", 
                na.strings = "")
rna.df = rna.df %>% dplyr::select(c("V1", colnames(rna.df)[colnames(rna.df) %in% clinical.df$sampleId]))

#### !!!!!!! same patients in all df
clinical.df = clinical.df[clinical.df$sampleId %in% colnames(rna.df), ]
mut.df = mut.df[mut.df$sampleId %in% colnames(rna.df), ]



######## SIGNATURES
DNA_REPAIR.df = read.delim("/home/l.leek/data/signatures/HALLMARK_DNA_REPAIR.txt",
                    col.names = "gene")
G2M_CHECKPOINT.df = read.delim("/home/l.leek/data/signatures/HALLMARK_G2M_CHECKPOINT.txt",
                     col.names = "gene")
KRAS_DOWN.df = read.delim("/home/l.leek/data/signatures/HALLMARK_KRAS_DOWN.txt",
                     col.names = "gene")
KRAS_UP.df = read.delim("/home/l.leek/data/signatures/HALLMARK_KRAS_UP.txt",
                     col.names = "gene")

DNA_REPAIR.df = rna.df[rna.df$V1 %in% DNA_REPAIR.df$gene, ]
DNA_REPAIR.m = as.matrix(DNA_REPAIR.df[ ,-c("V1")])
clin_group = as.factor(clinical.df$mut)#blauw is CB
#mut.col = brewer.pal(length(unique(clin_group)), "Set1")[clin_group]
mut.col = c("salmon", "skyblue")[clin_group] #WT is blue, red is MUT
mut.col = mut.col[!is.na(mut.col)]
heatmap(DNA_REPAIR.m, Colv = TRUE, Rowv = NA,
        labRow = DNA_REPAIR.df$V1, labCol = colnames(DNA_REPAIR.m),
        ColSideColors = mut.col )
#legend(x="bottomleft", legend=unique(clinical.df$mut),fill= unique(mut.col) , ncol = 3)

G2M_CHECKPOINT.df = rna.df[rna.df$V1 %in% G2M_CHECKPOINT.df$gene, ]
G2M_CHECKPOINT.m = as.matrix(G2M_CHECKPOINT.df[ ,-c("V1")])
clin_group = as.factor(clinical.df$mut)#blauw is CB
#mut.col = brewer.pal(length(unique(clin_group)), "Set1")[clin_group]
mut.col = c("salmon", "skyblue")[clin_group] #WT is blue, red is MUT
mut.col = mut.col[!is.na(mut.col)]
heatmap(G2M_CHECKPOINT.m, Colv = TRUE, Rowv = NA,
        labRow = G2M_CHECKPOINT.df$V1, labCol = colnames(G2M_CHECKPOINT.m),
        ColSideColors = mut.col )
#legend(x="bottomleft", legend=unique(clinical.df$mut),fill= unique(mut.col) , ncol = 3)


KRAS_DOWN.df = rna.df[rna.df$V1 %in% KRAS_DOWN.df$gene, ]
KRAS_DOWN.m = as.matrix(KRAS_DOWN.df[ ,-c("V1")])
clin_group = as.factor(clinical.df$mut)#blauw is CB
#mut.col = brewer.pal(length(unique(clin_group)), "Set1")[clin_group]
mut.col = c("salmon", "skyblue")[clin_group] #WT is blue, red is MUT
mut.col = mut.col[!is.na(mut.col)]
heatmap(KRAS_DOWN.m, Colv = TRUE, Rowv = NA,
        labRow = KRAS_DOWN.df$V1, labCol = colnames(KRAS_DOWN.m),
        ColSideColors = mut.col )
#legend(x="bottomleft", legend=unique(clinical.df$mut),fill= unique(mut.col) , ncol = 3)


KRAS_UP.df = rna.df[rna.df$V1 %in% KRAS_UP.df$gene, ]
KRAS_UP.m = as.matrix(KRAS_UP.df[ ,-c("V1")])
clin_group = as.factor(clinical.df$mut)#blauw is CB
#mut.col = brewer.pal(length(unique(clin_group)), "Set1")[clin_group]
mut.col = c("salmon", "skyblue")[clin_group] #WT is blue, red is MUT
mut.col = mut.col[!is.na(mut.col)]
heatmap(KRAS_UP.m, Colv = TRUE, Rowv = NA,
        labRow = KRAS_UP.df$V1, labCol = colnames(KRAS_UP.m),
        ColSideColors = mut.col )
#legend(x="bottomleft", legend=unique(clinical.df$mut),fill= unique(mut.col) , ncol = 3)



#################################
#SUM
sign.df = t(data.frame(DNA_REPAIR = colSums(DNA_REPAIR.m), 
             G2M_CHECKPOINT = colSums(G2M_CHECKPOINT.m),
             KRAS_DOWN = colSums(KRAS_DOWN.m),
             KRAS_UP = colSums(KRAS_UP.m)))
#MEAN
sign.df = t(data.frame(DNA_REPAIR = colMeans(DNA_REPAIR.m), 
                       G2M_CHECKPOINT = colMeans(G2M_CHECKPOINT.m),
                       KRAS_DOWN = colMeans(KRAS_DOWN.m),
                       KRAS_UP = colMeans(KRAS_UP.m)))
#AVERAGE
sign.df = t(data.frame(DNA_REPAIR = colMeans(DNA_REPAIR.m, na.rm = T), 
                       G2M_CHECKPOINT = colMeans(G2M_CHECKPOINT.m, na.rm = T),
                       KRAS_DOWN = colMeans(KRAS_DOWN.m, na.rm = T),
                       KRAS_UP = colMeans(KRAS_UP.m, na.rm = T)))

clin_group = as.factor(clinical.df$mut)#blauw is CB
mut.col = c("salmon", "skyblue")[clin_group] #WT is blue, red is MUT
heatmap(sign.df, Colv = TRUE, Rowv = TRUE,
        labRow = rownames(sign.df), labCol = colnames(sign.df),
        ColSideColors = mut.col )
#legend(x="bottomleft", legend=unique(clinical.df$mut),fill= unique(mut.col) , ncol = 3)

#write.table(mut.df, file = "trifluridine_WGS_KRAS.tsv", sep = "\t", row.names = F)
#tmp = fread("/home/l.leek/data/trifluridine_WGS_KRAS.tsv")
#rsync -avhzP -e 'ssh -o "ProxyCommand ssh -A l.leek@rhpc.nki.nl -W %h:%p"' l.leek@darwin:/home/l.leek/data/trifluridine_WGS_KRAS.tsv ~/Desktop/

