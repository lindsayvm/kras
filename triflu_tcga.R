library(data.table)
library("org.Hs.eg.db") 
library(edgeR)

#https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html
#https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf


##### Load RNAseq
fn = list.files(path = "/home/l.leek/data/triflu/tcga_RNA",
                pattern = ".htseq.counts$",
                full.names = TRUE,  
                recursive = TRUE)

tmp = gsub(".*\\/", "", fn)
names = gsub(".htseq.counts","", tmp)

# read
ls = lapply(fn, function(x) {
  df = fread(x)
  return(df)})
#Retrieve IDs
names(ls) = gsub(".*\\/|.driver.catalog.tsv$", "", fn)

#Concatenate
rna.df = ls[[1]][,1]
for(i in names(ls)){
  print(i)
  tmp.df = ls[[i]]
  colnames(tmp.df)[length(tmp.df)] = gsub(".htseq.counts","",names(ls[i]))
  rna.df = merge(rna.df, tmp.df, by = "V1")
}




###############
#EdgeR scaling/normalising

genes = rna.df$V1
rna.df = rna.df[ , -c(1)]

d = DGEList(counts = rna.df)

#recalculate the library sizes of the DGEList object AFTER the filtering
#RNA-seq measure relative expression rather than absolute expression. 
d$samples$lib.size = colSums(d$counts)

#Note that the “size factor” from DESeq is not equal to the “norm factor” in the edgeR. In edgeR, the library size and additional normalization scaling factors are separated. 
#normalizes the library sizes by finding a set of scaling factors for the library sizes that minimizes the log-fold changes between the samples for most genes.
d = calcNormFactors(d)

rna.df = as.data.frame(d$counts)



#Curate ENSEMBL name
rownames(rna.df) = genes
dim(rna.df)
# Unique ENSEMBL
length(unique(gsub("\\..*$","",rownames(rna.df))))
rownames(rna.df) = gsub("\\..*$","",rownames(rna.df))
length(unique(rownames(rna.df))) == length(rownames(rna.df))

########### DONT
# ENSEMBL is unique but there are multipel ENSEMBLs to 1 gene SYMBOL
# ===> therefore only do ENSEMBL --> gene at visualization
# https://www.researchgate.net/post/How-to-deal-with-multiple-ensemble-IDs-mapping-to-one-gene-symbol-in-a-RNA-Seq-dataset
# NOT merge the expression values of genes that have the same symbol. Gene symbols seem to shuffle around more haphazardly than Ensembl IDs or Entrez IDs, and so it is pretty common for completely different genes to be labeled with the same gene symbol depending on which gene annotation builds are being used.

#tmp = mapIds(org.Hs.eg.db, keys = rownames(rna.df), 
#                 keytype = "ENSEMBL", column="SYMBOL")
#tmp
#length(unique(tmp)) == nrow(rna.df)


###########SIGNATURES

fn = list.files(path = "/home/l.leek/data/signatures",
                pattern = "^HALLMARK",
                full.names = TRUE,  
                recursive = TRUE)
tmp = gsub(".*\\/", "", fn)
names = gsub(".txt","", tmp)

# read as list
ls = lapply(fn, function(x) {
  df = read.delim(x,
                  col.names = "gene",
                  skip = 1)
  df$gene = mapIds(org.Hs.eg.db, keys = df$gene, 
                      keytype = "SYMBOL", column="ENSEMBL")
  return(df)})
#Retrieve IDs
names(ls) = names


ls[[1]]
ls[["HALLMARK_ADIPOGENESIS"]]


tmp.df = rna.df[rownames(rna.df) %in% ls[[1]]$gene, ]
tmp.m = as.matrix(tmp.df[ ,!names(tmp.df) %in% c("mut" ,"gene")])  
finalSum.df = as.data.frame(colSums(tmp.m))
colnames(finalSum.df) = names(ls)[1]
for (i in 2:length(fn)){
  print(i)
  tmp.df = rna.df[rownames(rna.df) %in% ls[[i]]$gene, ]
  tmp.m = as.matrix(tmp.df[ ,!names(tmp.df) %in% c("mut" ,"gene")])  
  sum.df = as.data.frame(colSums(tmp.m))
  finalSum.df = cbind(finalSum.df, sum.df)
  colnames(finalSum.df)[length(finalSum.df)] = names(ls)[i]
}
finalSum.m = as.matrix(finalSum.df)

heatmap(finalSum.m)


###DUMMY ??????
mut = rep("WT", nrow(finalSum.m))
mut[0:(0.5*nrow(finalSum.m))] = "G12"
mut
mut_group = as.factor(mut)
mut.col = c("salmon", "skyblue")[mut_group] #WT is blue, red is MUT
heatmap(finalSum.m, Rowv = NA,
        labCol = colnames(finalSum.m), labRow = NA,
        RowSideColors = mut.col )


