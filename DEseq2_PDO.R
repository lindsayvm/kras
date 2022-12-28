


#' ###########################################################################
#' ###########################################################################
#' Process: DEG  
#' ###########################################################################
#' ###########################################################################

#DEseq object from raw data
dds = DESeqDataSetFromMatrix(countData = rna.m,
                             colData = metadata,
                             design = ~ pdo + Treatment)
#deseq2 normalization
dds = DESeq(dds)

cutoff = 1
keep = rowSums(counts(dds)) >= cutoff
dds = dds[keep,]

#comparison
dds.g12_ftd = dds[,dds$Group %in% c("WT-Control", "G12-Control", "WT-FTD", "G12-FTD")]
dds.g12_ftd$Group = droplevels(dds.g12_ftd$Group)
dds.g12_ftd$Genotype = as.factor(dds.g12_ftd$Genotype)
dds.g12_ftd$Treatment = as.factor(dds.g12_ftd$Treatment)
design(dds.g12_ftd) = ~Genotype*Treatment + Genotype:Treatment


# Get table of results
res = DESeq2::results(dds.g12_ftd)

res$Geneid = rownames(res) 
res$Geneid = gsub("\\..*","",res$Geneid)
res = as.data.frame(res)
res = res %>% 
  dplyr::left_join(grch38, by = c("Geneid" = "symbol")) %>% 
  dplyr::filter(biotype == "protein_coding") %>% 
  dplyr::filter(!is.na(padj))
#
EnhancedVolcano(res,
                lab = res$Geneid,
                x = 'log2FoldChange',
                y = 'padj',
                title = '',
                pCutoff = 0.01,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 4.0,#label size
                col=c('lightgrey', '#006994', 'darkgrey', "#ffdb58"),
                colAlpha = 1,
                legendLabels=c('Non-significant','Log2 FC','P-value',
                               'P-value & Log2 FC'),
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 5.0) 

#PCA
dds.vst = DESeq2::vst(dds)
DESeq2::plotPCA(dds.vst, "pdo")


