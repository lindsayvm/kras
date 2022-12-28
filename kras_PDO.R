#Libraries
library(GSVA)
library(edgeR)
library(dplyr)
library(DESeq2)
library(limma)
library(ggsci)
library(ggpubr)
library(ggthemes)
library(tidyverse)
library(patchwork)
library(EnhancedVolcano)
library(biomaRt)  
# install.packages("devtools")
# devtools::install_github("stephenturner/annotables")
library(annotables)

setwd(dir = "/home/l.leek/src/kras/")

#' ###########################################################################
#' ###########################################################################
#' Input:  
#' ###########################################################################
#' ###########################################################################

# Import gene rna.m
fn = "/DATA/share/Voesties/data/organoid-data/data/kras/gcf_analysis/genecounts.txt"
rna.m = read_tsv(fn) %>% 
  filter(!is.na(ensembl_gene_id)) %>% 
  janitor::clean_names() %>% 
  dplyr::select(external_gene_id, x1_sns26:x85_tum65_5_fu_1_u_m) %>% 
  distinct(external_gene_id, .keep_all = T) %>% 
  column_to_rownames("external_gene_id") 

#Strip names
colnames(rna.m) = gsub("x[0-9][0-9]_|x[0-9]_","",colnames(rna.m))
colnames(rna.m) = gsub("_ftd_1_u_m","_ftd",colnames(rna.m))
colnames(rna.m) = gsub("_5_fu_1_u_m","_fu",colnames(rna.m))
colnames(rna.m) = toupper(paste0("KRAS_",colnames(rna.m)))
#Add genotype
#TUM42, TUM50, TUM65 = KRAS wt
#SNS26, TUM3, TUM10, TUM52 = KRAS G12x
colnames(rna.m)[str_detect(string = colnames(rna.m), pattern = "TUM42|TUM50|TUM65")] = paste0(colnames(rna.m)[str_detect(string = colnames(rna.m), pattern = "TUM42|TUM50|TUM65")],"_WT")
colnames(rna.m)[str_detect(string = colnames(rna.m), pattern = "SNS26|TUM3|TUM10|TUM52")] = paste0(colnames(rna.m)[str_detect(string = colnames(rna.m), pattern = "SNS26|TUM3|TUM10|TUM52")],"_G12")


# Make metadata
metadata = data.frame(Sample = colnames(rna.m),
                      Genotype = c(rep("G12", 9),
                                   rep("WT", 6),
                                   rep("G12", 3),
                                   rep("WT", 3)),
                      pdo = c(rep("SNS26",3),
                              rep("TUM3",3),
                              rep("TUM10",3),
                              rep("TUM42",3),
                              rep("TUM50",3),
                              rep("TUM52",3),
                              rep("TUM65",3)),
                      Treatment = rep(c("Untreated", "FTD", "FU"), 7)) %>% 
  mutate(Group = paste(Genotype, Treatment, sep = "-")) %>% 
  mutate(sample_id = Sample) %>% 
  column_to_rownames("sample_id")

#' ###########################################################################
#' ###########################################################################
#' Process: Normalize 
#' ###########################################################################
#' ###########################################################################

#edger object
# Normalized the data with Limma Voom
dge0 = DGEList(counts = as.matrix(rna.m),
               group = metadata$Group)

#make design
group = dge0$samples$group
design = model.matrix(~0 + group)

keep = filterByExpr(dge0, design)
dge = dge0[keep,]

dge = calcNormFactors(dge)

dge.voom = voom(counts = dge, 
                design = design, 
                plot = T)

# Write table to disk
dge.voom$E %>% 
  as.data.frame() %>% 
  rownames_to_column("Symbol") %>% 
  write_csv("rna.m-voom-normalized.csv")



#' ###########################################################################
#' ###########################################################################
#' Hallmark signature RNA
#' ###########################################################################
#' ###########################################################################

# Import hallmarks
hallmark = GSEABase::getGmt("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/h.all.v7.5.1.symbols.gmt")

# Run ssGSEA analysis
ssgsea.hallmark = gsva(expr = dge.voom$E, 
                       gset.idx.list = hallmark, 
                       verbose = T, 
                       parallel.sz = 12, 
                       method = "ssgsea")

# Convert to dataframe
ssgsea.hallmark.df = ssgsea.hallmark %>% 
  as.data.frame() %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") 

# Write to disk
ssgsea.hallmark.df %>% 
  left_join(metadata) %>%
  relocate(Sample, Genotype, Treatment, Group) %>% 
  write_csv("hallmark-signatures-ssgsea.csv")

# Convert to zscore and write to disk
ssgsea.hallmark.df %>%
  column_to_rownames("Sample") %>% 
  scale() %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  left_join(metadata) %>%
  relocate(Sample, Genotype, Treatment, Group) %>% 
  write_csv("hallmark-signatures-ssgsea-zscore.csv")


sig.df = ssgsea.hallmark.df %>%
  column_to_rownames("Sample") %>% 
  scale() %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  left_join(metadata) %>%
  relocate(Sample, Genotype, Treatment, Group) 
sig.df = sig.df %>% mutate(Group = factor(sig.df$Group , 
                                            levels=c( "WT-FTD", "G12-FTD",
                                                      "WT-Untreated", "G12-Untreated",
                                                      "WT-FU","G12-FU" )))




#' ###########################################################################
#' ###########################################################################
#' Output:  
#' ###########################################################################
#' ###########################################################################

#because you want to know about the interaction of genotype and treatment
lm(genotype:treatment (1-subject))


my_barplot = function(df, VAR){
  tmp.df = df[,c(VAR,"Group", "Treatment")]
  colnames(tmp.df) = c("interest","Group", "Treatment")
  anno_df = compare_means(interest ~ Group, 
                          group.by = "Treatment", 
                          data = tmp.df, 
                          method = "t.test",
                          p.adjust.method = "holm") %>%
    mutate(y_pos = 2.3)
  
  p = df %>% ggplot(aes(x = Group, y = .data[[VAR]], fill = Group)) +
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(position=position_dodge()) + 
    # geom_point(aes(color=Group)) +
    xlab("") +
    ylab(paste("ssGSEA", gsub("HALLMARK_|_"," ",VAR), "(zscore)")) + 
    ggsignif::geom_signif(
      data=anno_df, 
      inherit.aes = FALSE,
      aes(xmin=group1, xmax=group2, annotations=p.adj, y_position=y_pos), 
      manual=TRUE) +
    theme_bw() +
    
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 8),
          axis.title.x=element_blank(), 
          legend.position = "none")  +
    ylim(-2,2.5)
  
  return(p)
  
}



#Pathways related to cell cycle
vars = c("HALLMARK_DNA_REPAIR",
         "HALLMARK_KRAS_SIGNALING_UP",
         "HALLMARK_KRAS_SIGNALING_DN",
         "HALLMARK_P53_PATHWAY",
         "HALLMARK_MITOTIC_SPINDLE",
         "HALLMARK_G2M_CHECKPOINT")
p  = list()
#Plot
for(i in 1:length(vars)){
  print(i)
  p[[i]] = my_barplot(sig.df, vars[i])
}

(p[[1]]  | p[[2]]  | p[[3]] ) /
  (p[[4]]  | p[[5]]  | p[[6]] )

#Plot all
p  = list()
vars = colnames(sig.df)[str_detect(colnames(sig.df),"HALLMARK_")]
for(i in 1:length(vars)){
  print(i)
  p[[i]] = my_barplot(sig.df, vars[i])
}

p = do.call(marrangeGrob, args = list(grobs = p, ncol=4, nrow=2))
ggsave("/home/l.leek/src/kras/KRAS_PDO_mSigDBhallmark.pdf", p, width=12, height=8.5)
while (!is.null(dev.list())) dev.off()

