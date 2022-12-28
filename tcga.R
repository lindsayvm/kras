library(data.table)
library(GenomicDataCommons)

rna.fn = list.files(path = "data/triflu/RNAseq",
                    pattern = "$",
                    full.names = TRUE,  
                    recursive = TRUE)

fread("data/triflu/RNAseq/")

GenomicDataCommons::status()
stopifnot(GenomicDataCommons::status()$status=="OK")

RNA_manifest = files() %>%
  filter( cases.project.project_id == 'TCGA-COAD') %>% 
  filter( type == 'gene_expression' ) %>%
  filter( analysis.workflow_type == 'HTSeq - Counts')  %>%
  manifest()
head(RNA_manifest)

WES_manifest = files() %>%
  filter( cases.project.project_id == 'TCGA-COAD') %>% 
  filter( type == 'gene_expression' ) %>%
  filter( analysis.workflow_type == 'HTSeq - Counts')  %>%
  manifest()
head(RNA_manifest)


fnames = lapply(RNA_manifest$id[1:20], gdcdata)

case_ids = cases() %>% results(size=10) %>% ids()
clindat = gdc_clinical(case_ids)
names(clindat)
