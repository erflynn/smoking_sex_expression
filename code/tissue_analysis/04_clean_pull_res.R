# Code for pulling the results


library(tidyverse)
load("data/all_studies_v2.RData") # all_studies


fix_colnames_ds <- function(ds){
  print(ds$study)
  if (str_detect(ds$study, "E-MTAB")){
    if (nrow(ds$gene_smok)==0){
      ds$gene_smok <- ds$probe_smok %>% 
        mutate(gene=probes,
               chromosome="") %>%
        rename(adj.p=adj.P.Val,
               logFC.l=CI.L,
               p=P.Value) 
      
    } 
    if (nrow(ds$gene_int)==0){
      ds$gene_int <- ds$probe_int %>% 
        mutate(gene=probes,
               chromosome="") %>%
        rename(adj.p=adj.P.Val,
               logFC.l=CI.L,
               p=P.Value)
    }  
    if (nrow(ds$gene_sex)==0){
      ds$gene_sex <- ds$probe_sex %>% 
        mutate(gene=probes,
               chromosome="") %>%
        rename(adj.p=adj.P.Val,
               logFC.l=CI.L,
               p=P.Value)
      
    }
  }

  if ("geneSymbol" %in% colnames(ds$gene_smok)){
      ds$gene_smok <- ds$gene_smok %>%
        dplyr::rename(gene=geneSymbol, p=P.Value,
                      probes=ID) %>%
        mutate(adj.p=p.adjust(p, "fdr")) %>%
        mutate(logFC.l=logFC-SD*1.96)
  }

  if (!all(is.na(ds$gene_int))){
    
    if ("geneSymbol" %in% colnames(ds$gene_int)){
        ds$gene_int <- ds$gene_int %>%
          dplyr::rename(gene=geneSymbol, p=P.Value,
                        probes=ID) %>%
          mutate(adj.p=p.adjust(p, "fdr")) %>%
          mutate(logFC.l=logFC-SD*1.96) 
    }

  }

  if (!all(is.na(ds$gene_sex))){
   
    if ("geneSymbol" %in% colnames(ds$gene_sex)){
          ds$gene_sex <- ds$gene_sex %>%
            dplyr::rename(gene=geneSymbol, p=P.Value,
                          probes=ID) %>%
            mutate(adj.p=p.adjust(p, "fdr")) %>%
            mutate(logFC.l=logFC-SD*1.96) 
      
    }
  
  }
    return(ds)
}
names(all_studies) <-  lapply(all_studies, function(x) x$study)

all_studies2 <- lapply(all_studies, fix_colnames_ds)
names(all_studies2) <- names(all_studies)

save(all_studies2, file="data/all_studies_v2_clean.RData") 


# all_studies[["GSE20189"]]$gene_int %>% 
#   filter(gene %in% c("NMOX1", "TAGLN", "SYNE1", "LGR6", "EPS8",
#                      "RGS9", "MYOF", "ZNF618", "RARRES3", "TCFL2",
#                      "CHD1L", "TMEM106C", "PLIN5", "TGFBR3", "CSF1R",
#                      "PPP2R2B", "PLCD4", "LOXL3", "BFSP1"))
# all_studies[["GSE46699"]]$gene_int %>% 
#   filter(gene %in% c("NMOX1", "TAGLN", "SYNE1", "LGR6", "EPS8",
#                      "RGS9", "MYOF", "ZNF618", "RARRES3", "TCFL2",
#                      "CHD1L", "TMEM106C", "PLIN5", "TGFBR3", "CSF1R",
#                      "PPP2R2B", "PLCD4", "LOXL3", "BFSP1"))
# table info:
# study, prop_var_smok, prop_var_sex, prop_var_int, # de probes, # de genes, # de probes int, # de genes int

table_ds <- function(ds){
  num_tot_probes <- ds$probe_smok %>% nrow()
  num_tot_genes <- ds$gene_smok %>% nrow()
  num_de_probes <- ds$probe_smok %>%
    filter(abs(logFC) >= 0.3 & adj.P.Val < 0.05 ) %>%
    nrow()

  num_de_genes <- ds$gene_smok %>%
    filter(abs(logFC) >= 0.3 & adj.p < 0.05 ) %>%
    nrow()
  
  if (!is.na(ds$probe_sex)){
    num_de_probes_sex <- ds$probe_sex %>%
      filter(abs(logFC) >= 0.3 & adj.P.Val < 0.05 ) %>%
      nrow()
    
    num_de_genes_sex <- ds$gene_sex %>%
      filter(abs(logFC) >= 0.3 & adj.p < 0.05 ) %>%
      nrow()
  } else {
    num_de_probes_sex <- NA
    num_de_genes_sex <- NA
  }
  
  
  
  if (!is.na(ds$probe_int)){
    num_de_probes_int <- ds$probe_int %>%
      filter(abs(logFC) >= 0.3 & adj.P.Val < 0.05 ) %>%
      nrow()

    num_de_genes_int <- ds$gene_int %>%
      filter(abs(logFC) >= 0.3 & adj.p < 0.05 ) %>%
      nrow()
  } else {
    num_de_probes_int <- NA
    num_de_genes_int <- NA
  }
  sex_n = NA
  sex_n_g = NA
  int_n_g = NA
  int_n = NA
  smok_n = NA
  smok_n_g = NA
  if (!is.na(ds$var)){
    var_smok =  ds$var[["smoking"]]
    if ( !is.na(ds$rand)){
      rand_smok = ds$rand %>% filter(covariate=="smoking") 
      smok_n = rand_smok %>% nrow()
      smok_n_g = rand_smok %>% filter(variance > var_smok) %>% nrow()
    }
    if ("sex" %in% names(ds$var)){
      var_sex = ds$var[["sex"]]
      var_int= ds$var[["sex:smoking"]]
      if ( !is.na(ds$rand)){
        rand_sex = ds$rand %>% filter(covariate=="sex") 
        sex_n = rand_sex %>% nrow()
        sex_n_g = rand_sex %>% filter(variance > var_sex) %>% nrow()
        rand_int = ds$rand %>% filter(covariate=="sex:smoking") 
        int_n = rand_int %>% nrow()
        int_n_g = rand_int %>% filter(variance > var_int) %>% nrow()
      }
    } else {
      var_sex =NA
      var_int= NA
    }
  } else {
    var_smok = NA
    var_sex = NA
    var_int = NA
  }

  return(list(ds$study, var_smok, var_sex, var_int, 
              smok_n, smok_n_g, sex_n, sex_n_g, int_n, int_n_g,
              num_tot_probes, num_tot_genes, num_de_probes, 
              num_de_genes, num_de_probes_sex, num_de_genes_sex,
              num_de_probes_int, num_de_genes_int))
}


# DE GENES

gene_lists <- function(ds){
  print(ds$study)
  sig_gene <- ds$gene_smok %>%
    filter(abs(logFC) >= 0.3 & adj.p < 0.05 )
  
  if (!is.na(ds$gene_int)){
    sig_int <- ds$gene_int %>%
      filter(abs(logFC) >= 0.3 & adj.p < 0.05 )
  } else {
    sig_int <- NA
  }
  
  if (!is.na(ds$gene_sex)){
    sig_sex <- ds$gene_sex %>%
      filter(abs(logFC) >= 0.3 & adj.p < 0.05 )
  } else {
    sig_sex <- NA
  }
 return(list("sig_gene"=sig_gene, "sig_sex"=sig_sex, "sig_int"=sig_int))
}

my_gene_l <- lapply(all_studies2, gene_lists)
names(my_gene_l) <- names(all_studies2)


my_tab <- do.call(rbind, lapply(all_studies2, table_ds))
data_df <- data.frame(my_tab)
head(data_df)
colnames(data_df) <- c("study",  "var_smok", "var_sex", "var_int", 
                       "smok_n", "smok_n_g", "sex_n", "sex_n_g", "int_n", "int_n_g",
                       "nprobes", "ngenes", "de_probes_smok", 
                       "de_genes_smok", "de_probes_sex", "de_genes_sex",
                       "de_probes_int", "de_genes_int")




data_df2 <- data_df %>% 
  mutate(across(-study, ~as.numeric(as.character(.)))) %>%
  as_tibble() %>%
  mutate(study=unlist(study))

data_df2 %>% select(study, contains("de_genes")) %>% View()


list_studies <- read_csv("data/table1b_updated.csv")

list_studies <- list_studies %>%
  filter(study!="GSE44456") %>%
  mutate(tissue=factor(tissue, levels=c("airway epithelium", "trachea epithelium",
                                        "nasal epithelium", "oral cavity",
                                        "buccal mucosa", "sputum", "alveolar macrophages", 
                                        "lung", 
                                        "blood - b cells", "blood - pbmcs", "blood - whole",
                                        "liver", "kidney", 
                                        "brain - prefrontal cortex")))


data_df3 <- data_df2 %>%
  mutate(across(contains("var"), ~signif(., 3))) %>%
  inner_join(list_studies %>% dplyr::select(study, tissue)) %>%
  dplyr::select(study, tissue, everything()) %>%
  arrange(tissue, study) 
data_df4 <- data_df3 %>% mutate(var_smok_p=smok_n_g/smok_n,
                    var_sex_p=sex_n_g/sex_n,
                    var_int_p=int_n_g/int_n) %>%
  select(-contains("_n")) %>%
  select(study,tissue, contains("var"), everything())

data_df4 %>% write_csv("data/supp_tables/small_study_summary_v2.csv")




### TODO save these objects



