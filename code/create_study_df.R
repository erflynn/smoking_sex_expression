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
head(my_gene_l)
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

# ADD AE HERE!

# --- ae data --- #
load("data/results/ae_int.RData") 
load("data/results/ae_smok.RData")
probe_int_ae <- smok_ci_int_p
gene_int_ae <- comb_ma_dat %>%
  ungroup() %>%
  dplyr::select(-chromosome)
probe_smok_ae <- smok_ci_smok 
gene_smok_ae <- ma_smok %>%
  ungroup() %>%
  dplyr::select(-chromosome)
probe_sex_ae <- smok_ci_sex 
gene_sex_ae <- ma_sex %>%
  ungroup() %>%
  dplyr::select(-chromosome)

# overlap + correlation of ES?
ae_de_probe <- probe_smok_ae %>% 
  filter(abs(logFC) >= 0.3 & adj.P.Val < 0.05 ) # 2625

ae_de_gene <- gene_smok_ae %>%
  filter(abs(logFC) >= 0.3 & adj.p < 0.05 ) # 932

ae_de_probe_i <- probe_int_ae %>% 
  filter(abs(logFC) >= 0.3 & adj.P.Val < 0.05 )  # 1
ae_de_gene_i <- gene_int_ae %>%
  filter(abs(logFC) >= 0.3 & adj.p < 0.05 ) # 30

ae_de_probe_s <- probe_sex_ae %>% 
  filter(abs(logFC) >= 0.3 & adj.P.Val < 0.05 )  # 128
ae_de_gene_s <- gene_sex_ae %>%
  filter(abs(logFC) >= 0.3 & adj.p < 0.05 ) # 48


# where is the grouped ae data from
ae_l <-  list("Grouped AE", "airway epithelium", 0.0215, 0.0211, 0.00666, 54675,20183,2625,932,1,30)
names(ae_l) <- colnames(data_df3)
ae_l2 <- rbind(data_df3, ae_l) %>% arrange(tissue, study)



df <- ae_l2 %>% left_join(list_studies %>% 
                            dplyr::select(study, samples, females, smokers)) %>%
  mutate(frac_smokers=smokers/samples) %>%
  filter(!is.na(samples))
cor.test(df$samples, df$de_genes)
cor.test(df$frac_smokers, df$de_genes)

ae_l2 %>% left_join(list_studies %>% 
                      dplyr::select(study, samples, females, smokers)) %>%
  mutate(frac_smokers=smokers/samples) %>%
  mutate(`fraction female`=females/samples) %>%
  ggplot(aes(x=frac_smokers, y=de_genes, col=`fraction female`))+
  geom_point(aes(size=samples), alpha=0.9)+
  theme_bw()+
  scale_fill_continuous(low = "purple", high = "green", limit = c(0,1))+
  xlab("fraction smokers")+
  ylab("number of DE genes")
ggsave("figures/paper_figs/sample_frac.png")  


ae_l2 %>% left_join(list_studies %>% 
                      dplyr::select(study, samples, females)) %>%
  mutate(frac_females=females/samples) %>%
  ggplot(aes(x=frac_females, y=de_genes, size=samples))+
  geom_point()+
  theme_bw()+
  xlab("fraction female")+
  ylab("number of DE genes")

l_b <- read_csv("data/supp_tables/s4b_breakdown_smok_sex.csv")
ggplot(l_b %>% pivot_longer(contains("male"), names_to="group", values_to="count") %>%
         mutate(group=str_replace_all(group,"_NS", " non-smokers"),
                group=str_replace_all(group, "_S", " smokers")), 
       aes(x=study_acc, y=count,fill=group))+
  geom_bar(stat="identity")+
  theme_bw()+
  xlab("")+
  theme(axis.text.x = element_text(angle = 90, hjust=1))
ggsave("figures/paper_figs/sex_breakdown_sm.png")
# fraction male vs fraction female

#data_df4 %>% filter(var_smok_p < 0.05) %>% nrow() # 16
#data_df4 %>% filter(var_sex_p < 0.05) %>% nrow() #15
#data_df4 %>% filter(var_int_p < 0.05) %>% nrow()  #3


ae_l2$study <- factor(ae_l2$study, levels=ae_l2$study)
data_df4 %>% 
  dplyr::select(study, tissue, contains("var")) %>%
  pivot_longer(contains("var"), names_to="term", values_to="variance") %>%
  mutate(term=case_when(
    term=="var_int" ~ "smoking*sex",
    term=="var_sex" ~ "sex",
    term=="var_smok" ~ "smoking"
  )) %>%
  filter(!is.na(term)) %>%
  ggplot(aes(x=1, y=variance, fill=term)) +
  geom_bar(stat="identity", position=position_dodge(0.8))+
  xlab("")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank())+
  coord_flip()+
  facet_grid(study~., switch="y")+
  theme(strip.text.y.left = element_text(angle = 0))
ggsave("figures/paper_figs/var_frac.png")

ggplot(ae_l2, aes(x=study, y=de_genes))+
  geom_bar(stat="identity")+
  theme_bw()

ggplot(ae_l2, aes(x=study, y=de_genes_int))+
  geom_bar(stat="identity")+
  theme_bw()

# add AE to the list?
data_df3 %>% filter(de_probes==0 | de_genes==0) # 4 have both, 9 have one of them


get_pair_overlap_auto <- function(ds1, study1, ds2, study2){
  # look at gene overlap
  ds_smok <- ds1 %>% 
    ungroup() %>%
    dplyr::select(-chromosome) %>%
    mutate(id=1:nrow(ds1)) %>% 
    separate_rows(gene, sep=",")  %>%
    filter(!gene %in% c(xy_chr_map$hgnc_symbol, "TTTY15", "KAL1", "ARSE", "ARSL","LOC389906"))
  ds_de_gene <- ds_smok %>%
    filter(abs(logFC) >= 0.3 & adj.p < 0.05 ) %>%
    pull(gene)
  ds2_de_gene <- ds2 %>%
    filter(abs(logFC) >= 0.3 & adj.p < 0.05 )  
  
  overlap <- ds_smok %>%
    filter(gene %in% ds_de_gene | gene %in% ds2_de_gene$gene) %>%
    inner_join(ds2, by="gene") # 702
  
  ggplot(overlap, aes(x=logFC.x, y=logFC.y))+
    geom_point(alpha=0.5)+
    theme_bw()
  if (nrow(overlap) < 10){
    cor <- NA
    corp <- NA
  } else {
    pcor <- cor.test(overlap$logFC.x, overlap$logFC.y, method="pearson") # 0.00042, p=0.99
    cor <- pcor$estimate[["cor"]]
    corp <- pcor$p.value
  }
  
  same_dir <- overlap %>% filter(logFC.x*logFC.y > 0) # 345/702
  rep_sig <- same_dir %>% filter(p.x < 0.05/nrow(overlap) & p.y < 0.05/nrow(overlap))
  
  return(list("study1"=study1, "study2"=study2, "tot"=nrow(overlap), "cor"=cor, "cor.p"=corp,
              "same_dir"=nrow(same_dir), "rep_sig"=nrow(rep_sig), 
              "genes"=paste(rep_sig$gene, collapse=";"),
              "gene_es"=paste(round(rep_sig$logFC.x,3), collapse=";")))
}


get_pair_overlap <- function(ds1, study1, ds2, study2){
  # look at gene overlap
  ds_smok <- ds1 %>% 
    ungroup() %>%
    dplyr::select(-chromosome) %>%
    mutate(id=1:nrow(ds1)) %>% 
    separate_rows(gene, sep=",") 
  ds_de_gene <- ds_smok %>%
    filter(abs(logFC) >= 0.3 & adj.p < 0.05 ) %>%
    pull(gene)
  ds2_de_gene <- ds2 %>%
    filter(abs(logFC) >= 0.3 & adj.p < 0.05 )  
  
  overlap <- ds_smok %>%
    filter(gene %in% ds_de_gene | gene %in% ds2_de_gene$gene) %>%
    inner_join(ds2, by="gene") # 702

  ggplot(overlap, aes(x=logFC.x, y=logFC.y))+
    geom_point(alpha=0.5)+
    theme_bw()
  if (nrow(overlap) < 10){
   cor <- NA
   corp <- NA
  } else {
    pcor <- cor.test(overlap$logFC.x, overlap$logFC.y, method="pearson") # 0.00042, p=0.99
    cor <- pcor$estimate[["cor"]]
    corp <- pcor$p.value
  }

  same_dir <- overlap %>% filter(logFC.x*logFC.y > 0) # 345/702
  rep_sig <- same_dir %>% filter(p.x < 0.05/nrow(overlap) & p.y < 0.05/nrow(overlap))
  
  return(list("study1"=study1, "study2"=study2, "tot"=nrow(overlap), "cor"=cor, "cor.p"=corp,
              "same_dir"=nrow(same_dir), "rep_sig"=nrow(rep_sig), 
              "genes"=paste(rep_sig$gene, collapse=";"),
              "gene_es"=paste(round(rep_sig$logFC.x,3), collapse=";")))
}



res_smok <- lapply(all_studies2, function(x) get_pair_overlap(x$gene_smok,
                                                      x$study, 
                                                      gene_smok_ae,
                                                      "AE"))
ae_overlap <- data.table::rbindlist(res_smok)
ae_overlap_df <- ae_overlap %>%
  dplyr::select(-study2) %>%
  mutate(across(c(-study1,-genes), ~as.numeric(as.character(.)))) %>%
  mutate(study=unlist(study1)) %>%
  mutate(across(contains("cor"), ~signif(., 3))) %>%
  left_join(list_studies %>% dplyr::select(study, tissue, platform)) %>%
  dplyr::select(study, tissue, platform, everything()) %>%
  arrange(tissue, study)
ae_overlap_df %>% 
  #mutate(same_dir=round(same_dir/tot, digits=3)) %>%
  #dplyr::select(-cor.p) %>%
  write_csv("data/supp_tables/ae_overlap_results_smok2.csv")



res_sex <- lapply(all_studies2[!(names(all_studies2) %in% c("GSE16149", "GSE18723", "GSE87072"))], 
                  function(x) {
                    print(x$study);
                    get_pair_overlap(x$gene_sex, x$study, gene_sex_ae, "AE")})
ae_overlap_sex <- data.table::rbindlist(res_sex)
ae_overlap_df_sex <- ae_overlap_sex %>%
  dplyr::select(-study2) %>%
  mutate(across(c(-study1,-genes), ~as.numeric(as.character(.)))) %>%
  mutate(study=unlist(study1)) %>%
  mutate(across(contains("cor"), ~signif(., 3))) %>%
  left_join(list_studies %>% dplyr::select(study, tissue, platform)) %>%
  dplyr::select(study, tissue, platform, everything()) %>%
  arrange(tissue, study)
ae_overlap_df_sex %>% 
  #mutate(same_dir=round(same_dir/tot, digits=3)) %>%
  #dplyr::select(-cor.p) %>%
  write_csv("data/supp_tables/ae_overlap_results_sex2.csv")


res_sex_auto <- lapply(all_studies2[!(names(all_studies2) %in% c("GSE16149", "GSE18723", "GSE87072"))], 
                  function(x) {
                    print(x$study);
                    get_pair_overlap_auto(x$gene_sex, x$study, gene_sex_ae, "AE")})
ae_overlap_sex_auto <- data.table::rbindlist(res_sex_auto)
ae_overlap_df_sex_auto <- ae_overlap_sex_auto %>%
  dplyr::select(-study2) %>%
  mutate(across(c(-study1,-genes), ~as.numeric(as.character(.)))) %>%
  mutate(study=unlist(study1)) %>%
  mutate(across(contains("cor"), ~signif(., 3))) %>%
  left_join(list_studies %>% dplyr::select(study, tissue, platform)) %>%
  dplyr::select(study, tissue, platform, everything()) %>%
  arrange(tissue, study)

count_ae_overlap2 <- ae_overlap_df %>% 
  dplyr::select(study, tissue, genes) %>%
  as_tibble() %>%
  filter(genes!="") %>%
  separate_rows(genes, sep=";") %>%
  group_by(genes) %>%
  filter(study != "Grouped AE") %>%
  dplyr::summarize(tissue=paste(unique(tissue), collapse=";"),
            study=paste(study, collapse=";"),
            n=n()) %>%
  arrange(desc(n)) %>%
  filter(n>1)

count_ae_overlap2 %>% rename(gene=genes, nstudies=n) %>%
  left_join(gene_smok_ae, by="gene") %>% 
  left_join(chr_l %>% filter(chromosome_name %in% c(1:22, "X", "Y")), 
            by=c("gene"="hgnc_symbol")) %>%
  dplyr::select( -src, -p, -logFC.l, -logFC.u) %>% 
  write_csv("data/supp_tables/ae_overlap2v2.csv")

count_ae_overlap2_sex <- ae_overlap_df_sex %>% 
  dplyr::select(study, tissue, genes) %>%
  as_tibble() %>%
  filter(genes!="") %>%
  filter(study != "Grouped AE") %>%
  separate_rows(genes, sep=";") %>%
  group_by(genes) %>%
  dplyr::summarize(tissue=paste(unique(tissue), collapse=";"),
            study=paste(study, collapse=";"),
            n=n()) %>%
  arrange(desc(n)) %>%
  filter(n>1)
head(chr_l)
count_ae_overlap2_sex %>% rename(gene=genes, nstudies=n) %>%
  left_join(chr_l %>% filter(chromosome_name %in% c(1:22, "X", "Y")), 
            by=c("gene"="hgnc_symbol")) %>%
  left_join(gene_sex_ae, by="gene") %>% 
  dplyr::select( -src, -p, -logFC.l, -logFC.u) %>% 
  write_csv("data/supp_tables/ae_overlap2v2_sex.csv")
count_ae_overlap2 %>%
  filter(str_detect(tissue, "airway"), str_detect(tissue, "blood"))

count_ae_overlap2 %>%
  filter(str_detect(tissue, "airway"), str_detect(tissue, "brain"))
count_ae_overlap2 %>%
  filter(str_detect(tissue, "airway"), str_detect(tissue, "kidney"))
# make a correlation plot

all_pairs <- combn(names(all_studies2),2)
all_pairs_sex <- combn(setdiff(names(all_studies2),c("GSE16149", "GSE18723", "GSE87072")),2)

all_pair_overlap <- apply(all_pairs, 2, function(x) {
  study1 <- x[1]
  study2 <- x[2]
  print(sprintf("%s,%s", study1, study2))
  get_pair_overlap(all_studies2[[study1]]$gene_smok, study1, 
                   all_studies2[[study2]]$gene_smok, study2)
})
all_pair_overlap_sex <- apply(all_pairs_sex, 2, function(x) {
  study1 <- x[1]
  study2 <- x[2]
  print(sprintf("%s,%s", study1, study2))
  get_pair_overlap(all_studies2[[study1]]$gene_sex, study1, 
                   all_studies2[[study2]]$gene_sex, study2)
})

all_pair_overlap_sex_auto <- apply(all_pairs_sex, 2, function(x) {
  study1 <- x[1]
  study2 <- x[2]
  print(sprintf("%s,%s", study1, study2))
  get_pair_overlap_auto(all_studies2[[study1]]$gene_sex, study1, 
                   all_studies2[[study2]]$gene_sex, study2)
})



all_pair_df <- data.table::rbindlist(all_pair_overlap)
all_pair_df_sex <- data.table::rbindlist(all_pair_overlap_sex)
all_pair_df_sex_auto <- data.table::rbindlist(all_pair_overlap_sex_auto)

full_pair <- ae_overlap %>% 
  dplyr::select(-study2) %>%
  dplyr::rename(study2=study1)  %>%
  mutate(study1="Grouped AE") %>% 
  dplyr::select(study1, everything()) %>%
  bind_rows(all_pair_df)

full_pair_sex <- ae_overlap_sex %>% 
  dplyr::select(-study2) %>%
  dplyr::rename(study2=study1)  %>%
  mutate(study1="Grouped AE") %>% 
  dplyr::select(study1, everything()) %>%
  bind_rows(all_pair_df_sex)

full_pair_sex_auto <- ae_overlap_sex_auto %>% 
  dplyr::select(-study2) %>%
  dplyr::rename(study2=study1)  %>%
  mutate(study1="Grouped AE") %>% 
  dplyr::select(study1, everything()) %>%
  bind_rows(all_pair_df_sex_auto)


all_pair_overlap_g <- full_pair %>% filter(genes!="") %>% 
  dplyr::select(study1, study2, genes, gene_es) %>%
  as_tibble() %>%
  separate_rows(genes, sep=";") %>%
  group_by(study1, study2) %>%
  mutate(n=1:n()) %>%
  group_by(study1, study2, n) %>%
  mutate(gene_effect=as.numeric(str_split(gene_es, ";")[[1]][[n]])) %>%
  ungroup() %>%
  dplyr::select(-gene_es, -n) %>%
  dplyr::rename(gene=genes) %>%
  mutate(direction=ifelse(gene_effect < 0, "down", "up")) %>%
  dplyr::select(-gene_effect)


all_pair_overlap_g_sex <- full_pair_sex %>% filter(genes!="") %>% 
  dplyr::select(study1, study2, genes, gene_es) %>%
  as_tibble() %>%
  separate_rows(genes, sep=";") %>%
  group_by(study1, study2) %>%
  mutate(n=1:n()) %>%
  group_by(study1, study2, n) %>%
  mutate(gene_effect=as.numeric(str_split(gene_es, ";")[[1]][[n]])) %>%
  ungroup() %>%
  dplyr::select(-gene_es, -n) %>%
  rename(gene=genes) %>%
  mutate(direction=ifelse(gene_effect < 0, "down", "up")) %>%
  dplyr::select(-gene_effect)

all_pair_overlap_g_sex_auto <- full_pair_sex_auto %>% filter(genes!="") %>% 
  dplyr::select(study1, study2, genes, gene_es) %>%
  as_tibble() %>%
  separate_rows(genes, sep=";") %>%
  group_by(study1, study2) %>%
  mutate(n=1:n()) %>%
  group_by(study1, study2, n) %>%
  mutate(gene_effect=as.numeric(str_split(gene_es, ";")[[1]][[n]])) %>%
  ungroup() %>%
  dplyr::select(-gene_es, -n) %>%
  rename(gene=genes) %>%
  mutate(direction=ifelse(gene_effect < 0, "down", "up")) %>%
  dplyr::select(-gene_effect)


load("ref/xy_chr_map.RData")
chr_l1 <- grab_chr(all_pair_overlap_g %>% distinct(gene))

chr_l <- grab_chr(all_pair_overlap_g_sex %>% distinct(gene))

all_pair_overlap_g_2 <- all_pair_overlap_g %>% 
  mutate(gene=ifelse(gene=="FAM46C", "TENT5C", gene)) %>%
  left_join(chr_l1 %>% filter(chromosome_name %in% c(1:22, "X", "Y")),
            by=c("gene"="hgnc_symbol")) %>%
  mutate(chromosome_name=case_when(
    !is.na(chromosome_name) ~ chromosome_name,
    gene=="FKSG49" ~ "5",
    gene=="MARC2" ~ "1", # MTARC2
    gene=="C5orf56" ~ "5" # IRF1-AS1
  ))

  
all_pair_overlap_g_sex2 <- all_pair_overlap_g_sex %>% 
  mutate(gene=ifelse(gene=="SEPT6", "SEPTIN6",gene )) %>%
  left_join(chr_l %>% filter(chromosome_name %in% c(1:22, "X", "Y")),
            by=c("gene"="hgnc_symbol"))  %>%
  mutate(chromosome_name=case_when(
    !is.na(chromosome_name) ~ chromosome_name,
    gene=="TTTY15" ~ "Y",
    gene=="LOC389906" ~ "X",
    gene=="HIST1H2AJ" ~ "6", # H2AC14
    gene=="KAL1" ~ "X", # ANOS1
    gene=="ARSE" ~ "X" #ARSL
  )) 



all_pair_overlap_g2 <- all_pair_overlap_g_2 %>% 
  dplyr::select(study1, study2, direction, gene, chromosome_name) %>%
  filter(study1!="Grouped AE", study2!="Grouped AE") %>%
  pivot_longer(c(study1, study2), values_to="study") %>%
  dplyr::select(-name) %>%
  left_join(list_studies %>% dplyr::select(study, tissue), by="study") %>%
  mutate(tissue=as.character(tissue)) %>%
  mutate(tissue=ifelse(is.na(tissue), "airway epithelium", tissue)) %>%
  distinct(chromosome_name, gene, direction, study, tissue) %>%
  group_by(chromosome_name, gene, direction) %>%
  dplyr::summarize(studies=paste(unique(study), collapse=";"),
            tissues=paste(unique(tissue), collapse=";"),
            n=n()) %>%
  arrange(desc(n))

all_pair_overlap_g2_sex <- all_pair_overlap_g_sex2 %>% 
  dplyr::select(study1, study2, direction, gene, chromosome_name) %>%
  filter(study1!="Grouped AE", study2!="Grouped AE") %>%
  pivot_longer(c(study1, study2), values_to="study") %>%
  dplyr::select(-name) %>%
  left_join(list_studies %>% dplyr::select(study, tissue), by="study") %>%
  mutate(tissue=as.character(tissue)) %>%
  mutate(tissue=ifelse(is.na(tissue), "airway epithelium", tissue)) %>%
  distinct(chromosome_name, gene, direction, study, tissue) %>%
  group_by(chromosome_name, gene, direction) %>%
  dplyr::summarize(studies=paste(unique(study), collapse=";"),
            tissues=paste(unique(tissue), collapse=";"),
            n=n()) %>%
  arrange(desc(n))


 
all_pair_overlap_g2 %>% filter(n>=3, chromosome_name %in% c("X", "Y"))  # 5 XY
all_pair_overlap_g2 %>% filter(n >=3, !chromosome_name %in% c("X", "Y"))  # 68 autosomal
all_pair_overlap_g2 %>% group_by(gene) %>% dplyr::select(-n) %>% count() %>% filter(n>1)
all_pair_overlap_g2 %>%
  write_csv("data/supp_tables/overlapping_genes2_smok.csv")  

all_pair_overlap_g2_sex %>% filter(n>=3, chromosome_name %in% c("X", "Y")) # 41 XY
all_pair_overlap_g2_sex %>% filter(n>=3, !chromosome_name %in% c("X", "Y")) # 4 autosomal
all_pair_overlap_g2_sex %>%
  write_csv("data/supp_tables/overlapping_genes2_sex.csv")  


# TODO:
# - add AE
# - reorder by tissue

self_cor <- data.table::rbindlist(lapply(c(names(all_studies2),"Grouped AE"), function(x)
  list("study1"=x, "study2"=x, "cor"=1)))

self_cor_sex <- data.table::rbindlist(lapply(c(setdiff(names(all_studies2),c("GSE16149", "GSE18723", "GSE87072")),"Grouped AE"), function(x)
  list("study1"=x, "study2"=x, "cor"=1)))

autosomal_sex_cts <- all_pair_overlap_g_sex2 %>%
  filter(!chromosome_name %in% c("X", "Y")) %>%
  group_by(study1, study2) %>%
  count()

sex_chr_cts <- all_pair_overlap_g_sex2 %>%
  filter(chromosome_name %in% c("X", "Y")) %>%
  group_by(study1, study2) %>%
  count()


full_pair2 <- full_pair %>%
  bind_rows(full_pair %>% 
              mutate(study3=study1, study1=study2, study2=study3) %>%
              dplyr::select(-study3) %>%
              dplyr::select(study1, study2, everything())) %>%
  mutate(cor=ifelse(tot < 30, NA, cor )) %>%
  bind_rows(self_cor) 

full_pair2_sex <- full_pair_sex %>%
  bind_rows(full_pair_sex %>% 
              mutate(study3=study1, study1=study2, study2=study3) %>%
              dplyr::select(-study3) %>%
              dplyr::select(study1, study2, everything())) %>%
  mutate(cor=ifelse(tot < 30, NA, cor )) %>%
  bind_rows(self_cor)  ### TYPO

full_pair2_sex_auto <- full_pair_sex_auto %>%
  bind_rows(full_pair_sex_auto %>% 
              mutate(study3=study1, study1=study2, study2=study3) %>%
              dplyr::select(-study3) %>%
              dplyr::select(study1, study2, everything())) %>%
  mutate(cor=ifelse(tot < 30, NA, cor )) %>%
  bind_rows(self_cor) 


counts2 <- full_pair2 %>% 
  left_join(list_studies %>% dplyr::select(study, tissue), 
            by=c("study1"="study")) %>%
  mutate(tissue=as.character(tissue)) %>%
  mutate(tissue=ifelse(study1=="Grouped AE", "airway epithelium", tissue)) %>%
  filter(study1!="GSE44456", study2!="GSE44456") %>%
  mutate(tissue=factor(tissue, levels=c("airway epithelium", "trachea epithelium",
                                        "nasal epithelium", "oral cavity",
                                        "buccal mucosa", "sputum", "alveolar macrophages", 
                                        "lung", 
                                        "blood - b cells", "blood - pbmcs", "blood - whole",
                                        "liver", "kidney", 
                                        "brain - prefrontal cortex"))) %>% 
  arrange(tissue,study1)
  
  
counts2$study1=factor(counts2$study1, levels=unique(counts2$study1))
counts2$study2=factor(counts2$study2, levels=unique(counts2$study1))



counts2_sex <- full_pair2_sex %>% 
  left_join(list_studies %>% dplyr::select(study, tissue), 
            by=c("study1"="study")) %>%
  mutate(tissue=as.character(tissue)) %>%
  mutate(tissue=ifelse(study1=="Grouped AE", "airway epithelium", tissue)) %>%
  filter(study1!="GSE44456", study2!="GSE44456") %>%
  mutate(tissue=factor(tissue, levels=c("airway epithelium", "trachea epithelium",
                                        "nasal epithelium", "oral cavity",
                                        "buccal mucosa", "sputum", "alveolar macrophages", 
                                        "lung", 
                                        "blood - b cells", "blood - pbmcs", "blood - whole",
                                        "liver", "kidney", 
                                        "brain - prefrontal cortex"))) %>% 
  arrange(tissue,study1)


counts2_sex$study1=factor(counts2_sex$study1, levels=unique(counts2_sex$study1))
counts2_sex$study2=factor(counts2_sex$study2, levels=unique(counts2_sex$study1))


counts2_sex_auto <- autosomal_sex_cts %>% 
  left_join(list_studies %>% dplyr::select(study, tissue), 
            by=c("study1"="study")) %>%
  mutate(tissue=as.character(tissue)) %>%
  mutate(tissue=ifelse(study1=="Grouped AE", "airway epithelium", tissue)) %>%
  filter(study1!="GSE44456", study2!="GSE44456") %>%
  mutate(tissue=factor(tissue, levels=c("airway epithelium", "trachea epithelium",
                                        "nasal epithelium", "oral cavity",
                                        "buccal mucosa", "sputum", "alveolar macrophages", 
                                        "lung", 
                                        "blood - b cells", "blood - pbmcs", "blood - whole",
                                        "liver", "kidney", 
                                        "brain - prefrontal cortex"))) %>% 
  arrange(tissue,study1)


counts2_sex_auto$study1=factor(counts2_sex_auto$study1, levels=unique(counts2_sex_auto$study1))
counts2_sex_auto$study2=factor(counts2_sex_auto$study2, levels=unique(counts2_sex_auto$study1))


#library('RColorBrewer')
#set3 <- brewer.pal(12, "Paired")
#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A"
# "#FFFF99" "#B15928"
library('LaCroixColoR')
paired2 <- lacroix_palette(type = "paired")
paired3 <- c(paired2[[1]], paired2[[2]], paired2[[11]], paired2[[12]], 
paired2[[5]], paired2[[6]], paired2[[3]], paired2[[4]],
paired2[[10]],paired2[[13]], paired2[[14]],
paired2[[7]], paired2[[8]], paired2[[9]])
ggplot(counts2)+
  geom_bar(mapping = aes(x = study1, y = 1, fill = tissue), 
           stat = "identity", 
           width = 1)+
  scale_fill_manual(values=paired3)
ggsave("figures/smok_bar.png")

counts2.2 <- counts2 %>% 
  filter(as.numeric(study1) < as.numeric(study2))
# make it only the top half

 ggplot(counts2.2,aes(x=study1,y=study2)) +
  geom_tile(aes(fill=rep_sig)) +
  geom_text(aes(label=rep_sig), size=3)+
  scale_fill_gradient2(low = "lightblue", high = "blue", limit = c(0,70), 
                       space = "Lab", 
                       name="Number of \noverlapping genes") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust=-0.1, vjust=0.5))+
   scale_x_discrete(position = "top") +
  coord_fixed()+
  ylab("")+
  xlab("")
ggsave("figures/paper_figs/gene_overlap2.png")

counts2_sex <- counts2_sex %>% 
  select(-tissue) %>%
  
  filter(study1!="GSE44456", study2!="GSE44456",
         !study1 %in% c("GSE87072", "GSE18723", "GSE16149"),
         !study2 %in% c("GSE87072", "GSE18723", "GSE16149")
  ) %>%
  inner_join(list_studies %>% dplyr::select(study, tissue), 
            by=c("study1"="study")) %>%
  mutate(tissue=as.character(tissue)) %>%
  mutate(tissue=ifelse(study1=="Grouped AE", "airway epithelium", tissue)) %>%
  mutate(tissue=factor(tissue, levels=c("airway epithelium", "trachea epithelium",
                                        "nasal epithelium", "oral cavity",
                                        "buccal mucosa", "sputum", "alveolar macrophages", 
                                        "lung", "blood - b cell",
                                         "blood - pbmcs", "blood - whole",
                                        "liver", "kidney", 
                                        "brain - prefrontal cortex"))) %>%
  arrange(tissue, study1)
counts2_sex$study1=factor(counts2_sex$study1, levels=unique(counts2_sex$study1))
counts2_sex$study2=factor(counts2_sex$study2, levels=unique(counts2_sex$study1))



counts2.2_sex <- counts2_sex %>% 
  filter(as.numeric(study1) < as.numeric(study2))
ggplot(counts2.2_sex)+
  geom_bar(mapping = aes(x = study1, y = 1, fill = tissue), 
           stat = "identity", 
           width = 1)+
  scale_fill_manual(values=paired3[c(1:8,10:13)])
ggsave("figures/sex_bar.png")


ggplot(counts2.2_sex,aes(x=study1,y=study2)) +
  geom_tile(aes(fill=rep_sig)) +
  geom_text(aes(label=rep_sig), size=3)+
  scale_fill_gradient2(low = "lightblue", high = "blue", limit = c(0,70), 
                       space = "Lab", 
                       name="Number of \noverlapping genes") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust=-0.1, vjust=0.5))+
  scale_x_discrete(position = "top") +  
  coord_fixed()+
  ylab("")+
  xlab("")
ggsave("figures/paper_figs/gene_overlap2_sex.png")
# ... how many are autosomal?

counts2.2_sex_auto <- counts2_sex_auto %>% 
  filter(as.numeric(study1) < as.numeric(study2)) 
ggplot(counts2.2_sex_auto,aes(x=study1,y=study2)) +
  geom_tile(aes(fill=n)) +
  geom_text(aes(label=n), size=3)+
  scale_fill_gradient2(low = "lightblue", high = "blue", limit = c(0,70), 
                       space = "Lab", 
                       name="Number of \noverlapping genes") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust=-0.1, vjust=0.5))+
  scale_x_discrete(position = "top") +  
  coord_fixed()+
  ylab("")+
  xlab("")
# "GSE20681", "E-MTAB-5278" overlap on 4
# GSE21862    GSE20681 on 3


#  filter(study1 < study2)
 # Create a ggheatmap

long_reord2 <- full_pair2 %>% 
  left_join(list_studies %>% dplyr::select(study, tissue), 
            by=c("study1"="study")) %>%
  mutate(tissue=as.character(tissue)) %>%
  mutate(tissue=ifelse(study1=="Grouped AE", "airway epithelium", tissue)) %>%
  filter(study1!="GSE44456", study2!="GSE44456") %>%
  mutate(tissue=factor(tissue, levels=c("airway epithelium", "trachea epithelium",
                                        "nasal epithelium", "oral cavity",
                                        "buccal mucosa", "sputum", "alveolar macrophages", 
                                        "lung", 
                                        "blood - b cells", "blood - pbmcs", "blood - whole",
                                        "liver", "kidney", 
                                        "brain - prefrontal cortex"))) %>%
  arrange(tissue, study1)
long_reord2$study1=factor(long_reord2$study1, levels=unique(long_reord2$study1))
long_reord2$study2=factor(long_reord2$study2, levels=unique(long_reord2$study1))
ggplot(long_reord2, aes(study1, study2, fill = cor))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", high = "orange", mid = "white", 
                       midpoint = 0, limit = c(-0.5,1), space = "Lab", 
                       name="Pearson\nCorrelation",
                       na.value="#DCDCDC") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  coord_fixed()+
  ylab("")+
  xlab("")
ggsave("figures/paper_figs/correlation_plot2.png")
# add airway epithelium


long_reord2_sex <- full_pair2_sex %>% 
  left_join(list_studies %>% dplyr::select(study, tissue), 
            by=c("study1"="study")) %>%
  mutate(tissue=as.character(tissue)) %>%
  mutate(tissue=ifelse(study1=="Grouped AE", "airway epithelium", tissue)) %>%
  filter(study1!="GSE44456", study2!="GSE44456",
         !study1 %in% c("GSE87072", "GSE18723", "GSE16149"),
         !study2 %in% c("GSE87072", "GSE18723", "GSE16149")
         ) %>%
  mutate(tissue=factor(tissue, levels=c("airway epithelium", "trachea epithelium",
                                        "nasal epithelium", "oral cavity",
                                        "buccal mucosa", "sputum", "alveolar macrophages", 
                                        "lung", 
                                        "blood - b cells", "blood - pbmcs", "blood - whole",
                                        "liver", "kidney", 
                                        "brain - prefrontal cortex"))) %>%
  arrange(tissue, study1)
long_reord2_sex$study1=factor(long_reord2_sex$study1, levels=unique(long_reord2_sex$study1))
long_reord2_sex$study2=factor(long_reord2_sex$study2, levels=unique(long_reord2_sex$study1))
ggplot(long_reord2_sex, aes(study1, study2, fill = cor))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", high = "orange", mid = "white", 
                       midpoint = 0, limit = c(-0.5,1), space = "Lab", 
                       name="Pearson\nCorrelation",
                       na.value="#DCDCDC") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  coord_fixed()+
  ylab("")+
  xlab("")
ggsave("figures/paper_figs/correlation_plot2_sex.png")



long_reord2_sex_a <- full_pair2_sex_auto %>% 
  left_join(list_studies %>% dplyr::select(study, tissue), 
            by=c("study1"="study")) %>%
  mutate(tissue=as.character(tissue)) %>%
  mutate(tissue=ifelse(study1=="Grouped AE", "airway epithelium", tissue)) %>%
  filter(study1!="GSE44456", study2!="GSE44456",
         !study1 %in% c("GSE87072", "GSE18723", "GSE16149"),
         !study2 %in% c("GSE87072", "GSE18723", "GSE16149")
  ) %>%
  mutate(tissue=factor(tissue, levels=c("airway epithelium", "trachea epithelium",
                                        "nasal epithelium", "oral cavity",
                                        "buccal mucosa", "sputum", "alveolar macrophages", 
                                        "lung", 
                                        "blood - b cells", "blood - pbmcs", "blood - whole",
                                        "liver", "kidney", 
                                        "brain - prefrontal cortex"))) %>%
  arrange(tissue, study1)
long_reord2_sex_a$study1=factor(long_reord2_sex_a$study1, levels=unique(long_reord2_sex_a$study1))
long_reord2_sex_a$study2=factor(long_reord2_sex_a$study2, levels=unique(long_reord2_sex_a$study1))
ggplot(long_reord2_sex_a, aes(study1, study2, fill = cor))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", high = "orange", mid = "white", 
                       midpoint = 0, limit = c(-0.5,1), space = "Lab", 
                       name="Pearson\nCorrelation",
                       na.value="#DCDCDC") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  coord_fixed()+
  ylab("")+
  xlab("")
ggsave("figures/paper_figs/correlation_plot2_sex_auto.png")

# plot fraction same direction? or n????

# meta-analyze?

prep_meta_study <- function(ds, study){
  ds %>%
    ungroup() %>%
    dplyr::select(-chromosome) %>%
    mutate(SD=(logFC-logFC.l)/1.96) %>%
    mutate(ID=study) %>%
    dplyr::select(gene, logFC, SD, ID)
}

run_meta <- function(df){
  df<- df %>% filter(!is.na(logFC))
  
  return(meta::metagen(df %>% pull(logFC), # treatment estimate
              df %>% pull(SD), # standard error
              studlab=df %>% pull(ID),
              comb.fixed = FALSE,
              comb.random = TRUE,
              method.tau = "DL", # method for between study variance
              hakn = FALSE,
              prediction = FALSE,
              sm = "MD") )
}
meta_study <- function(df){
  ma <- run_meta(df)
  # ma <- meta::metagen(df %>% pull(logFC), # treatment estimate
  #                           df %>% pull(SD), # standard error
  #                           studlab=df %>% pull(ID),
  #                           comb.fixed = FALSE,
  #                           comb.random = TRUE,
  #                           method.tau = "DL", # method for between study variance
  #                           hakn = FALSE,
  #                           prediction = FALSE,
  #                           sm = "MD") 
  return(list("gene"=unique(df$gene),
              "logFC"=ma$TE.random,
              "logFC.l"=ma$lower.random,
              "logFC.u"=ma$upper.random,
              "p"=ma$pval.random,
              "n"=unique(df$n)))
}
valid_studies <- c("GSE7895", "GSE27002", "GSE21862", "E-MTAB-5279")

meta_in <- data.table::rbindlist(lapply(all_studies2, function(x) 
  prep_meta_study(x$gene_smok, x$study)))

meta_in_disc <- meta_in %>% filter(!ID %in% valid_studies)
meta_in_valid <- meta_in %>% filter(ID %in% valid_studies)

meta_in2 <- meta_in_disc %>% 
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup() %>%
  arrange(gene)
meta_in3 <- meta_in2 %>% filter(n >=15) %>% group_split(gene)
names(meta_in3) <- unique(meta_in2 %>% filter(n>=15) %>% pull(gene)) # 14367

ma_vals <- lapply(meta_in3, meta_study)

# 5min
mult_study <- data.frame(apply(data.table::rbindlist(ma_vals) ,
                                 c(1,2), unlist)) %>%
  mutate(across(c(contains("logFC"), p), 
                ~as.numeric(as.character(.)) )) %>%
  mutate(n=as.numeric(as.character(n))) %>% 
  arrange(p)
mult_study$adj.p <- p.adjust(mult_study$p, method="fdr")
save(mult_study, file="data/results/meta_across_sm2.RData")

mult_study %>% filter(adj.p < 0.05) %>% nrow() # 48
sig_genes <- mult_study %>% filter(adj.p < 0.05 & abs(logFC)>=0.3) # 7
valid_g <- do.call(rbind, 
        lapply(valid_studies,
       function(x) all_studies2[[x]]$gene_smok %>% 
         filter(gene %in% as.character(sig_genes$gene)) %>%
         mutate(study=x) %>%
         ungroup() %>%
         dplyr::select(gene, logFC, p, adj.p, study)))
valid_g %>% filter(p < 0.05/nrow(valid_g))

#gene   logFC        p    adj.p study      
# NQO1   1.99  4.83e-20 6.00e-16 GSE7895    
# CYP1B1 2.61  7.71e- 5 8.47e- 2 GSE27002   
# LRRN3  1.58  1.88e- 6 1.16e- 2 E-MTAB-5279
# AHRR   0.878 7.84e- 4 1.26e- 1 E-MTAB-5279

valid_g %>% group_by(gene) %>% arrange(gene, p) %>%
  filter(p < 0.05) %>%
  dplyr::select(-adj.p) %>%
  write_csv("data/supp_tables/meta_smok_valid.csv")
sig_genes %>%
  write_csv("data/supp_tables/meta_smok.csv")


 

# sex
meta_in_sex <- data.table::rbindlist(lapply(all_studies2[! names(all_studies2) %in% c("GSE16149", "GSE18723", "GSE87072")], function(x) 
  prep_meta_study(x$gene_sex, x$study)))

meta_in_disc_s <- meta_in_sex %>% filter(!ID %in% valid_studies)
meta_in_valid_s <- meta_in_sex %>% filter(ID %in% valid_studies)

meta_in2_s <- meta_in_disc_s %>% 
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup() %>%
  arrange(gene)
summary(meta_in2_s$n)
meta_in3_s <- meta_in2_s %>% filter(n >=10) %>% group_split(gene)
names(meta_in3_s) <- unique(meta_in2_s %>% filter(n>=10) %>% pull(gene)) # 14367

ma_vals_s <- lapply(meta_in3_s, meta_study)

# 5min
mult_study_s <- data.frame(apply(data.table::rbindlist(ma_vals_s) ,
                               c(1,2), unlist)) %>%
  mutate(across(c(contains("logFC"), p), 
                ~as.numeric(as.character(.)) )) %>%
  mutate(n=as.numeric(as.character(n))) %>% 
  arrange(p)
mult_study_s$adj.p <- p.adjust(mult_study_s$p, method="fdr")
save(mult_study_s, file="data/results/meta_across_sm2_sex.RData")

mult_study_s %>% filter(adj.p < 0.05) %>% nrow() # 60
mult_study_s %>% filter(adj.p < 0.05 & abs(logFC)>=0.3) #22
mult_study_s %>% 
  filter(adj.p < 0.05) %>%
  left_join(xy_chr_map, by=c("gene"="hgnc_symbol")) %>%
  filter(!chromosome_name %in% c("X", "Y"))
sig_genes_sex <- mult_study_s %>% filter(adj.p < 0.05 & abs(logFC)>=0.3)





plot_gene(meta_in2_s, "DRAM2")


# forestPlot
plot_gene <- function(my_meta, my_gene){
  df <- my_meta %>% filter(gene==my_gene)
  ma <- run_meta(df)
  meta::forest.meta(ma)
}

plot_gene(meta_in2, "GZMA")

plot_gene(meta_in2, "LRRN3")


# 5. meta grouped by tissue
get_meta_out <- function(list_g){
  ma_vals <- lapply(list_g, meta_study)
  
  # 5min
  mult_study <- data.frame(apply(data.table::rbindlist(ma_vals) ,
                                 c(1,2), unlist)) %>%
    mutate(across(c(contains("logFC"), p), 
                  ~as.numeric(as.character(.)) )) %>%
    mutate(n=as.numeric(as.character(n))) %>% 
    arrange(p)
  mult_study$adj.p <- p.adjust(mult_study$p, method="fdr")
  return(mult_study)
}

# AE/BE
# TODO - remove validations studies
ae_studies <-list_studies %>% 
  filter(!study %in% valid_studies) %>%
  filter(str_detect(tissue,"airway") | str_detect(tissue, "trachea")) %>%
  pull(study)
ae_meta <- meta_in %>% filter(ID %in% ae_studies) %>%
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup() %>%
  filter(n>=4) %>%
  group_split(gene)
ae_meta_out <- get_meta_out(ae_meta)
ae_sig <- ae_meta_out %>% filter(abs(logFC) >= 0.3 & adj.p < 0.05)
ae_sig # 18

# blood or blood component
blood_studies <-list_studies %>% 
  filter(!study %in% valid_studies) %>%
  filter(str_detect(tissue, "blood"), !str_detect(tissue, "b cell") ) %>%
  pull(study) # PBMC, whole blood
blood_meta <- meta_in %>% 
  filter(ID %in% blood_studies) %>%
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup() %>%
  filter(n>=5) %>%
  group_split(gene)
blood_meta_out <- get_meta_out(blood_meta)
blood_meta_out %>% filter(adj.p < 0.05) # 15
blood_sig <- blood_meta_out %>% filter(abs(logFC) >= 0.3 & adj.p < 0.05)
blood_sig
plot_gene(meta_in %>% filter(ID %in% blood_studies), "FAM81B")


# compare - 314 genes
meta2 <- ae_meta_out %>% 
  left_join(blood_meta_out, by="gene") %>%
  filter(adj.p.x < 0.1 | adj.p.y < 0.1) 
meta2 %>% 
  filter(logFC.x*logFC.y > 0 & p.x < 0.05 & 
           p.y < 0.05) %>%
  dplyr::select(gene, logFC.x, logFC.y, p.x, p.y) # 22
meta2 %>% filter(logFC.x*logFC.y < 0) %>% nrow()

# meta-analyze by sex
ae_meta_sex <- meta_in_sex %>% 
  filter(ID %in% ae_studies) %>%
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup() %>%
  filter(n>=4) %>%
  group_split(gene)
ae_meta_out_s <- get_meta_out(ae_meta_sex)
ae_sig_sex <- ae_meta_out_s %>% filter(abs(logFC) >= 0.3 & adj.p < 0.05)
ae_sig_sex

blood_meta_sex <- meta_in_sex %>% 
  filter(ID %in% blood_studies) %>%
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup() %>%
  filter(n>=5) %>%
  group_split(gene)
blood_meta_out_s <- get_meta_out(blood_meta_sex)
blood_meta_out_s %>% filter(adj.p < 0.05) # 15
blood_sig_sex <- blood_meta_out_s %>% filter(abs(logFC) >= 0.3 & adj.p < 0.05)

# save
save(ae_meta_out, blood_meta_out, ae_meta_out_s, blood_meta_out_s, file="data/tissue_spec_ma.RData")


# add chromosomes to all the sig genes
sig_df <- do.call(rbind, 
                  list(sig_genes %>% mutate(term="smoking", ds="all"),
  sig_genes_sex  %>% mutate(term="sex", ds="all"),
  blood_sig  %>% mutate(term="smoking", ds="blood"),
  blood_sig_sex  %>% mutate(term="sex", ds="blood"),
  ae_sig %>% mutate(term="smok", ds="ae"),
  ae_sig_sex %>% mutate(term="sex", ds="ae")))
chr_map <- grab_chr(sig_df) 
sig_df2 <- sig_df %>% 
  left_join(chr_map %>% filter(chromosome_name %in% c(1:22, "X", "Y")), by=c("gene"="hgnc_symbol")) %>%
  mutate(chromosome_name=case_when(
    !is.na(chromosome_name) ~ chromosome_name,
    gene == "TTTY15" ~ "Y",
    gene=="MARCH2" ~ "19"
  ))

# - write out output

sig_meta_comb <- sig_df2 %>% 
  as_tibble() %>%
  dplyr::select(ds, term, everything()) %>%
  mutate(across(c(contains("logFC"), contains("p")), ~signif(., 3)))
sig_meta_comb %>% write_csv("data/supp_tables/sig_meta_comb.csv")


ae_de_gene_s %>% 
  filter(adj.p < 0.05 & abs(logFC) >= 0.3) %>% 
  left_join(chr_ae %>% filter(chromosome_name %in% c(1:22, "X", "Y")),
            by=c("gene"="hgnc_symbol")) %>% 
  write_csv("data/supp_tables/ae_genes_sex.csv")

# look at validation

valid_g <- do.call(rbind, 
                   lapply(valid_studies,
                          function(x) all_studies2[[x]]$gene_smok %>% 
                            filter(gene %in% as.character(sig_genes$gene)) %>%
                            mutate(study=x) %>%
                            ungroup() %>%
                            dplyr::select(gene, logFC, p, adj.p, study)))
valid_g2 <- valid_g %>% 
  group_by(gene) %>% 
  arrange(gene, p) %>%
  filter(p < 0.05/nrow(valid_g)) %>%
  dplyr::select(-adj.p)
sig_genes %>% inner_join(valid_g2, by="gene"  ) %>%
  dplyr::select(gene, logFC.x, logFC.y) %>%
  filter(logFC.x*logFC.y > 0)
valid_g2 %>% write_csv("data/supp_tables/meta_smok_valid.csv")

valid_g_sex <- do.call(rbind, 
                   lapply(valid_studies,
                          function(x) all_studies2[[x]]$gene_sex %>% 
                            filter(gene %in% as.character(sig_genes_sex$gene)) %>%
                            mutate(study=x) %>%
                            ungroup() %>%
                            dplyr::select(gene, logFC, p, adj.p, study)))
valid_g_sex2 <- valid_g_sex %>% 
  filter(p < 0.05/nrow(valid_g_sex)) %>%
  arrange(gene, p) %>%
  dplyr::select(-adj.p)

valid_g_sex2 %>% 
  left_join(sig_genes_sex, by="gene") %>%
  filter(logFC.x*logFC.y < 0) # none

ae_valid <- all_studies2[["GSE7895"]]$gene_smok %>% 
  filter(gene %in% as.character(ae_sig$gene)) %>%
  ungroup() %>%
  dplyr::select(gene, logFC, p) 
ae_valid2 <- ae_valid %>% 
  filter(p < 0.05/nrow(ae_valid)) %>% 
  arrange(gene, p) %>%
  mutate(study="GSE7895")# 21
ae_valid2 %>% 
  left_join(ae_sig, by="gene") %>%
  filter(logFC.x*logFC.y < 0) # none

ae_valid_sex <- all_studies2[["GSE7895"]]$gene_sex %>% 
  filter(gene %in% as.character(ae_sig_sex$gene)) %>%
  ungroup() %>%
  dplyr::select(gene, logFC, p)
ae_valid_sex2 <- ae_valid_sex %>%
  filter(p < 0.05/nrow(ae_valid_sex)) %>% 
  arrange(gene, p) %>%
  mutate(study="GSE7895")

ae_valid_sex2 %>%
  left_join(ae_sig_sex, by="gene") %>%
  filter(logFC.x*logFC.y < 0) # none

blood_valid <- do.call(rbind, 
                       lapply(c("GSE21862", "E-MTAB-5279"),
                              function(x) all_studies2[[x]]$gene_smok %>% 
                                filter(gene %in% as.character(blood_sig$gene)) %>%
                                mutate(study=x) %>%
                                ungroup() %>%
                                dplyr::select(gene, logFC, p, study)))
blood_valid2 <- blood_valid %>% 
  filter(p < 0.05/nrow(blood_valid)) %>%
  arrange(gene, p) # 3

blood_valid2 %>%
  left_join(blood_sig, by="gene") %>%
  filter(logFC.x*logFC.y < 0) # none


blood_valid_sex <- do.call(rbind, 
                           lapply(c("GSE21862", "E-MTAB-5279"),
                                  function(x) all_studies2[[x]]$gene_sex %>% 
                                    filter(gene %in% as.character(blood_sig_sex$gene)) %>%
                                    mutate(study=x) %>%
                                    ungroup() %>%
                                    dplyr::select(gene, logFC, p, study)))
blood_valid_sex2 <- blood_valid_sex %>% 
  filter(p < 0.05/nrow(blood_valid_sex)) %>%
  arrange(gene, p) 

blood_valid_sex2 %>%
  left_join(blood_sig_sex, by="gene") %>%
  filter(logFC.x*logFC.y < 0) # none

all_smok <- sig_genes %>% 
  left_join(valid_g2 %>% group_by(gene) %>%
              dplyr::summarize(study=paste(study, collapse=";"),
                               p.valid=min(p),
                               logFC.valid=max(abs(logFC))*(logFC/(abs(logFC)))),
            by="gene")

all_sex <- sig_genes_sex %>% 
  left_join(valid_g_sex2 %>% 
              group_by(gene) %>%
              dplyr::summarize(study=paste(study, collapse=";"),
                               p.valid=min(p),
                               logFC.valid=max(abs(logFC))*(max(logFC)/(abs(max(logFC))))),
            by="gene")



blood_smok <- blood_sig %>% 
  left_join(blood_valid2 %>% group_by(gene) %>%
              dplyr::summarize(study=paste(study, collapse=";"),
                               p.valid=min(p),
                               logFC.valid=max(abs(logFC))*(logFC/(abs(logFC)))),
            by="gene")

blood_sex <- blood_sig_sex %>% 
  left_join(blood_valid_sex2 %>% 
              group_by(gene) %>%
              dplyr::summarize(study=paste(study, collapse=";"),
                               p.valid=min(p),
                               logFC.valid=max(abs(logFC))*(max(logFC)/(abs(max(logFC))))),
            by="gene")

ae_smok <- ae_sig %>% 
  left_join(ae_valid2 %>% group_by(gene) %>%
              dplyr::summarize(study=paste(study, collapse=";"),
                               p.valid=min(p),
                               logFC.valid=max(abs(logFC))*(logFC/(abs(logFC)))),
            by="gene")

ae_sex <- ae_sig_sex %>% 
  left_join(ae_valid_sex2 %>% 
              group_by(gene) %>%
              dplyr::summarize(study=paste(study, collapse=";"),
                               p.valid=min(p),
                               logFC.valid=max(abs(logFC))*(max(logFC)/(abs(max(logFC))))),
            by="gene")



sig_df_v <- do.call(rbind, 
                  list(all_smok %>% mutate(term="smoking", ds="all"),
                       all_sex  %>% mutate(term="sex", ds="all"),
                       blood_smok  %>% mutate(term="smoking", ds="blood"),
                       blood_sex  %>% mutate(term="sex", ds="blood"),
                       ae_smok %>% mutate(term="smok", ds="ae"),
                       ae_sex %>% mutate(term="sex", ds="ae")))
sig_df2_v <- sig_df_v %>% 
  left_join(chr_map %>% filter(chromosome_name %in% c(1:22, "X", "Y")), by=c("gene"="hgnc_symbol")) %>%
  mutate(chromosome_name=case_when(
    !is.na(chromosome_name) ~ chromosome_name,
    gene == "TTTY15" ~ "Y",
    gene=="MARCH2" ~ "19"
  ))

# - write out output

sig_meta_comb_v <- sig_df2_v %>% 
  as_tibble() %>%
  rename(valid_study=study) %>%
  dplyr::select(ds, term, gene, chromosome_name, everything()) %>%
  mutate(across(c(contains("logFC"), contains("p")), ~signif(., 3)))
sig_meta_comb_v %>% write_csv("data/supp_tables/sig_meta_comb_valid.csv")
sig_meta_comb_v %>% 
  group_by(ds, term) %>%
  dplyr::summarize(n=n(), n_valid=sum(!is.na(logFC.valid)))
# ds    term        n n_valid
# <chr> <chr>   <int>   <int>
#   1 ae    sex         6       4
# 2 ae    smok       66      21
# 3 all   sex        22      20
# 4 all   smoking     7       4
# 5 blood sex        32      13
# 6 blood smoking    19       3

bs_auto <- sig_meta_comb_v %>% 
  filter(ds=="blood", term=="sex") %>%
  mutate(auto=(!chromosome_name %in% c("X", "Y"))) %>%
  dplyr::select(-term, -ds, -logFC.l, -logFC.u, -chromosome_name, -n, -p)  
table(bs_auto$auto)
# overlap?
intersect(blood_sig$gene, sig_genes$gene) # GZMH
intersect(ae_sig$gene, sig_genes$gene) # NQO1
intersect(blood_sig$gene, ae_sig$gene) # AKR1C3 -- opp dir
intersect(ae_sig$gene, sig_genes$gene) # NQO1
intersect(blood_sig_sex$gene, ae_sig_sex$gene) # 4



# heatmap for gene data  #

blood_plot_dat <- meta_in %>% filter(gene %in% blood_sig$gene, 
                                     ID %in% blood_studies ) %>%
  dplyr::select(-SD) %>%
  bind_rows(blood_sig %>% 
              dplyr::select(gene, logFC) %>% 
              mutate(ID="pooled")) 

blood_plot_dat_w <- blood_plot_dat %>% 
  pivot_wider(names_from="ID", values_from="logFC") %>%
  filter(!is.na(pooled)) 
blood_plot_dat2 <- blood_plot_dat_w %>% 
  pivot_longer(-gene, names_to="ID", values_to="logFC") %>%
  mutate(gene=factor(gene, levels=blood_plot_dat_w$gene[order(blood_plot_dat_w$pooled)]))


ggplot(blood_plot_dat2, aes(y=ID, x=gene, fill=logFC))+
  geom_tile()+
  scale_fill_gradient2(low = "turquoise", high = "gold4", mid = "white", 
                       midpoint = 0,  na.value="#DCDCDC", 
                       limit=c(-1.75, 0.75))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))+
  xlab("")+
  ylab("")
ggsave("figures/paper_figs/blood_meta_heatmap2.png")


ae_plot_dat <- meta_in %>% 
  filter(gene %in% ae_sig$gene, ID %in% ae_studies ) %>%
  dplyr::select(-SD) %>%
  bind_rows(ae_sig %>% 
              dplyr::select(gene, logFC) %>% 
              mutate(ID="pooled"))
ae_plot_dat_w <- ae_plot_dat %>% 
  pivot_wider(names_from="ID", values_from="logFC") %>%
  filter(!is.na(pooled)) 
ae_plot_dat2 <- ae_plot_dat_w %>% 
  pivot_longer(-gene, names_to="ID", values_to="logFC") %>%
  mutate(gene=factor(gene, levels=ae_plot_dat_w$gene[order(ae_plot_dat_w$pooled)]))

ggplot(ae_plot_dat2, aes(y=ID, x=gene, fill=logFC))+
  geom_tile()+
  scale_fill_gradient2(low = "turquoise", high = "gold4", mid = "white", 
                       midpoint = 0,  na.value="#DCDCDC")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))+
  xlab("")+
  ylab("")
ggsave("figures/paper_figs/ae_meta_heatmap2.png")


# by sex
blood_plot_dat_s <- meta_in_sex %>% filter(gene %in% blood_sig_sex$gene, 
                                     ID %in% blood_studies ) %>%
  dplyr::select(-SD) %>%
  bind_rows(blood_sig_sex %>% 
              dplyr::select(gene, logFC) %>% 
              mutate(ID="pooled")) 

blood_plot_dat_s_w <- blood_plot_dat_s %>% 
  pivot_wider(names_from="ID", values_from="logFC") %>%
  filter(!is.na(pooled)) 
blood_plot_dat_s_2 <- blood_plot_dat_s_w %>% 
  pivot_longer(-gene, names_to="ID", values_to="logFC") %>%
  mutate(gene=factor(gene, levels=blood_plot_dat_s_w$gene[order(blood_plot_dat_s_w$pooled)]))


ggplot(blood_plot_dat_s_2 %>% 
         mutate(logFC=ifelse(logFC < -6, logFC+4, logFC)), aes(y=ID, x=gene, fill=logFC))+
  geom_tile()+
  scale_fill_gradient2(low = "turquoise", high = "gold4", mid = "white", 
                       midpoint = 0,  na.value="#DCDCDC")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  xlab("")+
  ylab("")
ggsave("figures/paper_figs/blood_meta_heatmap2_sex.png")


ae_plot_dat_s <- meta_in_sex %>% 
  filter(gene %in% ae_sig_sex$gene, ID %in% ae_studies ) %>%
  dplyr::select(-SD) %>%
  bind_rows(ae_sig_sex %>% 
              dplyr::select(gene, logFC) %>% 
              mutate(ID="pooled"))
ae_plot_dat_w_s <- ae_plot_dat_s %>% 
  pivot_wider(names_from="ID", values_from="logFC") %>%
  filter(!is.na(pooled)) 
ae_plot_dat2_s <- ae_plot_dat_w_s %>% 
  pivot_longer(-gene, names_to="ID", values_to="logFC") %>%
  mutate(gene=factor(gene, levels=ae_plot_dat_w_s$gene[order(ae_plot_dat_w_s$pooled)]))

ggplot(ae_plot_dat2_s, aes(y=ID, x=gene, fill=logFC))+
  geom_tile()+
  scale_fill_gradient2(low = "turquoise", high = "gold4", mid = "white", 
                       midpoint = 0,  na.value="#DCDCDC")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  xlab("")+
  ylab("")
ggsave("figures/paper_figs/ae_meta_heatmap2_sex.png")


plot_dat_s <- meta_in %>% 
  filter(gene %in% sig_genes$gene ) %>%
  dplyr::select(-SD) %>%
  bind_rows(sig_genes %>% 
              dplyr::select(gene, logFC) %>% 
              mutate(ID="pooled"))
plot_dat_w_s <- plot_dat_s %>% 
  pivot_wider(names_from="ID", values_from="logFC") %>%
  filter(!is.na(pooled)) 
plot_dat_w_s2 <- plot_dat_w_s %>% 
  pivot_longer(-gene, names_to="ID", values_to="logFC") %>%
  mutate(gene=factor(gene, 
                     levels=plot_dat_w_s$gene[order(plot_dat_w_s$pooled)])) 
plot_dat <- plot_dat_w_s2 %>% 
  rename(study=ID) %>%
  left_join(list_studies %>% 
            dplyr::select(study, tissue, samples) ) %>%
  mutate(study=as.character(study), tissue=as.character(tissue)) %>%
  filter(study!="GSE44456") %>%
  mutate(tissue=ifelse(study=="pooled", "all", tissue)) %>%
  mutate(tissue=factor(tissue, levels=c("all", "airway epithelium", "trachea epithelium",
                                        "nasal epithelium", "oral cavity",
                                        "buccal mucosa", "sputum", "alveolar macrophages", 
                                        "lung", 
                                        "blood - b cells", "blood - pbmcs", "blood - whole",
                                        "liver", "kidney", 
                                        "brain - prefrontal cortex"))) %>%
  arrange(tissue, study)
plot_dat$study <- factor(plot_dat$study, levels=unique(plot_dat$study)  )

ggplot(plot_dat, aes(y=study, x=gene, fill=logFC))+
  geom_tile()+
  scale_fill_gradient2(low = "turquoise", high = "gold4", mid = "white", 
                       midpoint = 0,  na.value="#DCDCDC")+
  theme_bw()+
  #theme(axis.text.x = element_text(angle = 90, hjust=1))+
  xlab("")+
  ylab("")
ggsave("figures/paper_figs/heatmap_all.png")
color_bar(plot_dat)
ggsave("figures/paper_figs/heatmap_all_colorbar.png")


plot_dat_s <- meta_in_sex %>% 
  filter(gene %in% sig_genes_sex$gene ) %>%
  dplyr::select(-SD) %>%
  bind_rows(sig_genes_sex %>% 
              dplyr::select(gene, logFC) %>% 
              mutate(ID="pooled"))
plot_dat_w_s <- plot_dat_s %>% 
  pivot_wider(names_from="ID", values_from="logFC") %>%
  filter(!is.na(pooled)) 
plot_dat_w_s2 <- plot_dat_w_s %>% 
  pivot_longer(-gene, names_to="ID", values_to="logFC") %>%
  mutate(gene=factor(gene, 
                     levels=plot_dat_w_s$gene[order(plot_dat_w_s$pooled)])) 
plot_dat <- plot_dat_w_s2 %>% 
  rename(study=ID) %>%
  left_join(list_studies %>% 
              dplyr::select(study, tissue, samples) ) %>%
  mutate(study=as.character(study), tissue=as.character(tissue)) %>%
  filter(study!="GSE44456") %>%
  mutate(tissue=ifelse(study=="pooled", "all", tissue)) %>%
  mutate(tissue=factor(tissue, levels=c("all", "airway epithelium", "trachea epithelium",
                                        "nasal epithelium", "oral cavity",
                                        "buccal mucosa", "sputum", "alveolar macrophages", 
                                        "lung", 
                                        "blood - b cells", "blood - pbmcs", "blood - whole",
                                        "liver", "kidney", 
                                        "brain - prefrontal cortex"))) %>%
  arrange(tissue, study)
plot_dat$study <- factor(plot_dat$study, levels=unique(plot_dat$study)  )

ggplot(plot_dat, aes(y=study, x=gene, fill=logFC))+
  geom_tile()+
  scale_fill_gradient2(low = "turquoise", high = "gold4", mid = "white", 
                       midpoint = 0,  na.value="#DCDCDC")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust=1))+
  xlab("")+
  ylab("")
  
ggsave("figures/paper_figs/heatmap_all_sex.png")
color_bar(plot_dat)
ggsave("figures/paper_figs/heatmap_all_sex_colorbar.png")


### FOREST PLOTS
my_gene <- "AHRR"
plot_gene2 <- function(my_gene){
gene_df <- meta_in %>% filter(gene %in% my_gene ) %>%
  dplyr::rename(study=ID) %>%
  bind_rows(mult_study %>% 
              mutate(SD=(logFC-logFC.l)/1.96) %>%
              dplyr::select(gene, logFC, SD) %>% 
              filter(gene %in% my_gene) %>%
              mutate(study="pooled"))  %>%
  as_tibble() %>% left_join(list_studies %>% 
                              dplyr::select(study, tissue, samples) ) %>%
  mutate(study=as.character(study), tissue=as.character(tissue)) %>%
  filter(study!="GSE44456") %>%
  mutate(tissue=ifelse(study=="pooled", "all", tissue)) %>%
  mutate(tissue=factor(tissue, levels=c("all", "airway epithelium", "trachea epithelium",
                                        "nasal epithelium", "oral cavity",
                                        "buccal mucosa", "sputum", "alveolar macrophages", 
                                        "lung", 
                                        "blood - b cells", "blood - pbmcs", "blood - whole",
                                        "liver", "kidney", 
                                        "brain - prefrontal cortex"))) %>%
  arrange(tissue, study)
gene_df$study <- factor(gene_df$study, levels=unique(gene_df$study)  )
ggplot(gene_df, aes(x=logFC, y=study))+
  geom_vline(xintercept=0, lty=2, col="darkgray") +
  geom_point(aes(size=samples), alpha=0.6)+
  geom_errorbarh(aes(xmin=logFC-1.96*SD, xmax=logFC+1.96*SD))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab("")+
  facet_grid(.~gene, scales="free")
}
plot_gene2(sig_genes$gene)
ggsave("figures/all_forest.png")

ggplot(gene_df)+
  geom_bar(mapping = aes(x = study, y = 1, fill = tissue), 
           stat = "identity", 
           width = 1)+
  scale_fill_manual(values=c("black", paired3))
ggsave("figures/smok_bar.png")

plot_gene2("AHRR")
ggsave("figures/ahrr_forest_plot2.png")  

plot_gene2("LRRN3")
ggsave("figures/lrrn3_forest_plot2.png")  

plot_gene2("CYP1B1")
ggsave("figures/cy1b1_forest_plot2.png")  

plot_gene2("NQO1")
ggsave("figures/nqo1_forest_plot2.png") 

plot_gene2("CCL4")
ggsave("figures/nccl4_forest_plot2.png")

plot_gene2("GZMH")
ggsave("figures/gzmh_forest_plot2.png") 

#facet_grid(tissue~., scales="free")
plot_gene(meta_in, "AHRR")


# ....... #
# [x] 0. fix IDs

fu_int <- data_df4 %>% 
  #filter(de_probes_int > 0 | de_genes_int > 0 ) %>%
  dplyr::select(study, tissue, contains("int")) %>%
  filter(!study %in% c("GSE4302", "GSE19027", "GSE8987", "GSE42743", "GSE16149",
                       "GSE2125", "GSE103174", "GSE56768", "GSE21862", "GSE87072", "GSE18723"))
# remove studies that are too small


# 1. results for sex-related studies

# // TODO cleanup and write these out again?
get_rep <- function(ds1, study1, ds2, study2){
  # look at gene overlap
  ds_smok <- ds1 %>% 
    ungroup() %>%
    dplyr::select(-chromosome) %>%
    mutate(id=1:nrow(ds1)) %>% 
    separate_rows(gene, sep=",") 
  ds_de_gene <- ds_smok %>%
    filter(abs(logFC) >= 0.3 & adj.p < 0.05 ) %>%
    pull(gene)
  
  overlap <- ds_smok %>%
    filter(gene %in% ds_de_gene) %>%
    inner_join(ds2, by="gene") # 702
  
  ggplot(overlap, aes(x=logFC.x, y=logFC.y))+
    geom_point(alpha=0.5)+
    theme_bw()
  if (nrow(overlap) < 10){
    cor <- NA
    corp <- NA
  } else {
    pcor <- cor.test(overlap$logFC.x, overlap$logFC.y, method="pearson") # 0.00042, p=0.99
    cor <- pcor$estimate[["cor"]]
    corp <- pcor$p.value
  }
  
  same_dir <- overlap %>% filter(logFC.x*logFC.y > 0) # 345/702
  rep_sig <- same_dir %>% filter(p.y < 0.05)#/nrow(overlap))
  
  return(list("study1"=study1, "study2"=study2, "tot"=nrow(overlap), "cor"=cor, "cor.p"=corp,
              "same_dir"=nrow(same_dir), "rep_sig"=nrow(rep_sig), 
              "genes"=paste(rep_sig$gene, collapse=";"),
              "gene_es"=paste(round(rep_sig$logFC.y,3), collapse=";")))
}

res_int <- lapply(all_studies2[fu_int %>% pull(study)], 
                  function(x) get_rep( gene_int_ae %>% mutate(chromosome=""),
                                       "AE", x$gene_int, x$study))
int_overlap <- data.frame(do.call(rbind, res_int))

chatz_de %>% filter(gene=="BFSP1")
int_overlap %>% 
  filter(study2 %in% 
           c("E-MTAB-3604", "GSE17913", "E-MTAB-5278", "GSE30272")) 
head(ae_de_gene_i)


replicated_int <- ae_de_gene_i %>% 
  dplyr::select(-src, -logFC.l, -logFC.u) %>%
  inner_join( all_dfs_int2[["E-MTAB-5278"]] %>%
               mutate(study="E-MTAB-5278") %>% 
  bind_rows(all_dfs_int2[["GSE30272"]] %>% 
              mutate(study="GSE30272")) %>%
  dplyr::select(study, gene, logFC, p), by="gene") %>%
  filter(p.y < 0.05 & logFC.x*logFC.y > 0) %>%
  rename(valid_study=study, logFC.valid=logFC.y, p.valid=p.y,
         logFC.disc=logFC.x, adj.p.disc=adj.p) %>%
  dplyr::select(valid_study, gene, logFC.disc, logFC.valid,
                adj.p.disc, p.valid)
replicated_int %>%
  mutate(across(contains("\\."), ~signif(., 3))) %>%
  write_csv("data/supp_tables/replicated.csv")

g1 <- get_gene("E-MTAB-5278", "SLC25A37")

g1_probes <- probe_gene %>% 
  filter(gene=="SLC25A37") %>%
  pull(probes)
g1_e <- colMeans(expDat5[g1_probes,] )
g1_e2 <- tibble(sample_acc=names(g1_e), expr=unlist(g1_e)) %>%
  left_join(pDat5.1 %>% dplyr::select(geo_accession, sex, smoking) %>%
              rename(sample_acc=geo_accession)) %>%
  mutate(study_acc="Grouped AE")
g1_2 <- g1 %>% rename(smoking=smok)

g1_2 %>% 
  bind_rows(g1_e2 %>% dplyr::select(colnames(g1_2))) %>%
  filter(!is.na(study_acc), !is.na(smoking), !is.na(sex)) %>%
  mutate(smoking=ifelse(smoking=="S", "smokers", "non-smokers"),
         sex=sprintf("%ss", sex) ,
         study_acc=factor(study_acc, c("Grouped AE", "E-MTAB-5278"))) %>%
  ggplot(aes(x=sex, col=smoking, y=expr))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(alpha=0.5, position=position_jitterdodge(jitter.width = 0.2))+
  facet_wrap(.~study_acc, scales="free")+
  theme_bw()+
  xlab("")+
  scale_color_manual(values=c("turquoise", "gold4")) +
  ylab("log expression of SLC25A37")
ggsave("figures/paper_figs/SLC25A37.png")



g1 <- get_gene("E-MTAB-5278", "OPN3")
g1_probes <- probe_gene %>% 
  filter(gene=="OPN3") %>%
  pull(probes)
g1_e <- colMeans(expDat5[g1_probes,] )
g1_e2 <- tibble(sample_acc=names(g1_e), expr=unlist(g1_e)) %>%
  left_join(pDat5.1 %>% dplyr::select(geo_accession, sex, smoking) %>%
              rename(sample_acc=geo_accession)) %>%
  mutate(study_acc="Grouped AE")
g1_2 <- g1 %>% rename(smoking=smok)

g1_2 %>% 
  bind_rows(g1_e2 %>% dplyr::select(colnames(g1_2))) %>%
  filter(!is.na(study_acc), !is.na(smoking), !is.na(sex)) %>%
  mutate(smoking=ifelse(smoking=="S", "smokers", "non-smokers"),
         sex=sprintf("%ss", sex) ,
         study_acc=factor(study_acc, c("Grouped AE", "E-MTAB-5278"))) %>%
  ggplot(aes(x=sex, col=smoking, y=expr))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(alpha=0.5, position=position_jitterdodge(jitter.width = 0.2))+
  facet_wrap(.~study_acc, scales="free")+
  theme_bw()+
  xlab("")+
  scale_color_manual(values=c("turquoise", "gold4")) +
  ylab("log expression of OPN3")
ggsave("figures/paper_figs/OPN3.png")


all_dfs_int2[["GSE30272"]] %>%
  filter(gene %in% c("RALGDS", "KCNJ1", "MS4A7", "MUC5B", "LINC00342"))


plot_gene_int <- function(study, my_gene, orig=F){
  g1 <- get_gene(study, my_gene, orig)
  g1_probes <- probe_gene %>% 
    filter(gene==my_gene) %>%
    pull(probes)
  g1_e <- colMeans(expDat5[g1_probes,] )
  g1_e2 <- tibble(sample_acc=names(g1_e), expr=unlist(g1_e)) %>%
    left_join(pDat5.1 %>% dplyr::select(geo_accession, sex, smoking) %>%
                rename(sample_acc=geo_accession)) %>%
    mutate(study_acc="Grouped AE")
  g1_2 <- g1 %>% rename(smoking=smok)
  
  g1_2 %>% 
    bind_rows(g1_e2 %>% dplyr::select(colnames(g1_2))) %>%
    filter(!is.na(study_acc), !is.na(smoking), !is.na(sex), 
           sex!="unlabeled", smoking %in% c("NS", "S")) %>%
    mutate(smoking=ifelse(smoking=="S", "smokers", "non-smokers"),
           sex=sprintf("%ss", sex) ,
           study_acc=factor(study_acc, c("Grouped AE", study))) %>%
    ggplot(aes(x=sex, col=smoking, y=expr))+
    geom_boxplot(outlier.shape=NA)+
    geom_point(alpha=0.5, position=position_jitterdodge(jitter.width = 0.2))+
    facet_wrap(.~study_acc, scales="free")+
    theme_bw()+
    xlab("")+
    scale_color_manual(values=c("turquoise", "gold4")) +
    ylab(sprintf("log expression of %s", my_gene))
  
}
# SLC25A37 decreases in males 
# OPN3 - decreases in females but not males

plot_gene_int("GSE30272", "LINC00342") # nope

ggsave("figures/paper_figs/LINC00342.png")


plot_gene_int("GSE30272", "RALGDS") # increases in females, decreases in males
ggsave("figures/paper_figs/RALGDS.png")

plot_gene_int("GSE30272", "KCNJ1") # decreases more in females
ggsave("figures/paper_figs/KCNJ1.png")

plot_gene_int("GSE30272", "MS4A7") # decreases in females not males
ggsave("figures/paper_figs/MS4A7.png")

plot_gene_int("GSE30272", "MUC5B") # nope
ggsave("figures/paper_figs/MUC5B.png")

plot_gene_int("GSE7895", "CAPN9", orig=T)
ggsave("figures/paper_figs/capn9.png")

plot_gene_int("GSE7895", "LLGL2", orig=T)
ggsave("figures/paper_figs/llgl2.png")

int_overlap_df <- apply(int_overlap , c(1,2), unlist) %>%
  as_tibble() %>%
  mutate(across(contains("cor"), ~signif(as.numeric(.), 3))) %>%
  dplyr::select(-study1) 

all_dfs_int2[["GSE20681"]] %>% filter(gene=="MOCOS")
all_dfs_int2[["Grouped AE"]] %>% filter(gene=="MOCOS")
int_overlap_df

int_overlap_df %>%
  write_csv("data/supp_tables/ae_int_overlap_results_smok2.csv")

#
# ---  plot each of these --- #




# [x] pairwise correlation?
# - is correlation associated with platform
# 3. divide into disc/valid for meta --> re-run
# 4. meta sex-related

meta_int <- data.table::rbindlist(lapply(all_studies[fu_int %>% pull(study)], function(x) 
  prep_meta_study(x$gene_int, x$study)))
meta_int2 <- meta_int %>%
  filter(gene %in% ae_de_gene_i$gene) %>%
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup()
table(meta_int2$n)
meta_int3 <- meta_int2 %>% filter(n >=5) %>% group_split(gene)


# 5min
ma_vals_int <- lapply(meta_int3, meta_study)

mult_study_int <- data.frame(apply(data.table::rbindlist(ma_vals_int) ,
                               c(1,2), unlist)) %>%
  mutate(across(c(contains("logFC"), p), 
                ~as.numeric(as.character(.)) )) %>%
  mutate(n=as.numeric(as.character(n))) %>% 
  arrange(p)
mult_study_int$adj.p <- p.adjust(mult_study_int$p, method="fdr")
mult_study_int %>% filter(adj.p < 0.05)

ae_int <- data.table::rbindlist(lapply(all_studies[fu_int %>% filter(str_detect(tissue, "airway")) %>% pull(study)], function(x) 
  prep_meta_study(x$gene_int, x$study)))
ae_int2 <- ae_int %>%
  filter(gene %in% ae_de_gene_i$gene) %>%
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup() %>%
  group_split(gene)
  

ae_vals_int <- lapply(ae_int2, meta_study)
ae_study_int <- data.frame(apply(data.table::rbindlist(ae_vals_int) ,
                                   c(1,2), unlist)) %>%
  mutate(across(c(contains("logFC"), p), 
                ~as.numeric(as.character(.)) )) %>%
  mutate(n=as.numeric(as.character(n))) %>% 
  arrange(p)
ae_study_int$adj.p <- p.adjust(ae_study_int$p, method="fdr")
ae_study_int %>% filter(adj.p < 0.05)

chatz_de <- read_csv("ref/chatz_de.csv", skip=1) 
colnames(chatz_de) <- c("probe", "gene", "descript", "rank_m", "rank_f", "p_rank", "es_m", "es_f")
chatz_de$gene
blood_int <- data.table::rbindlist(lapply(all_studies[fu_int %>% filter(str_detect(tissue, "blood")) %>% pull(study)], function(x) 
  prep_meta_study(x$gene_int, x$study)))
blood_int2 <- blood_int %>%
  filter(gene %in% chatz_de$gene) %>%
  group_by(gene) %>% 
  mutate(n=n())  %>% 
  ungroup() %>%
  group_split(gene)


blood_vals_int <- lapply(blood_int2, meta_study)
blood_study_int <- data.frame(apply(data.table::rbindlist(blood_vals_int) ,
                                 c(1,2), unlist)) %>%
  mutate(across(c(contains("logFC"), p), 
                ~as.numeric(as.character(.)) )) %>%
  mutate(n=as.numeric(as.character(n))) %>% 
  arrange(p)
blood_study_int$adj.p <- p.adjust(blood_study_int$p, method="fdr")
blood_study_int %>% filter(adj.p < 0.05)

blood_int %>% filter(gene=="BFSP1")
ae_int %>% filter(gene=="BFSP1")
# grab a single gene

# --- AKR1C3 --- #

# -->
#  sample | expr | smoking_status 
get_gene <- function(ds, gene_name, orig=F){
  
  if (orig) {
    study_f <- sprintf("data/small_studies/%s.RData", ds)
    if (!file.exists(study_f)){
      return(NA)
    }
    phe_f <- sprintf("data/pdata_filt/%s.csv", tolower(ds))
    if (!file.exists(phe_f)){
      
      return(NA)
    }
    load(study_f)
    df <- gse$originalData[[1]]
    phe = read_csv(phe_f)
    if ("sex" %in% colnames(phe)){
      phe <- phe %>% dplyr::select(-sex)
    }
    phe <- phe %>% rename(sex=sex_lab)
  } else {
    study_f <- sprintf("data/small_gses2/%s.RData", ds)
    if (!file.exists(study_f)){
      return(NA)
    }
    phe_f <- sprintf("data/pdata_filt2_v2/%s.csv", tolower(ds))
    if (!file.exists(phe_f)){
      
      return(NA)
    }
    
    load(study_f)
    df <- my_gse
    
    phe = read_csv(phe_f) %>%
      mutate(study_acc=ds)
    if (!"sample_acc" %in% colnames(phe)){
      phe <- phe %>% rename(sample_acc=gsm) 
    }
  }
  phe <- phe %>% dplyr::select(study_acc, sample_acc, sex, smok)
  
  if (!"keys" %in% names(df)){
    df$keys <- rownames(df$expr)
  }
  vals = df$expr[which(df$keys==gene_name),] # FIX
  vals <- colMeans(vals)
  
  if (is.null(dim(vals))) {
    val_df = tibble("sample_acc"=names(vals),
                    "expr"=unlist(vals))
  } else {
    val_df = data.frame(t(vals))
    val_df$sample_acc = colnames(vals)
    
  }
  val_phe = val_df %>% left_join(phe)
  
  return(val_phe)
}

get_gene("E-MTAB-5278", "AKR1C3")
my_studies2 <- c("GSE14633", "GSE5056", "GSE27002",  "E-MTAB-5278", "E-MTAB-5279" ,
                 "GSE20681" , "GSE23323" , "GSE23515" , "GSE30272"  , "GSE32504"  ,  "E-MTAB-3604")

# BFSP1
load("data/ae_full_exp.RData")
ae <- expDat5["206746_at",]
ae_dat <- tibble(sample_acc=names(ae),
 expr=unlist(ae)) %>%
left_join(pDat5.1 %>% 
            dplyr::select(sex, geo_accession, smoking) %>%
            rename(smok=smoking), by=c("sample_acc"="geo_accession"))

bfsp1 <- get_gene("GSE13896", "BFSP1", orig=T) %>%
  bind_rows(get_gene("GSE14633", "BFSP1", orig=F)) 

bfsp1 %>% bind_rows(ae_dat %>% 
                      mutate(study_acc="Grouped AE") %>%
                      dplyr::select(colnames(bfsp1))) %>%
  filter(sex %in% c("male", "female"), !is.na("study_acc")) %>%
  ggplot(  
       aes(x=smok, col=sex, y=expr))+
  geom_boxplot()+
  geom_point(alpha=0.5, position=position_jitterdodge(dodge.width=0.7, jitter.width = 0.2))+
  theme_bw()+
  facet_grid(.~study_acc)

chatz_de %>% filter(gene=="BFSP1") 
# ES = ratio current to never smokers
# this does not match


my_studies <- c('GSE103174',
                'GSE13896','GSE16149','GSE17913','GSE18723','GSE19027',
                'GSE20189','GSE2125','GSE21862','GSE31210','GSE32539',
                'GSE42057','GSE42743','GSE4302','GSE44456','GSE46699',
                'GSE55962','GSE56768',
                'GSE7895','GSE87072', 'GSE8987','GSE994')

# AKR1C3
list_val_phe = lapply(my_studies, function(x) get_gene(x, "AKR1C3", orig=T))
#save(list_val_phe, file="data/akr1c3_vals.RData")
list_val_phe2 = lapply(my_studies2, function(x) get_gene(x, "AKR1C3"))
list2 <- list_val_phe2[sapply(list_val_phe2, function(x) "expr" %in% colnames(x))]
list1 <- list_val_phe[sapply(list_val_phe, function(x) "expr" %in% colnames(x))]


df <- do.call(rbind, c(list1, list2)) %>%
  filter(!is.na(smok), !is.na(sex), !is.na(study_acc)) 
  
  ggplot(  df, aes(x=smok, col=sex, y=expr))+
  geom_boxplot()+
  geom_point(alpha=0.5, 
             position=position_jitterdodge(dodge.width=0.7, jitter.width = 0.2))+
  theme_bw()+
  facet_grid(.~study_acc)


load("data/akr1c3_vals.RData")
df2 <- df %>% filter(smok %in% c("NS", "S"), 
                     sex %in% c("male", "female"))
df3 <- df2 %>% 
  left_join(list_studies %>% 
              dplyr::select(study, tissue,
                            samples), by=c("study_acc"="study"))
df4 <- df3 %>% filter(str_detect(tissue, "airway") | 
                 str_detect(tissue, "trachea") |
                 str_detect(tissue, "blood")) %>%
  mutate(smok=ifelse(smok=="S", "smoker", "non-smoker")) %>%
  filter(!str_detect(tissue, "b cell")) %>%
  #mutate(tissue2=ifelse(str_detect(tissue, "epithelium"), 
  #                      "airway", "blood")) %>%
  rename(smoking=smok) %>%
  rename(`sample size`=samples,
         study=study_acc) %>%
  mutate(tissue=factor(tissue, 
                       levels=c("airway epithelium", "trachea epithelium",
                                "blood - pbmcs", "blood - whole"))) %>%
  arrange(tissue)
df4$study <- factor(df4$study, levels=unique(df4$study))
color_bar(df4)  
ggsave("figures/color_bar_akr1c3.png")
ggplot(df4, aes(x=study, y=expr, col=smoking))+
  geom_violin()+
  stat_summary(fun=mean, geom="point", alpha=0.5, 
               position=position_dodge(0.9),
               aes(size=`sample size`)) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar",
               position=position_dodge(0.9),
               width=0.3)+
  #geom_point(alpha=0.5, position=position_jitterdodge(dodge.width = 0.9))+
  #geom_boxplot(width=0.2, position=position_dodge(0.9), outlier.shape=NA, coef = 0)+
  #facet_grid(.~tissue2, scales="free")+
  theme_bw()+
  xlab("")+
  theme(panel.grid.minor = element_blank())+
  scale_color_manual(values=c('turquoise', "gold4"))+
  ylab("log2 expression of AKR1C3")+
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
# add tissue labels
ggsave("figures/paper_figs/akr1c3_comb.png")

#df3 %>% ggplot(aes(x=study_acc, y=expr, col=smok))+
#  geom_boxplot()+
#  facet_grid(.~tissue, scales="free")+
#  theme_bw()+
#  xlab("")+
#  ylab("expression of AKR1C3")+
#  theme(axis.text.x = element_text(angle = 90, hjust=1))

