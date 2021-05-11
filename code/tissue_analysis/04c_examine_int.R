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

#####


# meta_int <- data.table::rbindlist(lapply(all_studies[fu_int %>% pull(study)], function(x) 
#   prep_meta_study(x$gene_int, x$study)))
# meta_int2 <- meta_int %>%
#   filter(gene %in% ae_de_gene_i$gene) %>%
#   group_by(gene) %>% 
#   mutate(n=n())  %>% 
#   ungroup()
# table(meta_int2$n)
# meta_int3 <- meta_int2 %>% filter(n >=5) %>% group_split(gene)
# 
# 
# # 5min
# ma_vals_int <- lapply(meta_int3, meta_study)
# 
# mult_study_int <- data.frame(apply(data.table::rbindlist(ma_vals_int) ,
#                                    c(1,2), unlist)) %>%
#   mutate(across(c(contains("logFC"), p), 
#                 ~as.numeric(as.character(.)) )) %>%
#   mutate(n=as.numeric(as.character(n))) %>% 
#   arrange(p)
# mult_study_int$adj.p <- p.adjust(mult_study_int$p, method="fdr")
# mult_study_int %>% filter(adj.p < 0.05)
# 
# ae_int <- data.table::rbindlist(lapply(all_studies[fu_int %>% filter(str_detect(tissue, "airway")) %>% pull(study)], function(x) 
#   prep_meta_study(x$gene_int, x$study)))
# ae_int2 <- ae_int %>%
#   filter(gene %in% ae_de_gene_i$gene) %>%
#   group_by(gene) %>% 
#   mutate(n=n())  %>% 
#   ungroup() %>%
#   group_split(gene)
# 
# 
# ae_vals_int <- lapply(ae_int2, meta_study)
# ae_study_int <- data.frame(apply(data.table::rbindlist(ae_vals_int) ,
#                                  c(1,2), unlist)) %>%
#   mutate(across(c(contains("logFC"), p), 
#                 ~as.numeric(as.character(.)) )) %>%
#   mutate(n=as.numeric(as.character(n))) %>% 
#   arrange(p)
# ae_study_int$adj.p <- p.adjust(ae_study_int$p, method="fdr")
# ae_study_int %>% filter(adj.p < 0.05)
# 
# chatz_de <- read_csv("ref/chatz_de.csv", skip=1) 
# colnames(chatz_de) <- c("probe", "gene", "descript", "rank_m", "rank_f", "p_rank", "es_m", "es_f")
# chatz_de$gene
# blood_int <- data.table::rbindlist(lapply(all_studies[fu_int %>% filter(str_detect(tissue, "blood")) %>% pull(study)], function(x) 
#   prep_meta_study(x$gene_int, x$study)))
# blood_int2 <- blood_int %>%
#   filter(gene %in% chatz_de$gene) %>%
#   group_by(gene) %>% 
#   mutate(n=n())  %>% 
#   ungroup() %>%
#   group_split(gene)
# 
# 
# blood_vals_int <- lapply(blood_int2, meta_study)
# blood_study_int <- data.frame(apply(data.table::rbindlist(blood_vals_int) ,
#                                     c(1,2), unlist)) %>%
#   mutate(across(c(contains("logFC"), p), 
#                 ~as.numeric(as.character(.)) )) %>%
#   mutate(n=as.numeric(as.character(n))) %>% 
#   arrange(p)
# blood_study_int$adj.p <- p.adjust(blood_study_int$p, method="fdr")
# blood_study_int %>% filter(adj.p < 0.05)
# 
# blood_int %>% filter(gene=="BFSP1")
# ae_int %>% filter(gene=="BFSP1")