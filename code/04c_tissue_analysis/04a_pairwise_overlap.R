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