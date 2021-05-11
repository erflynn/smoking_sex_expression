
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