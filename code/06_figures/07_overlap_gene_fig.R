
# code for examining tissue specificity
# calculates tau
# plots figures with this information

library('tidyverse')
library('gridExtra')

# set up colors #
library(RColorBrewer)
dark2 <- brewer.pal(8, "Dark2")
set2 <- brewer.pal(8, "Set2")
alluv_col <- c(dark2[1], set2[1], set2[4], set2[3], dark2[3], dark2[8])
load("data/fig_color_scale.RData") # --> paired3

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# load general data 
load("data/all_studies_v2_clean.RData") # --> all_studies2
tiss_levels =c("airway epithelium", "trachea epithelium",
               "nasal epithelium", "oral cavity",
               "buccal mucosa", "sputum", "alveolar macrophages", 
               "lung", 
               "blood - b cells", "blood - pbmcs", "blood - whole",
               "liver", "kidney", 
               "brain - prefrontal cortex")
tiss_to_col <- c(paired3, "black", "gray")
names(tiss_to_col) <- c(tiss_levels, "all", "other")

list_studies <- read_csv("data/table1b_updated.csv")

list_studies2 <- list_studies %>% 
  filter(study!="GSE44456") %>%
  dplyr::select(study, tissue) %>%
  mutate(across(c(tissue, study), as.character)) %>%
  bind_rows(tibble("study"="Grouped AE", "tissue"="airway epithelium"))

# ------- FUNCTIONS ------ #

group_tissues <- function(df){
  df %>%
    mutate(tissue=case_when(
      str_detect(tissue, "blood") ~ "blood - whole", 
      str_detect(tissue, "airway") | str_detect(tissue, "trachea") ~ "airway epithelium",
      str_detect(tissue, "oral") | str_detect(tissue, "nasal") | 
        str_detect(tissue, "buccal") | 
        tissue=="sputum" ~ "oral cavity",
      tissue %in% c("liver", "brain - prefrontal cortex", "kidney") ~ "other")) 
}

calc_tau <- function(df, my_gene){
  filt_df <- df %>% 
    filter(gene==my_gene) %>% 
    filter(!is.na(logFC))
  
  tiss_df0 <- filt_df %>% 
    group_tissues() %>%
    filter(tissue!="other")
  # make sure at least two per group
  tiss_df0.1 <- tiss_df0 %>% semi_join(tiss_df0 %>% 
                                         group_by(tissue) %>% 
                                         count() %>% filter(n>1))
  tiss_df0.2 <- tiss_df0.1 %>%
    group_by(tissue) %>%
    summarize(medianLF=abs(median(logFC)), .groups="drop_last") # TODO take into acct direction?
  
  maxFC <- max(tiss_df0.2$medianLF) # max - 1?
  if (nrow(tiss_df0.2) < 2){
    return(NA)
  }
  tiss_df2 <- tiss_df0.2 %>% mutate(hat_xi=medianLF/maxFC) 
  return(sum(1-tiss_df2$hat_xi)/(nrow(tiss_df2)-1))
}



# function to get the list of genes
get_gene_list <- function(genes, ae_genes) {
  genes1 <- genes %>% 
    filter(n>=3)
  ae_overlap1 <- ae_genes %>%
    filter(nstudies>=3)
  
  ae2 <- ae_genes %>% filter(nstudies < 3) %>% pull(gene)
  other2 <- genes %>% filter(n < 3) %>% pull(gene)
  to_add <- intersect(ae2, other2)
  comb_l <- c(union(genes1$gene, ae_overlap1$gene), to_add) # 79 genes
  return(comb_l)
}


# get overlapping gene df in order
get_olap_df <- function(my_dfs, my_genes){
    olap_genes <- do.call(rbind, 
                        lapply(names(my_dfs), function(x) 
                        {my_dfs[[x]] %>% 
                            ungroup() %>%
                            filter(gene %in% my_genes) %>%
                            select(gene, logFC, adj.p) %>%
                            mutate(study=x)})) #%>%
      # as_tibble() %>%
      # group_by(study) %>%
      # mutate(min_study=min(logFC, na.rm=T),
      #        max_study=max(logFC, na.rm=T)) %>%
      # ungroup() %>%
      # group_by(study, gene) %>%
      # mutate(logFC=case_when(
      #   logFC > 0 ~ logFC/max_study,
      #   logFC < 0 ~ logFC/abs(min_study),
      #   logFC == 0 ~ 0
      # )) %>%
      # select(-min_study, -max_study) %>%
      # ungroup()

  
    olap_matrix <- olap_genes  %>% 
      select(-adj.p) %>%
      pivot_wider(names_from=study, values_from=logFC, values_fill=0)
    olap_matrix2 <- as.matrix(olap_matrix %>% select(-gene))

    rownames(olap_matrix2) <- olap_matrix$gene
    dist_olap <- dist(x = olap_matrix2)
    hc <- hclust(d = dist_olap)
    new_order <- olap_matrix$gene[hc$order]
    
    dist_olap_t <- dist(x = t(olap_matrix2))
    hc_t <- hclust(d = dist_olap_t)
    new_order2 <- colnames(olap_matrix2)[hc_t$order]
    
    olap_genes2 <- olap_genes %>%
      select(-adj.p) %>%
      pivot_wider(names_from="study", values_from="logFC", values_fill=NA) %>%
      pivot_longer(-gene, names_to="study", values_to="logFC") %>%
      left_join(list_studies2 %>% select(study, tissue), by="study") %>%
      filter(!is.na(tissue)) %>%
      arrange(tissue, study) %>%
      mutate(gene=factor(gene, levels=new_order)) %>%
      mutate(study=factor(study, levels=new_order2)) 
    return(olap_genes2)
}



plot_olap_fig <- function(olap_genes, dn, up, my_height, tau_df){
  olap_genes2 <- olap_genes %>% group_tissues()
  longest_gene <- unique(olap_genes2$gene)[which.max(sapply(unique(olap_genes2$gene), 
                                            function(x) str_length(as.character(x))))]
  p0 <- olap_genes2 %>%
    ggplot(aes(x=study, y=gene))+
    geom_tile(aes(fill=logFC))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))+
    scale_fill_gradient2(low = dn, high=up, mid = "white", 
                                          midpoint = 0,  na.value="#DCDCDC")+
    #scale_fill_gradient2(low = "turquoise", high = "gold4", mid = "white", 
    #                     midpoint = 0,  na.value="#DCDCDC")+
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    coord_flip()
  
  
  p1 <- ggplot(olap_genes2)+
    geom_bar(mapping = aes(x = study, y = longest_gene, fill = tissue), 
             stat = "identity", 
             width = 1)+
    scale_fill_manual(values=
                        tiss_to_col[olap_genes2 %>% 
                                      distinct(tissue, study) %>% 
                                      mutate(tissue=as.character(tissue)) %>%
                                      pull(tissue)])+
    theme_minimal()+
    theme(panel.grid.major=element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x= element_text(angle = 90, hjust=1, vjust=0.5, color="white"),
          axis.ticks.y = element_blank(),
          axis.title.y=element_blank(),
          axis.title.x=element_blank())+
    coord_flip()
  

  p3 <- ggplot(tau_df %>% 
           mutate(gene=factor(gene, levels=levels(olap_genes2$gene))), 
           aes(y=tau, x=gene, group=1))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          #axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(color=c("black", "white", "black" ,"white", "black")))+
    #geom_segment(aes(xend=gene), yend=0, alpha=0.5) +   
    geom_line()+
    geom_point( aes(color=tau)) +
    #expand_limits(y=0) +
    scale_y_continuous(
      limits=c(0,1),
      breaks=c(0,0.25, 0.5, 0.75,1),
      labels=c(0, 0.25, 0.5, "E-MTAB-3604", 1))+
    theme(axis.title.x=element_blank(), 
          axis.title.y=element_blank())
  p0.legend <- get_legend(p0) 
  p0.clean <- p0+theme(legend.position="none")
  p1.clean <-  p1+theme(legend.position="none")
  p3.clean <-  p3+theme(legend.position="none")
  p1.legend <-get_legend(p1) 
  p3.legend <- get_legend(p3) 
  
  legend_plot <- arrangeGrob(p0.legend, p1.legend, p3.legend, 
                             layout_matrix = rbind(c(1,3), c(2,2)) )
  grid.arrange(p1.clean, p0.clean, p3.clean, legend_plot,
               heights=c(5, 1), widths=c(my_height, 100, 20), 
               layout_matrix=rbind(c(1,2, 4), c(5, 3, 6)))
}




##### get data for plotting ######

# 0. load data for both
ae_overlap <- read_csv("data/supp_tables/ae_overlap2v2.csv")
smok_genes <- read_csv("data/supp_tables/overlapping_genes2_smok.csv")
ae_overlap_s <- read_csv("data/supp_tables/ae_overlap2v2_sex.csv")
sex_genes <- read_csv("data/supp_tables/overlapping_genes2_sex.csv") 

load("data/results/ae_smok.RData")
gene_smok_ae <- ma_smok %>%
  ungroup() %>%
  dplyr::select(-chromosome)
gene_sex_ae <- ma_sex %>%
  ungroup() %>%
  dplyr::select(-chromosome)


# 1. smoking
smok_dfs <- lapply(all_studies2, function(x) x$gene_smok)
names(smok_dfs) <- names(all_studies2)
smok_dfs2 <- c(smok_dfs[!names(smok_dfs) %in% c("GSE44456")], list("Grouped AE"=gene_smok_ae))
comb_l <- get_gene_list(smok_genes, ae_overlap)
og <- get_olap_df(smok_dfs2, comb_l)
tau_smok <- sapply(unique(og$gene), function(gene){ calc_tau(og, gene)})
smok_df <- tibble("gene"=unique(og$gene),"tau"=tau_smok)

# write out table
comb_tab <- function(overlap_df, gene_df, tau_df){
  ae_overlap_reform <- overlap_df %>% 
    mutate(direction=ifelse(logFC > 0, "up", "down"),
           studies=sprintf("Grouped AE;%s", study),
           tissues=case_when(
             str_detect(tissue, "airway") ~ tissue,
             TRUE ~ sprintf("airway epithelium;%s", tissue)
           ),
           n=nstudies+1) %>%
    select(-chromosome_name) %>%
    left_join(gene_df %>% select(chromosome_name, gene)) %>%
    select(colnames(gene_df)) 
  comb_genes <- gene_df %>% 
    anti_join(ae_overlap_reform, by="gene") %>%
    bind_rows(ae_overlap_reform) %>%
    arrange(desc(n)) 
  comb_genes %>% filter(n>2) %>% left_join(tau_df)
}

my_tab_smok <- comb_tab(ae_overlap, smok_genes, smok_df) %>%
  mutate(chromosome_name=ifelse(gene=="MT1M", 16, chromosome_name)) %>%
  select(chromosome_name, gene, direction, tau, n, everything()) %>%
  rename(`in smokers`=direction)
my_tab_smok %>% write_csv("data/supp_tables/gene_spec_smok.csv")

plot_olap_fig(og, dn="turquoise", up = "gold4", 3, smok_df)
#ggsave("figures/overlapping_genes.pdf") # 10.5

# 2. sex
all_dfs_s <- lapply(all_studies2 , function(x) x$gene_sex)
names(all_dfs_s) <- names(all_studies2)
missing <- names(which(is.na(all_dfs_s)))
all_dfs2_s <- c(all_dfs_s[!names(all_dfs_s) %in% c("GSE44456", missing)], 
                "Grouped AE"=list(gene_sex_ae))
comb_sex <- get_gene_list(sex_genes, ae_overlap_s)
olap_genes_s <- get_olap_df(all_dfs2_s, comb_sex)


tau_sex <- sapply(unique(olap_genes_s$gene), function(gene){calc_tau(olap_genes_s, gene)})

df_sex <- tibble("gene"=unique(olap_genes_s$gene), "tau"=tau_sex) 

my_tab_sex <- comb_tab(ae_overlap_s, sex_genes, df_sex) %>%
  select(chromosome_name, gene, direction, tau, n, everything()) %>%
  mutate(chromosome_name=case_when(
    !is.na(chromosome_name) ~ chromosome_name,
    gene=="STS" ~ "X",
    gene=="OFD1" ~ "X"),
  direction=ifelse(direction=="down", "females", "males")) %>%
  rename(`higher in`=direction)
my_tab_sex %>% write_csv("data/supp_tables/gene_spec_sex.csv")
plot_olap_fig(olap_genes_s, dn= dark2[1], up = dark2[3], 5, df_sex)
#ggsave("figures/overlapping_genes_sex.pdf")  

my_tab_sex %>% filter(!chromosome_name %in% c("X", "Y") )


# 3. compare across the two
t.test(df_sex$tau, smok_df$tau) # p= 4.916e-10

df_sex %>% 
  mutate(term="sex") %>% 
  bind_rows(smok_df %>% mutate(term="smoking")) %>%
  ggplot(aes(x=tau, fill=term))+
  geom_histogram(position="identity", alpha=0.6)+
  geom_density(alpha=0.1, aes(color=term))+
  #geom_histogram(data=smok_df,fill = "red", alpha = 0.2) +
  #geom_histogram(data=df_sex,fill = "blue", alpha = 0.2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab("tissue specificity (tau)")+
  scale_fill_manual(values=c("purple", "gold4"))+
  scale_color_manual(values=c("purple", "gold4"))+
  ylab("number of genes")
ggsave("figures/tissue_spec_dist.pdf")


df_sex %>% 
  mutate(term="sex") %>% 
  bind_rows(smok_df %>% mutate(term="smoking")) %>%
  ggplot(aes(y=tau, x=term, color=term))+
  geom_violin(alpha=0.5)+
  geom_boxplot(width=0.2)+
  geom_point(alpha=0.4, position=position_jitter(width=0.1))+
  scale_color_manual(values=c("purple", "gold4"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")+
  xlab("")+
  ylab("tissue specificity (tau)")
ggsave("figures/tissue_spec_violin.pdf")
# [ ] add in tissue specific
# 


#### significant specificity??
# this does not work well

# calc_tau_sm <- function(tiss_df0.1){
#   tiss_df0.2 <- tiss_df0.1 %>%
#     mutate(absLogFC=abs(logFC)) %>%
#     group_by(tissue) %>%
#     summarize(medianLF=abs(median(logFC)), .groups="drop_last") # TODO take into acct direction?
#   
#   maxFC <- max(tiss_df0.2$medianLF) # max - 1?
#   if (nrow(tiss_df0.2) < 2){
#     return(NA)
#   }
#   tiss_df2 <- tiss_df0.2 %>% mutate(hat_xi=medianLF/maxFC) 
#   return(sum(1-tiss_df2$hat_xi)/(nrow(tiss_df2)-1))
# }
# 
# shuffle_tau <- function(df, my_gene){
#   filt_df <- df %>% 
#     filter(gene==my_gene) %>% 
#     filter(!is.na(logFC))
#   
#   tiss_df0 <- filt_df %>% 
#     mutate(tissue=case_when(
#       str_detect(tissue, "blood") ~ "blood", 
#       str_detect(tissue, "airway") | str_detect(tissue, "trachea") ~ "airway",
#       str_detect(tissue, "oral") | str_detect(tissue, "nasal") | 
#         str_detect(tissue, "buccal") | 
#         tissue=="sputum" ~ "oral airway",
#       TRUE ~ tissue))  
#   # make sure at least two per group
#   tiss_df0.1 <- tiss_df0 %>% semi_join(tiss_df0 %>% 
#                                          group_by(tissue) %>% 
#                                          count() %>% filter(n>1))
#   # now we shuffle
#   shuff_tau <- sapply(1:500, function(i){
#     tiss_shuff <- tiss_df0.1 %>% 
#       mutate(tissue=sample(tiss_df0.1$tissue, nrow(tiss_df0.1), replace=F))
#     calc_tau_sm(tiss_shuff)
#   })
#   return(shuff_tau)
# }
# 
# shuffs <- shuffle_tau(og, "AKR1B10")
# sum(shuffs > smok_df %>% filter(gene=="AKR1B10") %>% pull(tau)) # 49/500
