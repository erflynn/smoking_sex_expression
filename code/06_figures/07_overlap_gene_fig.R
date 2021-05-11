
library(RColorBrewer)
dark2 <- brewer.pal(8, "Dark2")
set2 <- brewer.pal(8, "Set2")
alluv_col <- c(dark2[1], set2[1], set2[4], set2[3], dark2[3], dark2[8])

# make a figure that shows the overlapping genes
library('tidyverse')
library('gridExtra')
load("data/all_studies_v2_clean.RData") # --> all_studies2
load("data/fig_color_scale.RData") # --> paired3

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


### functions for plotting ###


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
          axis.title.y=element_blank())
  
  
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
    ylab("")+
    xlab("")+
    theme(panel.grid.major=element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y=element_text(colour="white"),
          axis.ticks.y = element_blank(),
          axis.title.y=element_blank(),
          legend.position = "none")
  
  
  tmp <- ggplot_gtable(ggplot_build(p0))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  
  p0.clean <- p0+theme(legend.position="none")
  
  
  
  p3 <- ggplot(tau_df %>% 
           mutate(gene=factor(gene, levels=levels(olap_genes2$gene))), 
           aes(y=tau, x=gene, group=1))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position="none",
          axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, 
                                     color=c("white")))+
    #geom_segment(aes(xend=gene), yend=0, alpha=0.5) +   
    geom_line()+
    geom_point( aes(color=tau)) +
    expand_limits(y=0) +
    scale_y_reverse(
      breaks=c(0,0.5,1),
      labels=c(0, "E-MTAB-3604", 1))+
    coord_flip()+
    theme(axis.title.x=element_blank(), 
          axis.title.y=element_blank())
  
  
  grid.arrange(p3, p0.clean, p1,
               heights=c(100,my_height), widths=c(1,4.5), 
               layout_matrix=rbind(c(1,2), c(4, 3)))
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

plot_olap_fig(og, dn="turquoise", up = "gold4", 5, smok_df)
ggsave("figures/overlapping_genes.pdf")


# 2. sex
all_dfs_s <- lapply(all_studies2 , function(x) x$gene_sex)
names(all_dfs_s) <- names(all_studies2)
missing <- names(which(is.na(all_dfs_s)))
all_dfs2_s <- c(all_dfs_s[!names(all_dfs_s) %in% c("GSE44456", missing)], 
                "Grouped AE"=list(gene_sex_ae))
comb_sex <- get_gene_list(sex_genes, ae_overlap_s)
olap_genes_s <- get_olap_df(all_dfs2_s, comb_sex)

plot_olap_fig(olap_genes_s, dn= dark2[1], up = dark2[3], 7)
ggsave("figures/overlapping_genes_sex.pdf")  

sex_genes %>% filter(gene %in% comb_sex) %>% filter(!chromosome_name %in% c("X", "Y"))

# significant autosomal genes from all?



# plot sex [x]
# plot autosomal [ ]

# change colors for plotting [x] --> function
# figure out how to add legend
# add color bar on L side with gene membership?

# group together tissues:
#  airway
#  airway-adjacent
#  blood
#  other



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
tau_smok <- sapply(unique(og$gene), function(gene){
  calc_tau(og, gene)
})

# 1 = specific, 0 = ubiquitous

smok_df <- tibble("gene"=unique(og$gene),
       "tau"=tau_smok)
smok_df %>% arrange(desc(tau))
smok_df %>% arrange(tau) 

plot_olap_fig(og %>% 
                mutate(gene=factor(gene, 
                                   levels=smok_df %>% arrange(tau) %>% pull(gene))) %>%
                mutate(tissue=case_when(
                  str_detect(tissue, "blood") ~ "blood - whole", 
                  str_detect(tissue, "airway") | str_detect(tissue, "trachea") ~ "airway epithelium",
                  str_detect(tissue, "oral") | str_detect(tissue, "nasal") | 
                    str_detect(tissue, "buccal") | 
                    tissue=="sputum" ~ "oral cavity",
                  tissue %in% c("liver", "brain - prefrontal cortex", "kidney") ~ "other",
                  TRUE ~ tissue))  %>%
                filter(!tissue %in% c("liver", "brain - prefrontal cortex", "kidney")))



plot(density(smok_df$tau))

ggplot(smok_df %>% mutate(gene=factor(gene, levels=levels(og$gene))), aes(y=tau, x=gene))+
 # geom_bar(stat="identity")+
  #geom_point(alpha=0.5)+
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, colour="white"))+
  geom_segment(aes(xend=gene), yend=0, alpha=0.5) +   
  geom_point( color="darkgray") +
  expand_limits(y=0) +
  scale_y_continuous(
    breaks=c(0,0.5,1),
    labels=c("A", "E-MTAB-3604", "B"))+
  scale_y_reverse()+
  coord_flip()



tau_sex <- sapply(unique(olap_genes_s$gene), function(gene){
  calc_tau(olap_genes_s, gene)
})

df_sex <- tibble("gene"=unique(olap_genes_s$gene),
             "tau"=tau_sex) 
df_sex %>% arrange(tau)
df_sex %>% arrange(desc(tau))

t.test(df_sex$tau, smok_df$tau)
plot(density(df_sex$tau))

# UPDATE THE TISSUE LABELS + REMOVE THE OTHER TISSUES
o_in <- olap_genes_s %>% 
  mutate(gene=factor(gene, 
                     levels=df_sex %>% arrange(tau) %>% pull(gene))) %>%
  mutate(tissue=case_when(
    str_detect(tissue, "blood") ~ "blood - whole", 
    str_detect(tissue, "airway") | str_detect(tissue, "trachea") ~ "airway epithelium",
    str_detect(tissue, "oral") | str_detect(tissue, "nasal") | 
      str_detect(tissue, "buccal") | 
      tissue=="sputum" ~ "oral cavity",
    TRUE ~ tissue))  %>%
  filter(!tissue %in% c( "liver", "brain - prefrontal cortex", "kidney"))
plot_olap_fig(o_in)


ggplot(df_sex %>% mutate(gene=factor(gene, levels=df_sex %>% arrange(tau) %>% pull(gene))), aes(y=gene, x=tau, group=1))+
  geom_bar(stat="identity")+
  #geom_line()+
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_reverse()


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
