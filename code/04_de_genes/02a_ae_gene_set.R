# 02a_ae_gene_set.R
#
# Incomplete scratch code for gene set analysis
# re-evaluate what gene sets/data before running.



library(tidyverse)
# GO analysis - GSVA
library("GSVA")
library(GSEABase)

load("data/ae_only_eset.RData") # --> ae_only, ae_only_meta

# how does this match the previous analysis?
go_gsva <- gsva(as.matrix(ae_only), go.gs, method="ssgsea")
# NOTE - this takes ten mins and not full GO dataset
save(go_gsva, file="data/go_gsva_ae.RData") 

# filter to remove low var genes
library(GSVAdata)
data(c2BroadSets)

gene_sds <- apply(ae_only, 1, sd)
summary(gene_sds)
ae_filt <- ae_only[gene_sds > median(gene_sds),]

# TAKES 10-15min
c2_gsva <- gsva(as.matrix(ae_filt), c2BroadSets, method="ssgsea")
save(c2_gsva, file="data/c2_gsva_ae.RData")



# adj.P < 0.05 
smok <- factor(ae_only_meta$smok)
sex <- factor(ae_only_meta$expr_sex)

#batch <- factor(blood_pDat2$batch)
design<- model.matrix(~smok+sex+sex*smok)
fit_g <- lmFit(c2_gsva, design)
fit_g <- eBayes(fit_g)

fix_rows <- function(df) {
  df$gs <- rownames(df)
  df %>% as_tibble()
}
gs_smok <- topTable(fit_g, coef="smokS", number=nrow(c2_gsva)) %>% fix_rows()
gs_sex <- topTable(fit_g, coef="sexmale", number=nrow(c2_gsva)) %>% fix_rows()
gs_int <- topTable(fit_g, coef="smokS:sexmale", number=nrow(c2_gsva)) %>% fix_rows()
gs_int  %>% filter(adj.P.Val < 0.05) %>% nrow() # 30
gs_sex  %>% filter(adj.P.Val < 0.05) %>% nrow() # 129
gs_smok  %>% filter(adj.P.Val < 0.05) %>% nrow() # 1165


# //todo: add in gene set size
plotTopGeneSets <- function(df, pcut=0.05){
  df %>% 
    filter(adj.P.Val < 0.05) %>%
    mutate(gs=factor(gs, 
                     levels=df %>% filter(adj.P.Val < 0.05) %>% pull(gs))) %>%
  ggplot(aes(x=logFC, y=gs, col=adj.P.Val))+
  geom_point()+
  theme_bw()+
  ylab("")
}

plotTopGeneSets(gs_int)+ggtitle("Top gene sets smoking*sex")
plotTopGeneSets(gs_sex)+ggtitle("Top gene sets: sex")
plotTopGeneSets(gs_smok)+ggtitle("Top gene sets: smoking")


# --- deprecated: wordcloud --- #
# remove_src <- function(x) {
#   y <- str_split(x, "_")[[1]]; 
#   return(paste(setdiff(y, y[[1]]), collapse=";"))
#   }
# gs_int2 <- gs_int %>% filter(adj.P.Val < 0.05) %>%
#   group_by(gs) %>%
#   mutate(gs2=remove_src(gs))
# words <- gs_int2 %>%
#   ungroup() %>%
#   dplyr::select(gs2) %>%
#   separate_rows(gs2, sep=";") %>%
#   group_by(gs2) %>%
#   count() %>%
#   arrange(desc(n)) %>%
#   filter(!gs2 %in% c("DN", "UP", "TO", "TARGETS",
#                      "VS", "PATHWAY", "RESPONSE", "48HR",
#                      "8", "ACROSS", "EUROPEAN"))
# library("wordcloud")
# wordcloud(words = words$gs2, freq = words$n, min.freq = 1,
#           max.words=200, random.order=FALSE, rot.per=0.35, 
#           colors=brewer.pal(8, "Dark2"))

