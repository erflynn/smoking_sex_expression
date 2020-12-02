# 01_run_lioness.R
#
# Code for runnning lioness
# TODO:
# - clean up the code for processing the resulting edges/genes
# - identify parameter cutoffs and assess sensitivity

library(lionessR)
library(igraph)
library(limma)
library(SummarizedExperiment)
library(tidyverse)
library(ggrepel)



# ---- helpful functions---- #

# make it a summarized EXP
makeSumExp <- function(exp_data, phe_data){
  rowData <- DataFrame(row.names = rownames(exp_data), 
                       gene = rownames(exp_data))
  colData <- phe_data %>% dplyr::rename("sample"="geo_accession")
  rownames(colData) <- colData$sample
  
  se <- SummarizedExperiment(assays = list(counts = as.matrix(exp_data)), colData = colData, rowData = rowData)
  return(se)
}


# grab 500 most variably expressed genes
varExp <- function(se, nsel=500) {
  cvar <- apply(assay(se), 1, sd)
  dat <- se[tail(order(cvar), nsel), ]
  return(dat)
}

varExpEither <- function(se1, se2, nsel=500){
  cvar1 <- apply(assay(se1), 1, sd)
  cvar2 <- apply(assay(se2), 1, sd)
  l1 <- tail(names(cvar1)[order(cvar1)], nsel)
  l2 <- tail(names(cvar2)[order(cvar2)], nsel)
  genes <- union(l1, l2)
  dat1 <- se1[genes,]
  dat2 <- se2[genes,]
  return(list("d1"=dat1, "d2"=dat2))
}

calcNetDiff <- function(dat, my_col){
  netyes <- cor(t(assay(dat)[, dat[[my_col]] == "Y"]))
  netno  <- cor(t(assay(dat)[, dat[[my_col]] =="N"]))
  return(netyes-netno)
} # //TODO - change to general y/n


# matrices --> edge lists
mat2elist <- function(netdiff, cor.cut=0.5){
  library(reshape2)
  nrows=nrow(netdiff)
  cormat2 <- rep(1:nrows, each=nrows)
  cormat1 <- rep(1:nrows,nrows)
  el <- cbind(cormat1, cormat2, c(netdiff))
  melted <- reshape2::melt(upper.tri(netdiff))
  melted <- melted[which(melted$value),]
  values <- netdiff[which(upper.tri(netdiff))]
  melted <- cbind(melted[,1:2], values)
  genes <- row.names(netdiff)
  melted[,1] <- genes[melted[,1]]
  melted[,2] <- genes[melted[,2]]
  row.names(melted) <- paste(melted[,1], melted[,2], sep="_")
  tosub <- melted
  tosel <- row.names(tosub[which(abs(tosub[,3])>cor.cut),])
  return(tosel)
}



deEdges <- function(se, corsub, my_col){
  # get 50 most DE edges -- limma on EDGES
  group <- factor(se[[my_col]]) 
  design <- model.matrix(~0+group)
  cont.matrix <- makeContrasts(yesvsno = (groupY - groupN), levels = design)   
  fit <- lmFit(corsub, design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2e <- eBayes(fit2)
  toptable <- topTable(fit2e, number=nrow(corsub), adjust="fdr")
  return(toptable)
}


deGenes <- function(se, edges, my_col){
  #topgeneslist <-unique(unlist(lapply(rownames(edges)[1:nrow()], function(x) str_split(x, "_")[[1]])) )
  group <- factor(se[[my_col]]) 
  design <- model.matrix(~0+group)
  cont.matrix <- makeContrasts(yesvsno = (groupY - groupN), levels = design)    # TODO generalize
  fit <- lmFit(assay(se), design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2e <- eBayes(fit2)
  topDE <- topTable(fit2e, number=nrow(assay(se)), adjust="fdr")
  #topDE <- topDE[which(row.names(topDE) %in% topgeneslist),]
  #topgenesDE <- tibble("gene"=row.names(topDE), "t"=topDE$t)
  return(topDE)
}

reformatEdge <- function(x) {
  y <- str_split(x, "_")[[1]];
  z <- y[order(y)];
  return(paste(z, collapse="_"))
}


# plot!!!

plotG <- function(genes, g, g_weight, my_col="t.s"){
  E(g)$weight <- as.numeric(g_weight)
  E(g)$color[E(g)$weight<0] <- "blue"
  E(g)$color[E(g)$weight>0] <- "red"
  E(g)$weight <- 1
  
  # coloring
  nodeorder <- tibble("gene"=V(g)$name, "order"=1:length(V(g)))
  nodes <- left_join(nodeorder, genes) %>%
    replace_na(replace=list("t.s"=0, "t.ns"=0))
  V(g)$weight <- as.numeric(nodes[[my_col]])
  
  # make a color palette
  mypalette4 <- colorRampPalette(c("blue","white","white","red"), space="Lab")(256) 
  breaks2a <- seq(min(V(g)$weight), 0, length.out=128)
  breaks2b <- seq(0.00001, max(V(g)$weight)+0.1,length.out=128)
  breaks4 <- c(breaks2a,breaks2b)
  
  # select bins for colors
  bincol <- rep(NA, length(V(g)))
  for(i in 1:length(V(g))){
    bincol[i] <- min(which(breaks4>V(g)$weight[i]))
  }
  bincol <- mypalette4[bincol]
  V(g)$weight <- 1
  
  # add colors to nodes
  V(g)$color <- bincol
  par(mar=c(0,0,0,0))
  plot(g, vertex.label.cex=0.7, vertex.size=10, 
       vertex.label.color = "black", vertex.label.font=3,
       edge.width=5*(abs(as.numeric(g_weight))))
}

# ---- end functions ---- #


load("data/ae_only_eset.RData") # --> ae_only, ae_only_meta

# reformat metadata fields and grab the subsets of the data
alt_meta <- ae_only_meta %>%
  dplyr::select(geo_accession, smok, expr_sex) %>%
  mutate(is_male=ifelse(expr_sex=="male", "Y", "N")) %>%
  mutate(is_smoker=ifelse(smok=="S", "Y", "N")) %>%
  dplyr::select(-smok, -expr_sex)

f_meta <- alt_meta %>%
  filter(is_male=="N") %>%
  dplyr::select(-is_male)

m_meta <- alt_meta %>%
  filter(is_male=="Y") %>%
  dplyr::select(-is_male)

s_meta <- alt_meta %>% 
  filter(is_smoker=="Y") %>%
  dplyr::select(-is_smoker)

ns_meta <- alt_meta %>% 
  filter(is_smoker=="N") %>%
  dplyr::select(-is_smoker)


# prepare the data and filter by variance
f_se <- makeSumExp(ae_only[,f_meta$geo_accession], f_meta )
m_se <- makeSumExp(ae_only[,m_meta$geo_accession], m_meta )

s_se <- makeSumExp(ae_only[,s_meta$geo_accession], s_meta )
ns_se <- makeSumExp(ae_only[,ns_meta$geo_accession], ns_meta )

f_se2 <- varExp(f_se)
m_se2 <- varExp(m_se)
fm_vars <- varExpEither(f_se, m_se)
f_se2 <- fm_vars[[1]]
m_se2 <- fm_vars[[2]]

s_se2 <- varExp(s_se)
ns_se2 <- varExp(ns_se)
sn_vars <- varExpEither(s_se, ns_se)
s_se2 <- sn_vars[[1]]
ns_se2 <- sn_vars[[2]]

f_ndf <- calcNetDiff(f_se2, "is_smoker")
m_ndf <- calcNetDiff(m_se2, "is_smoker")

s_ndf <- calcNetDiff(s_se2, "is_male")
ns_ndf <- calcNetDiff(ns_se2, "is_male")

f_tosel <- mat2elist(f_ndf)
m_tosel <- mat2elist(m_ndf)
fm_tosel <- union(f_tosel, m_tosel)

s_tosel <- mat2elist(s_ndf)
ns_tosel <- mat2elist(ns_ndf)
sn_tosel <- union(s_tosel, ns_tosel)

# create coexp network
f_cormat <- lioness(f_se2, netFun)
f_corsub <- assay(f_cormat[which(row.names(f_cormat) %in% fm_tosel), ])
m_cormat <- lioness(m_se2, netFun)
m_corsub <- assay(m_cormat[which(row.names(m_cormat) %in% fm_tosel), ])

s_cormat <- lioness(s_se2, netFun)
s_corsub <- assay(s_cormat[which(row.names(s_cormat) %in% sn_tosel), ])
ns_cormat <- lioness(ns_se2, netFun)
ns_corsub <- assay(ns_cormat[which(row.names(ns_cormat) %in% sn_tosel), ])


# grab DE edges and genes
f_edges <- deEdges(f_se, f_corsub, "is_smoker")
m_edges <- deEdges(m_se, m_corsub, "is_smoker")

s_edges <- deEdges(s_se, s_corsub, "is_male")
ns_edges <- deEdges(ns_se, ns_corsub, "is_male")

f_genes <- deGenes(f_se2, f_edges, "is_smoker")
m_genes <- deGenes(m_se2, m_edges, "is_smoker")

f_genes$gene <- rownames(f_genes)
m_genes$gene <- rownames(m_genes)

# put together the DE genes
comb_genes <- f_genes %>% 
  full_join(m_genes, by=c("gene")) %>%
  mutate(gene=as.numeric(gene)) %>%
  left_join(convert_genes, by=c("gene"="entrezgene_id")) %>%
  dplyr::rename("entrez"="gene", "gene"="hgnc_symbol")
sig_genes <- comb_genes %>% 
  dplyr::select(-contains("AveExpr"), -contains("P.Value"), -contains("B")) %>% 
  filter(adj.P.Val.x < 0.05 | adj.P.Val.y < 0.05) %>%
  as_tibble()

opp_sig <- sig_genes %>% filter(abs((logFC.x - logFC.y)/logFC.x )> 3)

# plot gene effect sizes
ggplot(comb_genes, aes(x=logFC.x, y=logFC.y))+
  geom_point(alpha=0.3)+
  geom_point(data=sig_genes, col="red", alpha=0.5)+
  geom_point(data=opp_sig, col="blue", alpha=0.5)+
  
  geom_vline(xintercept=0, col="gray", alpha=0.8)+
  geom_hline(yintercept=0, col="gray", alpha=0.8)+
  theme_bw()+
  xlab("log(S-NS) in females")+
  ylab("log(S-NS) in males")+
  geom_label_repel(data=sig_genes %>% 
                     filter(abs(logFC.x) > 20 |
                              abs(logFC.y) >20), 
                   aes(label=gene), size=2)+
  geom_label_repel(data=opp_sig, col="purple", 
                   aes(label=gene), size=2)



# data processing for edges/gennes
# TODO: make into functions and clean up
s_genes <- deGenes(s_se2, s_edges, "is_male")
ns_genes <- deGenes(ns_se2, ns_edges, "is_male")
s_genes$gene <- rownames(s_genes)
ns_genes$gene <- rownames(ns_genes)

s_genes_sig <- s_genes %>% filter(adj.P.Val < 0.1)
ns_genes_sig <- ns_genes %>% filter(adj.P.Val < 0.1)
setdiff(s_genes_sig$gene, ns_genes_sig$gene)


sns_genes <- full_join(
  s_genes %>% dplyr::rename(t.s=t), 
  ns_genes %>% dplyr::rename(t.ns=t), by=c("gene")) %>%
  mutate(gene=as.numeric(gene)) %>%
  left_join(convert_genes, by=c("gene"="entrezgene_id")) %>%
  dplyr::rename("entrez"="gene", "gene"="hgnc_symbol")

s_edges2 <- s_edges %>% dplyr::select(logFC) 
s_edges2$edge <- rownames(s_edges)
s_edges2$edge <- sapply(s_edges2$edge, reformatEdge)

ns_edges2 <- ns_edges %>% dplyr::select(logFC) 
ns_edges2$edge <- rownames(ns_edges)
ns_edges2$edge <- sapply(ns_edges2$edge, reformatEdge)

sns_edges <- full_join(
  s_edges2 %>% dplyr::rename(logFC.s=logFC), 
  ns_edges2 %>% dplyr::rename(logFC.ns=logFC), by=c("edge"))

sns_edges2 <- sns_edges %>% 
  separate(edge, into=c("e1", "e2"), sep="_") %>% 
  dplyr::select(e1, e2, everything()) %>%
  mutate(e2=as.numeric(e2)) %>%
  left_join(convert_genes %>% 
              dplyr::select(entrezgene_id, hgnc_symbol), by=c("e2"="entrezgene_id")) %>%
  dplyr::rename("entrez2"="e2", "e2"="hgnc_symbol") %>%
  mutate(e1=as.numeric(e1)) %>%
  left_join(convert_genes %>% 
              dplyr::select(entrezgene_id, hgnc_symbol), by=c("e1"="entrezgene_id")) %>%
  dplyr::rename("entrez1"="e1", "e1"="hgnc_symbol") 
  
# top 50 edges
sns_edges3 <- bind_rows(sns_edges2 %>% 
                          arrange(desc(abs(logFC.s))) %>% head(25), 
                        sns_edges2 %>% arrange(desc(abs(logFC.ns))) %>% head(25)) %>% 
  unique()


fm_genes <- full_join(
  f_genes %>% dplyr::rename(t.f=t), 
  m_genes %>% dplyr::rename(t.m=t), by=c("gene")) %>%
  mutate(gene=as.numeric(gene)) %>%
  left_join(convert_genes, by=c("gene"="entrezgene_id")) %>%
  dplyr::rename("entrez"="gene", "gene"="hgnc_symbol")

f_edges2 <- f_edges %>% dplyr::select(logFC) 
f_edges2$edge <- rownames(f_edges)
f_edges2$edge <- sapply(f_edges2$edge, reformatEdge)

m_edges2 <- m_edges %>% dplyr::select(logFC) 
m_edges2$edge <- rownames(m_edges)
m_edges2$edge <- sapply(m_edges2$edge, reformatEdge)

fm_edges <- full_join(
  f_edges2 %>% dplyr::rename(logFC.f=logFC), 
  m_edges2 %>% dplyr::rename(logFC.m=logFC), by=c("edge"))

fm_edges2 <- fm_edges %>% separate(edge, into=c("e1", "e2"), sep="_") %>% 
  dplyr::select(e1, e2, everything()) %>%
  mutate(e2=as.numeric(e2)) %>%
  left_join(convert_genes %>% dplyr::select(entrezgene_id, hgnc_symbol), by=c("e2"="entrezgene_id")) %>%
  dplyr::rename("entrez2"="e2", "e2"="hgnc_symbol") %>%
  mutate(e1=as.numeric(e1)) %>%
  left_join(convert_genes %>% dplyr::select(entrezgene_id, hgnc_symbol), by=c("e1"="entrezgene_id")) %>%
  dplyr::rename("entrez1"="e1", "e1"="hgnc_symbol") 


fm_edges3 <- bind_rows(fm_edges2 %>% 
                         arrange(desc(abs(logFC.f))) %>% head(25), 
                       fm_edges2 %>% arrange(desc(abs(logFC.m))) %>% head(25)) %>% unique()




# plot the network data
# make sure to set the same seeds
g1 <- graph.data.frame(sns_edges3 %>% dplyr::select(e1, e2), directed=FALSE)
set.seed(1)
plotG(sns_genes, g1, sns_edges3$logFC.s, my_col="t.s")
set.seed(1)
plotG(sns_genes, g1, sns_edges3$logFC.ns, my_col="t.ns")
g2 <- graph.data.frame(fm_edges3 %>% dplyr::select(e1, e2), directed=FALSE)
set.seed(1)
plotG(fm_genes, g2, fm_edges3$logFC.f, my_col="t.f")
set.seed(1)
plotG(fm_genes, g2, fm_edges3$logFC.m, my_col="t.m")

# plot the effect sizes of the edges
ggplot(fm_edges2, 
       aes(x=logFC.f, y=logFC.m))+
  geom_vline(xintercept=0, col="gray", alpha=0.8)+
  geom_hline(yintercept=0, col="gray", alpha=0.8)+
  geom_point(alpha=0.3)+theme_bw()+
  geom_point(data=fm_edges4, col="red")+
  geom_label_repel(data=fm_edges4 %>% 
                     unite(col="edge", c(e1, e2), sep="-"),size=3, aes(label=edge))+
  xlab("log(S-NS) in females")+
  ylab("log(S-NS) in males")

# look at the lists of genes
fm_edges3 %>% 
  filter(logFC.f < 0) %>% 
  dplyr::select(e1, e2) %>%
  pivot_longer(c(e1, e2), values_to=c("gene")) %>%
  group_by(gene)  %>%
  dplyr::count() %>% 
  dplyr::arrange(desc(n))

fm_edges3 %>% 
  filter(logFC.f > 0) %>% 
  dplyr::select(e1, e2) %>%
  pivot_longer(c(e1, e2), values_to=c("gene")) %>%
  group_by(gene)  %>%
  dplyr::count() %>% 
  dplyr::arrange(desc(n))

