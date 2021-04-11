library(limma)
source("code/00_utils.R")
load("../tissue_meta_analysis/liver_meta_analysis_old/data/processed/ds8g.RData")
head(datasetInfo8$expr)
head(datasetInfo8$pheno)
head(datasetInfo8$keys)

# run DE analysis
sex <- as.factor(datasetInfo8$pheno$sexLabels2) 
design <- model.matrix(~sex)
fit <- lmFit(datasetInfo8$expr[,rownames(datasetInfo8$pheno)], design)
fit <- eBayes(fit)
tt <- data.frame(topTable(fit, number=nrow(datasetInfo8$expr),  confint = TRUE))
tt$probe <- rownames(tt)
tt$gene <- datasetInfo8$keys[tt$probe]

# MA genes to probes
liver_dat <- tt %>% mutate(chromosome="") %>% rename(geneSymbol=gene, ID=probe) 

liver_dat_ma <- run_clean_ma(liver_dat)
nrow(liver_dat_ma) # 19484
sig_genes <- liver_dat_ma %>% select(-chromosome) %>% filter(adj.p < 0.05) # 2070 (10.6%)
sig_genes %>% head(207) %>% ungroup() %>% select(gene) %>% write_tsv("data/list_liver_genes_sig.txt")
# convert to STRING for STAMS

ggplot(gene_mapped, aes(x=Pvalue))+
  geom_density()+
  theme_bw()+
  ggtitle("Input p-values for liver genes")

map_to_stams(liver_dat_ma, "data/gene_mapped_liver.RData")

###### look at the results #######
load.Rdata("data/stams_liver/norm_res.RData", "norm_liver") 
chosen_liver=moduleChoose(norm_liver, top=0.01, plot=FALSE)
length(chosen_liver[[1]]) # --> 74
norm_liver$ordered.module.score.matrix %>% filter(p==0) %>% nrow() # 1230
norm_liver$ordered.module.score.matrix %>% filter(p<0.05) %>% nrow() # 6621

norm_liver$ordered.module.score.matrix %>% nrow() # 7433


top_mods_liver <- data.frame(do.call(rbind, lapply(names(chosen_liver$modules), function(x) {
  vals <- chosen_liver$modules[[x]];
  name_convert <-  paste(gene_mapped %>% 
    filter(STRING_id==x)  %>% pull(Gene), collapse=";");
  gene_mapped %>% 
    filter(STRING_id %in% vals) %>% 
    select(Gene) %>%
    mutate(parent_gene=name_convert)})))


counts_gene_liver <- top_mods_liver %>% 
  pivot_longer(c(Gene, parent_gene), names_to="type", values_to="gene") %>%
  group_by(gene) %>%
  count() %>%
  arrange(desc(n))


# look at GO enrichment
counts_gene_liver %>% select(gene) %>% write_tsv("data/list_liver_genes.txt")


library('readxl')
s1 = read_excel("~/Downloads/pone.0023506.s006.xls", 1, skip=1)
symbols1 <- s1 %>% 
  filter(!is.na(`Gene Symbol`)) %>% 
  pull(`Gene Symbol`)

symbols2 <- s1 %>% 
  filter(!is.na(`Other gene symbol`)) %>% 
  separate_rows(`Other gene symbol`, sep="\\|") %>%
  pull(`Other gene symbol`)
de.genes <- union(symbols1, symbols2)
length(de.genes) # 1089

sum(counts_gene_liver$gene %in% de.genes) # 38/193 = 0.197
# what is the probably at random?
# 1089/19484 = 5.6%
# randomly select 193 out of 19484 genes
# select 38 genes out of a sample of 193 from urn containing 19848 with 1089 correct
dhyper(x=38, m=1089, k=193, n=17625-1089) # significant


setdiff(counts_gene_liver$gene, de.genes)
counts_gene_liver %>%
  mutate(rep=(gene %in% de.genes)) %>%
  View()

sum(sig_genes %>% pull(gene) %in% de.genes) # 249/2070 = 0.120
dhyper(x=249, k=2070, m=1089, n=19484-1089)


## plot
proteins <- names(V(chosen_liver$subnetwork))
prot_df <- string_db$add_proteins_description(data.frame("STRING_id"=proteins) )
createNetworkFromIgraph(chosen_liver$subnetwork,"my_network_liver")

colnames(prot_df)
loadTableData(prot_df %>% select(STRING_id, preferred_name) %>%
                rename(name=STRING_id,
                       node=preferred_name), data.key.column ="name" )
gene_mapped$node_weight <- 1- gene_mapped$Pvalue

loadTableData(gene_mapped %>% select(STRING_id, node_weight) %>%
                rename(name=STRING_id,
                       node_weight=node_weight), data.key.column ="name" )

string_db_instance <- STRINGdb$new()

# function for plotting modules
plot_module <- function(my_module, subnetwork,gene_mapped){
  my_mat <- as.matrix(subnetwork[my_module ,my_module ])
  prot_df <- string_db_instance$add_proteins_description(
    data.frame("STRING_id"=my_module) )
  prot_df %>% dplyr::select(STRING_id, preferred_name)
  rownames(my_mat) <- prot_df$preferred_name
  colnames(my_mat) <- prot_df$preferred_name
  g <- graph_from_adjacency_matrix(my_mat, weighted=T, mode="upper", add.rownames="code")
  
  list_p2 <- gene_mapped$Pvalue
  names(list_p2) <- gene_mapped$Gene
  vert_w <- -log(list_p2[names(V(g) )])
  plot(g, vertex.size=vert_w*3, directed=F, layout=layout_with_lgl, multiple=F, 
       vertex.label.color="black", vertex.color="lightblue")
}

plot_module(chosen_liver$modules[[1]],chosen_liver$subnetwork, gene_mapped)
