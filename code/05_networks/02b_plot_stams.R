# code for plotting STAMS output
#  uses the normalized result file + mapped gene pvals
library(STAMS)
library(tidyverse)

stams_dir <- "data/stams_blood_int"
map_file <- "data/gene_mapped_blood_int.RData"
load(sprintf("%s/norm_res.RData", stams_dir)) # --> normalized_results
normalized_results$ordered.module.score.matrix %>% filter(p==0) %>% nrow() # 422/12817
# 3199/7319
normalized_results$GWPI[[1]] # there are 10 here
# 6943/13779

load(map_file) # --> gene_mapped
#ae_p <- read_csv("data/ae_interaction_pvals_hgnc.csv")

ae_p2 <- gene_mapped$Pvalue
names(ae_p2) <- gene_mapped$Gene
#ae_p2 <- ae_p$P.Value
#names(ae_p2) <- ae_p$gene
chosen=moduleChoose(normalized_results, top=0.01, plot=FALSE) # --> 138
string_db_instance <- STRINGdb$new()

# function for plotting modules
plot_module <- function(my_module){
  my_mat <- as.matrix(chosen$subnetwork[my_module ,my_module ])
  prot_df <- string_db_instance$add_proteins_description(
    data.frame("STRING_id"=my_module) )
  prot_df %>% dplyr::select(STRING_id, preferred_name)
  rownames(my_mat) <- prot_df$preferred_name
  colnames(my_mat) <- prot_df$preferred_name
  g <- graph_from_adjacency_matrix(my_mat, weighted=T, mode="upper", add.rownames="code")
  vert_w <- -log(ae_p2[names(V(g) )])
  plot(g, vertex.size=vert_w*3, directed=F, layout=layout_with_lgl, multiple=F, 
       vertex.label.color="black", vertex.color="lightblue")
}

plot_subnetwork <- function(){
  proteins <- names(V(chosen$subnetwork))
  my_mat <- as.matrix(chosen$subnetwork[proteins, proteins])
  prot_df <- string_db_instance$add_proteins_description(
    data.frame("STRING_id"=proteins) )
  prot_df %>% dplyr::select(STRING_id, preferred_name)
  rownames(my_mat) <- prot_df$preferred_name
  colnames(my_mat) <- prot_df$preferred_name
  g <- graph_from_adjacency_matrix(my_mat, weighted=T, mode="upper", add.rownames="code")
  vert_w <- -log(ae_p2[names(V(g) )])
  plot(g, vertex.size=vert_w*3, directed=F, layout=layout_with_lgl, multiple=F, 
       vertex.label.color="black", vertex.color="lightblue")
}
plot_subnetwork(1)

plot.new()
par(mfrow=c(1,1))
plot_module(chosen$modules[[1]])

plot_module(chosen$modules[[2]])

plot.new()
par(mfrow=c(1,2))
plot_module(chosen$modules[[3]])
plot_module(chosen$modules[[4]])

plot.new()
par(mfrow=c(1,2))
plot_module(chosen$modules[[5]])
plot_module(chosen$modules[[6]])

plot.new()
par(mfrow=c(1,2))
plot_module(chosen$modules[[7]])
plot_module(chosen$modules[[8]])

#descript = string_db_instance$add_proteins_description(
#  data.frame("STRING_id"=unique(unlist(chosen$modules[1:6]))))


stams_dir <- "data/stams_run4"
load(sprintf("%s/norm_res.RData", stams_dir)) # --> normalized_results
chosen=moduleChoose(normalized_results, top=0.01, plot=FALSE) # --> 138

library(RCy3)
# source	target	interaction	directed	symbol	value
cytoscapePing()
proteins <- names(V(chosen$subnetwork))
prot_df <- string_db_instance$add_proteins_description(data.frame("STRING_id"=proteins) )
createNetworkFromIgraph(chosen$subnetwork,"my_network_ae")

colnames(prot_df)
loadTableData(prot_df %>% select(STRING_id, preferred_name) %>%
                rename(name=STRING_id,
                       node=preferred_name), data.key.column ="name" )

plot_module(chosen$modules[[5]])

