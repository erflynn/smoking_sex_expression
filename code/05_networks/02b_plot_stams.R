
library(STAMS)
library(tidyverse)
load("data/norm_res.RData")
normalized_results$ordered.module.score.matrix %>% filter(p==0) %>% nrow() # 422/12817

ae_p <- read_csv("data/ae_interaction_pvals_hgnc.csv")
ae_p2 <- ae_p$P.Value
names(ae_p2) <- ae_p$gene
chosen=moduleChoose(normalized_results, top=0.01, plot=FALSE) 
string_db_instance <- STRINGdb$new()#version="9_1")
#string_db_instance$plot_network(chosen$modules[[1]]) # does not work

# clean_my_mat <- function(mat){
#   my_df <- tibble("v1"="a", "v2"="b", "e"=0 )
#   for (i in 1:(ncol(mat)-1)){
#     for (j in (i+1):ncol(mat)){
#       if(mat[i,j]!=0){
#         v1 <- colnames(mat)[i]
#         v2 <- colnames(mat)[j]
#         verts <- sort(c(v1, v2))
#         my_df <- my_df %>% bind_rows(tibble("v1"=verts[[1]], "v2"=verts[[2]], 
#                                                "e"=mat[i,j]))
#       }
#     }
#   }
#   my_df %>% filter(e!=0)
# }

plot_module <- function(my_module){
  my_mat <- as.matrix(chosen$subnetwork[my_module ,my_module ])
  prot_df <- string_db_instance$add_proteins_description(
    data.frame("STRING_id"=my_module) )
  prot_df %>% select(STRING_id, preferred_name)
  rownames(my_mat) <- prot_df$preferred_name
  colnames(my_mat) <- prot_df$preferred_name
  g <- graph_from_adjacency_matrix(my_mat, weighted=T, mode="upper", add.rownames="code")
  vert_w <- -log(ae_p2[names(V(g) )])
  plot(g, vertex.size=vert_w*3, directed=F, layout=layout_with_lgl, multiple=F, 
       vertex.label.color="black", vertex.color="lightblue")
}
plot.new()
par(mfrow=c(1,2))
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

descript = string_db_instance$add_proteins_description(
  data.frame("STRING_id"=unique(unlist(chosen$modules[1:6]))))

