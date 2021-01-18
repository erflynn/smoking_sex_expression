# TODO: figure out if we need an upper or lower bound

library(STAMS)
library(tidyverse)

# ---- map data to STRINGdb ---- #
## LOAD PVALUES
string_db <- STRINGdb$new() 
gene_values=read_csv("ae_interaction_pvals.csv")
gene_values <- gene_values %>% rename(Pvalue=P.Value, Gene=entrezgene_id)

# get rid of pvalue = 0 lines and replace with upper bound 
#gene_values$Pvalue=pmax(gene_values$Pvalue, 1/gene_values$nSims)
# skip this - nSims is from VEGAS

# remove pvalue = 1 - there are none here? 
#print("warning, removing genes with a pvalue = 1") 
print(length(which(gene_values$Pvalue ==1))) 

gene_values2=gene_values %>% 
  filter(Pvalue != 1) %>% 
  filter(!is.na(Gene)) #%>%
  #group_by(Gene) %>%
  #mutate(Pvalue=min(Pvalue)) %>%
  #ungroup()
stopifnot(length(unique(gene_values$Gene))==nrow(gene_values))
#STRINGdb$help("map") # -- to get info about STRING functions
gene_mapped = string_db$map(data.frame(gene_values2), "Gene", removeUnmappedRows = TRUE )
# for HGNC symbol: 11% did not map (out of 17966 )
# for Entrez: 13% did not map (out of 19236) -- this is ~700 more present tho
# // TODO: improve this!

save(gene_mapped, file="data/mapped_ae.RData")

# check mapping 
#  TODO: look more into this for ones that arent HGNC??
data3=string_db$add_proteins_description(gene_mapped) 
data4=data3[which(data3$Gene==toupper(data3$preferred_name)),] 
ensg_genes=data3[grep('ENSG',data3$preferred_name),] 
data4=rbind(data4,ensg_genes) 
data_mapped=data4[,-which(names(data4)=='annotation')]



# ---- example ---- #
mapped_data <- gene_mapped
# Already done in example data: map input P-value to weight:
mapped_data$weight=qnorm(1-mapped_data$Pvalue)
#data("mapped_data")
string_db_instance <- STRINGdb$new()#version="9_1")
mapped_data$weight=qnorm(1-mapped_data$Pvalue)
mapped_data <- mapped_data %>% head(200)
G=get_string_edges(mapped_data, "combined_score", string_db_instance)

# to use the default method to estimate lambda, run this line. Otherwise, set lambda to a value between 0 and 1
lambda=est_lambda(G, nperm=10) # for a real experiment, set nperm to 10,000 +

# TODO: save here!

# set up a cluster to run the searches in parallel. Each node ndeeds to know these four variables, and have STAMS and igraph loaded.
num_cores=detectCores()

r=0.1; search_r=1; min_size=5

cl <- makeCluster(num_cores)
clusterEvalQ(cl, library(igraph))
clusterEvalQ(cl, library(STAMS))
clusterExport(cl, c("G", "r", "lambda", "search_r"))


search_results=stams_search(cl, G,  min_size=5)
stopCluster(cl)

# TODO: save here!

# For a real experiment, increase the number of permutations here to 100,000 +
normalized_results=normalize_scores(G, search_results, 
                                    filename='normalized_scores.txt', perm=100)

# TODO: save here!

chosen=moduleChoose(normalized_results, top=0.01, plot=FALSE) 
# This method is from dmGWAS, the plot uses tcl/tk. To plot the results this way, set plot=TRUE



png=string_db_instance$get_png(chosen$modules[[1]], file="this.png") 
# This plot doesn't need tcl/tk or X
# string_db_instance$plot_network(chosen$modules[[1]], add_summary=FALSE) # This plot uses X

