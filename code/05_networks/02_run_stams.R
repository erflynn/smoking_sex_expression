library(STAMS)
library(miceadds)
load.Rdata("data/mapped_ae.RData", "mapped_data")
mapped_data$weight=qnorm(1-mapped_data$Pvalue)
#data("mapped_data")
string_db_instance <- STRINGdb$new()#version="9_1")
mapped_data$weight=qnorm(1-mapped_data$Pvalue)
#mapped_data <- mapped_data %>% head(200)
G=get_string_edges(mapped_data, "combined_score", string_db_instance)
print("Edges obtained")

# to use the default method to estimate lambda, run this line. Otherwise, set lambda to a value between 0 and 1
lambda=est_lambda(G, nperm=100) # for a real experiment, set nperm to 10,000 +

print("Lambda estimated")
save(lambda, file="data/stams_run/lambda.RData")

# TODO: save here!

# set up a cluster to run the searches in parallel. 
# Each node ndeeds to know these four variables, and have STAMS and igraph loaded.
num_cores=detectCores()
r=0.1; search_r=1; min_size=5

cl <- makeCluster(num_cores)
clusterEvalQ(cl, library(igraph))
clusterEvalQ(cl, library(STAMS))
clusterExport(cl, c("G", "r", "lambda", "search_r"))

search_results=stams_search(cl, G,  min_size=5)
stopCluster(cl)

print("Search complete")
# TODO: save here!
save(search_res, file="data/stams_run/search_res.RData")

# For a real experiment, increase the number of permutations here to 100,000 +
normalized_results=normalize_scores(G, search_results, 
                                    filename='normalized_scores.txt', perm=100)

print("Results normalized")
# TODO: save here!
save(normalized_results, file="data/stams_run/norm_res.RData")


chosen=moduleChoose(normalized_results, top=0.01, plot=FALSE) 
png=string_db_instance$get_png(chosen$modules[[1]], file="test.png") 

