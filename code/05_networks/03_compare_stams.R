# compare STAMS results from a normal-length vs long run
# ... using blood data

# do we need a minimum P-val? PROBABLY

library(STAMS)
library(tidyverse)
# 100,000
load("data/stams_blood_int/norm_res.RData") # --> normalized_results
norm_res1 <- normalized_results
chosen1=moduleChoose(norm_res1, top=0.01, plot=FALSE) # --> 73
length(chosen1[[1]]) # --> 73
norm_res1$ordered.module.score.matrix %>% filter(p==0) %>% nrow() 
# 3199/7319



# more permutations: 1 mil
load("data/stams_blood_int_long/norm_res.RData") # --> normalized_results
norm_res2 <- normalized_results
chosen2=moduleChoose(norm_res2, top=0.01, plot=FALSE) 
length(chosen2[[1]]) # --> 73
norm_res2$ordered.module.score.matrix %>% filter(p==0) %>% nrow() # 3211
nrow(norm_res2$ordered.module.score.matrix) # 7319
chosen2$modules[[4]]
chosen1$modules[[4]]


# what is the distribution of input p-values?
# ... perhaps it needs to be closer to uniform?
library(miceadds)
load.Rdata("data/gene_mapped_blood_int.RData", "mapped_data_int")
plot(density(mapped_data_int$Pvalue)) # approx uniform?
load.Rdata("data/gene_mapped_blood_sex.RData", "mapped_data_sex")
plot(density(mapped_data_sex$Pvalue)) # not really uniform
mapped_data_sex %>% filter(chromosome=="X") %>% nrow() # 316/8412 are X chrom
# NOTE: may need to set the 0 p-value to an actual number, this is probably a problem

load.Rdata("data/gene_mapped_blood_smok.RData", "mapped_data_smok")
plot(density(mapped_data_smok$Pvalue)) # not really uniform
load.Rdata("data/stams_blood_smok/norm_res.RData", "norm_res_smok") 
norm_res_smok$ordered.module.score.matrix %>% filter(p==0) %>% nrow() # 676 -- this seems much more reasonable...
norm_res_smok$ordered.module.score.matrix  %>% nrow() # 6886
chosen_smok=moduleChoose(norm_res_smok, top=0.01, plot=FALSE) # --> 69
length(chosen_smok$modules)

head(chosen_smok$modules[[1]])


load.Rdata("data/stams_blood_sex/norm_res.RData", "norm_res_sex") 
norm_res_sex$ordered.module.score.matrix %>% filter(p==0) %>% nrow() # 385 -- this seems much more reasonable...
norm_res_sex$ordered.module.score.matrix  %>% nrow() # 3160
chosen_sex_c=moduleChoose(norm_res_sex, top=0.01, plot=FALSE) # --> 32
length(chosen_sex_c$modules)

top_mods_sex <- data.frame(do.call(rbind, lapply(names(chosen_sex_c$modules), function(x) 
  {vals <- chosen_sex_c$modules[[x]];
  name_convert <-  mapped_data_sex %>% 
    filter(STRING_id==x)  %>% pull(Gene);
  mapped_data_sex %>% 
    filter(STRING_id %in% vals) %>% 
    select(Gene, chromosome) %>%
    mutate(parent_gene=name_convert) %>%
    left_join(mapped_data_sex %>% select(Gene, chromosome) %>% rename(parent_chr=chromosome), by=c("parent_gene"="Gene"))})))
top_mods_sex %>% group_by(parent_gene) %>% count() %>% arrange(desc(n))
top_mods_sex %>% group_by(Gene) %>% count() %>% arrange(desc(n))

top_mods_sex %>% group_by(parent_chr) %>% count() %>% arrange(desc(n)) 
top_mods_sex %>% group_by(chromosome) %>% count() %>% arrange(desc(n)) 
top_mods_sex %>% filter(parent_chr=="X" | chromosome=="X") # only 3/190

chosen_sex_c$modules[[1]]
head(mapped_data_sex)
#### 
# what fraction of X chromosome genes are included?


# results are the same!
