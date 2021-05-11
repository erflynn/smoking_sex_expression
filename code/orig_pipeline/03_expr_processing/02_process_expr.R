# 02_process_expr.R
# 11/22/2020
# Code for pre-processing the expression data
#  
# This code is specifically set up for Affy airway epithelium data
# which is GPL570 and GPL96. 
#
# It performs the following:
#  - deduplicate and identify duplicated columns for improving phenotyping
#  - run AffyImpute to move the GPL96 data to GPL570 space
#  - convert to genes
# 
# Note: this is run server-side.
# TODO:
#   - move saved objects to intermediate directory
#   - move convert to genes to utils, save probe_gene df

library(AffyImpute)
library(tidyverse)


# ------- 1. Run AFFYIMPUTE ------ #
# load the GPL96 data
load("data/eset_disc2.RData") # --> eset2, 22283

# note - the date needs to be RMA processed before running imputation
# (affyImpute can run this too but it breaks, also no pt)
hgu133plus2.matrix = affyImpute(exprs(eset2)) # 20089 x 77
# note this converts to ENTREZ genes
save(hgu133plus2.matrix, file="data/eset2_affyI.RData")

# remove duplicates
eset2_mat <- hgu133plus2.matrix[,!duplicated(t(hgu133plus2.matrix))] # 76

# put the data together, consider normalizing again?

# rownames correspond to ENTREZ IDs
load("data/eset_disc1.RData") # --> eset1 for GPL570 data
# 54675   796

# ----- 2. map plus2 data to genes ---- #
library(hgu133plus2.db)
x <- hgu133plus2ENTREZID
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])

probe_gene <- tibble(
  "probes"=names(xx),
  "gene" = unlist(xx)
)
# TODO save this object, slow


convertToGenes <- function(expData, key.df, gene.list){
  # TODO
  #  optimize, this is slow!! <-- possibly switch to MetaIntegrator function
  #  make sure the object contains expression data, keys
  #  change so that this is INPUT expr matrix + key mapping
  #  should work with both GEOQuery -AND- MetaIntegrator Objects
  
  list.keys <- key.df[key.df$gene %in% gene.list,]
  
  gene.to.probe <- split(key.df$probes,  key.df$gene) # this is slow... mb store for each platform
  expData2 <- do.call(cbind, lapply(1:length(gene.to.probe), function(x) {
    # get the gene and the probe
    g <- names(gene.to.probe)[x]
    p <- unlist(gene.to.probe[g])
    if (length(p)>1){
      expD <- expData[p,]
      df <- (data.frame(colMeans(expD, na.rm=TRUE)))
      return(df)
    }
    else {
      df <- data.frame(expData[p,])
      return(df)
    }})) ### ALSO SLOW...
  
  colnames(expData2) <- names(gene.to.probe)
  expData2.2 <- data.frame(t(expData2)) # columns are samples, rows are genes
  
  # create a data fram of NAs for missing genes
  missing.genes <- setdiff(gene.list, list.keys)
  missing.vec <- rep(NA, ncol(expData2.2))
  missing.df <- do.call(rbind, lapply(1:length(missing.genes), function(x) missing.vec))
  rownames(missing.df) <- missing.genes
  colnames(missing.df) <- colnames(expData2.2)
  
  # put together and reorder
  expDataPlusMiss <- rbind(expData2.2, missing.df )
  expData2.3 <- expDataPlusMiss[gene.list,] # REORDER so it matches other data
  
  return(expData2.3)
}

eset1_gene <- convertToGenes(exprs(eset1), probe_gene, rownames(eset2_mat))
save(eset1_gene, file="eset1_gene.RData") # 796
# save

# ---- 3. deduplicate and identify duplicates ---- #
eset1_mat <- eset1_gene[,!duplicated(t(eset1_gene))] # 595
stopifnot(rownames(eset1_mat)==rownames(eset2_mat))

# figure out which are duplicated...
eset1_gene2 <- data.frame(t(eset1_gene[1:2,order(eset1_gene[1,])]))
eset1_gene2$sample <- rownames(eset1_gene2)
dup_split <- eset1_gene2 %>% 
  as_tibble() %>% 
  group_by(`X3310`) %>% 
  summarize(n=n(), sample=paste(str_extract(sample, "^[0-9A-Za-z]+"), 
                                collapse=";")) %>%
  filter(n>=2) %>%
  arrange(desc(n), sample) %>%
  dplyr::select(-`X3310`)
  
dup_split %>% write_csv("data/duplicated_samples.csv")

# ----- 4. put data together and save ----- #
eset_comb <- cbind(eset1_mat, eset2_mat)
new_colnames <- str_extract(colnames(eset_comb), "^[0-9A-Za-z]+")
colnames(eset_comb) <- new_colnames
save(eset_comb, dup_split, file="data/eset_comb.RData")

