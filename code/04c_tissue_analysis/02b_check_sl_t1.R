# GSE31908, GSE40013, GSE17913, GSE2125, GSE36868 (no)

gse1 <- getGEO("GSE31908") # 20 missing sex labels
dim(exprs(gse1$`GSE31908-GPL570_series_matrix.txt.gz`)) # 17 
dim(exprs(gse1$`GSE31908-GPL96_series_matrix.txt.gz`)) # 50
dim(exprs(gse1$`GSE31908-GPL97_series_matrix.txt.gz`)) # 50
gse31908 # ... 40 ...
ggplot(sl %>% filter(study_acc=="GSE31908"),
       aes(x=p_male))+
  geom_histogram()

# gse2 <- getGEO("GSE40013") # sex labels are a mess
# exp2 <- exprs(gse2$GSE40013_series_matrix.txt.gz)
# unique(pData(gse2$GSE40013_series_matrix.txt.gz)$platform) # GPL6244
# f_df <- fData(gse2$GSE40013_series_matrix.txt.gz)
# massir_gpl6244 <- f_df %>% filter(str_detect(mrna_assignment, "chrY")) %>% pull(ID)
# 
# xist_gpl6244 <- f_df %>% filter(str_detect(gene_assignment, "XIST")) %>% pull(ID)
# rps_gpl6244 <- f_df %>% filter(str_detect(gene_assignment, "RPS4Y1")) %>% pull(ID)
# kdm_gpl6244 <- f_df %>% filter(str_detect(gene_assignment, "KDM5D")) %>% pull(ID)
# plot(exp2[as.character(xist_gpl6244),], 
#      exp2[as.character(rps_gpl6244),])
# #massir_lab2 <- massiRAcc(exp2, massir_gpl6244, plot=T)
# # this is a mess
# ggplot(sl %>% filter(study_acc=="GSE40013"),
#        aes(x=p_male))+
#   geom_histogram()

gse3 <- getGEO("GSE17913")
exp3 <- exprs(gse3$GSE17913_series_matrix.txt.gz)
unique(pData(gse3$GSE17913_series_matrix.txt.gz)$platform) # GPL570
f_df <- fData(gse3$GSE17913_series_matrix.txt.gz)
xist_gpl570 <- f_df %>% filter(str_detect(`Gene Symbol`, "XIST")) %>% pull(ID)
rps_gpl570 <- f_df %>% filter(str_detect(`Gene Symbol`, "RPS4Y1")) %>% pull(ID)
kdm_gpl570 <- f_df %>% filter(str_detect(`Gene Symbol`, "KDM5D")) %>% pull(ID)

ts <- tokerSexLab(exp3, f.genes =xist_gpl570, m.genes=c(rps_gpl570, kdm_gpl570))
my_s <- ts[gse17913$sample_acc]

plot(colMeans(exp3[xist_gpl570,]), colMeans(exp3[c(rps_gpl570, kdm_gpl570),]))

table(gse17913$sex_lab==my_s) # two of ours do not
table(gse17913$gender==my_s) # 7 of theirs do not


View(head(f_df))
gse4 <- getGEO("GSE2125")
unique(pData(gse4$GSE2125_series_matrix.txt.gz)$platform) # GPL570
exp4 <- exprs(gse4$GSE2125_series_matrix.txt.gz)

ts4 <- tokerSexLab(exp4, f.genes =xist_gpl570, m.genes=c(rps_gpl570, kdm_gpl570))
my_s4 <- ts4[gse2125$sample_acc]

plot(colMeans(exp4[xist_gpl570,]), colMeans(exp4[c(rps_gpl570, kdm_gpl570),]))
my_s4[is.na(gse2125$sex_lab)]
table(gse2125$sex_lab==my_s4) # all match
table(gse17913$gender==my_s) # 7 of theirs do not


# gse5 <- getGEO("GSE36868")
# exp5 <- exprs(gse5$GSE36868_series_matrix.txt.gz)
# dim(exp5) # 960
# str(gse5$GSE36868_series_matrix.txt.gz, 2) # GPL6883
# unique(pData(gse5$GSE36868_series_matrix.txt.gz)$platform)
# f_df <- fData(gse5$GSE36868_series_matrix.txt.gz)
# 
# f_df %>%filter(str_detect(Symbol, "RPS4Y1"))
# f_df %>% filter(str_detect(Synonyms, "XIST")) # appears to have no XIST transcripts...
# f_df %>% filter(Chromosome=="X") %>% pull(Symbol)
# massir_genes <- f_df %>% filter(Chromosome=="Y") %>% pull(ID)
# massir_lab <- massiRAcc(exp5, massir_genes, plot=T)
# ... very unclear separation


