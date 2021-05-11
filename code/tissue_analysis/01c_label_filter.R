library(tidyverse)
library(GEOmetadb)

######## NOW WHAT DOES THIS LOOK LIKE? ########
# 27 studies!!!
new_kept <- toupper(str_replace_all(list.files("data/pdata_filt2/"), ".csv", ""))

kept_r2 <- post_annot %>% filter(study_acc %in% new_kept) %>%
  select(study_acc, title, description, num_samples, tissue2, n, X13)
gse_kept <- kept_r2 %>% filter(str_detect(study_acc, "GSE")) %>% pull(study_acc)
paste(gse_kept , collapse='","') # 23 are GSEs
# GSE20250 failed
kept_r2 %>% filter(!str_detect(study_acc, "GSE")) %>% pull(study_acc) # 5


# steps:
# 0) get more metadata
#   - platform
con <- dbConnect(SQLite(), "../GEOmetadb.sqlite") # 11/8/2020

gse_gpl <- dbGetQuery(con, sprintf("SELECT gse, gpl FROM gse_gpl WHERE gse IN ('%s');",
                                   paste()))#, collapse="','")))
gse_gpl2 <- gse_gpl %>%
  filter(gpl!="GPL6804") # skip SNP array
list_gpl <- gse_gpl2 %>% distinct(gpl) %>% pull(gpl)
gpl_info <- dbGetQuery(con, sprintf("SELECT * FROM gpl WHERE gpl IN ('%s');",
                                    paste(list_gpl, collapse="','")))

# GPL21290, GPL6804
array_info <- gpl_info %>% 
  select(gpl, technology, manufacturer, data_row_count,title)

dl_info <- gpl_info %>% select(gpl, title, manufacturer, bioc_package, supplementary_file, web_link) 

gse_gpl2 %>% 
  group_by(gpl) %>%
  count() %>% arrange(desc(n)) %>%
  left_join(array_info) %>%
  select(-title)

# 1) try to download
# -- FAILED b/c RNA-seq -- #
# GSE101353 - SRP111742
# [GSE111819 - SRP135671  [possibly cut, diverticulitis], on RB]
# [GSE131391 -  SRP198757 [possibly cut, single cell]]
# GSE134174 - SRP214324, available on ds
# GSE134692 -  available on the ds

# GSE102556 -  	SRP115956, skip mouse (GPL13112)  available on ds
# GSE110907 -   	SRP133217, [is this miRNA only?] (looks like RAW avail?)     
# GSE47718 -  	SRP024274 (looks like RAW avail?) , on RB, on recount2
# GSE136262 -   	SRP219430, available on ds
# GSE155588 -  	SRP275616, available on ds

list_dl <- str_replace_all(list.files("data/small_ds2/"), ".RData", "")


# --- GOOD ---- #
# 12 

to_read <- c("GSE32504", "GSE15289", "GSE20681", "GSE30272", 
             "GSE14633","GSE16008","GSE27002", "GSE23323",
             "GSE23515", "GSE47415", "GSE19738", "GSE5056" )
gse_gpl %>% filter(gse %in% to_read) %>% left_join(array_info)
my_gses <- lapply(to_read, function(x){
  load(sprintf("data/small_ds2/%s.RData", x)); return(gse)})
lapply(my_gses, function(x) length(x$originalData)) # all 1
lapply(my_gses, function(x) x$originalData[[1]]$exp_comment)
lapply(my_gses, function(x) x$originalData[[1]]$key_comment)


# --- ArrayExpress ones --- #
unique(ae1$`Array Design REF`) # A-MEXP-2072 - Illumina HumanHT-12_V4_0_R2_15002873_B
unique(ae2$`Array Design REF`) # A-MEXP-691 - Illumina Human-6 v1 Expression BeadChip
unique(ae3$`Array Design REF`) # "A-AFFY-44", GPL570
unique(ae4$`Array Design REF`) # "A-AFFY-44", GPL570
unique(ae5$`Array Design REF`) # "A-AFFY-44", GPL570
unique(ae6$`Array Design REF`) # A-AGIL-9 - Agilent Human 1A Microarray (V2)

# get list of files to download
kept_r2 %>% filter(!str_detect(study_acc, "GSE")) %>% pull(study_acc) # 5

#"E-MTAB-6667" - need to download from raw, Illumina HumanHT
#"E-MTAB-5278" - processed available, GPL570
#"E-MTAB-3604" - processed available, GPL570
#"E-MTAB-5279" - processed available, GPL570
#"E-MEXP-1277" - processed available - agilent



######### 2) SEX LABEL ############
# get XY chromosome genes
mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
chr_map <- biomaRt::getBM(attributes=c("hgnc_symbol", "chromosome_name"), mart=mart)
xy_chr_map <- chr_map %>% filter(chromosome_name %in% c("X", "Y"))
save(xy_chr_map, file="ref/xy_chr_map.RData")
# try labeling w both

library('MetaIntegrator')

# ---- GSE32504, liver
my_gse <- getGEOData("GSE32504")
my_gse <- my_gse$originalData$GSE32504
save(my_gse, file="data/small_gses2/GSE32504.RData")
pdat <- read_csv("data/pdata_filt2/gse32504.csv")

map_to_title <- my_gse$pheno %>% select(title, geo_accession)
meta_ds <- read_excel("ref/gse32504_meta.xls", skip=4)
meta_ds %>% filter(!is.na(`Sample ID`)) %>% 
  rename(title=`Sample ID`) %>%
  inner_join(map_to_title) %>%
  inner_join(pdat, by=c("geo_accession"="gsm")) %>%
  select(geo_accession, title, Sex, gender, Smoking, smok) %>%
  group_by(Sex, gender, Smoking, smok) %>%
  count()
# looks like they swapped M and F

#my_gse <- my_gses[[1]][[1]][[1]]
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
# appears to be probe --> one gene
probe_gene %>% filter(str_detect(gene, ","))
probe_gene %>% filter(str_detect(gene, ";"))

xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 

sl <- massiRAcc(my_gse$expr, xy_genes %>% 
                  filter(chromosome_name=="Y") %>% 
                  pull(probe), plot=T)

# woot - good separation
sl_df <- tibble("gsm"=names(sl),
                "expr_sex"=unlist(sl))

# uhhh everything is wrong except 2??? wtf?
# ... but the paper says 79 females, 71 males --> keep ours
pdat2 <- pdat %>% 
  mutate(gender_fixed=ifelse(gender=="male", "female", "male")) %>%
  left_join(sl_df) %>% 
  filter(gender_fixed==expr_sex)  %>%
  rename(tissue=source_name_ch1,
         sex=expr_sex,
         metadata_sex=gender_fixed) %>%
  mutate(race_ethnicity="white") %>%
  select(-`ethnic background`, -gender) 


table(pdat$gender)
table(sl_df$expr_sex)

plot(t(my_gse$expr[probe_gene %>% 
                     filter(gene %in% c("XIST", "RPS4Y1")) %>%
                     pull(probe),] ), col=ifelse(sl=="male", "blue", "red"))
# do we trust the smoking data tho?
# ... check -- it checks out in other ways?

pdat2 %>% write_csv("data/pdata_filt2_v2/gse32504.csv")

#  "GSE15289" -- possibly skip
# TODO: check annot
my_gse <- getGEOData("GSE15289")
my_gse <- my_gse$originalData$GSE15289
my_gse$platform
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
# appears to be probe --> one gene
probe_gene %>% filter(str_detect(gene, ","))
probe_gene %>% filter(str_detect(gene, ";"))

xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 

sl <- massiRAcc(my_gse$expr, xy_genes %>% 
                  filter(chromosome_name=="Y") %>% 
                  pull(probe),  plot=T)

# not good separation, try toker
probe_gene %>% filter(!is.na(gene)) # only 7.6k genes... hmm


# "GSE20681"
gse <- getGEOData("GSE20681")
my_gse <- gse$originalData$GSE20681
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 

sl <- massiRAcc(my_gse$expr, xy_genes %>% 
                  filter(chromosome_name=="Y") %>% 
                  pull(probe),  plot=T)
# decent separation
sl_df <- tibble("gsm"=names(sl),
                "expr_sex"=unlist(sl))

pdat <- read_csv("data/pdata_filt2/gse20681.csv")
save(my_gse, file="data/small_gses2/GSE20681.RData")
# all match!
pdat %>% left_join(sl_df) %>% filter(sex==expr_sex)
pdat %>% 
  left_join(sl_df) %>% 
  rename(metadata_sex=sex, sex=expr_sex) %>%
  write_csv("data/pdata_filt2_v2/gse20681.csv")

# "GSE30272"
gse <- getGEOData("GSE30272")
my_gse <- gse$originalData$GSE30272
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 

sl <- massiRAcc(my_gse$expr, xy_genes %>% 
                  filter(chromosome_name=="Y") %>% 
                  pull(probe),  plot=T) # ok sep
pdat <- read_csv("data/pdata_filt2/gse30272.csv")
save(my_gse, file="data/small_gses2/GSE30272.RData")
sl_df <- tibble("gsm"=names(sl),
                "expr_sex"=unlist(sl))
pdat2 <- pdat %>% 
  mutate(metadata_sex=case_when(
    sex=="m" ~ "male",
    sex=="f" ~ "female"
  )) %>%
  left_join(sl_df) %>%
  select(-sex) %>%
  rename(sex=expr_sex, tissue=tisuse)  %>%
  mutate(race_ethnicity=case_when(
    race=="aa" ~ "black",
    race=="as" ~ "asian",
    race=="hisp" ~ "hispanic",
    race=="cauc" ~ "white")) %>%
  select(-race) %>%
  select(gse, gsm, smok, sex, metadata_sex, age, race_ethnicity, everything())

pdat2 %>%  
  select(gsm, metadata_sex, expr_sex) %>%
  filter(metadata_sex!=expr_sex) # all match
pdat2 %>% write_csv("data/pdata_filt2_v2/gse30272.csv")

# "GSE14633"
gse <- getGEOData("GSE14633")
my_gse <- gse$originalData$GSE14633
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 

sl <- massiRAcc(my_gse$expr, xy_genes %>% 
                  filter(chromosome_name=="Y") %>% 
                  pull(probe),  plot=T)  # good sep!
sl_df <- tibble("gsm"=names(sl),
                "expr_sex"=unlist(sl))
pdat <- read_csv("data/pdata_filt2/gse14633.csv")
save(my_gse, file="data/small_gses2/GSE14633.RData")

pdat %>% left_join(sl_df) %>% filter(sex!=expr_sex) # all match
pdat2 <- pdat %>% 
  left_join(sl_df) %>%
  rename(metadata_sex=sex,
         sex=expr_sex) 
pdat2 %>% write_csv("data/pdata_filt2_v2/gse14633.csv")


# "GSE16008" - bronchial and nasal epithelium
# GPL5175, Exon Array data
gse <- getGEOData("GSE16008")
my_gse <- gse$originalData$GSE16008
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 


pdat <- read_csv("data/pdata_filt2/gse16008.csv")
pdat %>% group_by(tissue) %>% count()
nasal_samples <- pdat %>% filter(tissue=="nasal epithelium") %>% pull(gsm)
sl_n <- massiRAcc(my_gse$expr[,nasal_samples], xy_genes %>% 
                    filter(chromosome_name=="Y") %>% 
                    pull(probe),  plot=T, threshold=4)  # looks like there are batch effects?

bronchial_samples <- pdat %>% filter(tissue=="bronchial epithelium") %>% pull(gsm)
sl_b <- massiRAcc(my_gse$expr[,bronchial_samples], xy_genes %>% 
                    filter(chromosome_name=="Y") %>% 
                    pull(probe),  plot=T, threshold=4)  # good sep

sl_df <- tibble("gsm"=names(sl_n),
                "expr_sex"=unlist(sl_n)) %>%
  bind_rows(tibble("gsm"=names(sl_b),
                   "expr_sex"=unlist(sl_b)))

# counts match for BE (minus 1 white f NS)
pdat %>% left_join(sl_df) %>% 
  filter(tissue=="bronchial epithelium")  %>%
  filter(sex!=expr_sex) # 8/26
lapply(pdat %>% left_join(sl_df) %>% 
         filter(tissue=="bronchial epithelium") %>%
         group_split(smok), function(x) fct_summ(x))

pdat %>% left_join(sl_df) %>% 
  filter(tissue=="nasal epithelium")  %>%
  filter(sex!=expr_sex) # 16/37

lapply(pdat %>% left_join(sl_df) %>% 
         filter(tissue=="nasal epithelium") %>%
         group_split(smok), function(x) fct_summ(x))


# some of these are from previous studies??

save(my_gse, file="data/small_gses2/GSE16008.RData")



# "GSE27002" 
gse <- getGEOData("GSE27002")
my_gse <- gse$originalData$GSE27002
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 

sl <- massiRAcc(my_gse$expr, xy_genes %>% 
                  filter(chromosome_name=="Y") %>% 
                  pull(probe),  plot=T)  # good sep
sl_df <- tibble("gsm"=names(sl),
                "expr_sex"=unlist(sl))
pdat <- read_csv("data/pdata_filt2/gse27002.csv")
save(my_gse, file="data/small_gses2/GSE27002.RData")
# no sex information tho
pdat %>% left_join(sl_df) %>% rename(sex=expr_sex) %>%
  write_csv("data/pdata_filt2_v2/gse27002.csv")


# "GSE23323",
gse <- getGEOData("GSE23323")
my_gse <- gse$originalData$GSE23323
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 

y_probes <- xy_genes %>% 
  filter(chromosome_name=="Y") %>% 
  pull(probe)

count_nas <- apply(my_gse$expr[y_probes,],1, function(x) sum(is.na(x)))
y_probes_present <- names(count_nas[count_nas==0])

sl <- massiRAcc(my_gse$expr, y_probes_present,  plot=T)  # good sep
sl_df <- tibble("gsm"=names(sl),
                "expr_sex"=unlist(sl))

pdat <- read_csv("data/pdata_filt2/gse23323.csv")
save(my_gse, file="data/small_gses2/GSE23323.RData")
# no sex information tho
pdat %>% left_join(sl_df) %>% rename(sex=expr_sex) %>%
  write_csv("data/pdata_filt2_v2/gse23323.csv")

# "GSE23515"
# Paul and Amundson
gse <- getGEOData("GSE23515")
my_gse <- gse$originalData$GSE23515
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 


y_probes <- xy_genes %>% 
  filter(chromosome_name=="Y") %>% 
  pull(probe)

count_nas <- apply(my_gse$expr[y_probes,],1, function(x) sum(is.na(x)))
y_probes_present <- names(count_nas[count_nas==0])


sl <- massiRAcc(my_gse$expr, y_probes_present,  plot=T)  # good sep
sl_df <- tibble("gsm"=names(sl),
                "expr_sex"=unlist(sl))

pdat <- read_csv("data/pdata_filt2/gse23515.csv")
save(my_gse, file="data/small_gses2/GSE23515.RData")
pdat %>% 
  left_join(sl_df) %>%
  filter(gender!=expr_sex) # all match
pdat %>% left_join(sl_df) %>%
  rename(metadata_sex=gender,
         sex=expr_sex) %>%
  write_csv("data/pdata_filt2_v2/gse23515.csv")

# "GSE47415" - I think this might be a duplicate
# Paul and Amundson
# gse <- getGEOData("GSE47415")
# my_gse <- gse$originalData$GSE47415
# probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
# xy_genes <- probe_gene %>% 
#   inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 
# 
# 
# y_probes <- xy_genes %>% 
#   filter(chromosome_name=="Y") %>% 
#   pull(probe)
# 
# count_nas <- apply(my_gse$expr[y_probes,],1, function(x) sum(is.na(x)))
# y_probes_present <- names(count_nas[count_nas==0])
# 
# 
# sl <- massiRAcc(my_gse$expr, y_probes_present,  plot=T)  # good sep
# sl_df <- tibble("gsm"=names(sl),
#                 "expr_sex"=unlist(sl))
# 
# pdat <- read_csv("data/pdata_filt2/gse47415.csv")
# 
# save(my_gse, file="data/small_gses2/GSE47415.RData")
# pdat %>% left_join(sl_df) %>%
#   rename(metadata_sex=gender,
#          sex=expr_sex) %>%
#   fct_summ()



# "GSE19738"
gse <- getGEOData("GSE19738")
my_gse <- gse$originalData$GSE19738
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 

pdat <- read_csv("data/pdata_filt2/gse19738.csv")
pdat %>% fct_summ()
y_probes <- xy_genes %>% 
  filter(chromosome_name=="Y") %>% 
  pull(probe)

expr2 <- my_gse$expr[,pdat$gsm]

count_nas <- apply(expr2[y_probes,],1, function(x) sum(is.na(x)))
y_probes_present <- names(count_nas[count_nas==0])
y_probes_present

sl <- massiRAcc(expr2, y_probes_present,  plot=T, threshold=4)  
# not great sep


save(my_gse, file="data/small_gses2/GSE19738.RData")

sl_df <- tibble("gsm"=names(sl),
                "expr_sex"=unlist(sl))

unlist(my_gse$pheno %>%
         filter(`treatment:ch1`=="none",
                `sample type:ch1`=="control") %>%
         filter(!str_detect( characteristics_ch1.3, "quit")) %>%
         mutate(title=as.character(title)) %>% pull(title) )

#fct_summ() # 20f, 14 m
pdat %>% fct_summ() # study says 23f, 12 m

# we have 16 and 9
pdat %>% 
  left_join(sl_df) %>%
  filter(gender!=expr_sex) # 11/25 do not match

#  write_csv("data/pdata_filt2_v2/gse19738.csv")

# "GSE5056" - note 2 per each
gse <- getGEOData("GSE5056")
my_gse <- gse$originalData$GSE5056
probe_gene <- tibble("probe"=names(my_gse$keys), "gene"=unlist(my_gse$keys))
xy_genes <- probe_gene %>% 
  inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 


y_probes <- xy_genes %>% 
  filter(chromosome_name=="Y") %>% 
  pull(probe)

count_nas <- apply(my_gse$expr[y_probes,],1, function(x) sum(is.na(x)))
y_probes_present <- names(count_nas[count_nas==0])

sl <- massiRAcc(my_gse$expr, y_probes_present,  plot=T, threshold=2)  
# ok sep at 2
sl_df <- tibble("gsm"=names(sl),
                "expr_sex"=unlist(sl))
pdat <- read_csv("data/pdata_filt2/gse5056.csv")
save(my_gse, file="data/small_gses2/GSE5056.RData")
pdat %>% 
  left_join(sl_df) %>%
  filter(sex!=expr_sex) # all match
pdat2 <- pdat %>% left_join(sl_df) %>%
  rename(metadata_sex=sex,
         sex=expr_sex,
         race_ethnicity=`ethnic group`)  
pdat2 %>% fct_summ()
pdat2 %>% write_csv("data/pdata_filt2_v2/gse5056.csv")

# download the array express ones
#"E-MTAB-6667" - need to download from raw, Illumina HumanHT
res <- sapply(1:16, function(i)
  sprintf("wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6667/E-MTAB-6667.raw.%s.zip", i))

res %>% write_tsv("list_cmd.txt")
#"E-MTAB-3604" - processed available, GPL570
# added a comma at the start
exp <- read_csv("data/array_exp/E-MTAB-3604/COPD_SPUTUM_GCRMANormalizedMatrix.txt")
# check parsing fails?
pdat <- read_csv("data/pdata_filt2/e-mtab-3604.csv")
full_p <- read_tsv("data/array_exp/E-MTAB-3604.sdrf.txt") 
pdat2 <- pdat %>% left_join(full_p %>% select(`Source Name`, `Array Data File`))

exp_sm <- exp[,c("X1", pdat2$`Array Data File`)] %>% rename(probes=X1)
head(exp_sm[,1:5])
exp_mat <- as.matrix(exp_sm %>% select(-probes))
rownames(exp_mat) <- exp_sm$probes

.GEM_log_check_v1 <- function(GEM){
  
  #pick the first column that is not all NAs
  ref_p <- which(sapply(1:ncol(GEM$expr),
                        function(i)
                          !all(is.na(GEM$expr[,i]))))[1]
  
  #this range is obserbed only in log-normal case [in this part context]
  if(range(GEM$expr[,ref_p],na.rm=T)[1]<0){
    return(TRUE)
  }
  #if it's all positives all the bets are off
  else{
    #The data should be normally distributed
    #once in log scale
    #Get qq-norm objects
    obj_real <- qqnorm(GEM$expr[,ref_p],plot.it=FALSE);
    obj_exp  <- qqnorm(exp(GEM$expr[,ref_p]),plot.it=FALSE);
    obj_log  <- qqnorm(log(abs(GEM$expr[,ref_p])+0.00001),plot.it=FALSE);
    
    #remove infinite
    obj_exp$x[which(is.infinite(obj_exp$x))] <- NA;
    obj_exp$y[which(is.infinite(obj_exp$y))] <- NA;  
    
    #look at correlations
    cor_real <- cor(obj_real$x,obj_real$y,use='pairwise.complete');
    cor_exp  <- cor(obj_exp$x, obj_exp$y, use='pairwise.complete');
    cor_log  <- cor(obj_log$x ,obj_log$y ,use='pairwise.complete');
    
    #compute R^2 difference and take ratio
    log_check_score <- log2(abs(cor_real**2 - cor_log**2)/abs(cor_real**2 - cor_exp**2))
    
    #dataset is in log scale
    if(log_check_score<=0){
      return(TRUE);
    }
    #dataset is not in log scale
    else{
      return(FALSE);
    }
  }
}

.GEM_log_check_v1(my_gse) # log scale
my_gse <- list("expr"=exp_mat, "platform"="GPL570")
save(my_gse, file="data/small_gses2/E-MTAB-3604.RData")



load("ref/gpl570_probe_gene.RData") # probe_gene
xy_genes <- probe_gene %>% inner_join(xy_chr_map, by=c("gene"="hgnc_symbol")) 
sl <- massiRAcc(exp_mat, xy_genes %>% filter(chromosome_name=="Y") %>%pull(probes),  plot=T)  # good sep
sl_df <- tibble("sample_acc"=names(sl),
                "expr_sex"=unlist(sl))

pdat3 <- sl_df %>% inner_join(pdat2, by=c("sample_acc"="Array Data File"))  %>%
  rename(metadata_sex=Sex,
         sex=expr_sex,
         age=Age,
         race_ethnicity=`ethnic group`,
         BMI=`body mass index`,
         pack_years=`Pack Year`) %>%
  mutate(race_ethnicity=ifelse(str_detect(race_ethnicity, "Black"), "black", "white")) %>%
  mutate(tissue=tolower(tissue))

pdat3 %>% 
  select(sample_acc, sex, metadata_sex, smok, tissue, age, race_ethnicity, pack_years, everything()) %>%
  write_csv("data/pdata_filt2_v2/e-mtab-3604.csv")

#  filter(expr_sex != Sex) # all match!

#"E-MTAB-5278" - processed available, GPL570
exp <- read_tsv("data/array_exp/E-MTAB-5278/H1_QASMC_Blood_Hs_MRNA_normMatrix.txt")
pdat <- read_csv("data/pdata_filt2/e-mtab-5278.csv")
full_p <- read_tsv("data/array_exp/E-MTAB-5278.sdrf.txt") 
pdat2 <- pdat %>% inner_join(full_p %>% select(`Source Name`, `Array Data File`)) %>% distinct()
length(unique(pdat2$`Source Name`)) # 118
length(unique(pdat2$`Array Data File`)) # 119
length(intersect(pdat2$`Array Data File`, colnames(exp))) # 116
setdiff(pdat2$`Array Data File`, colnames(exp)) # 3 missing
pdat2.2 <- pdat2 %>% filter(`Array Data File` %in% colnames(exp))
length(unique(pdat2.2$`Source Name`)) # 115 - one person has two samples
pdat2.2[duplicated(pdat2.2$`Source Name`),]
pdat2.2 %>% filter(`Source Name`=="DD3056")


exp_sm <- exp[,c("Hybridization Name", pdat2.2$`Array Data File`)] %>% 
  rename(gene=`Hybridization Name`) %>% filter(gene!="Reporter REF")
exp_mat <- as.matrix(apply(exp_sm %>% select(-gene), c(1,2), as.numeric))
rownames(exp_mat) <- exp_sm$gene
colnames(exp_mat) <- str_replace_all(colnames(exp_mat), "-", "\\.")
my_gse <- list("expr"=exp_mat)

my_gse <- list("expr"=exp_mat, platform="GPL570")
save(my_gse, file="data/small_gses2/E-MTAB-5278.RData")



sl <- massiRAcc(exp_mat, xy_genes %>% 
                  filter(gene %in% rownames(exp_mat)) %>%
                  filter(chromosome_name=="Y") %>% 
                  pull(gene),  plot=F)  # good sep

sl_df <- tibble("sample_acc"=names(sl),
                "expr_sex"=unlist(sl))

pdat3 <- sl_df %>% inner_join(pdat2.2 %>% 
                                mutate(`Array Data File`=str_replace_all(`Array Data File`, "-", "\\.")), 
                              by=c("sample_acc"="Array Data File"))

pdat4 <- pdat3 %>% filter(sex == expr_sex) %>%  # 4 mismatch
  rename(metadata_sex=sex,
         sex=expr_sex,
         id=`Source Name`,
         race_ethnicity=`ethnic group`) %>%
  mutate(race_ethnicity=ifelse(race_ethnicity=="African American", "black", "white")) 
pdat4 %>% write_csv("data/pdata_filt2_v2/e-mtab-5278.csv")

#"E-MTAB-5279" - processed available, GPL570
exp <- read_tsv("data/array_exp/E-MTAB-5279/H2_BLDSMK01_Blood_Hs_MRNA_normMatrix.txt")
pdat <- read_csv("data/pdata_filt2/e-mtab-5279.csv")
full_p <- read_tsv("data/array_exp/E-MTAB-5279.sdrf.txt") 

pdat2 <- pdat %>% inner_join(full_p %>% select(`Source Name`, `Array Data File`)) %>% distinct()
colnames(exp)
setdiff(pdat2$`Array Data File`, colnames(exp)) # 4 missing
pdat2.2 <- pdat2 %>% filter(`Array Data File` %in% colnames(exp))

exp_sm <- exp[,c("Hybridization Name", pdat2.2$`Array Data File`)] %>% 
  rename(gene=`Hybridization Name`) %>% filter(gene!="Reporter REF")
exp_mat <- as.matrix(apply(exp_sm %>% select(-gene), c(1,2), as.numeric))
rownames(exp_mat) <- exp_sm$gene

my_gse <- list("expr"=exp_mat, platform="GPL570")
save(my_gse, file="data/small_gses2/E-MTAB-5279.RData")



sl <- massiRAcc(exp_mat, xy_genes %>% 
                  filter(gene %in% rownames(exp_mat)) %>%
                  filter(chromosome_name=="Y") %>% 
                  pull(gene),  plot=T)  # good sep

sl_df <- tibble("sample_acc"=names(sl),
                "expr_sex"=unlist(sl))

pdat3 <- sl_df %>% inner_join(pdat2.2, 
                              by=c("sample_acc"="Array Data File"))

pdat3 %>% 
  rename(
    metadata_sex=sex,
    sex=expr_sex,
    id=`Source Name`,
    race_ethnicity=`ethnic group`
  ) %>%
  mutate(race_ethnicity="white") %>%
  write_csv("data/pdata_filt2_v2/e-mtab-5279.csv")


#"E-MEXP-1277" - processed available - agilent
load.Rdata("data/array_exp/E-MEXP-1277.eSet.r", "my_eset")
colnames(fData(my_eset)) # hgnc_symbol: "Composite.Element.Database.Entry.locus."
exp_mat <- eData(my_eset)
head(exp_mat[,1:5])
#my_eset@assayData$
# G, Gb, R, Rb
