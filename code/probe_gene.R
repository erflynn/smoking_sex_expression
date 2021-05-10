# create a probe/gene mapping for each platform

list_studies <- read_csv("data/supp_tables/s4_list_smok_studies_formatted.csv")
list_studies %>%
  filter(str_detect(platform, ";")) %>%
  select(study, platform, tissue)

# GSE8987 --> GPL571
# GSE32539 --> GPL6244 
l2 <- list_studies %>% 
  mutate(platform=case_when(
    study=="GSE8987" ~ "GPL571",
    study=="GSE32539" ~ "GPL6244",
    TRUE ~ platform
  )) %>%
  filter(study != "GSE55962")

l2 %>% group_by(platform) %>% count()
plats <- l2 %>% distinct(platform) %>% pull(platform)

# "GPL13667" "GPL570"   "GPL96"    "GPL571"   "GPL6104"  "GPL6244" 
# GPL6104 - Illumina
library('GEOmetadb')
con <- dbConnect(SQLite(), "../GEOmetadb.sqlite")
plat_d <- dbGetQuery(con, sprintf("SELECT * FROM gpl WHERE gpl IN ('%s');",
                        paste(plats, collapse="','")))
dbDisconnect(con)

plat_d %>% select(gpl, title, bioc_package)





probe_gene_map <- function(x, gpl){
  mapped_probes <- mappedkeys(x)
  xx <- as.list(x[mapped_probes])
  probe_gene <- tibble(
    "probes"=names(xx),
    "gene" = unlist(xx)
  )
  print(probe_gene %>% group_by(probes) %>% count() %>% filter(n>1) %>% nrow())
  save(probe_gene, file=sprintf("ref/%s_probe_gene.RData", gpl))
}
# --- GPL570, hgu133plus2 --- #
library(hgu133plus2.db)
probe_gene_map(hgu133plus2SYMBOL, "gpl570")

#  GPL13667                      hgu219
library(hgu219.db)
probe_gene_map(hgu219SYMBOL, "gpl13667")

#    GPL96                     hgu133a
library(hgu133a.db)
probe_gene_map(hgu133aSYMBOL, "gpl96")

#   GPL571                    hgu133a2
library(hgu133a2.db)
probe_gene_map(hgu133a2SYMBOL, "gpl571")

#  GPL6244 hugene10sttranscriptcluster
library(hugene10sttranscriptcluster.db)
probe_gene_map(hugene10sttranscriptclusterSYMBOL, "gpl6244")

#  GPL6104                        ??
gpl6104 <- GEOquery::getGEO("GPL6104")
probe_gene <- dataTable(gpl6104)@table[,c("ID", "ILMN_Gene")] %>%
  dplyr::rename(probes="ID", gene="ILMN_Gene")
print(probe_gene %>% group_by(probes) %>% count() %>% filter(n>1) %>% nrow())
save(probe_gene, file=sprintf("ref/%s_probe_gene.RData", "gpl6104"))
