
library('MetaIntegrator')

# code for downloading the small datasets
list_studies <- read_csv("data/supp_tables/s4_list_smok_studies_formatted.csv")
list_studies$study
paste(list_studies$study, collapse="','")
my_studies <- c('GSE103174',
                'GSE13896','GSE16149','GSE17913','GSE18723','GSE19027',
                'GSE20189','GSE2125','GSE21862','GSE31210','GSE32539',
                'GSE42057','GSE42743','GSE4302','GSE44456','GSE46699',
                'GSE55962','GSE56768',
                'GSE7895','GSE87072', 'GSE8987','GSE994')

my_studies <- c("GSE20189", "GSE42057")
res <- lapply(my_studies, function(study){
  print(study)
  gse <- getGEOData(study )
  save(gse, file=sprintf("data/small_studies/%s.RData", study))
})


