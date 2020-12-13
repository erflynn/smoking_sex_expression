
library(tidyverse)
clinical <- read_tsv("~/Downloads/clinical.cases_selection.2020-12-09.tar-1/clinical.tsv")
exposure <- read_tsv("~/Downloads/clinical.cases_selection.2020-12-09.tar-1/exposure.tsv")
fh <- read_tsv("~/Downloads/clinical.cases_selection.2020-12-09.tar-1/family_history.tsv") # NO ROWS

# ... hmm there are only 40 clinical rows and 20 exposure rows? what is going on here
# ok it's only selectinng the first page...

exposure2 <- exposure %>% mutate(across(everything(), ~ifelse(.=="'--", NA, .)))
clinical2 <- clinical %>% mutate(across(everything(), ~ifelse(.=="'--", NA, .)))
count_missing <- apply(exposure, 2, function(x) sum(x=="'--"))
exposure2 %>% 
  select(age_at_onset, cigarettes_per_day, 
         contains("smok"), contains("tobacco")) %>%
  mutate(across(c(contains("age"), contains("year"), 
                  contains("time"), contains("day")), as.numeric)) %>%
  summary()

clinical %>% select(gender, ethnicity, race, age_at_diagnosis)
