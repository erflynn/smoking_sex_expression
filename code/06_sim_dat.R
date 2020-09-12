# code for exploring simulations
#
# NOTE: seems like disc/validation is a POOR idea
#
# TODO:
#  - what do we expect to see for the gene mean distributions? (g)
#  - ranges of B, A, AB
#  - plot how things are affected (ranges of coef, comparison to each other, FDR cutoffs)
#  - does this look like real data?
# 
# UPDATE. tried a bunch of ranges but it did not seem to help! :/


require('limma')
require('tidyverse')


simGene <- function(betaA, betaB, betaInt, n){
  g = rnorm(1,0,0.5)
  
  # Generate n*4 trials
  A = c(rep(c(0), n*2), rep(c(1), n*2))
  B = rep(c(rep(c(0), n), rep(c(1), n)), 2) 
  e = rnorm(n*4, 0, sd=1) 
  
  # Generate your data using the regression equation
  y = g + betaA*A + betaB*B + betaInt*A*B + e
  
  # Join the variables in a data frame
  #tibble("A"=A, "B"=B, "Y"=Y, gene_id="gene_id")
  names(y) <- seq(1:n)
  return(y)
}  

fullSimDat <- function(betaA, betaB, betaInt, n, num_int, num_both, num_main, num_no){
  gex1 <- do.call(rbind, 
          c(lapply(seq(1,num_int), function(x) simGene(betaA, betaB, betaInt, n)),
            lapply(seq(1,num_both), function(x) simGene(betaA, betaB, 0, n)),
            lapply(seq(1,num_main), function(x) simGene(betaA, 0, 0, n)) ,
            lapply(seq(1,num_main), function(x) simGene(0, betaB, 0, n)) ,
            lapply(seq(1,num_no), function(x) simGene(0, 0, 0, n))) ) 
  colnames(gex1) <- seq(1,n*4) %>% map(~paste("x",.x, sep="")) %>% unlist()
  gex <- gex1 %>% as_tibble() %>%
    mutate(gene_id=seq(1,nrow(gex1)) %>% map(~paste("g",.x, sep="")) %>% unlist()) %>%
    mutate(gene_type=c(rep("intAB", num_int),
                       rep("both", num_both),
                       rep("mainA", num_main),
                       rep("mainB", num_main),
                       rep("no", num_no)))
  phe <- tibble("person_id"= colnames(gex1),
                   "A" = c(rep(c(0), n*2), rep(c(1), n*2)),
                   "B" = rep(c(rep(c(0), n), rep(c(1), n)), 2) )
  return(list("gex"=gex, "phe"=phe))
}



  

simSplit <- function(){
  dat <- fullSimDat(betaA, betaB, betaInt, ndisc/4, num_int, num_both, num_main, num_no)
  design <- model.matrix(~A+B+A*B, data=dat$phe)
  gex <- dat$gex %>% select(-gene_id, -gene_type)
  rownames(gex) <- dat$gex$gene_id
  fit <- lmFit(gex, design=design)
  fit <- eBayes(fit)
  int_top <- topTable(fit, coef="A:B", n=nrow(gex))
  int_top$gene_id <- rownames(int_top)
  
  int_top <- int_top %>% left_join(dat$gex %>% dplyr::select(gene_id, gene_type))
  
  
  # divide by covariate
  a0  <- dat$phe %>% filter(A==0)
  a1  <- dat$phe %>% filter(A==1)
  gex0 <- gex[,a0$person_id]
  rownames(gex0) <- rownames(gex)
  gex1 <- gex[,a1$person_id]
  rownames(gex1) <- rownames(gex)
  
  fit0 <- lmFit(gex0, design=model.matrix(~B, dat=a0))
  fit1 <- lmFit(gex1, design=model.matrix(~B, dat=a1))
  fit0 <- eBayes(fit0)
  fit0_top <- topTable(fit0, coef="B", n=nrow(gex)) 
  fit0_top$gene <- rownames(fit0_top)
  fit1 <- eBayes(fit1)
  fit1_top <- topTable(fit1, coef="B", n=nrow(gex))
  fit1_top$gene <- rownames(fit1_top)
  
  fit0_hit <- fit0_top %>% filter(adj.P.Val < 0.05)
  fit1_hit <- fit1_top %>% filter(adj.P.Val < 0.05)
  sig_genes <- union(fit0_hit$gene, fit1_hit$gene)
  
  dat2 <- fullSimDat(betaA, betaB, betaInt, nvalid/4, num_int, num_both, num_main, num_no)
  design2<- model.matrix(~A+B+A*B, data=dat2$phe)
  gex2 <- dat2$gex %>% select(-gene_id, -gene_type) %>% as.matrix()
  rownames(gex2) <- dat2$gex$gene_id
  fit2 <- lmFit(gex2[sig_genes,], design=design2)
  fit2 <- eBayes(fit2)
  int_top2 <- topTable(fit2, coef="A:B", n=length(sig_genes))
  int_top2$gene_id <- rownames(int_top2)
  int_top2 <- int_top2 %>% left_join(dat2$gex %>% dplyr::select(gene_id, gene_type))
  int_top2 %>% filter(adj.P.Val < 0.05) # FP = 1, FN=5, TP=5, TN=989
  split_res <- int_top2 %>% filter(adj.P.Val < 0.05)  %>% pull(gene_type)
  split_res
  tp = sum(split_res=="intAB")
  fp=length(split_res) - tp
  fn=num_int-tp
  return(tibble("tp"=tp, "fp"=fp, "fn"=fn))
}
# what if we had just used the whole dataset?
simFull <- function(){
  dat <- fullSimDat(betaA, betaB, betaInt, (ndisc+nvalid)/4,num_int, num_both, num_main, num_no)
  design <- model.matrix(~A+B+A*B, data=dat$phe)
  gex <- dat$gex %>% select(-gene_id, -gene_type)
  rownames(gex) <- dat$gex$gene_id
  fit <- lmFit(gex, design=design)
  fit <- eBayes(fit)
  int_top <- topTable(fit, coef="A:B", n=nrow(gex))
  int_top$gene_id <- rownames(int_top)
  
  int_top <- int_top %>% left_join(dat$gex %>% dplyr::select(gene_id, gene_type))
  full_res <- int_top %>% filter(adj.P.Val < 0.05)  %>% pull(gene_type)
  tp = sum(full_res=="intAB")
  fp=length(full_res) - tp
  fn=num_int-tp
  return(tibble("tp"=tp, "fp"=fp, "fn"=fn))
}

betaA = 1
betaB = -1
betaInt = -1
ndisc = 200
nvalid = 200
num_int =10
num_both=90
num_main=100
num_no=700

betaAs <- c(-1, -0.5, 0.5, 1)
betaBs <- c(-1, -0.5,0.5, 1)
betaInts <- c(-1, -0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5, 1)
cross_list <- cross3(betaAs, betaBs, betaInts)
simAll <- lapply(cross_list, function(x){
  betaA <- x[[1]]; betaB <- x[[2]]; betaInt <- x[[3]];
  simF <- do.call(rbind, lapply(1:10, function(y) simFull()))
  simS <- do.call(rbind, lapply(1:10, function(y) simSplit()))
  simF <- simF %>% mutate(type="full", betaA=betaA, betaB=betaB, betaInt=betaInt)
  simS <- simS %>% mutate(type="split",betaA=betaA, betaB=betaB, betaInt=betaInt)
  return(bind_rows(simF, simS))
  })
simAll2 <- do.call(rbind, simAll)
betaA <- cross_list[[1]][[1]]

simAll2 <- lapply(1:length(simAll), function(x){
  betaA <- cross_list[[x]][[1]]; betaB <- cross_list[[x]][[2]]; 
  betaInt <- cross_list[[x]][[3]];
  simF <- simAll[[x]][1:10,]
  simA <- simAll[[x]][11:20,]
  simF <- simF %>% mutate(type="full", betaA=betaA, betaB=betaB, betaInt=betaInt)
  simS <- simS %>% mutate(type="split",betaA=betaA, betaB=betaB, betaInt=betaInt)
  return(bind_rows(simF, simS))
})
simAll3 <- do.call(rbind, simAll2)
simAll4 <- simAll3 %>% 
  select(-fn) %>%
  pivot_longer(tp:fp, names_to="metric", values_to="count")  %>%
  group_by(betaA, betaB, betaInt, type, metric) %>%
  summarize(mu=mean(count), sd=sd(count)) %>%
  mutate(mu_l=mu-1.96*sd, mu_u=mu+1.96*sd) %>%
  ungroup()

ggplot(simAll4, aes(y=mu, x=betaInt))+
  #geom_point(aes(col=metric), alpha=0.5)+
  #geom_errorbar(aes(ymin=mu_l, ymax=mu_u, col=metric))+
  geom_smooth(aes(col=metric, lty=type))+
  facet_grid(betaA ~ betaB)+
  theme_bw()+
  ylab("number of tp or fp")+
  ylim(-2, 12)+
  xlab("interaction coefficient")


betaA <- 0.5
betaB <- 0.5
betaInt <- 0.1
simSamp <- do.call(rbind, lapply(c(20, 50, 100, 150, 200),
       function(x){
         ndisc <- x;
         nvalid <- x;
         simF <- do.call(rbind, lapply(1:100, function(y) simFull()))
         simF$nper_grp <- nvalid/2
         return(simF)
       }))

ggplot(simSamp %>% select(-fn) %>%
         pivot_longer(tp:fp, names_to="metric", values_to="count")  %>%
         group_by(nper_grp, metric) %>%
         summarize(mu=mean(count), sd=sd(count)) %>%
         mutate(mu_l=mu-1.96*sd, mu_u=mu+1.96*sd) %>%
         ungroup()
         , aes(x=nper_grp, y=mu, col=metric))+
  geom_errorbar(aes(ymin=mu_l, ymax=mu_u, col=metric))+
  geom_point(alpha=0.8)

ggplot(simSamp %>% 
         pivot_longer(tp:fp, names_to="metric", 
                      values_to="count") %>%
       unite(grp, c(nper_grp, metric), remove=FALSE),
       aes(x=nper_grp, y=count, group=grp, col=metric))+
  geom_boxplot()+
  xlab("number of samples per group")+
  theme_bw()

nvalid <- 100
ndisc <- 100
simN <- do.call(rbind, lapply(c(700, 1700, 2700, 4700, 9700),
                                 function(x){
                                   num_no <- x
                                   simF <- do.call(rbind, 
                                                   lapply(1:100, function(y) simFull()))
                                   simF$num_no <- x
                                   return(simF)
                                 }))

ggplot(simN %>%
         pivot_longer(tp:fp, names_to="metric", values_to="count") %>%
         mutate(num_genes=num_no+300) %>%
         unite(grp, c(num_genes, metric), remove=FALSE),
       aes(x=num_genes, y=count, group=grp, col=metric))+
  geom_boxplot()+
  xlab("number of genes")+
  theme_bw()

ggplot(simN %>% select(-fn) %>%
         pivot_longer(tp:fp, names_to="metric", values_to="count")  %>%
         group_by(num_no, metric) %>%
         summarize(mu=mean(count), sd=sd(count)) %>%
         mutate(mu_l=mu-1.96*sd, mu_u=mu+1.96*sd) %>%
         ungroup() %>%
         mutate(num_genes=num_no+300)
       , aes(x=num_genes, y=mu, col=metric))+
  geom_errorbar(aes(ymin=mu_l, ymax=mu_u, col=metric))+
  geom_point(alpha=0.8)+
  xlab("number of genes")+
  ylab("count")+
  theme_bw()


# plot of estimated number of genes beta, direction
# same dir
# opp dir
# beta-range
# fp, fn
