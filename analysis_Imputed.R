

######################## Loading packages #######################
library(sgof) 
library(xtable)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)

setwd("~/Library/CloudStorage/Dropbox/Nyabada project/Submission/DataAnalysis/OSCC_data_analysis")
source("Source/Preparing_data.R")

do.impute <-TRUE
group_type <- "original"


####################################################################
################# t-test and the p-values ###########################
####################################################################

regulation <- rep(0, nrow(mydata))

#--------------------- t-test  --------------------------
#The NA values are omitted
my.t.test <- function(x,y)
{
  if(median(x-y, na.rm = TRUE)>0)
    pval <- t.test(x,y, alternative="greater", mu=1, paired=TRUE,  na.action=na.omit)$p.value
  
  if(median(x-y, na.rm=TRUE)<0)
    pval <- t.test(x,y, alternative="less", mu=-1, paired=TRUE,  na.action=na.omit)$p.value
  
  pval
}

#------------------- The main function for t-test ------------------------------
pool.t.test <- function(mydata, impute)
{
  n <- nrow(mydata)
  
  #diseased people
  dis <- mydata[,6:23]
  
  #normal people
  nor <- mydata[,24:41]
  #dis <- mean_impute(dis)
  #nor <- mean_impute(nor)
  
  #---------Imputation
  
  if(impute){
    dis <- median_impute(dis)
    nor <- median_impute(nor)
  }
  
  pv <- c(1:n)
  for (i in 1:n)
  {
    pv[i] <-my.t.test(as.numeric(dis[i,]), as.numeric(nor[i,]))
    x <- as.numeric(dis[i,])
    y <- as.numeric(nor[i,])
    regulation[i] <- ifelse(median(x-y, na.rm=TRUE)>0,1,-1)
  }
  cbind(pv, regulation)
}

#do.impute <- FALSE # Median imputation not done

# t-test without median imputation
pval <- pool.t.test(mydata, do.impute)[, 1]
# t-tester direction diche, vector of length 522
regulation <- pool.t.test(mydata, do.impute)[, 2]

####################################################################
################# Group assignments ###########################
####################################################################

group.info <- mydata %>% mutate(strand=factor(strand, levels = c("-", "+"))) %>% group_by(Chr, strand, centro) %>%
  summarize(freq=n(), .groups = "drop") %>%
  group_by(Chr, strand) %>%
  mutate(min_freq=min(freq), centro_new=ifelse(min_freq == 1, 2-row_number(), 1)) %>%
  ungroup() %>%
  mutate(group_id = cumsum(centro_new)) %>% select(Chr, strand, centro, group_id)

# Old group
orig_group <- mydata %>% inner_join(group.info, by = c("Chr", "strand", "centro")) %>% pull(group_id)

# Updated on July 2024
# Groups are assigned randomly, but the group sizes remained fixed
random_group <- sample(orig_group)

# Larger group
group.info.large <- mydata %>% mutate(strand=factor(strand, levels = c("-", "+"))) %>% group_by(Chr, strand)%>%
  summarize(freq=n(), .groups = "drop") %>%
  group_by(Chr, strand) %>%
  mutate(min_freq=min(freq), centro_new=ifelse(min_freq == 1, 2-row_number(), 1)) %>%
  ungroup() %>%
  mutate(group_id = cumsum(centro_new)) %>% select(Chr, strand, group_id)

coarse_group <- mydata %>% inner_join(group.info.large, by = c("Chr", "strand")) %>% pull(group_id)

# Much smaller group assignment : updated on Aug 1, 2024
sgrid <- seq(0,1,length.out=101)
max_group_size <- 4
# new_group_dat <- mydata %>% inner_join(group.info, by = c("Chr", "strand", "centro")) %>%
#              group_by(group_id) %>% mutate(group_freq = n(), row_num = row_number(), 
#                                            num_split = floor((group_freq-1)/max_group_size) + 1) %>%
#                                 mutate(new_group_pos = group_id + sgrid[floor((row_num-1)/max_group_size)+1] )
# 
# new_group_dat <- mydata %>% inner_join(group.info, by = c("Chr", "strand", "centro")) %>%
#   group_by(group_id) %>% mutate(group_freq = n(), num_split = floor((group_freq-1)/max_group_size) + 1) 

further_cut <- function(x, n){ 
  if (n == 1) return(rep(1,length(x)))
  else cut(x, n, labels=FALSE)
}

new_group_dat <- mydata %>% inner_join(group.info, by = c("Chr", "strand", "centro"))  %>% group_by(group_id) %>% 
  mutate(new_group_pos = group_id + sgrid[further_cut(1:n(), floor((n()-1)/max_group_size) + 1)] )

fine_group <- new_group_dat %>% inner_join(data.frame("new_group_pos" = unique(new_group_dat$new_group_pos),
                         "new_group_id"  = rank(unique(new_group_dat$new_group_pos))), by="new_group_pos") %>%
         pull(new_group_id) # This is the new group assignment

group <- switch(group_type,  "original" = orig_group, "random" = random_group, "smaller" = fine_group, "larger" = coarse_group)


# # New group: This one is the group used to check if we split the group larger than 30, is the result difference,
# before submitting the paper to Nature
# new_group_dat <- mydata %>% inner_join(group.info, by = c("Chr", "strand", "centro")) %>% 
#              group_by(group_id) %>% mutate(group_freq = n(), row_num = row_number(),
#                                 new_group_pos = if_else(group_freq > 30, 
#                                                 if_else(row_num <= group_freq/2, group_id + 0.1, group_id + 0.6),
#                                                 group_id))
# 
# group <- new_group_dat %>% inner_join(data.frame("new_group_pos" = unique(new_group_dat$new_group_pos), 
#                          "new_group_id"  = rank(unique(new_group_dat$new_group_pos))), by="new_group_pos") %>%
#          pull(new_group_id) # This is the new group assignment

# Check: every group has less than 30 freq
#table(group)

#gives number of missing pairs for 522 genes
miss.val <- read.table("Source/miss.txt",header=TRUE)



#The number of groups
K <- max(group)

# The assihnment into groups
assign <- as.list(1:K)
for(i in 1:K){
  assign[[i]] <- which(group==i)
}

#################################################################
################# Auxilliary functions for #######################
################# calculating inclusion probability ##############
##################################################################

#___________________________ LSL ______________________________________
# g is the sorted p-values of a g
lsl <- function(g)
{
  g <- sort(g)
  ng <- length(g)
  lg <- (ng+1-1:ng)/(1-g)
  #What happens if ng=1, which means the cluster has length only 1?
  i=2
  while(lg[i]<lg[i-1])
  {
    if(i==ng)
      break
    i=i+1
  }
  
  min((floor(lg)+1)/ng, 1)
}

#__________________________ TST procedure _____________________________________________

tst <- function(g, al)
{
  al <- al/(1+al)
  #why this is 0.05?
  r <- BH(g, al)$Rejections
  ng <- length(g)
  (ng-r)/ng
}

###################################################################
###################### Main function ##############################
###################################################################

# Method: lsl or tst
# chr and mydata comes from the sourced R files
#---------------------  GBH ----------------------------------------------------
gbh <- function(al,pval,group,assign, mydata, method, chr)
{
  n <- length(pval)
  wp <- 1:n
  k <- length(assign)
  s <- 0
  for(i in 1:k)
  {
    ind <- assign[[i]]
    ng <- length(ind)
    if (method=="tst")
    {pig <- tst(pval[ind], al)} else { pig <- lsl(pval[ind])}
    s <- ng*pig+s
    wp[ind] <- pig*pval[ind]/(1-pig)
  }
  pio <- s/n
  ial <- al
  if (method=="tst")
    al <- al/(1+al)
  #pi-hat is pio, aw alpha-w, wp is p-w
  aw <- al/(1-pio)
  swp <- sort(wp)
  vec <- (c(1:n)*aw)/n
  
  i=1
  while(swp[i]<=vec[i])
    i=i+1
  i=i-1
  if (i==0) return(data.frame(mydata[0,c(1:3,5)], pval= round(pval[0],4),
                              missing = miss.val[0,1], regulation = regulation[0],
                              group_id = group[0], type=NULL))
  cut <- swp[i]
  name <- which(wp<=cut)
  class.prob <- name #Probability of the classes where the significant genes belong
  for(i in 1:length(name))
    if(method=="tst")
    {
      class.prob[i] <- tst(pval[assign[[group[name[i]]]]], ial)
    } else {class.prob[i] <- lsl(pval[assign[[group[name[i]]]]])}
  
  data.frame(mydata[name,c(1:3,5)], pval= round(pval[name],4),
             missing = miss.val[name,1], regulation = regulation[name],
             group_id = group[name],
             type=toupper(method))
}

##################################################################
################    Application ##################################
##################################################################

# Change ALPHA here
alpha.fix <- 0.05

#------------------- BH--------------------------------------------
bh <- p.adjust(pval, method="BH")
name <- which(bh < alpha.fix)
if(length(name)==0) print("Error: no tests are significant for BH")

my.bh <- data.frame(mydata[name,c(1:3,5)], pval= round(pval[name],4),
                    missing = miss.val[name,1], regulation = regulation[name],
                    group_id = group[name],
                    type="BH")



#------------------- LSL GBH ----------------------------------------------------------
#al: alpha
# pval: pvalue
# group: read kora hoechilo, asign ei filei defined hoeche

lsl.gbh <- gbh(alpha.fix, pval, group, assign, mydata, "lsl", chr)

# I don't know why chr is there, probably will run without this chr.
# probably does not need chr
# jodi chr dite chas tahole mydata$Chr etar elementgulo (23ta unique elements)
# etake input hisabe use koris.

#-------------------- TST GBH ----------------------------------------------------------
# al is the level
tst.gbh <- gbh(alpha.fix, pval, group, assign, mydata, "tst", chr)


#-------------------- SABHA Method  ----------------------------------------------------------

source('Source/All_q_est_functions.R')

tau = 0.5; eps = 0.1 # parameters for SABHA # Robust to the choice of eps
ADMM_params = c(10^2, 10^3, 2, 5000, 1e-3) # alpha_ADMM,beta,eta,max_iters,converge_thr


SABHA_method = function(pvals, qhat, alpha, tau){
  # Use the original, or estimated q as input
  pvals[pvals>tau] = Inf
  khat=max(c(0,which(sort(qhat*pvals)<=alpha*(1:length(pvals))/length(pvals))))
  which(qhat*pvals<=alpha*khat/length(pvals))
}

qhat = Solve_q_block(pval,tau,eps,group,ADMM_params)
SABHA_Res = rep(0,length(pval))
SABHA_Res[SABHA_method(pval,qhat,alpha.fix,tau)] = 1
name <- which(SABHA_Res==1)
if(length(name)==0) print("Error: no tests are significant for SABHA")

sabha <- data.frame(mydata[name,c(1:3,5)], pval= round(pval[name],4),
                    missing = miss.val[name,1], regulation = regulation[name],
                    group_id = group[name],
                    type="SABHA")


discoveries <- full_join(my.bh, tst.gbh, by = names(my.bh)[-length(names(my.bh))], suffix = c("_BH", "_TST")) %>%
  full_join(lsl.gbh, by = names(my.bh)[-length(names(my.bh))]) %>%
  full_join(sabha, by = names(my.bh)[-length(names(my.bh))], suffix = c("_LSL", "_SABHA")) %>%
  unite(signif_tests, starts_with("type"), na.rm = TRUE, sep=", ") %>%
  mutate(n_sig = stringi::stri_count_words(signif_tests) )

df_cleaned <- discoveries %>%
  # Create the new columns based on whether the test is in signif_tests
  mutate(
    signif_list = str_split(signif_tests, ",\\s*"),
    BH = map_chr(signif_list, ~ ifelse("BH" %in% .x, "x", "")),
    TST = map_chr(signif_list, ~ ifelse("TST" %in% .x, "x", "")),
    LSL = map_chr(signif_list, ~ ifelse("LSL" %in% .x, "x", "")),
    SABHA = map_chr(signif_list, ~ ifelse("SABHA" %in% .x, "x", "")),
    Arm = ifelse(centro == 0, "p", "q")
  ) %>%
  select(-signif_tests, -group_id, -n_sig, -signif_list, -centro) %>%
  relocate(BH, TST, LSL, SABHA, .after = last_col())

df_cleaned <- df_cleaned %>%
  arrange(pval) %>%
  mutate(Row = row_number())%>%
  select(Row, everything())
# View cleaned result
print(df_cleaned)