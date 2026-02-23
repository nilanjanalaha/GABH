###############################################################################
# Project: OSCC Data Analysis
# Content: Grouping Logic and Discovery Calculation Functions
###############################################################################

library(tidyverse)
library(dplyr)
library(tidyr)

# 1. Coarse group: Group by Chromosome only
level_2_group <- function(mydata, group.info) {
  group.info.large <- mydata %>%
    group_by(Chr) %>%
    summarize(freq = n(), .groups = "drop") %>%
    group_by(Chr) %>%
    mutate(min_freq = min(freq),
           centro_new = ifelse(min_freq == 1, 2 - row_number(), 1)) %>%
    ungroup() %>%
    mutate(group_id = cumsum(centro_new)) %>%
    select(Chr, group_id)
  
  mydata %>%
    inner_join(group.info.large, by = c("Chr")) %>%
    pull(group_id)
}

# 2. Semi-coarse group: Group by Chromosome and Strand
level_1_group <- function(mydata, group.info) {
  group.info.large <- mydata %>% 
    mutate(strand = factor(strand, levels = c("-", "+"))) %>% 
    group_by(Chr, strand) %>%
    summarize(freq = n(), .groups = "drop") %>%
    group_by(Chr, strand) %>%
    mutate(min_freq = min(freq), 
           centro_new = ifelse(min_freq == 1, 2 - row_number(), 1)) %>%
    ungroup() %>%
    mutate(group_id = cumsum(centro_new)) %>% 
    select(Chr, strand, group_id)
  
  mydata %>% 
    inner_join(group.info.large, by = c("Chr", "strand")) %>% 
    pull(group_id)
}

# 3. Fine-grained group: Sub-dividing groups by size k
level_k_groups <- function(mydata, k, group.info) {
  sgrid <- seq(0, 1, length.out = 101)
  max_group_size <- k
  
  further_cut <- function(x, n) { 
    if (n == 1) return(rep(1, length(x)))
    else cut(x, n, labels = FALSE)
  }
  
  new_group_dat <- mydata %>% 
    inner_join(group.info, by = c("Chr", "strand", "centro")) %>% 
    group_by(group_id) %>% 
    mutate(new_group_pos = group_id + sgrid[further_cut(1:n(), floor((n() - 1) / max_group_size) + 1)])
  
  fine_group <- new_group_dat %>% 
    inner_join(data.frame("new_group_pos" = unique(new_group_dat$new_group_pos),
                          "new_group_id" = rank(unique(new_group_dat$new_group_pos))), 
               by = "new_group_pos") %>%
    pull(new_group_id)
  
  fine_group
}

# 4. Wrapper to select group based on k
give_group <- function(k, mydata, group.info) {
  # Note: 'orig_group' must exist in the environment where this is called
  group <- switch(
    as.character(k),
    "-2" = level_2_group(mydata, group.info),
    "-1" = level_1_group(mydata, group.info),
    "0"  = orig_group, 
    level_k_groups(mydata, k, group.info)
  )
  group
}

# 5. Helper to get the total number of groups
n_group <- function(k, mydata, group.info) {
  group <- give_group(k, mydata, group.info)
  max(group)
}

# 6. Discovery calculations (BH, TST, LSL, SABHA)
n_discoveries <- function(k, mydata, miss.val, pval, chr, group.info, alpha.fix, regulation) {
  group <- give_group(k, mydata, group.info)
  K <- max(group)
  
  assign <- as.list(1:K)
  for(i in 1:K) { assign[[i]] <- which(group == i) }
  
  # --- BH ---
  bh_adj <- p.adjust(pval, method = "BH")
  BH.no <- sum(bh_adj < alpha.fix)
  
  # --- LSL & TST GBH ---
  # These functions (gbh) must be loaded from analysis.R
  lsl.gbh <- gbh(alpha.fix, pval, group, assign, mydata, "lsl", chr)
  lsl.no <- nrow(lsl.gbh)
  
  tst.gbh <- gbh(alpha.fix, pval, group, assign, mydata, "tst", chr)
  tst.no <- nrow(tst.gbh)
  
  # --- SABHA ---
  # These functions (Solve_q_block, SABHA_method) must be loaded from All_q_est_functions.R
  tau <- 0.5; eps <- 0.1
  ADMM_params <- c(10^2, 10^3, 2, 5000, 1e-3)
  qhat <- Solve_q_block(pval, tau, eps, group, ADMM_params)
  SABHA_Res <- SABHA_method(pval, qhat, alpha.fix, tau)
  sabha.no <- length(SABHA_Res)
  
  return(c(BH = BH.no, TST = tst.no, LSL = lsl.no, SABHA = sabha.no))
}