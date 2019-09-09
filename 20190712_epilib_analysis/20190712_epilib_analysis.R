#Written for the analysis of the follow-up CRE library in the episomal MPRA


#D3_BC: DNA
#R0A_BC and R0B_BC: RNA at 0 µM Forsk replicate A and B
#R4A_BC and R4B_BC: RNA at 4 µM Forsk replicate A and B


#Establish workspace------------------------------------------------------------

options(stringsAsFactors = F)

#import necessary libraries

library(tidyverse)
library(cowplot)

#Load index and bcmap files-----------------------------------------------------

bc_DNA <- read_tsv('BCreads_txts/D3_BC.txt')
bc_R0A <- read_tsv('BCreads_txts/R0A_BC.txt')
bc_R0B <- read_tsv('BCreads_txts/R0B_BC.txt')
bc_R4A <- read_tsv('BCreads_txts/R4A_BC.txt')
bc_R4B <- read_tsv('BCreads_txts/R4B_BC.txt')

#Load barcode mapping table, sequences (most_common) are rcomp due to sequencing
#format

Follow_up_map <- read_tsv('../BCMap/CRE_15K/barcode_statistics.txt') %>%
  select(-num_unique_constructs, -num_reads, -num_reads_most_common)

#Join reads to bcmap------------------------------------------------------------

#Determine normalized BC reads per million. Join BC reads to BC mapping, keeping
#the reads only appearing in BC mapping and replacing na with 0 reads. 

bc_map_join_bc <- function(df1, df2) {
  df2 <- df2 %>%
    mutate(normalized = as.numeric((num_reads * 1000000) / (sum(num_reads))))
  keep_bc <- left_join(df1, df2, by = 'barcode') %>%
    mutate(normalized = if_else(is.na(normalized), 
                                0, 
                                normalized)) %>%
    mutate(num_reads = if_else(is.na(num_reads), 
                               0, 
                               num_reads))
  return(keep_bc)
}

bc_join_DNA <- bc_map_join_bc(Follow_up_map, bc_DNA)
bc_join_R0A <- bc_map_join_bc(Follow_up_map, bc_R0A)
bc_join_R0B <- bc_map_join_bc(Follow_up_map, bc_R0B)
bc_join_R4A <- bc_map_join_bc(Follow_up_map, bc_R4A)
bc_join_R4B <- bc_map_join_bc(Follow_up_map, bc_R4B)

#Median BC expression-----------------------------------------------------------

#Retain BCs with > 6 reads in DNA sample, left-join RNA to DNA BCs. Take ratio 
#of RNA/DNA normalized reads per million

dna7_join_rna_rep <- function(df1, df2) {
  filter_DNA <- filter(df1, num_reads > 6)
  DNA_RNA_join <- left_join(filter_DNA, df2,
                            by = c("barcode", "name", "most_common"), 
                            suffix = c('_DNA', '_RNA')) %>%
    mutate(ratio = normalized_RNA/normalized_DNA)
  return(DNA_RNA_join)
}

bc_DNA_RNA_0A <- dna7_join_rna_rep(bc_join_DNA, bc_join_R0A)
bc_DNA_RNA_0B <- dna7_join_rna_rep(bc_join_DNA, bc_join_R0B)
bc_DNA_RNA_4A <- dna7_join_rna_rep(bc_join_DNA, bc_join_R4A)
bc_DNA_RNA_4B <- dna7_join_rna_rep(bc_join_DNA, bc_join_R4B)

#Count barcodes per variant per DNA and RNA sample, set minimum of 8 BC's per 
#variant in DNA sample, take median BC RNA/DNA per variant, then per variant 
#determine the median absolute deviation of all barcode ratios. Only look at 
#variants with greater than 0 median expression

ratio_bc_med_var <- function(df) {
  bc_count_DNA <- df %>%
    group_by(name, most_common) %>%
    summarize(barcodes_DNA = n()) %>%
    filter(barcodes_DNA > 7)
  bc_count_RNA <- df %>%
    group_by(name, most_common) %>%
    filter(num_reads_RNA != 0) %>%
    summarize(barcodes_RNA = n())
  bc_DNA_RNA <- left_join(bc_count_DNA, bc_count_RNA, 
                          by = c('name', 'most_common')) %>%
    ungroup()
  bc_min_8_df <- left_join(bc_DNA_RNA, df, 
                           by = c('name', 'most_common')) %>%
    ungroup()
  med_ratio <- bc_min_8_df %>%
    group_by(name, most_common) %>%
    summarize(med_ratio = median(ratio)) %>%
    filter(med_ratio > 0)
  mad_ratio <- bc_min_8_df %>%
    group_by(name, most_common) %>%
    summarize(mad = mad(ratio, constant = 1))
  med_mad <- left_join(med_ratio, mad_ratio, 
                       by = c('name', 'most_common')) %>%
    mutate(mad_over_med = as.double(mad/med_ratio))
  bc_med <- inner_join(med_mad, bc_DNA_RNA, 
                       by = c('name', 'most_common')) %>%
    ungroup()
  return(bc_med)
}

med_ratio_R0A <- ratio_bc_med_var(bc_DNA_RNA_0A)
med_ratio_R0B <- ratio_bc_med_var(bc_DNA_RNA_0B)
med_ratio_R4A <- ratio_bc_med_var(bc_DNA_RNA_4A)
med_ratio_R4B <- ratio_bc_med_var(bc_DNA_RNA_4B)

#Combine biological replicates across concentrations

var_conc_rep_med <- function(df0A, df0B, df4A, df4B) {
  join_0 <- inner_join(df0A, df0B, 
                       by = c("name", "most_common", "barcodes_DNA"), 
                       suffix = c("_0A", "_0B"))
  join_4 <- inner_join(df4A, df4B, 
                         by = c("name", "most_common", "barcodes_DNA"), 
                         suffix = c("_4A", "_4B"))
  join_0_4 <- inner_join(join_0, join_4, 
                           by = c("name", "most_common", "barcodes_DNA")) %>%
    ungroup()
  print('processed dfs in order of samples: 0A, 0B, 4A, 4B')
  return(join_0_4)
}

med_rep_0_4_A_B <- var_conc_rep_med(med_ratio_R0A, med_ratio_R0B,
                                    med_ratio_R4A, med_ratio_R4B)

#Separate into subpools---------------------------------------------------------

#The scramble library consists of 10 bp scrambles tiling the 93 generated 
#backgrounds and 3 original backgrounds in 5 bp intervals. The distance measures
#the distance of the scramble to the minimal promoter. Similarity measures the 
#similarity of the scramble to the original background on a per position basis. 
#Two datasets used for the different backgrounds sets.

scramble_newback <- med_rep_0_4_A_B %>%
  filter(grepl('scramble', name)) %>%
  filter(grepl('background', name)) %>%
  separate(name, 
           into = c("subpool", "fluff1", "dist", "fluff2", "similarity", 
                    "fluff3", "background"),
           sep = "_", convert = TRUE) %>%
  select(-subpool, -fluff1, -fluff2, -fluff3) %>%
  mutate(dist = dist + 64)

scramble_oldback <- med_rep_0_4_A_B %>%
  filter(grepl('scramble', name)) %>%
  filter(grepl('oldback', name)) %>%
  separate(name, 
           into = c("subpool", "fluff1", "dist", "fluff2", "similarity", 
                    "fluff3", "background"),
           sep = "_", convert = TRUE) %>%
  select(-subpool, -fluff1, -fluff2, -fluff3) %>%
  mutate(dist = dist + 64)

#The sixsite library corresponds to the CRE Number Library across a wide range
#of gc contents. This library contains 6 equally spaced sites with 17 bp CRE 
#Spacing. Per variant, each site is one of: the consensus CRE, or no CRE. Here 
#site 1, 2, 3, 4, 5, and 6 equate to -191, -166, -141, -116, -91 and -66 site 
#distances to the downstream promoter. Separation lists identity of CRE per 
#site, the total concensus CRE sites per variants. Two datasets used for the 
#different backgrounds sets.

sixsite_newback <- med_rep_0_4_A_B %>%
  filter(grepl('sixsite', name))  %>%
  filter(grepl('background', name)) %>%
  separate(name, 
           into = c("subpool", "site1", "site2", "site3", "site4", "site5",
                    "site6", "fluff1", "background", "fluff2", "gc"),
           sep = "_", convert = TRUE) %>%
  select(-subpool, -fluff1, -fluff2) %>%
  mutate(consensus = str_detect(site1, "consensus") + 
           str_detect(site2, "consensus") + 
           str_detect(site3, "consensus") + 
           str_detect(site4, "consensus") + 
           str_detect(site5, "consensus") + 
           str_detect(site6, "consensus"))

sixsite_oldback <- med_rep_0_4_A_B %>%
  filter(grepl('sixsite', name))  %>%
  filter(grepl('oldback', name)) %>%
  separate(name, 
           into = c("subpool", "site1", "site2", "site3", "site4", "site5",
                    "site6", "fluff1", "background", "fluff2", "gc"),
           sep = "_", convert = TRUE) %>%
  select(-subpool, -fluff1, -fluff2) %>%
  mutate(consensus = str_detect(site1, "consensus") + 
           str_detect(site2, "consensus") + 
           str_detect(site3, "consensus") + 
           str_detect(site4, "consensus") + 
           str_detect(site5, "consensus") + 
           str_detect(site6, "consensus"))

#The twosite library corresponds to further spacings tested with the CRE Spacing
#and Distance Library. This library contains 2 consensus CRE sites with flanks 
#(ATTGACGTCAGC) that vary in CRE Spacing from one another by 1-13. Both CREs are
#then placed at variable distances to the 3' end of the backgrounds. Separation 
#lists the CRE spacing between sites and CRE distance. Distances measured from 
#the end of the background to the CRE proximal to the promoter. Added 64 bp to 
#measure to the minimal promoter. Only the 3 original backgrounds were used 
#here.

twosite <- med_rep_0_4_A_B %>%
  filter(grepl('twosite_', name)) %>%
  separate(name, 
           into = c("subpool", "fluff1", "spacing", "fluff2", "dist", 
                    "fluff3", "background"),
           sep = "_", convert = TRUE) %>%
  select(-subpool, -fluff1, -fluff2, -fluff3) %>%
  mutate(dist = dist + 64)


#Background-normalize subpools--------------------------------------------------

#Per each subpool separate out background values for each condition and left 
#join to original dataset. Normalize expression of each variant to its 
#background in that condition. Determine average expression and average 
#background-normalized expression across biological replicates.

newback <- sixsite_newback %>%
  filter(consensus == 0) %>%
  select(background, med_ratio_0A, med_ratio_0B, med_ratio_4A, 
         med_ratio_4B) %>%
  rename(med_ratio_0A_back = med_ratio_0A) %>%
  rename(med_ratio_0B_back = med_ratio_0B) %>%
  rename(med_ratio_4A_back = med_ratio_4A) %>%
  rename(med_ratio_4B_back = med_ratio_4B)

oldback <- sixsite_oldback %>%
  filter(consensus == 0) %>%
  select(background, med_ratio_0A, med_ratio_0B, med_ratio_4A, 
         med_ratio_4B) %>%
  rename(med_ratio_0A_back = med_ratio_0A) %>%
  rename(med_ratio_0B_back = med_ratio_0B) %>%
  rename(med_ratio_4A_back = med_ratio_4A) %>%
  rename(med_ratio_4B_back = med_ratio_4B)

back_norm <- function(backgrounds, df1) {
  back_join_norm <- left_join(df1, backgrounds, 
                              by = 'background') %>%
    mutate(ave_med_ratio_0 = (med_ratio_0A + med_ratio_0B)/2) %>%
    mutate(ave_med_ratio_4 = (med_ratio_4A + med_ratio_4B)/2) %>%
    mutate(med_ratio_0A_norm = med_ratio_0A/med_ratio_0A_back) %>%
    mutate(med_ratio_0B_norm = med_ratio_0B/med_ratio_0B_back) %>%
    mutate(med_ratio_4A_norm = med_ratio_4A/med_ratio_4A_back) %>%
    mutate(med_ratio_4B_norm = med_ratio_4B/med_ratio_4B_back) %>%
    mutate(ave_med_ratio_norm_0 = (med_ratio_0A_norm + med_ratio_0B_norm)/2) %>%
    mutate(ave_med_ratio_norm_4 = (med_ratio_4A_norm + med_ratio_4B_norm)/2)
}

scramble_newback_norm <- back_norm(newback, scramble_newback)
scramble_oldback_norm <- back_norm(oldback, scramble_oldback)
sixsite_newback_norm <- back_norm(newback, sixsite_newback)
sixsite_oldback_norm <- back_norm(oldback, sixsite_oldback)
twosite_norm <- back_norm(oldback, twosite)

#Export dfs---------------------------------------------------------------------

#Export dfs to perform main analysis in 20170921_epi_lib

med_rep_0_4_A_B %>% write.table("med_rep_follow_up.txt", sep = '\t', row.names =
                                  FALSE)

scramble_newback_norm %>%
  write.table(
    "scramble_newback_norm.txt", 
    sep = '\t', row.names = FALSE)

scramble_oldback_norm %>%
  write.table(
    "scramble_oldback_norm.txt", 
    sep = '\t', row.names = FALSE)

sixsite_newback_norm %>%
  write.table(
    "sixsite_newback_norm.txt", 
    sep = '\t', row.names = FALSE)

sixsite_oldback_norm %>%
  write.table(
    "sixsite_oldback_norm.txt", 
    sep = '\t', row.names = FALSE)

twosite_norm %>%
  write.table(
    "twosite_norm.txt", 
    sep = '\t', row.names = FALSE)



#Determining generated background motif activity--------------------------------

#Visualize effects of scrambles across different backgrounds. Plot expression of
#each scramble normalized to the backgrounds in both uninduced and induced
#conditions

p_newback_scramble_4 <- scramble_newback_norm %>%
  ggplot(aes(background, ave_med_ratio_norm_4, group = background)) +
  geom_boxplot(size = 0.3, outlier.size = 0.5) + 
  geom_hline(yintercept = 1, color = 'red') +
  ylab('Average normalized\nexpression 4 µM forskolin') +
  figurefont_theme

ggsave('../plots/p_newback_scramble_4.png', p_newback_scramble_4,
       width = 5, height = 2)

p_newback_scramble_0 <- scramble_newback_norm %>%
  ggplot(aes(background, ave_med_ratio_norm_0, group = background)) +
  geom_boxplot(size = 0.3, outlier.size = 0.5) + 
  geom_hline(yintercept = 1, color = 'red') +
  ylab('Average normalized\nexpression 0 µM forskolin') +
  figurefont_theme

ggsave('../plots/p_newback_scramble_0.png', p_newback_scramble_0,
       width = 5, height = 2)

p_oldback_scramble <- scramble_oldback_norm %>%
  ggplot(aes(background, ave_med_ratio_norm_4, group = background)) +
  geom_boxplot() + 
  geom_hline(yintercept = 1, color = 'red')

#output BC dfs of scrambled background expressions and background alone 
#expressions for t.test analysis of signigicant differences

bc_scramble_sep <- function(df) {
  df <- df %>%
    filter(grepl('scramble', name)) %>%
    filter(grepl('background', name)) %>%
    mutate(fluff = name) %>%
    separate(fluff, 
             into = c("subpool", "fluff1", "dist", "fluff2", "similarity", 
                      "fluff3", "background"),
             sep = "_", convert = TRUE) %>%
    select(barcode, name, background, ratio)
}

scramble_0A <- bc_scramble_sep(bc_DNA_RNA_0A) %>%
  mutate(sample = '0A')
scramble_0B <- bc_scramble_sep(bc_DNA_RNA_0B) %>%
  mutate(sample = '0B')
scramble_4A <- bc_scramble_sep(bc_DNA_RNA_4A) %>%
  mutate(sample = '4A')
scramble_4B <- bc_scramble_sep(bc_DNA_RNA_4B) %>%
  mutate(sample = '4B')
scrambles <- rbind(scramble_0A, scramble_0B, scramble_4A, scramble_4B)

bc_backgrounds_sep <- function(df) {
  df <- df %>%
    filter(grepl(
      'sixsite_nosite_nosite_nosite_nosite_nosite_nosite',
      name)) %>%
    filter(grepl('background', name)) %>%
    mutate(fluff = name) %>%
    separate(fluff, 
             into = c("subpool", "site1", "site2", "site3", "site4", "site5",
                      "site6", "fluff1", "background", "fluff2", "gc"),
             sep = "_", convert = TRUE) %>%
    select(barcode, name, background, ratio)
}

backgrounds_0A <- bc_backgrounds_sep(bc_DNA_RNA_0A) %>%
  mutate(sample = '0A')
backgrounds_0B <- bc_backgrounds_sep(bc_DNA_RNA_0B) %>%
  mutate(sample = '0B')
backgrounds_4A <- bc_backgrounds_sep(bc_DNA_RNA_4A) %>%
  mutate(sample = '4A')
backgrounds_4B <- bc_backgrounds_sep(bc_DNA_RNA_4B) %>%
  mutate(sample = '4B')
backgrounds <- rbind(backgrounds_0A, backgrounds_0B, 
                     backgrounds_4A, backgrounds_4B)

#Do same process for old backgrounds

bc_scramble_sep_old <- function(df) {
  df <- df %>%
    filter(grepl('scramble', name)) %>%
    filter(grepl('oldback', name)) %>%
    mutate(fluff = name) %>%
    separate(fluff, 
             into = c("subpool", "fluff1", "dist", "fluff2", "similarity", 
                      "fluff3", "background"),
             sep = "_", convert = TRUE) %>%
    select(barcode, name, background, ratio)
}

scramble_0A_old <- bc_scramble_sep_old(bc_DNA_RNA_0A) %>%
  mutate(sample = '0A')
scramble_0B_old <- bc_scramble_sep_old(bc_DNA_RNA_0B) %>%
  mutate(sample = '0B')
scramble_4A_old <- bc_scramble_sep_old(bc_DNA_RNA_4A) %>%
  mutate(sample = '4A')
scramble_4B_old <- bc_scramble_sep_old(bc_DNA_RNA_4B) %>%
  mutate(sample = '4B')
scrambles_old <- rbind(scramble_0A_old, scramble_0B_old, 
                       scramble_4A_old, scramble_4B_old)

bc_backgrounds_sep_old <- function(df) {
  df <- df %>%
    filter(grepl(
      'sixsite_nosite_nosite_nosite_nosite_nosite_nosite',
      name)) %>%
    filter(grepl('oldback', name)) %>%
    mutate(fluff = name) %>%
    separate(fluff, 
             into = c("subpool", "site1", "site2", "site3", "site4", "site5",
                      "site6", "fluff1", "background", "fluff2", "gc"),
             sep = "_", convert = TRUE) %>%
    select(barcode, name, background, ratio)
}

backgrounds_0A_old <- bc_backgrounds_sep_old(bc_DNA_RNA_0A) %>%
  mutate(sample = '0A')
backgrounds_0B_old <- bc_backgrounds_sep_old(bc_DNA_RNA_0B) %>%
  mutate(sample = '0B')
backgrounds_4A_old <- bc_backgrounds_sep_old(bc_DNA_RNA_4A) %>%
  mutate(sample = '4A')
backgrounds_4B_old <- bc_backgrounds_sep_old(bc_DNA_RNA_4B) %>%
  mutate(sample = '4B')
backgrounds_old <- rbind(backgrounds_0A_old, backgrounds_0B_old, 
                         backgrounds_4A_old, backgrounds_4B_old)

#Determine significant scrambles based on the difference in the distribution of 
#barcode expression ratios between scramble and background

#this part was written by Kim

library(dplyr)
library(tidyr)
library(purrr)

options(stringsAsFactors = F)
options(scipen = 10000)

# assign category to make t-test easier

backgrounds$category <- 'background'
scrambles$category <- 'scramble'

backgrounds_old$category <- 'background'
scrambles_old$category <- 'scramble'

#For each scramble we want to grab all barcodes (for a given sample) and compare
#the mean expression for the barcoded scrambles to the mean expression for the 
#barcoded background.

ttest_cre_custom <- function(df, df_bg) {
  bg_name <- df$background[1]
  sample_name <- df$sample[1]
  # bind unscrambled barcodes to df
  df_with_bg <- bind_rows(df,
                          filter(df_bg, 
                                 background == bg_name,
                                 sample == sample_name)) %>% 
    mutate(category_fctr = factor(category))
  
  result <- tryCatch(
    {
      t.test(ratio ~ category_fctr, df_with_bg)
    }, warning = function(cond) {
      return(NA)
    }, error = function(cond) {
      return(NA)
    }
  )
  return(result)
}

# test <- filter(scrambles, name == 'scramble_dist_60_similarity_0_background_78')

scramble_ttests <- scrambles %>% 
  group_by(name, sample) %>% 
  do(ttest = ttest_cre_custom(., df_bg = backgrounds)) %>% 
  broom::tidy(ttest) %>% 
  ungroup()

scramble_ttests_old <- scrambles_old %>%
  group_by(name, sample) %>% 
  do(ttest = ttest_cre_custom(., df_bg = backgrounds_old)) %>% 
  broom::tidy(ttest) %>% 
  ungroup()

scramble_ttests <- scramble_ttests %>% 
  select(name, sample, mean_diff = estimate, mean_bg = estimate1,
         mean_scramble = estimate2, tstat = statistic, p.value,
         conf.low, conf.high)

scramble_ttests_old <- scramble_ttests_old %>% 
  select(name, sample, mean_diff = estimate, mean_bg = estimate1,
         mean_scramble = estimate2, tstat = statistic, p.value,
         conf.low, conf.high)

scramble_ttests$p.value.fdr <- p.adjust(scramble_ttests$p.value, method = 'fdr')

scramble_ttests_old$p.value.fdr <- p.adjust(scramble_ttests_old$p.value, 
                                            method = 'fdr')

#Filter for adjusted p-value < 0.05. This corresponds to a false discovery rate 
#of 5%, meaning that for all significant results 5% will be false positive.

scramble_signif <- filter(scramble_ttests, p.value.fdr <= 0.05) %>%
  mutate(background = str_sub(name, nchar(name)-1, nchar(name))) %>%
  mutate(background = as.integer(background))

bgs_signif <- tibble(background = c(unique(scramble_signif$background)))

removed_bgs_sixsite_newback_norm <- anti_join(sixsite_newback_norm, 
                                              bgs_signif, by = 'background')

removed_bgs_sixsite_newback_norm %>%
  write.table(
    "removed_bgs_sixsite_newback_norm.txt", 
    sep = '\t', row.names = FALSE)


#Check if scramble is significant for old backgrounds. Returned one scramble for 
#background 52 but I checked this using the single CRE distance library and 
#there is no obvious dip in this portion of the background

scramble_signif_old <- filter(scramble_ttests, p.value.fdr <= 0.05) %>%
  mutate(background = str_sub(name, nchar(name)-1, nchar(name))) %>%
  mutate(background = as.integer(background))
  






