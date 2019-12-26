#Written for the analysis of a range of inductions in the episomal library
#DNA_BC: DNA
#R0A_BC and R0B_BC: RNA at 0 µM Forsk replicate A and B
#R2#A_BC and R2#B_BC: RNA at 2^# µM Forsk replicate A and B, includes negative 
#and positive #'s

#Tested concentrations: 0, 2^-5, 2^-4, 2^-3, 2^-2, 2^-1, 2^0, 2^2 µM forskolin

#All figures made with this dataset require sections "Load index and bcmap 
#files" through "Median BC expression".

#Most figures require section "Background-normalize expression" and it is
#recommended to run this section before plotting any figure

#Most figures require both genomic and episomal MPRAs, the genomic MPRA data
#processing is performed in 20171129_genlib.R and imported here in section
#"Import genomic MPRA and combine with episomal". Other episomal MPRAs used in
#supplemental figures are performed in 20170320_epilib_analysis_0_25.R, 
#20170631_epilib_analysis_0_64.R, and 20190712_eplib_analysis.R and are 
#imported per main and supplemental figure section.

#Beyond figure 1 and associated supplemental figures, the section "Separate into
#sublibraries" is required for figure generation.

#Establish workspace------------------------------------------------------------

options(stringsAsFactors = F)
set.seed(1234)

#import necessary libraries

library(tidyverse)
library(lemon)
library(viridis)
library(cowplot)
library(caTools)
library(broom)
library(modelr)
library(reshape2)
library(GGally)

#General figure customizations

cbPalette7 <- c('#440154FF', '#39568CFF', '#287D8EFF', '#20A387FF', '#73D055FF',
                '#B8DE29FF', '#FDE725FF')

cbPalette7_grad_light <- c('white', '#FDE725FF', '#B8DE29FF', '#55C667FF', 
                           '#1F968BFF', '#39568CFF', '#482677FF')

spacing_5_20_cbpalette <- c('gray20', 'deepskyblue2', 'orangered3', 'sandybrown')

cbPalette_cont_blue <- c('gray100', 'lightskyblue', 'dodgerblue1', 'mediumblue', 'navy')

background_cbPalette <- c('#E4CA2C', '#E46E2C', '#E42C46')

figurefont_theme <- theme(text = element_text(size = 8)) +
  theme(axis.title = element_text(size = 8)) +
  theme(legend.title = element_text(size = 8)) +
  theme(axis.text = element_text(size = 8))

#General functions

var_log10 <- function(df) {
  log_ratio_df <- df %>% 
    mutate_if(is.double, funs(log10(.)))
  return(log_ratio_df)
}

#3bp moving window function used in most CRE Spacing and Distance library 
#figures

moveavg_dist3 <- function(df) {
  df <- df %>%
    mutate(ave_3 = runmean(ave_ratio_norm, 3, alg = 'R', endrule = 'NA'))
}

#Load index and bcmap files-----------------------------------------------------

#Load indexed BC reads

bc_DNA <- read_tsv('BCreads_txts/20170921_DNA_BC.txt')
bc_R0A <- read_tsv('BCreads_txts/20170921_R0A_BC.txt')
bc_R0B <- read_tsv('BCreads_txts/20170921_R0B_BC.txt')
bc_R2_5A <- read_tsv('BCreads_txts/20170921_R2-5A_BC.txt')
bc_R2_5B <- read_tsv('BCreads_txts/20170921_R2-5B_BC.txt')
bc_R2_4A <- read_tsv('BCreads_txts/20170921_R2-4A_BC.txt')
bc_R2_4B <- read_tsv('BCreads_txts/20170921_R2-4B_BC.txt')
bc_R2_3A <- read_tsv('BCreads_txts/20170921_R2-3A_BC.txt')
bc_R2_3B <- read_tsv('BCreads_txts/20170921_R2-3B_BC.txt')
bc_R2_2A <- read_tsv('BCreads_txts/20170921_R2-2A_BC.txt')
bc_R2_2B <- read_tsv('BCreads_txts/20170921_R2-2B_BC.txt')
bc_R2_1A <- read_tsv('BCreads_txts/20170921_R2-1A_BC.txt')
bc_R2_1B <- read_tsv('BCreads_txts/20170921_R2-1B_BC.txt')
bc_R20A <- read_tsv('BCreads_txts/20170921_R20A_BC.txt')
bc_R20B <- read_tsv('BCreads_txts/20170921_R20B_BC.txt')
bc_R22A <- read_tsv('BCreads_txts/20170921_R22A_BC.txt')
bc_R22B <- read_tsv('BCreads_txts/20170921_R22B_BC.txt')

#Load barcode mapping table, sequences (most_common) are rcomp due to sequencing
#format. Pick out controls, subpool 3, and subpool 5 in the bcmap that were used
#in this assay

SP3_SP5_map <- read_tsv('../BCMap/uniqueSP2345.txt', 
                        col_names = c(
                          'fluff', 'barcode', 'name', 'most_common'), 
                        skip = 1) %>%
  select(-fluff) %>%
  mutate(subpool = ifelse(startsWith(name, 'subpool'), 
                          substr(name, 1, 8), 
                          'control')) %>%
  filter(subpool != 'subpool2' & subpool != 'subpool4')

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

bc_join_DNA <- bc_map_join_bc(SP3_SP5_map, bc_DNA)
bc_join_R0A <- bc_map_join_bc(SP3_SP5_map, bc_R0A)
bc_join_R0B <- bc_map_join_bc(SP3_SP5_map, bc_R0B)
bc_join_R2_5A <- bc_map_join_bc(SP3_SP5_map, bc_R2_5A)
bc_join_R2_5B <- bc_map_join_bc(SP3_SP5_map, bc_R2_5B)
bc_join_R2_4A <- bc_map_join_bc(SP3_SP5_map, bc_R2_4A)
bc_join_R2_4B <- bc_map_join_bc(SP3_SP5_map, bc_R2_4B)
bc_join_R2_3A <- bc_map_join_bc(SP3_SP5_map, bc_R2_3A)
bc_join_R2_3B <- bc_map_join_bc(SP3_SP5_map, bc_R2_3B)
bc_join_R2_2A <- bc_map_join_bc(SP3_SP5_map, bc_R2_2A)
bc_join_R2_2B <- bc_map_join_bc(SP3_SP5_map, bc_R2_2B)
bc_join_R2_1A <- bc_map_join_bc(SP3_SP5_map, bc_R2_1A)
bc_join_R2_1B <- bc_map_join_bc(SP3_SP5_map, bc_R2_1B)
bc_join_R20A <- bc_map_join_bc(SP3_SP5_map, bc_R20A)
bc_join_R20B <- bc_map_join_bc(SP3_SP5_map, bc_R20B)
bc_join_R22A <- bc_map_join_bc(SP3_SP5_map, bc_R22A)
bc_join_R22B <- bc_map_join_bc(SP3_SP5_map, bc_R22B)

#Median BC expression-----------------------------------------------------------

#Retain BCs with > 6 reads in DNA sample, left-join RNA to DNA BCs. Take ratio 
#of RNA/DNA normalized reads per million

dna7_join_rna_rep <- function(df1, df2) {
  filter_DNA <- filter(df1, num_reads > 6)
  DNA_RNA_join <- left_join(filter_DNA, df2,
                            by = c("barcode", "name", "subpool", "most_common"), 
                            suffix = c('_DNA', '_RNA')) %>%
    mutate(ratio = normalized_RNA/normalized_DNA)
  return(DNA_RNA_join)
}

bc_DNA_RNA_0A <- dna7_join_rna_rep(bc_join_DNA, bc_join_R0A)
bc_DNA_RNA_0B <- dna7_join_rna_rep(bc_join_DNA, bc_join_R0B)
bc_DNA_RNA_2_5A <- dna7_join_rna_rep(bc_join_DNA, bc_join_R2_5A)
bc_DNA_RNA_2_5B <- dna7_join_rna_rep(bc_join_DNA, bc_join_R2_5B)
bc_DNA_RNA_2_4A <- dna7_join_rna_rep(bc_join_DNA, bc_join_R2_4A)
bc_DNA_RNA_2_4B <- dna7_join_rna_rep(bc_join_DNA, bc_join_R2_4B)
bc_DNA_RNA_2_3A <- dna7_join_rna_rep(bc_join_DNA, bc_join_R2_3A)
bc_DNA_RNA_2_3B <- dna7_join_rna_rep(bc_join_DNA, bc_join_R2_3B)
bc_DNA_RNA_2_2A <- dna7_join_rna_rep(bc_join_DNA, bc_join_R2_2A)
bc_DNA_RNA_2_2B <- dna7_join_rna_rep(bc_join_DNA, bc_join_R2_2B)
bc_DNA_RNA_2_1A <- dna7_join_rna_rep(bc_join_DNA, bc_join_R2_1A)
bc_DNA_RNA_2_1B <- dna7_join_rna_rep(bc_join_DNA, bc_join_R2_1B)
bc_DNA_RNA_20A <- dna7_join_rna_rep(bc_join_DNA, bc_join_R20A)
bc_DNA_RNA_20B <- dna7_join_rna_rep(bc_join_DNA, bc_join_R20B)
bc_DNA_RNA_22A <- dna7_join_rna_rep(bc_join_DNA, bc_join_R22A)
bc_DNA_RNA_22B <- dna7_join_rna_rep(bc_join_DNA, bc_join_R22B)


#Count barcodes per variant per DNA and RNA sample, set minimum of 8 BC's per 
#variant in DNA sample, take median BC RNA/DNA per variant, then per variant 
#determine the median absolute deviation of all barcode ratios. Only look at 
#variants with greater than 0 median expression

ratio_bc_med_var <- function(df) {
  bc_count_DNA <- df %>%
    group_by(subpool, name, most_common) %>%
    summarize(barcodes_DNA = n()) %>%
    filter(barcodes_DNA > 7)
  bc_count_RNA <- df %>%
    group_by(subpool, name, most_common) %>%
    filter(num_reads_RNA != 0) %>%
    summarize(barcodes_RNA = n())
  bc_DNA_RNA <- left_join(bc_count_DNA, bc_count_RNA, 
                          by = c('subpool', 'name', 'most_common')) %>%
    ungroup()
  bc_min_8_df <- left_join(bc_DNA_RNA, df, 
                           by = c('subpool', 'name', 'most_common')) %>%
    ungroup()
  med_ratio <- bc_min_8_df %>%
    group_by(subpool, name, most_common) %>%
    summarize(med_ratio = median(ratio)) %>%
    filter(med_ratio > 0)
  mad_ratio <- bc_min_8_df %>%
    group_by(subpool, name, most_common) %>%
    summarize(mad = mad(ratio, constant = 1))
  med_mad <- left_join(med_ratio, mad_ratio, 
                        by = c('subpool', 'name', 'most_common')) %>%
    mutate(mad_over_med = as.double(mad/med_ratio))
  bc_med <- inner_join(med_mad, bc_DNA_RNA, 
                       by = c('subpool', 'name', 'most_common')) %>%
    ungroup()
  return(bc_med)
}

med_ratio_R0A <- ratio_bc_med_var(bc_DNA_RNA_0A)
med_ratio_R0B <- ratio_bc_med_var(bc_DNA_RNA_0B)
med_ratio_R2_5A <- ratio_bc_med_var(bc_DNA_RNA_2_5A)
med_ratio_R2_5B <- ratio_bc_med_var(bc_DNA_RNA_2_5B)
med_ratio_R2_4A <- ratio_bc_med_var(bc_DNA_RNA_2_4A)
med_ratio_R2_4B <- ratio_bc_med_var(bc_DNA_RNA_2_4B)
med_ratio_R2_3A <- ratio_bc_med_var(bc_DNA_RNA_2_3A)
med_ratio_R2_3B <- ratio_bc_med_var(bc_DNA_RNA_2_3B)
med_ratio_R2_2A <- ratio_bc_med_var(bc_DNA_RNA_2_2A)
med_ratio_R2_2B <- ratio_bc_med_var(bc_DNA_RNA_2_2B)
med_ratio_R2_1A <- ratio_bc_med_var(bc_DNA_RNA_2_1A)
med_ratio_R2_1B <- ratio_bc_med_var(bc_DNA_RNA_2_1B)
med_ratio_R20A <- ratio_bc_med_var(bc_DNA_RNA_20A)
med_ratio_R20B <- ratio_bc_med_var(bc_DNA_RNA_20B)
med_ratio_R22A <- ratio_bc_med_var(bc_DNA_RNA_22A)
med_ratio_R22B <- ratio_bc_med_var(bc_DNA_RNA_22B)

#Combine biological replicates across concentrations

var_conc_rep_med <- function(df0A, df0B, df2_5A, df2_5B, df2_4A, df2_4B, df2_3A, 
                             df2_3B, df2_2A, df2_2B, df2_1A, df2_1B, df20A, df20B, 
                             df22A, df22B) {
  join_0 <- inner_join(df0A, df0B, 
                       by = c("name", "subpool", "most_common", "barcodes_DNA"), 
                       suffix = c("_0A", "_0B"))
  join_2_5 <- inner_join(df2_5A, df2_5B, 
                         by = c("name", "subpool", "most_common", 
                                "barcodes_DNA"), 
                         suffix = c("_2_5A", "_2_5B"))
  join_2_4 <- inner_join(df2_4A, df2_4B, 
                         by = c("name", "subpool", "most_common", 
                                "barcodes_DNA"), 
                         suffix = c("_2_4A", "_2_4B"))
  join_2_3 <- inner_join(df2_3A, df2_3B, 
                         by = c("name", "subpool", "most_common", 
                                "barcodes_DNA"), 
                         suffix = c("_2_3A", "_2_3B"))
  join_2_2 <- inner_join(df2_2A, df2_2B, 
                         by = c("name", "subpool", "most_common", 
                                "barcodes_DNA"), 
                         suffix = c("_2_2A", "_2_2B"))
  join_2_1 <- inner_join(df2_1A, df2_1B, 
                         by = c("name", "subpool", "most_common", 
                                "barcodes_DNA"), 
                         suffix = c("_2_1A", "_2_1B"))
  join_20 <- inner_join(df20A, df20B, 
                        by = c("name", "subpool", "most_common", 
                               "barcodes_DNA"), 
                        suffix = c("_20A", "_20B"))
  join_22 <- inner_join(df22A, df22B, 
                        by = c("name", "subpool", "most_common", 
                               "barcodes_DNA"), 
                        suffix = c("_22A", "_22B"))
  join_0_2_5 <- inner_join(join_0, join_2_5, 
                           by = c("name", "subpool", "most_common", 
                                  "barcodes_DNA"))
  join_0_2_4 <- inner_join(join_0_2_5, join_2_4, 
                           by = c("name", "subpool", "most_common", 
                                  "barcodes_DNA"))
  join_0_2_3 <- inner_join(join_0_2_4, join_2_3, 
                           by = c("name", "subpool", "most_common", 
                                  "barcodes_DNA"))
  join_0_2_2 <- inner_join(join_0_2_3, join_2_2, 
                           by = c("name", "subpool", "most_common", 
                                  "barcodes_DNA"))
  join_0_2_1 <- inner_join(join_0_2_2, join_2_1, 
                           by = c("name", "subpool", "most_common", 
                                  "barcodes_DNA"))
  join_0_20 <- inner_join(join_0_2_1, join_20, 
                          by = c("name", "subpool", "most_common", 
                                 "barcodes_DNA"))
  join_0_22 <- inner_join(join_0_20, join_22, 
                          by = c("name", "subpool", "most_common", 
                                 "barcodes_DNA")) %>%
    ungroup()
  print('processed dfs in order of samples: 0A, 0B, 2_5A, 2_5B, 2_4A, 2_4B, 
        2_3A, 2_3B, 2_2A, 2_2B, 2_1A, 2_1B, 20A, 20B, 22A, 22B')
  return(join_0_22)
}

med_rep_0_22_A_B <- var_conc_rep_med(med_ratio_R0A, med_ratio_R0B, 
                                     med_ratio_R2_5A, med_ratio_R2_5B,
                                     med_ratio_R2_4A, med_ratio_R2_4B,
                                     med_ratio_R2_3A, med_ratio_R2_3B,
                                     med_ratio_R2_2A, med_ratio_R2_2B,
                                     med_ratio_R2_1A, med_ratio_R2_1B,
                                     med_ratio_R20A, med_ratio_R20B,
                                     med_ratio_R22A, med_ratio_R22B)

#Background-normalize expression------------------------------------------------

#Normalize MPRA expression of all variants to variant without CRE sites 
#(backgrounds). Each variant normalized to the background used in its design.

back_norm <- function(df1) {
  gsub_0_22 <- df1 %>%
    ungroup () %>%
    filter(subpool != 'control') %>%
    mutate(
      name = gsub('Smith R. Vista chr9:83712599-83712766', 'back_41', name),
      name = gsub('Vista Chr5:88673410-88674494', 'back_52', name),
      name = gsub('scramble pGL4.29 Promega 1-63 \\+ 1-87', 'back_55', name)
    ) %>%
    mutate(background = name) %>%
    mutate(background = str_sub(background, 
                                nchar(background)-1, 
                                nchar(background)))
  backgrounds <- gsub_0_22 %>%
    filter(startsWith(name, 
                      'subpool5_no_site_no_site_no_site_no_site_no_site_no_site')) %>%
    select(background, med_ratio_0A, med_ratio_0B, med_ratio_2_5A, 
           med_ratio_2_5B, med_ratio_2_4A, med_ratio_2_4B, med_ratio_2_3A, 
           med_ratio_2_3B, med_ratio_2_2A, med_ratio_2_2B, med_ratio_2_1A, 
           med_ratio_2_1B, med_ratio_20A, med_ratio_20B, med_ratio_22A, 
           med_ratio_22B) %>%
    rename(ratio_0A_back = med_ratio_0A) %>%
    rename(ratio_0B_back = med_ratio_0B) %>%
    rename(ratio_2_5A_back = med_ratio_2_5A) %>%
    rename(ratio_2_5B_back = med_ratio_2_5B) %>%
    rename(ratio_2_4A_back = med_ratio_2_4A) %>%
    rename(ratio_2_4B_back = med_ratio_2_4B) %>%
    rename(ratio_2_3A_back = med_ratio_2_3A) %>%
    rename(ratio_2_3B_back = med_ratio_2_3B) %>%
    rename(ratio_2_2A_back = med_ratio_2_2A) %>%
    rename(ratio_2_2B_back = med_ratio_2_2B) %>%
    rename(ratio_2_1A_back = med_ratio_2_1A) %>%
    rename(ratio_2_1B_back = med_ratio_2_1B) %>%
    rename(ratio_20A_back = med_ratio_20A) %>%
    rename(ratio_20B_back = med_ratio_20B) %>%
    rename(ratio_22A_back = med_ratio_22A) %>%
    rename(ratio_22B_back = med_ratio_22B)
  back_join_norm <- left_join(gsub_0_22, backgrounds, by = 'background') %>%
    mutate(ratio_0A_norm = med_ratio_0A/ratio_0A_back) %>%
    mutate(ratio_0B_norm = med_ratio_0B/ratio_0B_back) %>%
    mutate(ratio_2_5A_norm = med_ratio_2_5A/ratio_2_5A_back) %>%
    mutate(ratio_2_5B_norm = med_ratio_2_5B/ratio_2_5B_back) %>%
    mutate(ratio_2_4A_norm = med_ratio_2_4A/ratio_2_4A_back) %>%
    mutate(ratio_2_4B_norm = med_ratio_2_4B/ratio_2_4B_back) %>%
    mutate(ratio_2_3A_norm = med_ratio_2_3A/ratio_2_3A_back) %>%
    mutate(ratio_2_3B_norm = med_ratio_2_3B/ratio_2_3B_back) %>%
    mutate(ratio_2_2A_norm = med_ratio_2_2A/ratio_2_2A_back) %>%
    mutate(ratio_2_2B_norm = med_ratio_2_2B/ratio_2_2B_back) %>%
    mutate(ratio_2_1A_norm = med_ratio_2_1A/ratio_2_1A_back) %>%
    mutate(ratio_2_1B_norm = med_ratio_2_1B/ratio_2_1B_back) %>%
    mutate(ratio_20A_norm = med_ratio_20A/ratio_20A_back) %>%
    mutate(ratio_20B_norm = med_ratio_20B/ratio_20B_back) %>%
    mutate(ratio_22A_norm = med_ratio_22A/ratio_22A_back) %>%
    mutate(ratio_22B_norm = med_ratio_22B/ratio_22B_back) %>%
    mutate(ave_ratio_0_norm = (ratio_0A_norm + ratio_0B_norm)/2) %>%
    mutate(ave_ratio_2_5_norm = (ratio_2_5A_norm + ratio_2_5B_norm)/2) %>%
    mutate(ave_ratio_2_4_norm = (ratio_2_4A_norm + ratio_2_4B_norm)/2) %>%
    mutate(ave_ratio_2_3_norm = (ratio_2_3A_norm + ratio_2_3B_norm)/2) %>%
    mutate(ave_ratio_2_2_norm = (ratio_2_2A_norm + ratio_2_2B_norm)/2) %>%
    mutate(ave_ratio_2_1_norm = (ratio_2_1A_norm + ratio_2_1B_norm)/2) %>%
    mutate(ave_ratio_20_norm = (ratio_20A_norm + ratio_20B_norm)/2) %>%
    mutate(ave_ratio_22_norm = (ratio_22A_norm + ratio_22B_norm)/2) %>%
    mutate(induction = ave_ratio_22_norm/ave_ratio_0_norm)
  return(back_join_norm)
}

#Add back to df CRE control that can be normalized to a background for plotting

epi_back_norm_pc_spGl4 <- med_rep_0_22_A_B %>%
  filter(name == 'pGL4.29 Promega 1-63 + 1-87') %>%
  mutate(name = str_c(name, '_scramble pGL4.29 Promega 1-63 + 1-87')) %>%
  mutate(subpool = 'subpool3') %>%
  rbind(med_rep_0_22_A_B) %>%
  back_norm()

#Make untidy df with conc as a variable and normalized expression as single 
#columns. Conc is indicated in log2 and 0 µM forskolin concentration is 
#represented as 2^-7 for plotting on log scales.

var_conc_exp <- function(df) {
  df_0 <- df %>%
    mutate(ave_barcode_0 = (barcodes_RNA_0A + barcodes_RNA_0B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_0, 
           ratio_0A_norm, ratio_0B_norm, ave_ratio_0_norm) %>%
    mutate(conc = -7) %>%
    rename(ratio_A_norm = ratio_0A_norm) %>%
    rename(ratio_B_norm = ratio_0B_norm) %>%
    rename(ave_ratio_norm = ave_ratio_0_norm) %>%
    rename(ave_barcode = ave_barcode_0)
  df_2_5 <- df %>%
    mutate(ave_barcode_2_5 = (barcodes_RNA_2_5A + barcodes_RNA_2_5B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_2_5, 
           ratio_2_5A_norm, ratio_2_5B_norm, ave_ratio_2_5_norm) %>%
    mutate(conc = -5) %>%
    rename(ratio_A_norm = ratio_2_5A_norm) %>%
    rename(ratio_B_norm = ratio_2_5B_norm) %>%
    rename(ave_ratio_norm = ave_ratio_2_5_norm) %>%
    rename(ave_barcode = ave_barcode_2_5)
  df_2_4 <- df %>%
    mutate(ave_barcode_2_4 = (barcodes_RNA_2_4A + barcodes_RNA_2_4B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_2_4, 
           ratio_2_4A_norm, ratio_2_4B_norm, ave_ratio_2_4_norm) %>%
    mutate(conc = -4) %>%
    rename(ratio_A_norm = ratio_2_4A_norm) %>%
    rename(ratio_B_norm = ratio_2_4B_norm) %>%
    rename(ave_ratio_norm = ave_ratio_2_4_norm) %>%
    rename(ave_barcode = ave_barcode_2_4)
  df_2_3 <- df %>%
    mutate(ave_barcode_2_3 = (barcodes_RNA_2_3A + barcodes_RNA_2_3B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_2_3, 
           ratio_2_3A_norm, ratio_2_3B_norm, ave_ratio_2_3_norm) %>%
    mutate(conc = -3) %>%
    rename(ratio_A_norm = ratio_2_3A_norm) %>%
    rename(ratio_B_norm = ratio_2_3B_norm) %>%
    rename(ave_ratio_norm = ave_ratio_2_3_norm) %>%
    rename(ave_barcode = ave_barcode_2_3)
  df_2_2 <- df %>%
    mutate(ave_barcode_2_2 = (barcodes_RNA_2_2A + barcodes_RNA_2_2B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_2_2, 
           ratio_2_2A_norm, ratio_2_2B_norm, ave_ratio_2_2_norm) %>%
    mutate(conc = -2) %>%
    rename(ratio_A_norm = ratio_2_2A_norm) %>%
    rename(ratio_B_norm = ratio_2_2B_norm) %>%
    rename(ave_ratio_norm = ave_ratio_2_2_norm) %>%
    rename(ave_barcode = ave_barcode_2_2)
  df_2_1 <- df %>%
    mutate(ave_barcode_2_1 = (barcodes_RNA_2_1A + barcodes_RNA_2_1B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_2_1, 
           ratio_2_1A_norm, ratio_2_1B_norm, ave_ratio_2_1_norm) %>%
    mutate(conc = -1) %>%
    rename(ratio_A_norm = ratio_2_1A_norm) %>%
    rename(ratio_B_norm = ratio_2_1B_norm) %>%
    rename(ave_ratio_norm = ave_ratio_2_1_norm) %>%
    rename(ave_barcode = ave_barcode_2_1)
  df_20 <- df %>%
    mutate(ave_barcode_20 = (barcodes_RNA_20A + barcodes_RNA_20B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_20, 
           ratio_20A_norm, ratio_20B_norm, ave_ratio_20_norm) %>%
    mutate(conc = 0) %>%
    rename(ratio_A_norm = ratio_20A_norm) %>%
    rename(ratio_B_norm = ratio_20B_norm) %>%
    rename(ave_ratio_norm = ave_ratio_20_norm) %>%
    rename(ave_barcode = ave_barcode_20)
  df_22 <- df %>%
    mutate(ave_barcode_22 = (barcodes_RNA_22A + barcodes_RNA_22B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_22, 
           ratio_22A_norm, ratio_22B_norm, ave_ratio_22_norm) %>%
    mutate(conc = 2) %>%
    rename(ratio_A_norm = ratio_22A_norm) %>%
    rename(ratio_B_norm = ratio_22B_norm) %>%
    rename(ave_ratio_norm = ave_ratio_22_norm) %>%
    rename(ave_barcode = ave_barcode_22)
  df_0_22 <- rbind(df_0, df_2_5, df_2_4, df_2_3, df_2_2, df_2_1, df_20, df_22)
  return(df_0_22)
}

#Episomal MPRA expression across tested forskolin concentrations

epi_back_norm_conc <- epi_back_norm_pc_spGl4 %>%
  var_conc_exp()

#Figure 1D----------------------------------------------------------------------

#Plot normalized variant expression curves across forskolin concentrations. 
#Expression curves for backgrounds and control overlayed

p_titr_pc_back <- epi_back_norm_conc %>%
  ggplot(aes(conc, ave_ratio_norm)) +
  geom_line(aes(group = name), alpha = 0.1) +
  geom_point(data = filter(epi_back_norm_conc, 
                           startsWith(name, 
                                      'subpool5_no_site_no_site_no_site_no_site_no_site_no_site')),
             color = 'darkgoldenrod1', shape = 19, stroke = 0.75) +
  geom_line(data = filter(epi_back_norm_conc, 
                          startsWith(name, 
                                     'subpool5_no_site_no_site_no_site_no_site_no_site_no_site')),
            color = 'darkgoldenrod1', size = 1) +
  geom_point(data = filter(epi_back_norm_conc, 
                           startsWith(name, 
                                      'pGL4.29 Promega 1-63 + 1-87')),
             color = 'firebrick2', shape = 19, stroke = 0.75) +
  geom_line(data = filter(epi_back_norm_conc, 
                          startsWith(name, 
                                     'pGL4.29 Promega 1-63 + 1-87')),
            color = 'firebrick2', size = 1) +
  ylab('Average normalized\nexpression (a.u.)') +
  annotation_logticks(sides = 'b', short = unit(0.05, 'cm'), 
                      mid = unit(0.1, 'cm'), long = unit(0.15, 'cm')) +
  scale_x_continuous(breaks = (-7:2), 'log2 forskolin (µM)') +
  figurefont_theme

ggsave('../plots/p_titr_pc_back.pdf', p_titr_pc_back, width = 2.8, height = 2,
       units = 'in')

ggsave('../plots/p_titr_pc_back.png', p_titr_pc_back, width = 2.8, height = 2,
       units = 'in')

#Import genomic MPRA and combine with episomal----------------------------------

#Genomic MPRA expression determination and data processing performed in
#"20171129_genlib_analysis/20171129_genlib.R". Processed df with expression
#values is exported from there and imported here to plot together

gen_rep_1_2 <- read_tsv('../20171129_genlib_analysis/rep_1_2.txt')

#Combine genomic and episomal dfs, only comparing to expression at 2^2 µM 
#forskolin in the episomal dataset. Here med_ratio_br# and ave_med_ratio refers 
#to genomic expression, either across biological replicates or averaged.
#Episomal MPRA columns are represented by the annotation <column>_22 as 
#representing the forskolin concentration the sample was incubated with (2^2)

gen_epi <- epi_back_norm_pc_spGl4 %>%
  mutate(ave_ratio_22 = (med_ratio_22A + med_ratio_22B)/2) %>%
  select(subpool, name, most_common, barcodes_DNA, med_ratio_22A,
         barcodes_RNA_22A, ratio_22A_norm, med_ratio_22B, barcodes_RNA_22B, 
         ratio_22B_norm, ave_ratio_22, ave_ratio_22_norm) %>%
  inner_join(gen_rep_1_2, by = c('subpool', 'name', 'most_common'))

#Make untidy df where MPRA format is a variable according to average variant
#expressions and average # barcodes per variant in the RNA samples, also create 
#a background column

MPRA_ave <- gen_epi %>%
  ungroup() %>%
  filter(subpool != 'control') %>%
  mutate(
    name = gsub('Smith R. Vista chr9:83712599-83712766', 'back_41', name),
    name = gsub('Vista Chr5:88673410-88674494', 'back_52', name),
    name = gsub('scramble pGL4.29 Promega 1-63 \\+ 1-87', 'back_55', name)
  ) %>%
  mutate(background = name) %>%
  mutate(background = str_sub(background, 
                              nchar(background)-1,
                              nchar(background))) %>%
  select(subpool, name, most_common, background, barcodes_RNA_br1, barcodes_RNA_br2, 
         ave_med_ratio, ave_med_ratio_norm, barcodes_RNA_22A, barcodes_RNA_22B, 
         ave_ratio_22, ave_ratio_22_norm) %>%
  mutate(genomic = (barcodes_RNA_br1 + barcodes_RNA_br2)/2) %>%
  mutate(episomal = (barcodes_RNA_22A + barcodes_RNA_22B)/2) %>%
  select(-barcodes_RNA_br1, -barcodes_RNA_br2, -barcodes_RNA_22A, 
         -barcodes_RNA_22B) %>%
  gather(genomic, episomal, key = 'MPRA', value = 'barcodes') %>%
  mutate(genomic = ave_med_ratio) %>%
  mutate(episomal = ave_ratio_22) %>%
  gather(genomic, episomal, key = 'MPRA2', value = 'ave_ratio') %>%
  filter((MPRA == 'genomic' & MPRA2 == 'genomic') | (MPRA == 'episomal' & MPRA2 == 'episomal')) %>%
  select(-MPRA2) %>%
  mutate() %>%
  mutate(genomic = ave_med_ratio_norm) %>%
  mutate(episomal = ave_ratio_22_norm) %>%
  gather(genomic, episomal, key = 'MPRA3', value = 'ave_ratio_norm') %>%
  filter((MPRA == 'genomic' & MPRA3 == 'genomic') | (MPRA == 'episomal' & MPRA3 == 'episomal')) %>%
  select(-MPRA3)

#Figure 1E and F----------------------------------------------------------------

#Figure 1E

#Plot replicability with backgrounds in orange and positive control in red

p_fig1_epi_med_rep <- gen_epi %>%
  ggplot(aes(med_ratio_22A, med_ratio_22B)) +
  geom_point(alpha = 0.1, size = 0.75) +
  geom_point(data = filter(gen_epi, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             fill = 'orange', shape = 21, size = 1.75) + 
  geom_point(data = filter(gen_epi, 
                           name == 'pGL4.29 Promega 1-63 + 1-87_back_55'), 
             fill = 'red', shape = 21, size = 1.75) +
  annotation_logticks(scaled = TRUE) +
  xlab("Expression (a.u.) replicate 1") +
  ylab("Expression (a.u.) replicate 2") +
  scale_x_log10(limits = c(0.03, 20), breaks = c(0.1, 1, 10)) + 
  scale_y_log10(limits = c(0.03, 20), breaks = c(0.1, 1, 10)) +
  theme(strip.background = element_rect(colour="black", fill="white"),
        axis.line.y = element_line(), panel.spacing.x=unit(1, "lines")) +
  figurefont_theme

p_fig1_gen_med_rep <- gen_epi %>%
  ggplot(aes(med_ratio_br1, med_ratio_br2)) +
  geom_point(alpha = 0.1, size = 1) +
  geom_point(data = filter(gen_epi, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             fill = 'orange', shape = 21, size = 1.75) + 
  geom_point(data = filter(gen_epi, name == 'pGL4.29 Promega 1-63 + 1-87_back_55'), 
             fill = 'red', shape = 21, size = 1.75) +
  annotation_logticks(scaled = TRUE) +
  xlab("Expression (a.u.) replicate 1") +
  ylab("Expression (a.u.) replicate 2") +
  scale_x_log10(limits = c(0.001, 15), breaks = c(0.01, 0.1, 1, 10)) + 
  scale_y_log10(limits = c(0.001, 15), breaks = c(0.01, 0.1, 1, 10)) +
  background_grid(major = 'xy', minor = 'none') + 
  theme(strip.background = element_rect(colour="black", fill="white"),
        axis.line.y = element_line(), panel.spacing.x=unit(1, "lines")) 

ggsave('../plots/p_fig1_epi_med_rep.png', p_fig1_epi_med_rep,
       width = 2, height = 2, units = 'in')

ggsave('../plots/p_fig1_gen_med_rep.png', p_fig1_gen_med_rep, 
       width = 2, height = 2, units = 'in')

log10_gen_epi <- var_log10(gen_epi)

pearsons_epi <- tibble(
  sample = c('all', 'subpool3', 'subpool5'),
  pearsons = c(cor(log10_gen_epi$med_ratio_22A, log10_gen_epi$med_ratio_22B, 
                   use = "pairwise.complete.obs", method = "pearson"),
               cor(filter(log10_gen_epi, subpool == 'subpool3')$med_ratio_22A,
                   filter(log10_gen_epi, subpool == 'subpool3')$med_ratio_22B,
                   use = "pairwise.complete.obs", method = "pearson"),
               cor(filter(log10_gen_epi, subpool == 'subpool5')$med_ratio_22A,
                   filter(log10_gen_epi, subpool == 'subpool5')$med_ratio_22B,
                   use = "pairwise.complete.obs", method = "pearson")))

pearsons_gen <- tibble(
  sample = c('all', 'subpool3', 'subpool5'),
  pearsons = c(cor(log10_gen_epi$med_ratio_br1, log10_gen_epi$med_ratio_br2, 
                   use = "pairwise.complete.obs", method = "pearson"),
               cor(filter(log10_gen_epi, subpool == 'subpool3')$med_ratio_br1,
                   filter(log10_gen_epi, subpool == 'subpool3')$med_ratio_br2,
                   use = "pairwise.complete.obs", method = "pearson"),
               cor(filter(log10_gen_epi, subpool == 'subpool5')$med_ratio_br1,
                   filter(log10_gen_epi, subpool == 'subpool5')$med_ratio_br2,
                   use = "pairwise.complete.obs", method = "pearson")))

#Figure 1F

p_gen_epi_rep <- gen_epi %>%
  ggplot(aes(ave_med_ratio, ave_ratio_22)) +
  geom_point(alpha = 0.1, size = 0.75) +
  geom_point(data = filter(gen_epi, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             fill = 'orange', shape = 21, size = 1.75) + 
  geom_point(data = filter(gen_epi, name == 'pGL4.29 Promega 1-63 + 1-87_back_55'), 
             fill = 'red', shape = 21, size = 1.75) +
  annotation_logticks(scaled = TRUE) +
  xlab("Average genomic expression (a.u.)") +
  ylab("Average episomal expression (a.u.)") +
  scale_x_log10(limits = c(0.01, 20), breaks = c(0.1, 1, 10)) + 
  scale_y_log10(limits = c(0.05, 20), breaks = c(0.1, 1, 10)) +
  theme(strip.background = element_rect(colour="black", fill="white"),
        axis.line.y = element_line(), panel.spacing.x=unit(1, "lines")) +
  figurefont_theme

ggsave('../plots/p_gen_epi_rep.png', p_gen_epi_rep,
       width = 2, height = 1.9, units = 'in')

log10_gen_epi <- var_log10(gen_epi)

pearsons_epi_gen <- tibble(
  pearsons = c(cor(log10_gen_epi$ave_med_ratio, log10_gen_epi$ave_ratio_22, 
                   use = "pairwise.complete.obs", method = "pearson")))

#Read in follow-up library containing barcoded replicates-----------------------

med_rep_followup <- read_tsv(
  '../20190712_epilib_analysis/med_rep_follow_up.txt')

scramble_newback_norm <- read_tsv(
  '../20190712_epilib_analysis/scramble_newback_norm.txt')

scramble_oldback_norm <- read_tsv(
  '../20190712_epilib_analysis/scramble_oldback_norm.txt')

sixsite_newback_norm <- read_tsv(
  '../20190712_epilib_analysis/sixsite_newback_norm.txt')

sixsite_oldback_norm <- read_tsv(
  '../20190712_epilib_analysis/sixsite_oldback_norm.txt')

twosite_norm <- read_tsv(
  '../20190712_epilib_analysis/twosite_norm.txt')

removed_bgs_sixsite_newback_norm <- read_tsv(
  '../20190712_epilib_analysis/removed_bgs_sixsite_newback_norm.txt')

#Figure 1G and H----------------------------------------------------------------

#compare differently barcoded replicates between old and new library

#Figure 1G - compare average expression between technical replicates between
#separately barcoded experiments

barcode_rep <- med_rep_0_22_A_B %>%
  select(subpool, name, most_common, barcodes_DNA, med_ratio_0A, 
         barcodes_RNA_0A, med_ratio_0B, barcodes_RNA_0B, med_ratio_22A, 
         barcodes_RNA_22A, med_ratio_22B, barcodes_RNA_22B) %>%
  inner_join(med_rep_followup, by = c('most_common'), 
             suffix = c('_old', '_new')) %>%
  mutate(ave_old_0 = (med_ratio_0A_old + med_ratio_0B_old)/2) %>%
  mutate(ave_old_4 = (med_ratio_22A + med_ratio_22B)/2) %>%
  mutate(ave_new_0 = (med_ratio_0A_new + med_ratio_0B_new)/2) %>%
  mutate(ave_new_4 = (med_ratio_4A + med_ratio_4B)/2) %>%
  mutate(ind_old = ave_old_4/ave_old_0) %>%
  mutate(ind_new = ave_new_4/ave_new_0)

p_old_new_4 <- barcode_rep %>%
  ggplot(aes(ave_old_4, ave_new_4)) +
  geom_point(alpha = 0.1, size = 1) +
  geom_point(data = filter(barcode_rep, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name_old)), 
             fill = 'orange', shape = 21, size = 1.75) + 
  geom_point(data = filter(barcode_rep, name_old == 'pGL4.29 Promega 1-63 + 1-87'), 
             fill = 'red', shape = 21, size = 1.75) +
  annotation_logticks(scaled = TRUE) +
  xlab("Barcode replicate 1\nexpression (a.u.)") +
  ylab("Barcode replicate 2\nexpression (a.u.)") +
  scale_x_log10(limits = c(0.06, 16), breaks = c(0.1, 1, 10)) + 
  scale_y_log10(limits = c(0.06, 16), breaks = c(0.1, 1, 10))  +
  theme(strip.background = element_rect(colour="black", fill="white"),
        axis.line.y = element_line(), panel.spacing.x=unit(1, "lines")) +
  figurefont_theme

ggsave('../plots/p_old_new_4.png', p_old_new_4,
       width = 2.2, height = 2, units = 'in')

log10_barcode_rep <- var_log10(barcode_rep)

pearsons_barcode_rep <- tibble(
  pearsons = c(cor(log10_barcode_rep$ave_old_4, log10_barcode_rep$ave_new_4, 
                   use = "pairwise.complete.obs", method = "pearson")))

#Figure 1H - plot average barcodes per variant between 2 technical replicates
#and compare between the library members between separately-barcoded variants in
#the episomal MPRA and to the genomic MPRA

bc_epi_RNA <- gen_epi %>%
  mutate(ave_barcodes = (barcodes_RNA_22A + barcodes_RNA_22B)/2) %>%
  mutate(sample = 'RNA') %>%
  mutate(group = 'episomal') %>%
  select(ave_barcodes, sample, group)

bc_epi_DNA <- gen_epi %>%
  mutate(ave_barcodes = barcodes_DNA) %>%
  mutate(sample = 'DNA') %>%
  mutate(group = 'episomal') %>%
  select(ave_barcodes, sample, group)

bc_gen_RNA <- gen_epi %>%
  mutate(ave_barcodes = (barcodes_RNA_br1 + barcodes_RNA_br2)/2) %>%
  mutate(sample = 'RNA') %>%
  mutate(group = 'genomic') %>%
  select(ave_barcodes, sample, group)

bc_gen_DNA <- gen_epi %>%
  mutate(ave_barcodes = (barcodes_DNA_br1 + barcodes_DNA_br2)/2) %>%
  mutate(sample = 'DNA') %>%
  mutate(group = 'genomic') %>%
  select(ave_barcodes, sample, group)

bc_barcode_rep_RNA <- barcode_rep %>%
  mutate(ave_barcodes = (barcodes_RNA_4A + barcodes_RNA_4B)/2) %>%
  mutate(sample = 'RNA') %>%
  mutate(group = 'episomal rep') %>%
  select(ave_barcodes, sample, group)

bc_barcode_rep_DNA <- barcode_rep %>%
  mutate(ave_barcodes = barcodes_DNA_new) %>%
  mutate(sample = 'DNA') %>%
  mutate(group = 'episomal rep') %>%
  select(ave_barcodes, sample, group)

bc_comb <- rbind(bc_epi_RNA, bc_epi_DNA, bc_gen_RNA, bc_gen_DNA, 
                 bc_barcode_rep_RNA, bc_barcode_rep_DNA)

p_bc_distr <- bc_comb %>%
  ggplot(aes(group, ave_barcodes, fill = sample)) +
  scale_fill_manual('Sample', values = c('#B8DDA3', '#A3B8DD')) +
  geom_violin(position = 'dodge', size = 0.25) +
  geom_boxplot(position = 'dodge', alpha = 0.5, outlier.alpha = 0, size = 0.25,
               width = 0.9) +
  scale_y_log10('Average barcodes per variant', breaks = c(10, 100)) +
  annotation_logticks(sides = 'l') +
  xlab('') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  figurefont_theme

ggsave('../plots/p_bc_distr.pdf', p_bc_distr, 
       width = 3, height = 2, units = 'in')


#Supplemental Figure 1C---------------------------------------------------------

#Import plate reader luminescence readings

orientation_lum <- read_csv("../plate_reader/20160607_controls_R.csv")

#Make df plottable with library orientation, control, and replicate as variables

orientation_lum_sep <- orientation_lum %>%
  melt(id = 'forskolin') %>%
  rename(luminescence = value) %>%
  separate(variable, into = c('orientation', 'control', 'fluff', 'replicate'),
           sep = "_", convert = TRUE) %>%
  select(-fluff) %>%
  mutate(orientation = ifelse(startsWith(orientation, 'rev'),
                              '- strand',
                              '+ strand')) %>%
  mutate(control = ifelse(startsWith(control, 'neg'),
                          'background',
                          'CRE control'))

p_orientation_lum_sep <- orientation_lum_sep %>%
  group_by(control, orientation) %>%
  ggplot(aes(forskolin, luminescence, color = control)) +
  scale_color_manual(values = c('darkgoldenrod1', 'firebrick2'), name = '') +
  geom_point(aes(fill = orientation), shape = 21, size = 1) +
  geom_smooth(aes(linetype = orientation)) +
  scale_fill_manual(values = c('black', 'grey')) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100), name = 'forskolin (µM)') +
  annotation_logticks(sides = 'b') +
  ylab('luminescence (a.u.)') +
  figurefont_theme +
  theme(legend.position = 'top')

ggsave('../plots/p_orientation_lum_sep.pdf', p_orientation_lum_sep,
       width = 3.5, height = 2.25, units = 'in')

#Supplemental Figure 1D---------------------------------------------------------

#Import plate reader luminescence readings

titration_luc <- read_csv("../plate_reader/170726_epi_gen_R.csv") %>%
  mutate(RLU_epi = luciferase_epi/renilla_epi) %>%
  mutate(forskolin = log2(forskolin))

p_gen_titration_luc <- titration_luc %>%
  ggplot(aes(forskolin, luciferase_gen)) +
  geom_point() +
  geom_smooth() +
  ylab('bulk luminescence (a.u.)') +
  scale_x_continuous(breaks = (-2:5), 'log2 forskolin (µM)') +
  annotation_logticks(sides = 'b') +
  figurefont_theme

p_epi_titration_luc <- titration_luc %>%
  ggplot(aes(forskolin, RLU_epi)) +
  geom_point() +
  geom_smooth() +
  ylab('bulk relative\nluminescence (a.u.)') +
  scale_x_continuous(breaks = (-2:5), 'log2 forskolin (µM)') +
  annotation_logticks(sides = 'b') +
  figurefont_theme

ggsave('../plots/p_gen_titration_luc.pdf', p_gen_titration_luc, width = 3.75, 
       height = 2.1, units = 'in')

ggsave('../plots/p_epi_titration_luc.pdf', p_epi_titration_luc, width = 3.75, 
       height = 2.1, units = 'in')

#Supplemental figure 1E---------------------------------------------------------

#Show background-normalized expression per variant between replicates and across
#forskolin titrations. Backgrounds normalize at 1 and are shown in orange and 
#the positive CRE control is shown in red. Concentrations are shown in log2 
#increments and the -7 value represents 0 µM forskolin. Pearson correlations are
#determined between replicates per each concentration.

p_var_rep_all <- ggplot(epi_back_norm_conc, 
                        aes(ratio_A_norm, ratio_B_norm)) +
  facet_rep_wrap(~ conc, nrow = 2, ncol = 4) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_point(data = filter(epi_back_norm_conc, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             fill = 'orange', shape = 21, size = 1.25) + 
  geom_point(data = filter(epi_back_norm_conc, name == 'pGL4.29 Promega 1-63 + 1-87_back_55'), 
             fill = 'red', shape = 21, size = 1.25) +
  annotation_logticks(scaled = TRUE) +
  xlab("Normalized expression (a.u.) replicate 1") +
  ylab("Normalized expression (a.u.) replicate 2") +
  scale_x_log10(limits = c(0.1, 200), breaks = c(0.1, 1, 10, 100)) + 
  scale_y_log10(limits = c(0.1, 200), breaks = c(0.1, 1, 10, 100)) +
  background_grid(major = 'xy', minor = 'none') + 
  panel_border(colour = 'black') +
  theme(strip.background = element_rect(colour="black", fill="white"),
        axis.line.y = element_line(), panel.spacing.x=unit(1, "lines")) +
  figurefont_theme

ggsave('../plots/var_rep_all.png', p_var_rep_all, width = 5.5, height = 3.1, 
       units = 'in')

log10_epi_back_norm_pc_spGl4 <- var_log10(epi_back_norm_pc_spGl4)

pearsons_epi_conc <- tibble(
  conc = c(0, 2^-5, 2^-4, 2^-3, 2^-2, 2^-1, 2^0, 2^2),
  pearsons = c(round(cor(log10_epi_back_norm_pc_spGl4$ratio_0A_norm, 
                         log10_epi_back_norm_pc_spGl4$ratio_0B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_epi_back_norm_pc_spGl4$ratio_2_5A_norm, 
                         log10_epi_back_norm_pc_spGl4$ratio_2_5B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_epi_back_norm_pc_spGl4$ratio_2_4A_norm, 
                         log10_epi_back_norm_pc_spGl4$ratio_2_4B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_epi_back_norm_pc_spGl4$ratio_2_3A_norm, 
                         log10_epi_back_norm_pc_spGl4$ratio_2_3B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_epi_back_norm_pc_spGl4$ratio_2_2A_norm, 
                         log10_epi_back_norm_pc_spGl4$ratio_2_2B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_epi_back_norm_pc_spGl4$ratio_2_1A_norm, 
                         log10_epi_back_norm_pc_spGl4$ratio_2_1B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_epi_back_norm_pc_spGl4$ratio_20A_norm, 
                         log10_epi_back_norm_pc_spGl4$ratio_20B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_epi_back_norm_pc_spGl4$ratio_22A_norm, 
                         log10_epi_back_norm_pc_spGl4$ratio_22B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3))
)

#Supplemental figure 1F---------------------------------------------------------

#Compare expression between variants with similar CRE architectures but
#different background sequences. Comparison in both the episomal and genomic
#MPRA

back_epi <- epi_back_norm_pc_spGl4 %>%
  select(name, background, med_ratio_22A, med_ratio_22B) %>%
  mutate(ave_ratio_22 = (med_ratio_22A + med_ratio_22B)/2) %>%
  select(-med_ratio_22A, -med_ratio_22B) %>%
  mutate(name = str_sub(name, 0, nchar(name)-8)) %>%
  spread(key = background, value = c(ave_ratio_22)) %>%
  var_log10() %>%
  na.omit() %>%
  select(-name)

back_gen <- gen_rep_1_2 %>%
  select(name, background, ave_med_ratio) %>%
  mutate(name = str_sub(name, 0, nchar(name)-8)) %>%
  spread(key = background, value = c(ave_med_ratio)) %>%
  var_log10() %>%
  na.omit() %>%
  select(-name)

my_points <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.1, size = 0.5) +
    scale_x_continuous(limits = c(-2, 1.3), breaks = c(-2:1)) + 
    scale_y_continuous(limits = c(-2, 1.3), breaks = c(-2:1)) +
    annotation_logticks(sides = 'bl')
}

my_density <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_density(kernel = 'gaussian') +
    scale_x_continuous(limits = c(-2, 1.3), breaks = c(-2:1)) +
    scale_y_continuous(limits = c(-0.5, 2)) +
    annotation_logticks(sides = 'b')
}

p_back_epi <- ggpairs(back_epi, 
                       columnLabels = c('background 41', 
                                        'background 52',
                                        'background 55'),
                       lower = list(continuous = my_points),
                       diag = list(continuous = my_density),
                       upper = list(continuous = wrap("cor", size = 2.5))) +
  panel_border(colour = 'black') +
  theme(strip.background = element_rect(colour="black", fill="white"),
        panel.spacing.x=unit(0.5, "lines")) +
  theme(panel.grid.major = element_blank()) +
  figurefont_theme

p_back_gen <- ggpairs(back_gen, 
                       columnLabels = c('background 41', 
                                        'background 52',
                                        'background 55'),
                       lower = list(continuous = my_points),
                       diag = list(continuous = my_density),
                       upper = list(continuous = wrap("cor", size = 2.5))) +
  panel_border(colour = 'black') +
  theme(strip.background = element_rect(colour="black", fill="white"),
        panel.spacing.x=unit(1, "lines")) +
  theme(panel.grid.major = element_blank()) +
  figurefont_theme

ggsave('../plots/p_back_epi.png', p_back_epi,
       width = 3.1, height = 2.96, units = 'in')

ggsave('../plots/p_back_gen.png', p_back_gen,
       width = 3.1, height = 2.96, units = 'in')

#Separate into sublibraries-----------------------------------------------------

#Subpool 3 corresponds to the CRE Spacing and Distance Library. This library
#contains 2 consensus CRE sites with flanks (ATTGACGTCAGC) that vary in 
#CRE Spacing from one another by 0 (no inner flanks), 5, 10, 15, 20, and 
#70 bp (all but 0 originally appear as spacing - 4 bp due to 4 bp of sequence 
#flanking both CREs). Both CREs are then placed at variable distances to the 3' 
#end of the backgrounds. Separation lists the CRE spacing between sites and CRE 
#distance. Distances measured from the end of the background to the CRE proximal
#to the promoter. Added 2 bp to all distances here to measure to start of CRE 
#without the flanks and then added 64 bp to measure to the minimal promoter. 
#Added 4 to all spacings but 0 to measure difference between start of CREs 
#without flanks.

subpool3 <- function(df) {
  df <- df %>%
    filter(subpool == "subpool3") %>%
    filter(name != 'pGL4.29 Promega 1-63 + 1-87_back_55') %>%
    ungroup() %>%
    select(-subpool) %>%
    mutate(name = gsub('2BS ', '', name), 
           name = gsub(' bp spacing ', '_', name)) %>%
    separate(name, 
             into = c("subpool", "spacing", "fluff1", "fluff2", "dist", 
                      "fluff3", "fluff4"),
             sep = "_", convert = TRUE) %>%
    select(-subpool, -fluff1, -fluff2, -fluff3, -fluff4) %>%
    mutate(dist = as.integer(dist + 2 + 64))%>%
    mutate(spacing = 
             ifelse(spacing != as.integer(0), 
                    as.integer(spacing + 4), as.integer(spacing)))
}

s3_gen_epi <- MPRA_ave %>%
  subpool3()

s3_epi_back_norm_conc <- med_rep_0_22_A_B %>%
  back_norm() %>%
  var_conc_exp() %>%
  subpool3()

#Subpool 5 corresponds to the CRE Number and Affinity Library. This library
#contains 6 equally spaced sites with 17 bp CRE Spacing. Per variant, each site
#is one of: the consensus CRE, the weak CRE, or no CRE. Both the weak and 
#consensus CREs have the same flanking sequence. Here site 1, 2, 3, 4, 5, and 6 
#equate to -191, -166, -141, -116, -91 and -66 site distances to the downstream
#promoter. Separation lists identity of CRE per site, the total CRE sites per
#variants (consensus + weak), and if variants contain only consenus CREs, only
#weak CREs, mixed consensus and weak, or no CREs.

subpool5 <- function(df) {
  df <- df %>%
    filter(subpool == "subpool5") %>%
    ungroup() %>%
    select(-subpool) %>%
    mutate(name = gsub('no_site', 'nosite', name)) %>%
    separate(name, into = c("subpool", "site1", "site2", "site3", "site4", 
                            "site5", "site6", "fluff", "fluff2"), sep = "_") %>%
    select(-subpool, -fluff, -fluff2) %>%
    mutate(consensus = str_detect(site1, "consensus") + 
             str_detect(site2, "consensus") + 
             str_detect(site3, "consensus") + 
             str_detect(site4, "consensus") + 
             str_detect(site5, "consensus") + 
             str_detect(site6, "consensus")) %>%
    mutate(weak = str_detect(site1, "weak") +
             str_detect(site2, "weak") +
             str_detect(site3, "weak") +
             str_detect(site4, "weak") +
             str_detect(site5, "weak") +
             str_detect(site6, "weak")) %>%
    mutate(nosite = str_detect(site1, "nosite") +
             str_detect(site2, "nosite") +
             str_detect(site3, "nosite") +
             str_detect(site4, "nosite") +
             str_detect(site5, "nosite") +
             str_detect(site6, "nosite")) %>%
    mutate(total_sites = consensus + weak) %>%
    mutate(site_combo = 
             ifelse(weak == 0 & consensus > 0, 
                    'consensus', 'mixed')) %>%
    mutate(site_combo = 
             ifelse(consensus == 0 & weak > 0, 
                    'weak', site_combo)) %>%
    mutate(site_combo = 
             ifelse(consensus == 0 & weak == 0, 
                    'none', site_combo)) 
}

s5_gen_epi <- MPRA_ave %>%
  subpool5()

s5_epi_back_norm_conc <- med_rep_0_22_A_B %>%
  back_norm() %>%
  var_conc_exp() %>%
  subpool5()

#Figure 2----------------------------------------------------------------------

#Figure 2A

#Episomal MPRA, effect of CRE distance on expression across forskolin 
#concentrations. Using Background 55, 10 bp CRE spacing as an example

s3_tidy_moveavg3 <- s3_epi_back_norm_conc %>%
  select(background, spacing, dist, conc, ave_ratio_norm) %>%
  group_by(background, spacing, conc) %>%
  arrange(dist, .by_group = TRUE) %>%
  moveavg_dist3()

p_induction_spacing_10_back_55 <- s3_tidy_moveavg3 %>%
  filter(background == '55' & spacing == 10) %>%
  ggplot(aes(dist, ave_ratio_norm)) +
  geom_point(aes(color = conc), size = 1, alpha = 0.5) +
  geom_line(aes(dist, ave_3, color = conc, group = conc), size = 0.5) +
  geom_line(data = filter(s3_tidy_moveavg3, conc == -7 & background == '55' & spacing == 10), 
            aes(dist, ave_3), color = '#FDE725FF', size = 0.75) +
  scale_x_reverse("Distance to minimal promoter (bp)",
                     breaks = c(seq(from = 60, to = 190, by = 10))) +
  scale_y_continuous("Average normalized\nexpression (a.u.)") +
  background_grid(major = 'x', minor = 'none', colour.major = 'grey70') + 
  scale_color_viridis(name = 'log2 forskolin (µM)', begin = 1, end = 0,
                      breaks = c(-7, -5, -4, -3, -2, -1, 0, 1, 2)) + 
  guides(color = guide_colorbar(frame.colour = 'black', ticks.colour = 'black')) +
  panel_border(colour = 'black') +
  figurefont_theme 

ggsave('../plots/p_induction_spacing_10_back_55.pdf', 
       p_induction_spacing_10_back_55, width = 5.5, height = 1.75, 
       units = 'in')

#Figure 2B

#Both MPRAs, effect of CRE distance on expression across spacings, using 
#background 55 as an example

s3_tidy_moveavg3_MPRA <- s3_gen_epi %>%
  select(background, spacing, dist, MPRA, ave_ratio_norm) %>%
  group_by(background, spacing, MPRA) %>%
  arrange(dist, .by_group = TRUE) %>%
  moveavg_dist3()

p_MPRA_dist_spacing_10_back_55 <- s3_tidy_moveavg3_MPRA %>%
  filter(background == 55, spacing == 10) %>%
  ggplot(aes(dist, ave_ratio_norm)) +
  annotate("rect", xmin = 67, xmax = 96, ymin = 0, ymax = 10.25, alpha = 0.6,
           color = 'black', fill = 'grey92', size = 0.4) +
  annotate("rect", xmin = 147, xmax = 176, ymin = 0, ymax = 10.25, alpha = 0.6,
           color = 'black', fill = 'grey92', size = 0.4) +
  geom_point(aes(color = MPRA), size = 1, alpha = 0.75) +
  geom_line(aes(dist, ave_3, color = MPRA, group = MPRA), 
            size = 0.5) +
  scale_color_manual(values = c('#3CBB75FF', 'gray35'), 
                     name = 'MPRA') +
  scale_x_reverse("Distance to minimal promoter (bp)",
                     breaks = c(seq(from = 60, to = 190, by = 10))) +
  scale_y_continuous("Average normalized\nexpression (a.u.)") +
  background_grid(major = 'x', minor = 'none', colour.major = 'grey70') + 
  panel_border(colour = 'black') +
  figurefont_theme +
  theme(strip.background = element_rect(colour="black", fill="white"))
  
ggsave('../plots/p_MPRA_dist_spacing_10_back_55.pdf', 
       p_MPRA_dist_spacing_10_back_55, 
       width = 5.4, height = 1.75, units = 'in')

#Take median expression of variants with a certain background, MPRA format,
#and CRE spacing across CRE distances 67-96 and 147-176 bp. Then divide median
#expression 147-176 bp over median expression 67-96 bp.

med_range_dist <- function(df) {
  low <- df %>%
    filter(dist >= 67 & dist <= 96) %>%
    mutate(range = '67-96')
  high <- df %>%
    filter(dist >= 147 & dist <= 176) %>%
    mutate(range = '147-176')
  range <- rbind(low, high) %>%
    group_by(range, spacing, background, MPRA) %>%
    summarize(median_range = median(ave_ratio_norm)) %>%
    ungroup()
  return(range)
}

fold_change_range <- function(df) {
  med_low <- df %>%
    filter(range == '67-96') %>%
    mutate(median_6796 = median_range) %>%
    select(-range, -median_range)
  med_high <- df %>%
    filter(range == '147-176') %>%
    mutate(median_147176 = median_range) %>%
    select(-range, -median_range)
  range_change <- inner_join(med_low, med_high, 
                            by = c('spacing', 'background', 'MPRA')) %>%
    mutate(range_change = median_6796/median_147176) %>%
    filter(spacing != 0, spacing != 70) %>%
    group_by(background, MPRA) %>%
    summarize(med_fold_change = median(range_change))
  return(range_change)
}

#Because so few variants were retained in this library with background 41 in the
#genomic MPRA, those in the episomal MPRA are added after the join for analysis 

epi_41_rbind <- function(df1, df2) {
  epi41 <- df1 %>%
    filter(background == 41 & conc == 2) %>%
    select(subpool, name, most_common, background, ave_ratio_norm) %>%
    subpool3() %>%
    mutate(MPRA = 'episomal')
  combinedMPRA <- df2 %>%
    select(spacing, dist, most_common, background, MPRA, ave_ratio_norm) %>%
    filter(!(MPRA == 'episomal' &  background == 41)) %>%
    rbind(epi41)
  return(combinedMPRA)
}

s3_gen_epi_rbind41 <- epi_41_rbind(epi_back_norm_conc, s3_gen_epi)
  
s3_gen_epi_med_range <- med_range_dist(s3_gen_epi)
s3_gen_epi_med_range_41 <- med_range_dist(s3_gen_epi_rbind41)

#determine median fold-change between distance ranges

s3_gen_epi_med_change_dist <- fold_change_range(s3_gen_epi_med_range) %>%
  filter(!(background == 41 & MPRA == 'genomic'))

s3_gen_epi_med_change_dist_41 <- fold_change_range(s3_gen_epi_med_range) %>%
  filter(background == 41 & MPRA == 'episomal')

#plot per background
  
p_s3_gen_epi_med_change_dist_55 <- s3_gen_epi_med_range %>%
  filter(spacing != 0 & spacing != 70 & background == '55') %>%
  mutate(range = factor(range, levels = c('147-176', '67-96'))) %>%
  ggplot(aes(range, median_range)) +
  facet_grid(. ~ MPRA) +
  geom_jitter(aes(color = as.factor(spacing)), 
              position=position_jitter(width=0.3, height=0), alpha = 0.75,
              size = 1, show.legend = TRUE) +
  geom_boxplot(outlier.shape=NA, size = 0.2, position = position_dodge(1),
               show.legend = FALSE, alpha = 0) +
  scale_color_manual(values = spacing_5_20_cbpalette, 
                     name = 'CRE spacing (bp)') +
  scale_y_continuous(limits = c(0, 13), 
                     breaks = c(seq(from = 0, to = 13, by = 3))) + 
  panel_border(colour = 'black') +
  ylab('Median normalized\nexpression across range') +
  xlab('CRE distance range (bp)') +
  figurefont_theme +
  theme(strip.background = element_rect(colour="black", fill="white"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(), legend.position = 'top')

p_s3_gen_epi_med_change_dist_52 <- s3_gen_epi_med_range %>%
  filter(spacing != 0 & spacing != 70 & background == '52') %>%
  mutate(range = factor(range, levels = c('147-176', '67-96'))) %>%
  ggplot(aes(range, median_range)) +
  facet_grid(. ~ MPRA) +
  geom_jitter(aes(color = as.factor(spacing)), 
              position=position_jitter(width=0.3, height=0), alpha = 0.75,
              size = 1, show.legend = TRUE) +
  geom_boxplot(outlier.shape=NA, size = 0.2, position = position_dodge(1),
               show.legend = FALSE, alpha = 0) +
  scale_color_manual(values = spacing_5_20_cbpalette, 
                     name = 'CRE spacing (bp)') +
  scale_y_continuous(limits = c(0, 13), 
                     breaks = c(seq(from = 0, to = 13, by = 3))) + 
  panel_border(colour = 'black') +
  ylab('Median normalize\nexpression across range') +
  xlab('CRE distance range (bp)') +
  figurefont_theme +
  theme(strip.background = element_rect(colour="black", fill="white"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(), legend.position = 'top')

p_s3_gen_epi_med_change_dist_41 <- s3_gen_epi_med_range_41  %>%
  filter(spacing != 0 & spacing != 70 & background =='41', MPRA == 'episomal') %>%
  mutate(range = factor(range, levels = c('147-176', '67-96'))) %>%
  ggplot(aes(range, median_range)) +
  geom_jitter(aes(color = as.factor(spacing)), 
              position=position_jitter(width=0.3, height=0), alpha = 0.75,
              size = 1, show.legend = TRUE) +
  geom_boxplot(outlier.shape=NA, size = 0.2, position = position_dodge(1),
               show.legend = FALSE, alpha = 0) +
  scale_color_manual(values = spacing_5_20_cbpalette, 
                     name = 'CRE spacing (bp)') +
  scale_y_continuous(limits = c(0, 13), 
                     breaks = c(seq(from = 0, to = 13, by = 3))) + 
  panel_border(colour = 'black') +
  ylab('Median normalize\nexpression across range') +
  xlab('CRE distance range (bp)') +
  figurefont_theme +
  theme(strip.background = element_rect(colour="black", fill="white"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(), legend.position = 'top')

ggsave('../plots/p_s3_gen_epi_med_change_dist_55.pdf', 
       p_s3_gen_epi_med_change_dist_55, width = 2.16, height = 2.3, unit = 'in')

ggsave('../plots/p_s3_gen_epi_med_change_dist_52.pdf', 
       p_s3_gen_epi_med_change_dist_52, width = 2.16, height = 2.3, unit = 'in')

ggsave('../plots/p_s3_gen_epi_med_change_dist_41.pdf', 
       p_s3_gen_epi_med_change_dist_41, width = 1.4, height = 2.3, unit = 'in')


#Wilcox test between distances for figure 2B

range_dist <- function(df) {
  low <- df %>%
    filter(dist >= 67 & dist <= 96) %>%
    mutate(range = '67-96')
  high <- df %>%
    filter(dist >= 147 & dist <= 176) %>%
    mutate(range = '147-176')
  range <- rbind(low, high)
}

s3_gen_epi_range <- s3_gen_epi_rbind41 %>%
  filter(spacing != 0 & spacing != 70) %>%
  select(spacing, dist, most_common, background, MPRA, ave_ratio_norm) %>%
  range_dist()

s3_gen_epi_range_67 <- s3_gen_epi_range %>%
  filter(range == '67-96')

s3_gen_epi_range_147 <- s3_gen_epi_range %>%
  filter(range == '147-176')


#Make dfs for wilcox test (ave_ratio_norm, range) per MPRA and background

#episomal MPRA

s3_gen_epi_range_67_epi_41 <- s3_gen_epi_range_67 %>%
  filter(MPRA == 'episomal', background == 41) %>%
  select(ave_ratio_norm, range)

s3_gen_epi_range_147_epi_41 <- s3_gen_epi_range_147 %>%
  filter(MPRA == 'episomal', background == 41) %>%
  select(ave_ratio_norm, range)

rbind_epi_41 <- rbind(s3_gen_epi_range_67_epi_41, s3_gen_epi_range_147_epi_41)


s3_gen_epi_range_67_epi_52 <- s3_gen_epi_range_67 %>%
  filter(MPRA == 'episomal', background == 52) %>%
  select(ave_ratio_norm, range)

s3_gen_epi_range_147_epi_52 <- s3_gen_epi_range_147 %>%
  filter(MPRA == 'episomal', background == 52) %>%
  select(ave_ratio_norm, range)

rbind_epi_52 <- rbind(s3_gen_epi_range_67_epi_52, s3_gen_epi_range_147_epi_52)


s3_gen_epi_range_67_epi_55 <- s3_gen_epi_range_67 %>%
  filter(MPRA == 'episomal', background == 55) %>%
  select(ave_ratio_norm, range)

s3_gen_epi_range_147_epi_55 <- s3_gen_epi_range_147 %>%
  filter(MPRA == 'episomal', background == 55) %>%
  select(ave_ratio_norm, range)

rbind_epi_55 <- rbind(s3_gen_epi_range_67_epi_55, s3_gen_epi_range_147_epi_55)

#genomic MPRA

s3_gen_epi_range <- s3_gen_epi %>%
  filter(spacing != 0 & spacing != 70) %>%
  select(spacing, dist, most_common, background, MPRA, ave_ratio_norm) %>%
  range_dist()

s3_gen_epi_range_67 <- s3_gen_epi_range %>%
  filter(range == '67-96')

s3_gen_epi_range_147 <- s3_gen_epi_range %>%
  filter(range == '147-176')

s3_gen_epi_range_67_gen_41 <- s3_gen_epi_range_67 %>%
  filter(MPRA == 'genomic', background == 41) %>%
  select(ave_ratio_norm, range)

s3_gen_epi_range_147_gen_41 <- s3_gen_epi_range_147 %>%
  filter(MPRA == 'genomic', background == 41) %>%
  select(ave_ratio_norm, range)

rbind_gen_41 <- rbind(s3_gen_epi_range_67_gen_41, s3_gen_epi_range_147_gen_41)


s3_gen_epi_range_67_gen_52 <- s3_gen_epi_range_67 %>%
  filter(MPRA == 'genomic', background == 52) %>%
  select(ave_ratio_norm, range)

s3_gen_epi_range_147_gen_52 <- s3_gen_epi_range_147 %>%
  filter(MPRA == 'genomic', background == 52) %>%
  select(ave_ratio_norm, range)

rbind_gen_52 <- rbind(s3_gen_epi_range_67_gen_52, s3_gen_epi_range_147_gen_52)


s3_gen_epi_range_67_gen_55 <- s3_gen_epi_range_67 %>%
  filter(MPRA == 'genomic', background == 55) %>%
  select(ave_ratio_norm, range)

s3_gen_epi_range_147_gen_55 <- s3_gen_epi_range_147 %>%
  filter(MPRA == 'genomic', background == 55) %>%
  select(ave_ratio_norm, range)

rbind_gen_55 <- rbind(s3_gen_epi_range_67_gen_55, s3_gen_epi_range_147_gen_55)

#wilcox tests

wilcox.test(ave_ratio_norm ~ range, data = rbind_epi_41)
wilcox.test(ave_ratio_norm ~ range, data = rbind_epi_52)
wilcox.test(ave_ratio_norm ~ range, data = rbind_epi_55)
wilcox.test(ave_ratio_norm ~ range, data = rbind_gen_52)
wilcox.test(ave_ratio_norm ~ range, data = rbind_gen_55)


#Supplemental Figure 2A---------------------------------------------------------

#Import background-normalized expression table following first MPRA performed.
#MPRA performed at 0 and 25 µM forskolin with library designs corresponding to 
#the CRE Distance, CRE Spacing and Distance, CRE Number and Affinity library, 
#and another library based on the CRE Spacing and Distance library but with 8 
#combinations of site affinities across the two CREs and less Distances tested.

epi_rep_0_25 <- read_tsv('../20170320_epilib_analysis_0_25/rep_0_25_med_norm.txt')

#look at CRE Distance Library (Subpool 2)

#Subpool 2 contains either a consensus CRE (TGACGTCA) or consensus surrounded 
#by 2 bp flanks on either side (ATTGACGTCAGC) as is used in the rest of the 
#subpools. Each site is placed on the background starting closest to the minP 
#and are then moved along the backgrounds at 1 bp increments. Separation lists 
#the type of CRE, CRE distance, and the background used. Added 64 bp to measure 
#distance of site to minimal promoter then added 2 bp to consensusflank so that 
#the start of the binding site is represented instead of the start of the flank.

subpool2 <- function(df) {
  df <- df %>%
    filter(subpool == "subpool2") %>%
    ungroup() %>%
    select(-subpool) %>%
    separate(name, into = c("fluff1", "site", "fluff2", "dist",
                            "fluff3", "fluff4"),
             sep = "_", convert = TRUE) %>% 
    select(-fluff1, -fluff2, -fluff3, -fluff4) %>%
    mutate(dist = dist + 64) %>%
    mutate(dist = ifelse(startsWith(site, 'consensusflank'), dist + 2, dist))
}

s2_epi_rep_0_25 <- subpool2(epi_rep_0_25)

s2_untidy_moveavg3 <- s2_epi_rep_0_25 %>%
  select(background, site, dist, ave_ratio_norm, conc) %>%
  group_by(background, site, conc) %>%
  arrange(dist, .by_group = TRUE) %>%
  moveavg_dist3()

p_subpool2_dist_0_25 <- s2_untidy_moveavg3 %>%
  filter(site == 'consensusflank') %>%
  ggplot(aes(dist, ave_ratio_norm, color = as.factor(conc))) + 
  facet_grid(background ~ .) + 
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c('gray20', 'firebrick3'), 
                     name = 'forskolin (µM)') +
  scale_x_reverse("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 66, to = 206, by = 20)) +
  theme(legend.position = 'right',
        strip.background = element_rect(colour="black", fill="white")) +
  panel_border(colour = 'black') + 
  ylab('Average normalized expression (a.u.)') +
  background_grid(major = 'x', minor = 'x', colour.major = 'grey90',
                  colour.minor = 'grey95') +
  figurefont_theme

ggsave('../plots/subpool2_dist_0_25.pdf', p_subpool2_dist_0_25, units = 'in',
       width = 5, height = 2.66)

#Supplemental figure 2B---------------------------------------------------------

#Plot average background-normalized expression per MPRA, CRE Spacing, and CRE
#background acros CRE Distances. 

#fit 3 bp moving window function

s3_tidy_moveavg3_MPRA <- s3_gen_epi_rbind41 %>%
  select(background, spacing, dist, MPRA, ave_ratio_norm) %>%
  group_by(background, spacing, MPRA) %>%
  arrange(dist, .by_group = TRUE) %>%
  moveavg_dist3()

p_space_dist_gen_epi <- s3_tidy_moveavg3_MPRA %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = MPRA)) + 
  geom_point(alpha = 0.5, size = 0.25) +
  geom_line(data = filter(s3_tidy_moveavg3_MPRA, 
                          !(MPRA == "genomic" & background == 41)),
            aes(y = ave_3), size = 0.4) +
  facet_grid(spacing ~ background) +
  scale_color_manual(values = c('#29AF7FFF', 'gray35')) +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  background_grid(major = 'x', minor = 'x', colour.major = 'grey90',
                  colour.minor = 'grey95') +
  scale_y_continuous(limits = c(0.3, 47)) +
  scale_x_reverse("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 70, to = 190, by = 20)) +
  figurefont_theme +
  theme(legend.position = 'top',
        strip.background = element_rect(colour="black", fill="white"))

ggsave('../plots/space_dist_gen_epi.pdf', p_space_dist_gen_epi, units = 'in',
       width = 6.4, height = 5)

#Figure 3 and Supplemental Figure 3---------------------------------------------

#Figure 3A

p_subpool3_spa_4_vchr9_5_10 <- s3_tidy_moveavg3_MPRA %>%
  filter(MPRA == 'episomal' & background == 41 & dist < 124 & (spacing == 5 | spacing == 10)) %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.65) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c('gray20', 'dodgerblue2'), name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  geom_vline(xintercept = c(78, 88), color = 'gray20', linetype = 2, 
             alpha = 0.5) +
  geom_vline(xintercept = c(83, 92), color = 'dodgerblue2', linetype = 2, 
             alpha = 0.5) +
  scale_y_continuous(limits = c(1, 15)) +
  background_grid(major = 'x', minor = 'none') +
  scale_x_reverse("Distance to minimal promoter from proximal CRE (bp)", 
                     breaks = seq(from = 66, to = 126, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

p_subpool3_spa_4_vchr9_5_15 <- s3_tidy_moveavg3_MPRA %>%
  filter(MPRA == 'episomal' & background == 41 & dist < 124 & (spacing == 5 | spacing == 15)) %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.65) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c('gray20', 'orangered3'), 
                     name = 'spacing (bp)') +
  scale_fill_manual(values = c('gray20', 'orangered3'), 
                    name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  geom_vline(xintercept = c(78, 88), color = 'gray20', linetype = 2, 
             alpha = 0.5) +
  scale_y_continuous(limits = c(1, 15)) +
  background_grid(major = 'x', minor = 'none') +
  scale_x_reverse("Distance to minimal promoter from proximal CRE (bp)", 
                  breaks = seq(from = 66, to = 126, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

#10 and 20 bp spacing do not have a data point at 66 bp CRE Distance, which 
#makes plotting similar-sized graphs difficult, add "data" at 66 with alpha = 1

dumbdata <- tibble(dist = 66, ave_ratio_norm = 0.5)

p_subpool3_spa_4_vchr9_10_20 <- s3_tidy_moveavg3_MPRA %>%
  filter(MPRA == 'episomal' & background == 41 & dist < 124 & (spacing == 10 | spacing == 20)) %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.65) +
  geom_point(data = dumbdata, alpha = 1, color = 'white') +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c('dodgerblue2', 'sandybrown'),
                     name = 'spacing (bp)') +
  scale_fill_manual(values = c('dodgerblue2', 'sandybrown'),
                    name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  geom_vline(xintercept = c(83, 92), color = 'dodgerblue2', linetype = 2, 
             alpha = 0.5) +
  scale_y_continuous(limits = c(1, 15)) +
  background_grid(major = 'x', minor = 'none') +
  scale_x_reverse("Distance to minimal promoter from proximal CRE (bp)", 
                  breaks = seq(from = 66, to = 126, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

ggsave('../plots/p_subpool3_med_spa_4_vchr9_5_10.pdf', 
       p_subpool3_spa_4_vchr9_5_10, height = 1.35, width = 4.8, units = 'in')

ggsave('../plots/p_subpool3_med_spa_4_vchr9_5_15.pdf', 
       p_subpool3_spa_4_vchr9_5_15, height = 1.35, width = 4.8, units = 'in')

ggsave('../plots/p_subpool3_med_spa_4_vchr9_10_20.pdf', 
       p_subpool3_spa_4_vchr9_10_20, height = 1.35, width = 4.8, units = 'in')

#Figure 3B

#Plot all on one plot based on distance from distal CRE

twosite_4_moveavg3_dist_distal_MPRA <- s3_gen_epi_rbind41 %>%
  mutate(distal_dist = dist + spacing + 8) %>%
  select(background, spacing, distal_dist, MPRA, ave_ratio_norm) %>%
  group_by(background, spacing, MPRA) %>%
  arrange(distal_dist, .by_group = TRUE) %>%
  moveavg_dist3()

p_subpool3_spa_4_vchr9_dist_distal <- twosite_4_moveavg3_dist_distal_MPRA %>%
  filter(MPRA == 'episomal' & background == 41 & distal_dist < 137 & spacing != 0 & spacing != 70) %>%
  ggplot(aes(x = distal_dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.65) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c('gray20', 'dodgerblue2', 
                                'orangered3', 'sandybrown'),
                     name = 'spacing (bp)') +
  scale_fill_manual(values = c('gray20', 'dodgerblue2', 
                               'orangered3', 'sandybrown'),
                    name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  geom_vline(xintercept = c(91, 101, 111), color = 'black', linetype = 2, 
             alpha = 0.5) +
  scale_y_continuous(limits = c(1, 15)) +
  background_grid(major = 'x', minor = 'none') +
  scale_x_reverse("Distance to minimal promoter from distal CRE (bp)", 
                  breaks = seq(from = 66, to = 140, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

ggsave('../plots/p_subpool3_spa_4_vchr9_dist_distal.pdf', 
       p_subpool3_spa_4_vchr9_dist_distal, height = 1.35, width = 4.8, 
       units = 'in')

#Supplemental Figure 3

#Overlay of CRE expression profiles following CRE Distance and between the CRE 
#Spacings 5 and 10, 5 and 15, and 10 and 20 bp. Plots similar to figure 3A are
#shown with back 52 and back 55 and across both MPRAs. Add in fake point for 10
#bp spacing to align graph axes.

dumbdata <- tibble(dist = 191, ave_ratio_norm = 3, spacing = 10)

p_subpool3_spa_spgl4_trans_5_10 <- s3_tidy_moveavg3_MPRA %>%
  filter(background == 55 & (spacing == 5 | spacing == 10) & MPRA == 'episomal') %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_vline(xintercept = c(94, 104.5), color = 'gray20', linetype = 2, 
             alpha = 0.5) +
  geom_vline(xintercept = c(99, 109), color = 'dodgerblue2', linetype = 2, 
             alpha = 0.5) +
  scale_color_manual(values = c('gray20', 'dodgerblue2'), name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_continuous(limits = c(1, 7)) +
  background_grid(major = 'x', minor = 'none') +
  scale_x_reverse("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

p_subpool3_spa_spgl4_trans_5_15 <- s3_tidy_moveavg3_MPRA %>%
  filter(background == 55 & (spacing == 5 | spacing == 15) & MPRA == 'episomal') %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_vline(xintercept = c(94, 104.5), color = 'gray20', linetype = 2, 
             alpha = 0.5) +
  scale_color_manual(values = c('gray20', 'orangered3'), 
                     name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_continuous(limits = c(1, 7)) +
  background_grid(major = 'x', minor = 'none') +
  scale_x_reverse("Distance to minimal promoter (bp)", 
                  breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

p_subpool3_spa_spgl4_trans_10_20 <- s3_tidy_moveavg3_MPRA %>%
  filter(background == 55 & (spacing == 10 | spacing == 20) & MPRA == 'episomal') %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_point(data = dumbdata, alpha = 0) +
  geom_vline(xintercept = c(99, 109), color = 'dodgerblue2', linetype = 2, 
             alpha = 0.5) +
  scale_color_manual(values = c('dodgerblue2', 'sandybrown'),
                     name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_continuous(limits = c(1, 7)) +
  background_grid(major = 'x', minor = 'none') +
  scale_x_reverse("Distance to minimal promoter (bp)", 
                  breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

ggsave('../plots/p_subpool3_spa_spgl4_trans_5_10.pdf', 
       p_subpool3_spa_spgl4_trans_5_10, 
       height = 1.25, width = 4.4, units = 'in')

ggsave('../plots/p_subpool3_spa_spgl4_trans_5_15.pdf', 
       p_subpool3_spa_spgl4_trans_5_15, 
       height = 1.25, width = 4.4, units = 'in')

ggsave('../plots/p_subpool3_spa_spgl4_trans_10_20.pdf', 
       p_subpool3_spa_spgl4_trans_10_20, 
       height = 1.25, width = 4.4, units = 'in')

p_subpool3_spa_spgl4_int_5_10 <- s3_tidy_moveavg3_MPRA %>%
  filter(background == 55 & (spacing == 5 | spacing == 10) & MPRA == 'genomic') %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c('gray20', 'dodgerblue2'), name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_continuous(limits = c(0.4, 11)) +
  background_grid(major = 'x', minor = 'none') +
  scale_x_reverse("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

p_subpool3_spa_spgl4_int_5_15 <- s3_tidy_moveavg3_MPRA %>%
  filter(background == 55 & (spacing == 5 | spacing == 15) & MPRA == 'genomic') %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c('gray20', 'orangered3'), 
                     name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_continuous(limits = c(0.4, 11)) +
  background_grid(major = 'x', minor = 'none') +
  scale_x_reverse("Distance to minimal promoter (bp)", 
                  breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

p_subpool3_spa_spgl4_int_10_20 <- s3_tidy_moveavg3_MPRA %>%
  filter(background == 55 & (spacing == 10 | spacing == 20) & MPRA == 'genomic') %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_point(data = dumbdata, alpha = 0) +
  scale_color_manual(values = c('dodgerblue2', 'sandybrown'),
                     name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_continuous(limits = c(0.4, 11)) +
  background_grid(major = 'x', minor = 'none') +
  scale_x_reverse("Distance to minimal promoter (bp)", 
                  breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

ggsave('../plots/p_subpool3_spa_spgl4_int_5_10.pdf', 
       p_subpool3_spa_spgl4_int_5_10, 
       height = 1.25, width = 4.4, units = 'in')

ggsave('../plots/p_subpool3_spa_spgl4_int_5_15.pdf', 
       p_subpool3_spa_spgl4_int_5_15, 
       height = 1.25, width = 4.4, units = 'in')

ggsave('../plots/p_subpool3_spa_spgl4_int_10_20.pdf', 
       p_subpool3_spa_spgl4_int_10_20, 
       height = 1.25, width = 4.4, units = 'in')


p_subpool3_spa_vchr5_trans_5_10 <- s3_tidy_moveavg3_MPRA %>%
  filter(background == 52 & (spacing == 5 | spacing == 10) & MPRA == 'episomal') %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c('gray20', 'dodgerblue2'), name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_continuous(limits = c(1, 32)) +
  background_grid(major = 'x', minor = 'none') +
  scale_x_reverse("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

p_subpool3_spa_vchr5_trans_5_15 <- s3_tidy_moveavg3_MPRA %>%
  filter(background == 52 & (spacing == 5 | spacing == 15) & MPRA == 'episomal') %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c('gray20', 'orangered3'), 
                     name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_continuous(limits = c(1, 32)) +
  background_grid(major = 'x', minor = 'none') +
  scale_x_reverse("Distance to minimal promoter (bp)", 
                  breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

p_subpool3_spa_vchr5_trans_10_20 <- s3_tidy_moveavg3_MPRA %>%
  filter(background == 52 & (spacing == 10 | spacing == 20) & MPRA == 'episomal') %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_point(data = dumbdata, alpha = 0) +
  scale_color_manual(values = c('dodgerblue2', 'sandybrown'),
                     name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_continuous(limits = c(1, 32)) +
  background_grid(major = 'x', minor = 'none') +
  scale_x_reverse("Distance to minimal promoter (bp)", 
                  breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

ggsave('../plots/p_subpool3_spa_vchr5_trans_5_10.pdf', 
       p_subpool3_spa_vchr5_trans_5_10, 
       height = 1.25, width = 4.4, units = 'in')

ggsave('../plots/p_subpool3_spa_vchr5_trans_5_15.pdf', 
       p_subpool3_spa_vchr5_trans_5_15, 
       height = 1.25, width = 4.4, units = 'in')

ggsave('../plots/p_subpool3_spa_vchr5_trans_10_20.pdf', 
       p_subpool3_spa_vchr5_trans_10_20, 
       height = 1.25, width = 4.4, units = 'in')


p_subpool3_spa_vchr5_int_5_10 <- s3_tidy_moveavg3_MPRA %>%
  filter(background == 52 & (spacing == 5 | spacing == 10) & MPRA == 'genomic') %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c('gray20', 'dodgerblue2'), name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_continuous(limits = c(1, 35)) +
  background_grid(major = 'x', minor = 'none') +
  scale_x_reverse("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

p_subpool3_spa_vchr5_int_5_15 <- s3_tidy_moveavg3_MPRA %>%
  filter(background == 52 & (spacing == 5 | spacing == 15) & MPRA == 'genomic') %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c('gray20', 'orangered3'), 
                     name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_continuous(limits = c(1, 35)) +
  background_grid(major = 'x', minor = 'none') +
  scale_x_reverse("Distance to minimal promoter (bp)", 
                  breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

p_subpool3_spa_vchr5_int_10_20 <- s3_tidy_moveavg3_MPRA %>%
  filter(background == 52 & (spacing == 10 | spacing == 20) & MPRA == 'genomic') %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_point(data = dumbdata, alpha = 0) +
  scale_color_manual(values = c('dodgerblue2', 'sandybrown'),
                     name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_continuous(limits = c(1, 35)) +
  background_grid(major = 'x', minor = 'none') +
  scale_x_reverse("Distance to minimal promoter (bp)", 
                  breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

ggsave('../plots/p_subpool3_spa_vchr5_int_5_10.pdf', 
       p_subpool3_spa_vchr5_int_5_10, 
       height = 1.25, width = 4.4, units = 'in')

ggsave('../plots/p_subpool3_spa_vchr5_int_5_15.pdf', 
       p_subpool3_spa_vchr5_int_5_15, 
       height = 1.25, width = 4.4, units = 'in')

ggsave('../plots/p_subpool3_spa_vchr5_int_10_20.pdf', 
       p_subpool3_spa_vchr5_int_10_20, 
       height = 1.25, width = 4.4, units = 'in')

#Figure 4-----------------------------------------------------------------------

#Figure 4B

#tileplot of all backgrounds and spacings, take average expression per distance
#per background across spacings and plot 3 bp moving average to summarize
#periodicities above plots

twosite_4_avespace_moveavg3_dist_distal <- twosite_norm %>%
  mutate(distal_dist = dist + spacing + 8) %>%
  mutate(ave_ratio_norm = ave_med_ratio_norm_4) %>%
  select(background, spacing, distal_dist, ave_ratio_norm) %>%
  group_by(background, distal_dist) %>%
  summarize(ave_ratio_norm = mean(ave_ratio_norm)) %>%
  group_by(background) %>%
  arrange(distal_dist, .by_group = TRUE) %>%
  moveavg_dist3()

p_twosite_ave3 <- twosite_4_avespace_moveavg3_dist_distal %>%
  ggplot(aes(distal_dist, ave_ratio_norm)) +
  facet_grid(background ~ ., scales = 'free') +
  geom_point(color = 'gray50', alpha = 0.75, size = 0.5) +
  geom_line(aes(distal_dist, ave_3), size = 0.5) +
  scale_x_reverse("Distance to minimal promoter from distal CRE (bp)",
                  breaks = c(seq(from = 60, to = 200, by = 10))) +
  ylab('Average expression across spacings (a.u.)') +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

ggsave('../plots/p_twosite_ave3.pdf', p_twosite_ave3,
       width = 5.6, height = 1.25, units = 'in')

#plot according to distal CRE distance

twosite_4_moveavg3_dist_distal <- twosite_norm %>%
  mutate(distal_dist = dist + spacing + 8) %>%
  mutate(ave_ratio_norm = ave_med_ratio_norm_4) %>%
  select(background, spacing, distal_dist, ave_ratio_norm) %>%
  group_by(background, spacing) %>%
  arrange(distal_dist, .by_group = TRUE) %>%
  moveavg_dist3()

p_twosite_tile_back41 <- twosite_4_moveavg3_dist_distal %>%
  filter(background == 41) %>%
  ggplot(aes(x = distal_dist, y = spacing, fill = ave_ratio_norm)) +
  facet_grid(background ~ .) +
  geom_tile() +
  scale_fill_gradientn(name = 'Average normalized\nexpression (a.u.)',
                       colors = cbPalette_cont_blue) +
  guides(fill = guide_colorbar(frame.colour = 'black', 
                               ticks.colour = 'black')) +
  scale_x_reverse("Distance to minimal promoter from distal CRE (bp)",
                  breaks = c(seq(from = 60, to = 200, by = 10))) +
  ylab("Spacing between CREs (bp)") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

p_twosite_tile_back52 <- twosite_4_moveavg3_dist_distal %>%
  filter(background == 52) %>%
  ggplot(aes(x = distal_dist, y = spacing, fill = ave_ratio_norm)) +
  facet_grid(background ~ .) +
  geom_tile() +
  scale_fill_gradientn(name = 'Average normalized\nexpression (a.u.)',
                       colors = cbPalette_cont_blue) +
  guides(fill = guide_colorbar(frame.colour = 'black', 
                               ticks.colour = 'black')) +
  scale_x_reverse("Distance to minimal promoter from distal CRE (bp)",
                  breaks = c(seq(from = 60, to = 200, by = 10))) +
  ylab("Spacing between CREs (bp)") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

p_twosite_tile_back55 <- twosite_4_moveavg3_dist_distal %>%
  filter(background == 55) %>%
  ggplot(aes(x = distal_dist, y = spacing, fill = ave_ratio_norm)) +
  facet_grid(background ~ .) +
  geom_tile() +
  scale_fill_gradientn(name = 'Average normalized\nexpression (a.u.)',
                       colors = cbPalette_cont_blue) +
  guides(fill = guide_colorbar(frame.colour = 'black', 
                               ticks.colour = 'black')) +
  scale_x_reverse("Distance to minimal promoter from distal CRE (bp)",
                  breaks = c(seq(from = 60, to = 200, by = 10))) +
  ylab("Spacing between CREs (bp)") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

ggsave('../plots/p_twosite_tile_back41.pdf', p_twosite_tile_back41,
       width = 7, height = 1.75, units = 'in')

ggsave('../plots/p_twosite_tile_back52.pdf', p_twosite_tile_back52,
       width = 7, height = 1.75, units = 'in')

ggsave('../plots/p_twosite_tile_back55.pdf', p_twosite_tile_back55,
       width = 7, height = 1.75, units = 'in')

#Figure 4C

#Look at mean expression across spacings and backgrounds

ave_exp_spacings <- twosite_norm %>%
  group_by(background, spacing) %>%
  summarize(ave_ratio_norm = mean(ave_med_ratio_norm_4))

p_twosite_spacing_medexp <- ave_exp_spacings %>%
  ggplot(aes(spacing, ave_ratio_norm, color = as.factor(background))) +
  geom_point(size = 1) +
  scale_x_continuous(breaks = c(1:13), name = 'CRE spacing (bp)') +
  ylab('Average expression\nacross distances (a.u.)') +
  scale_color_manual(values = background_cbPalette,
                     name = 'Background') +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

ggsave('../plots/p_twosite_spacing_medexp.pdf', 
       p_twosite_spacing_medexp, height = 1.75, width = 4, units = 'in')


#Figure 5, Supplemental Figure 4B and 5-----------------------------------------

#Figure 5A

p_s5_num_cons_num_weak_back55 <- s5_gen_epi %>%
  filter(background == '55') %>%
  ggplot(aes(as.factor(consensus), ave_ratio_norm)) +
  facet_grid(MPRA ~ .) +
  geom_boxplot(aes(fill = as.factor(weak)), outlier.size = 0.75, size = 0.3, 
               outlier.shape = 21, outlier.alpha = 1, 
               position = position_dodge(0.75)) +
  scale_y_log10() + 
  panel_border(colour = 'black') +
  annotation_logticks(sides = 'l') +
  scale_fill_manual(name = 'number of\nweak CREs', 
                    values = cbPalette7_grad_light) +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = 'top') +
  ylab('Average normalized expression (a.u.)') +
  xlab('Number of consensus CREs') +
  figurefont_theme

ggsave('../plots/p_s5_num_cons_num_weak_back55.pdf', 
       p_s5_num_cons_num_weak_back55, width = 4.5, height = 3.75, units = 'in')

#Determine median expression of variants with one consensus CRE in background 55

back55_cons_1_weak_1_epi <- s5_gen_epi %>%
  filter(background == '55' & consensus == 1 & weak == 1 & MPRA == 'episomal') %>%
  summarize(median(ave_ratio_norm))

back55_cons_1_weak_5_epi <- s5_gen_epi %>%
  filter(background == '55' & consensus == 1 & weak == 5 & MPRA == 'episomal') %>%
  summarize(median(ave_ratio_norm))

back55_cons_1_weak_5_epi/back55_cons_1_weak_1_epi

back55_cons_1_weak_1_gen <- s5_gen_epi %>%
  filter(background == '55' & consensus == 1 & weak == 1 & MPRA == 'genomic') %>%
  summarize(median(ave_ratio_norm))

back55_cons_1_weak_5_gen <- s5_gen_epi %>%
  filter(background == '55' & consensus == 1 & weak == 5 & MPRA == 'genomic') %>%
  summarize(median(ave_ratio_norm))

back55_cons_1_weak_5_gen/back55_cons_1_weak_1_gen

#Supplemental Figure 4B

p_s5_num_cons_num_weak_back41 <- s5_gen_epi %>%
  filter(background == 41) %>%
  ggplot(aes(as.factor(consensus), ave_ratio_norm)) +
  facet_grid(MPRA ~ .) +
  geom_boxplot(aes(fill = as.factor(weak)), outlier.size = 0.75, size = 0.3, 
               outlier.shape = 21, outlier.alpha = 1, 
               position = position_dodge(0.75)) +
  scale_y_log10() + 
  panel_border(colour = 'black') +
  annotation_logticks(sides = 'l') +
  scale_fill_manual(name = 'number of\nweak CREs', 
                    values = cbPalette7_grad_light) +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = 'top') +
  ylab('Average normalized expression (a.u.)') +
  xlab('Number of consensus CREs') +
  figurefont_theme

ggsave('../plots/p_s5_num_cons_num_weak_back41.pdf', 
       p_s5_num_cons_num_weak_back41, width = 4.5, height = 3.75, units = 'in')

p_s5_num_cons_num_weak_back52 <- s5_gen_epi %>%
  filter(background == 52) %>%
  ggplot(aes(as.factor(consensus), ave_ratio_norm)) +
  facet_grid(MPRA ~ .) +
  geom_boxplot(aes(fill = as.factor(weak)), outlier.size = 0.75, size = 0.3, 
               outlier.shape = 21, outlier.alpha = 1, 
               position = position_dodge(0.75)) +
  scale_y_log10() + 
  panel_border(colour = 'black') +
  annotation_logticks(sides = 'l') +
  scale_fill_manual(name = 'number of\nweak CREs', 
                    values = cbPalette7_grad_light) +
  theme(axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = 'top') +
  ylab('Average normalized expression (a.u.)') +
  xlab('Number of consensus CREs') +
  figurefont_theme

ggsave('../plots/p_s5_num_cons_num_weak_back52.pdf', 
       p_s5_num_cons_num_weak_back52, width = 4.5, height = 3.75, units = 'in')

#Prepare df for lm fits

#change order of site types to fit consensus and weak CRE expression relative to
#no site

set.seed(1234)

subpool5_ncw <- s5_gen_epi %>%
  mutate(site1 = gsub('nosite', 'anosite', site1)) %>%
  mutate(site2 = gsub('nosite', 'anosite', site2)) %>%
  mutate(site3 = gsub('nosite', 'anosite', site3)) %>%
  mutate(site4 = gsub('nosite', 'anosite', site4)) %>%
  mutate(site5 = gsub('nosite', 'anosite', site5)) %>%
  mutate(site6 = gsub('nosite', 'anosite', site6)) %>%
  var_log10()

subpool5_epi <- subpool5_ncw %>%
  filter(MPRA == 'episomal')

subpool5_gen <- subpool5_ncw %>%
  filter(MPRA == 'genomic')

#split into training and test data

s5_smp_size <- floor(0.8 * nrow(subpool5_epi))

s5_smp <- sample(seq_len(nrow(subpool5_epi)), 
                 size = s5_smp_size)

subpool5_ncw_test_epi <- subpool5_epi[s5_smp, ]
subpool5_ncw_train_epi <- subpool5_epi[-s5_smp, ]

subpool5_ncw_test_gen <- subpool5_gen[s5_smp, ]
subpool5_ncw_train_gen <- subpool5_gen[-s5_smp, ]


#Figure 5C

#Fit log-linear model to the different locations of CRE (site 1-6) and to the 
#different backgrounds, allowing weights per weak and consensus CRE per site

ind_site_ind_back <- function(df) {
  model <- lm(ave_ratio ~ background + site1 + site2 + site3 + site4 + site5 + site6, 
              data = df)
}

#Attach model predictions to original df

pred_resid <- function(df1, x) {
  df2 <- df1 %>%
    add_predictions(x)
  df3 <- df2 %>%
    add_residuals(x)
  return(df3)
  print('processed pred_resid(df1, df2) in order of (data, model)')
}

#Make model weights output plottable

ind_site_ind_back_weights <- function(df) {
  df_conf <-confint_tidy(df)
  df_tidy <- tidy(df)
  df_comb <- cbind(df_tidy, df_conf)
  sites <- df_comb %>%
    filter(str_detect(term, '^site')) %>%
    mutate(term = gsub('consensus', '_consensus', term)) %>%
    mutate(term = gsub('weak', '_weak', term)) %>%
    separate(term, into = c('variable', 'type'), sep = "_")
  background <- df_comb %>%
    filter(str_detect(term, '^background')) %>%
    mutate(term = gsub('background', 'background_', term)) %>%
    separate(term, into = c('variable', 'type'), sep = '_')
  weights <- rbind(sites, background)
  return(weights)
}

#train independent site independent background model

ind_site_ind_back_epi <- subpool5_ncw_train_epi %>%
  filter(MPRA == 'episomal') %>%
  ind_site_ind_back()

#test independent site independent background model

ind_site_ind_back_epi_test_p_r <- pred_resid(filter(subpool5_ncw_test_epi, 
                                                    MPRA == 'episomal'), 
                                              ind_site_ind_back_epi)

#make model fits graphs

ind_site_ind_back_weights_epi <- ind_site_ind_back_weights(ind_site_ind_back_epi)

ind_site_ind_back_anova_epi <- tidy(anova(ind_site_ind_back_epi)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq)

#Plot

background_viridis <- c('#440154FF', '#33638DFF', '#29AF7FFF')

p_ind_site_ind_back_epi <- ind_site_ind_back_epi_test_p_r %>%
  mutate(background = factor(background, levels = c('52', '55', '41')))  %>%
  ggplot(aes(ave_ratio, pred, color = background)) +
  geom_point(alpha = 0.25, size = 0.8, show.legend = FALSE) +
  scale_color_manual(values = background_viridis) +
  scale_x_continuous(name = 'Average log10 expression (a.u.)', breaks = c(-2:2),
                     limits = c(-2, 2)) + 
  scale_y_continuous(name = 'Predicted log10 expression (a.u.)', breaks = c(-2:2),
                     limits = c(-2, 2)) +
  annotation_logticks(sides = 'bl') +
  figurefont_theme

rss <- summarize(ind_site_ind_back_epi_test_p_r, sum((resid)^2))

tss <- summarize(ind_site_ind_back_epi_test_p_r, sum((ave_ratio - mean(ave_ratio))^2))

1-rss/tss

p_ind_site_ind_back_weights_epi <- ind_site_ind_back_weights_epi %>%
  mutate(type = factor(type, 
                       levels = c('52', '55', 'consensus', 'weak'))) %>%
  ggplot() + 
  geom_pointrange(aes(variable, estimate, ymin = conf.low, ymax = conf.high,
                    fill = type), fatten = 2, shape = 21, stroke = 0.5, size = 1) +
  geom_hline(yintercept = 0, size = 0.25) +
  scale_x_discrete(position = 'bottom') + 
  scale_fill_viridis(discrete = TRUE) + 
  theme(axis.ticks.x = element_blank(), legend.position = 'top',
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank()) +
  ylab('Weight') +
  figurefont_theme

p_ind_site_ind_back_anova_epi <- ind_site_ind_back_anova_epi %>%
  ggplot(aes(term, per_sumsq)) + 
  geom_bar(stat = 'identity') + 
  ylab('Proportion of\nvariance explained') +
  scale_y_continuous(limits = c(0, 0.26), 
                     breaks = seq(from = 0, to = 0.25, by = 0.05)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  figurefont_theme

ggsave('../plots/p_ind_site_ind_back_epi.png', p_ind_site_ind_back_epi, 
       width = 2.3, height = 2.2, units = 'in')

ggsave('../plots/p_ind_site_ind_back_weights_epi.pdf', 
       p_ind_site_ind_back_weights_epi,
       width = 2.75, height = 2.6, units = 'in')

ggsave('../plots/p_ind_site_ind_back_anova_epi.pdf', 
       p_ind_site_ind_back_anova_epi,
       width = 2.5, height = 2.4)

#Fit log-linear model to genomic data, join model predictions to df and 
#determine the proportion of variance explained

ind_site_ind_back_gen <- subpool5_ncw_train_gen %>%
  filter(MPRA == 'genomic') %>%
  ind_site_ind_back()

ind_site_ind_back_gen_test_p_r <- pred_resid(filter(subpool5_ncw_test_gen, MPRA == 'genomic'), 
                                        ind_site_ind_back_gen)

ind_site_ind_back_weights_gen <- ind_site_ind_back_weights(ind_site_ind_back_gen)

ind_site_ind_back_anova_gen <- tidy(anova(ind_site_ind_back_gen)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq)

#Plot

p_ind_site_ind_back_gen <-  ind_site_ind_back_gen_test_p_r %>%
  mutate(background = factor(background, levels = c('52', '55', '41')))  %>%
  ggplot(aes(ave_ratio, pred, color = background)) +
  geom_point(alpha = 0.25, size = 0.8, show.legend = FALSE) +
  scale_color_manual(values = background_viridis) +
  scale_x_continuous(name = 'Average log10 expression (a.u.)', breaks = c(-2:2),
                     limits = c(-2.1, 2)) + 
  scale_y_continuous(name = 'Predicted log10 expression (a.u.)', breaks = c(-2:2),
                     limits = c(-2.1, 2)) +
  annotation_logticks(sides = 'bl') +
  figurefont_theme

rss <- summarize(ind_site_ind_back_gen_test_p_r, sum((resid)^2))

tss <- summarize(ind_site_ind_back_gen_test_p_r, sum((ave_ratio - mean(ave_ratio))^2))

1-rss/tss

p_ind_site_ind_back_weights_gen <- ind_site_ind_back_weights_gen %>%
  mutate(type = factor(type, 
                       levels = c('52', '55', 'consensus', 'weak'))) %>%
  ggplot() + 
  geom_pointrange(aes(variable, estimate, ymin = conf.low, ymax = conf.high,
                      fill = type), fatten = 2, shape = 21, stroke = 0.5, size = 1) +
  geom_hline(yintercept = 0, size = 0.25) +
  scale_x_discrete(position = 'bottom') + 
  scale_fill_viridis(discrete = TRUE) + 
  theme(axis.ticks.x = element_blank(), legend.position = 'top',
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank()) +
  ylab('Weight') +
  figurefont_theme

p_ind_site_ind_back_anova_gen <- ind_site_ind_back_anova_gen %>%
  ggplot(aes(term, per_sumsq)) + 
  geom_bar(stat = 'identity') + 
  ylab('Proportion of\nvariance explained') +
  scale_y_continuous(limits = c(0, 0.26), 
                     breaks = seq(from = 0, to = 0.25, by = 0.05)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  figurefont_theme

ggsave('../plots/p_ind_site_ind_back_gen.png', p_ind_site_ind_back_gen, 
       width = 2.3, height = 2.2, units = 'in')

ggsave('../plots/p_ind_site_ind_back_weights_gen.pdf', 
       p_ind_site_ind_back_weights_gen,
       width = 2.75, height = 2.6, units = 'in')

ggsave('../plots/p_ind_site_ind_back_anova_gen.pdf', 
       p_ind_site_ind_back_anova_gen,
       width = 2.5, height = 2.4)


#Supplemental Figure 5A

#site number and bg (independent)

site_ind_back = function(df) {
  model <- lm(ave_ratio ~ background + consensus + weak, data = df)
}

site_ind_back_epi <- subpool5_ncw_train_epi %>%
  filter(MPRA == 'episomal') %>%
  site_ind_back()

site_ind_back_p_r_epi <- pred_resid(filter(subpool5_ncw_test_epi, 
                                           MPRA == 'episomal'), 
                                        site_ind_back_epi)

#Plot

lessthan1_2color <- c('red', 'black', 'black', 'black', 'black', 'black', 'black')

p_site_ind_back_epi <- ggplot(site_ind_back_p_r_epi, 
                                  aes(ave_ratio, pred, 
                                      color = as.factor(consensus))) +
  geom_point(alpha = 0.1, size = 0.75, show.legend = FALSE) +
  scale_color_manual(values = lessthan1_2color) +
  scale_x_continuous(name = 'Average log10 expression (a.u.)', breaks = c(-2:2),
                     limits = c(-2, 2)) + 
  scale_y_continuous(name = 'Predicted log10 expression (a.u.)', breaks = c(-2:2),
                     limits = c(-2, 2)) +
  annotation_logticks(sides = 'bl') +
  figurefont_theme

rss <- summarize(site_ind_back_p_r_epi, sum((resid)^2))

tss <- summarize(site_ind_back_p_r_epi, sum((ave_ratio - mean(ave_ratio))^2))

1-rss/tss

ggsave('../plots/p_site_ind_back_epi.png', p_site_ind_back_epi, 
       width = 2.3, height = 2.2, units = 'in')


#Supplemental figure 5B

#Variance explained episomal

#Determine variance explained per feature

model_back <- function(df) {
  model <- lm(ave_ratio ~ background, data = df)
}

model_back_epi <- subpool5_ncw_train_epi %>%
  filter(MPRA == 'episomal') %>%
  model_back()

model_back_epi_variance <- tidy(anova(model_back_epi)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq) %>%
  filter(term != 'Residuals')


model_cons <- function(df) {
  model <- lm(ave_ratio ~ consensus, data = df)
}

model_cons_epi <- subpool5_ncw_train_epi %>%
  filter(MPRA == 'episomal') %>%
  model_cons()

model_cons_epi_variance <- tidy(anova(model_cons_epi)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq) %>%
  filter(term != 'Residuals')


model_weak <- function(df) {
  model <- lm(ave_ratio ~ weak, data = df)
}

model_weak_epi <- subpool5_ncw_train_epi %>%
  filter(MPRA == 'episomal') %>%
  model_weak()

model_weak_epi_variance <- tidy(anova(model_weak_epi)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq) %>%
  filter(term != 'Residuals')



model_site1 <- function(df) {
  model <- lm(ave_ratio ~ site1, data = df)
}

model_site1_epi <- subpool5_ncw_train_epi %>%
  filter(MPRA == 'episomal') %>%
  model_site1()

model_site1_epi_variance <- tidy(anova(model_site1_epi)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq) %>%
  filter(term != 'Residuals')


model_site2 <- function(df) {
  model <- lm(ave_ratio ~ site2, data = df)
}

model_site2_epi <- subpool5_ncw_train_epi %>%
  filter(MPRA == 'episomal') %>%
  model_site2()

model_site2_epi_variance <- tidy(anova(model_site2_epi)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq) %>%
  filter(term != 'Residuals')


model_site3 <- function(df) {
  model <- lm(ave_ratio ~ site3, data = df)
}

model_site3_epi <- subpool5_ncw_train_epi %>%
  filter(MPRA == 'episomal') %>%
  model_site3()

model_site3_epi_variance <- tidy(anova(model_site3_epi)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq) %>%
  filter(term != 'Residuals')


model_site4 <- function(df) {
  model <- lm(ave_ratio ~ site4, data = df)
}

model_site4_epi <- subpool5_ncw_train_epi %>%
  filter(MPRA == 'episomal') %>%
  model_site4()

model_site4_epi_variance <- tidy(anova(model_site4_epi)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq) %>%
  filter(term != 'Residuals')


model_site5 <- function(df) {
  model <- lm(ave_ratio ~ site5, data = df)
}

model_site5_epi <- subpool5_ncw_train_epi %>%
  filter(MPRA == 'episomal') %>%
  model_site5()

model_site5_epi_variance <- tidy(anova(model_site5_epi)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq) %>%
  filter(term != 'Residuals')


model_site6 <- function(df) {
  model <- lm(ave_ratio ~ site6, data = df)
}

model_site6_epi <- subpool5_ncw_train_epi %>%
  filter(MPRA == 'episomal') %>%
  model_site6()

model_site6_epi_variance <- tidy(anova(model_site6_epi)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq) %>%
  filter(term != 'Residuals')

var_variance_epi <- rbind(model_back_epi_variance, model_cons_epi_variance, 
                          model_weak_epi_variance, model_site1_epi_variance, 
                          model_site2_epi_variance, model_site3_epi_variance,
                          model_site4_epi_variance, model_site5_epi_variance, 
                          model_site6_epi_variance)

p_var_variance_epi <- var_variance_epi %>%
  ggplot(aes(term, per_sumsq)) + 
  geom_bar(stat = 'identity') + 
  scale_y_continuous(limits = c(0, 0.73)) +
  ylab('Proportion of\nvariance explained') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  figurefont_theme

ggsave('../plots/p_var_variance_epi.pdf', 
       p_var_variance_epi,
       width = 2.75, height = 1.5)


#Variance explained genomic

#Determine vriance explained per feature

model_back_gen <- subpool5_ncw_train_gen %>%
  filter(MPRA == 'genomic') %>%
  model_back()

model_back_gen_variance <- tidy(anova(model_back_gen)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq) %>%
  filter(term != 'Residuals')


model_cons_gen <- subpool5_ncw_train_gen %>%
  filter(MPRA == 'genomic') %>%
  model_cons()

model_cons_gen_variance <- tidy(anova(model_cons_gen)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq) %>%
  filter(term != 'Residuals')


model_weak_gen <- subpool5_ncw_train_gen %>%
  filter(MPRA == 'genomic') %>%
  model_weak()

model_weak_gen_variance <- tidy(anova(model_weak_gen)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq) %>%
  filter(term != 'Residuals')


model_site1_gen <- subpool5_ncw_train_gen %>%
  filter(MPRA == 'genomic') %>%
  model_site1()

model_site1_gen_variance <- tidy(anova(model_site1_gen)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq) %>%
  filter(term != 'Residuals')


model_site2_gen <- subpool5_ncw_train_gen %>%
  filter(MPRA == 'genomic') %>%
  model_site2()

model_site2_gen_variance <- tidy(anova(model_site2_gen)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq) %>%
  filter(term != 'Residuals')


model_site3_gen <- subpool5_ncw_train_gen %>%
  filter(MPRA == 'genomic') %>%
  model_site3()

model_site3_gen_variance <- tidy(anova(model_site3_gen)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq) %>%
  filter(term != 'Residuals')


model_site4_gen <- subpool5_ncw_train_gen %>%
  filter(MPRA == 'genomic') %>%
  model_site4()

model_site4_gen_variance <- tidy(anova(model_site4_gen)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq) %>%
  filter(term != 'Residuals')


model_site5_gen <- subpool5_ncw_train_gen %>%
  filter(MPRA == 'genomic') %>%
  model_site5()

model_site5_gen_variance <- tidy(anova(model_site5_gen)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq) %>%
  filter(term != 'Residuals')


model_site6_gen <- subpool5_ncw_train_gen %>%
  filter(MPRA == 'genomic') %>%
  model_site6()

model_site6_gen_variance <- tidy(anova(model_site6_gen)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq) %>%
  filter(term != 'Residuals')

var_variance_gen <- rbind(model_back_gen_variance, model_cons_gen_variance, 
                          model_weak_gen_variance, model_site1_gen_variance, 
                          model_site2_gen_variance, model_site3_gen_variance,
                          model_site4_gen_variance, model_site5_gen_variance, 
                          model_site6_gen_variance)

p_var_variance_gen <- var_variance_gen %>%
  ggplot(aes(term, per_sumsq)) + 
  geom_bar(stat = 'identity') + 
  scale_y_continuous(limits = c(0, 0.73)) +
  ylab('Proportion of\nvariance explained') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  figurefont_theme

ggsave('../plots/p_var_variance_gen.pdf', 
       p_var_variance_gen,
       width = 2.75, height = 1.5)


#Supplemental Figure 5D, 5E and background analyses-----------------------------

#Supplemental Figure 5D

#Fitting model to new backgrounds

set.seed(1234)

sixsite_newback_norm_sub <- removed_bgs_sixsite_newback_norm %>%
  mutate(site1 = gsub('nosite', 'anosite', site1)) %>%
  mutate(site2 = gsub('nosite', 'anosite', site2)) %>%
  mutate(site3 = gsub('nosite', 'anosite', site3)) %>%
  mutate(site4 = gsub('nosite', 'anosite', site4)) %>%
  mutate(site5 = gsub('nosite', 'anosite', site5)) %>%
  mutate(site6 = gsub('nosite', 'anosite', site6)) %>%
  mutate(background = as.character(background)) %>%
  mutate(ave_med_ratio_4 = log10(ave_med_ratio_4))

sixsite_smp_size <- floor(0.8 * nrow(sixsite_newback_norm_sub))

sixsite_smp <- sample(seq_len(nrow(sixsite_newback_norm_sub)), 
                      size = sixsite_smp_size)

sixsite_test_epi <- sixsite_newback_norm_sub[sixsite_smp, ]
sixsite_train_epi <- sixsite_newback_norm_sub[-sixsite_smp, ]

new_sixsite_weights <- function(df) {
  df <- tidy(df)
  sites <- df %>%
    filter(str_detect(term, '^site')) %>%
    mutate(term = gsub('consensus', '_consensus', term)) %>%
    separate(term, into = c('variable', 'type'), sep = "_")
  background <- df %>%
    filter(str_detect(term, '^background')) %>%
    mutate(term = gsub('background', 'background_', term)) %>%
    separate(term, into = c('variable', 'type'), sep = '_')
  weights <- rbind(sites, background)
  return(weights)
}

#Fit log-linear model to episomal data, join model predictions to df and 
#determine the proportion of variance explained

ind_site_ind_back_new <- function(df) {
  model <- lm(ave_med_ratio_4 ~ background + site1 + site2 + site3 + site4 + site5 + site6, 
              data = df)
}

ind_site_ind_back_epi_new <- sixsite_train_epi %>%
  ind_site_ind_back_new()

ind_site_ind_back_p_r_epi_new <- pred_resid(sixsite_test_epi, 
                                            ind_site_ind_back_epi_new)

ind_site_ind_back_weights_epi_new <- new_sixsite_weights(ind_site_ind_back_epi_new)

ind_site_ind_back_anova_epi_new <- tidy(anova(ind_site_ind_back_epi_new)) %>%
  mutate(term = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq)

#Plot

lessthan1_2color <- c('red', 'black', 'black', 'black', 'black', 'black', 'black')

p_ind_site_ind_back_epi_new <- ggplot(ind_site_ind_back_p_r_epi_new, 
                                      aes(ave_med_ratio_4, pred, 
                                          color = as.factor(consensus))) +
  geom_point(alpha = 0.1, size = 0.75, show.legend = FALSE) +
  scale_color_manual(values = lessthan1_2color) +
  scale_x_continuous(name = 'Average log10 expression (a.u.)', breaks = c(-2:2),
                     limits = c(-1.6, 1.8)) + 
  scale_y_continuous(name = 'Predicted log10 expression (a.u.)', breaks = c(-2:2),
                     limits = c(-1.6, 1.8)) +
  annotation_logticks(sides = 'bl') +
  figurefont_theme

rss <- summarize(ind_site_ind_back_p_r_epi_new, sum((resid)^2))

tss <- summarize(ind_site_ind_back_p_r_epi_new, sum((ave_med_ratio_4 - mean(ave_med_ratio_4))^2))

1-rss/tss

p_ind_site_ind_back_anova_epi_new <- ind_site_ind_back_anova_epi_new %>%
  ggplot(aes(term, per_sumsq)) + 
  geom_bar(stat = 'identity') + 
  ylab('Proportion of\nvariance explained') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  figurefont_theme

ggsave('../plots/p_ind_site_ind_back_epi_new.png', p_ind_site_ind_back_epi_new, 
       width = 2.3, height = 2.2, units = 'in')

ggsave('../plots/p_ind_site_ind_back_anova_epi_new.pdf', 
       p_ind_site_ind_back_anova_epi_new,
       width = 2.5, height = 2.4)

#Supplemental Figure 5E

#Determine relationship between total GC and CG content of regulatory elements
#and expression. Perform analysis on variants with high expression and no
#variation in CRE combinations (consensus == 6). Need to add in final base of
#KpnI cut site at the beginning of the 3' to 5' variant sequence as can form a 
#CG dinucleotide with variant sequence.

sixsite_gc_cg <- removed_bgs_sixsite_newback_norm %>%
  mutate(most_common_kpnI = str_c('C', most_common)) %>%
  mutate(most_common_200 = str_c('CTCGAGGCTAGCGAGCTCAGGTACC', most_common, 
                                 'ACGCGTGCTCTACGACTATGCTCTT')) %>%
  mutate(cg = str_count(most_common_kpnI, pattern = 'CG')) %>%
  mutate(c = str_count(most_common_200, pattern = 'C')) %>%
  mutate(g = str_count(most_common_200, pattern = 'G')) %>%
  mutate(gc_200 = (str_count(most_common_200, pattern = 'C') + str_count(most_common_200, pattern = 'G'))/200) %>%
  mutate(obs_exp_CpG = (cg*200)/(c*g)) %>%
  mutate(A5 = str_count(most_common,
                        pattern = 'AAAAA') + str_count(most_common, 
                                                       pattern = 'TTTTT'))

#Look at pearsons correlation of CG dinucleotides and plot

p_back_gc_consensus <- sixsite_gc_cg %>%
  ggplot(aes(x = gc)) +
  facet_grid(~ consensus) +
  geom_point(aes(y = ave_med_ratio_4), alpha = 0.1, size = 0.75) +
  scale_y_log10() +
  panel_border(colour = 'black') +
  annotation_logticks(sides = 'l') +
  guides(color = guide_colorbar(frame.colour = 'black', 
                                ticks.colour = 'black')) +
  ylab('Average normalized\nexpression (a.u.)') +
  xlab('GC content') +
  figurefont_theme +
  theme(strip.background = element_rect(colour="black", fill="white"))

ggsave('../plots/p_back_gc_consensus.png', p_back_gc_consensus, 
       width = 5.5, height = 1.75, units = 'in')

p_back_cpg_consensus <- sixsite_gc_cg %>%
  filter(gc_200 > 0.50) %>%
  ggplot(aes(x = obs_exp_CpG)) +
  facet_grid(~ consensus) +
  geom_point(aes(y = ave_med_ratio_4), alpha = 0.1, size = 0.75) +
  scale_y_log10() +
  scale_x_continuous(breaks = c(0.5, 1)) +
  panel_border(colour = 'black') +
  annotation_logticks(sides = 'l') +
  guides(color = guide_colorbar(frame.colour = 'black', 
                                ticks.colour = 'black')) +
  ylab('Average normalized\nexpression (a.u.)') +
  xlab('Observed/Expected CpG') +
  figurefont_theme +
  theme(strip.background = element_rect(colour="black", fill="white"))

ggsave('../plots/p_back_cpg_consensus.png', p_back_cpg_consensus, 
       width = 5.5, height = 1.5, units = 'in')

#Extra analyses

#Count 3 bp flanking motif around sites

sixsite_newback_motif <- removed_bgs_sixsite_newback_norm %>%
  mutate(site6_extra = 'C') %>%
  mutate(site6_5 = str_sub(most_common, 1, 2)) %>%
  mutate(site6_5 = str_c(site6_extra, site6_5)) %>%
  select(-site6_extra) %>%
  mutate(site6_3 = str_sub(most_common, 11, 13)) %>%
  mutate(site5_5 = str_sub(most_common, 25, 27)) %>%
  mutate(site5_3 = str_sub(most_common, 36, 38)) %>%
  mutate(site4_5 = str_sub(most_common, 50, 52)) %>%
  mutate(site4_3 = str_sub(most_common, 61, 63)) %>%
  mutate(site3_5 = str_sub(most_common, 75, 77)) %>%
  mutate(site3_3 = str_sub(most_common, 86, 88)) %>%
  mutate(site2_5 = str_sub(most_common, 100, 102)) %>%
  mutate(site2_3 = str_sub(most_common, 111, 113)) %>%
  mutate(site1_5 = str_sub(most_common, 125, 127)) %>%
  mutate(site1_3 = str_sub(most_common, 136, 138))

#Look at local GC content around each site (4 bp on either side)

sixsite_newback_norm_localgc <- removed_bgs_sixsite_newback_norm %>%
  mutate(site6_extra = 'CC') %>%
  mutate(site6_window = str_sub(most_common, 1, 14)) %>%
  mutate(site6_window = str_c(site6_extra, site6_window)) %>%
  select(-site6_extra) %>%
  mutate(site5_window = str_sub(most_common, 24, 39)) %>%
  mutate(site4_window = str_sub(most_common, 49, 64)) %>%
  mutate(site3_window = str_sub(most_common, 74, 89)) %>%
  mutate(site2_window = str_sub(most_common, 99, 114)) %>%
  mutate(site1_window = str_sub(most_common, 124, 139)) %>%
  mutate(site6_gc = (str_count(site6_window, 'C') + str_count(site6_window, 'G'))/16)%>%
  mutate(site5_gc = (str_count(site5_window, 'C') + str_count(site5_window, 'G'))/16)%>%
  mutate(site4_gc = (str_count(site4_window, 'C') + str_count(site4_window, 'G'))/16)%>%
  mutate(site3_gc = (str_count(site3_window, 'C') + str_count(site3_window, 'G'))/16)%>%
  mutate(site2_gc = (str_count(site2_window, 'C') + str_count(site2_window, 'G'))/16)%>%
  mutate(site1_gc = (str_count(site1_window, 'C') + str_count(site1_window, 'G'))/16)

#Local GC content plots and anovas

library(export)

#site 6

aov_site6 <- tidy(aov(ave_med_ratio_norm_4 ~ as.factor(site6_gc), 
                      data = filter(sixsite_newback_norm_localgc, 
                                    consensus == 6))) %>%
  mutate(term = gsub('as\\.factor\\(site6\\_gc\\)', 'GC content', term))

p_site6_localgc <- sixsite_newback_norm_localgc %>%
  filter(consensus == 6) %>%
  ggplot(aes(site6_gc, ave_med_ratio_norm_4)) +
  geom_violin(aes(group = site6_gc), fill = 'gray93') +
  geom_boxplot(aes(group = site6_gc), alpha = 0.5, outlier.alpha = 0, 
               size = 0.25) +
  xlab('Local GC content') +
  ylab('Average normalized\nexpression (a.u.)') +
  figurefont_theme 

table2csv(aov_site6, '../tables/aov_site6.csv', digits = 16)

ggsave('../plots/p_site6_localgc.pdf', p_site6_localgc, width = 2.5, height = 2,
       units = 'in')

#site 5

aov_site5 <- tidy(aov(ave_med_ratio_norm_4 ~ as.factor(site5_gc), 
                      data = filter(sixsite_newback_norm_localgc, 
                                    consensus == 6))) %>%
  mutate(term = gsub('as\\.factor\\(site5\\_gc\\)', 'GC content', term))

p_site5_localgc <- sixsite_newback_norm_localgc %>%
  filter(consensus == 6) %>%
  ggplot(aes(site5_gc, ave_med_ratio_norm_4)) +
  geom_violin(aes(group = site5_gc), fill = 'gray93') +
  geom_boxplot(aes(group = site5_gc), alpha = 0.5, outlier.alpha = 0, 
               size = 0.25) +
  xlab('Local GC content') +
  ylab('Average normalized\nexpression (a.u.)') +
  figurefont_theme 

table2csv(aov_site5, '../tables/aov_site5.csv', digits = 16)

ggsave('../plots/p_site5_localgc.pdf', p_site5_localgc, width = 2.5, height = 2,
       units = 'in')

#site 4

aov_site4 <- tidy(aov(ave_med_ratio_norm_4 ~ as.factor(site4_gc), 
                      data = filter(sixsite_newback_norm_localgc, 
                                    consensus == 6))) %>%
  mutate(term = gsub('as\\.factor\\(site4\\_gc\\)', 'GC content', term))

p_site4_localgc <- sixsite_newback_norm_localgc %>%
  filter(consensus == 6) %>%
  ggplot(aes(site4_gc, ave_med_ratio_norm_4)) +
  geom_violin(aes(group = site4_gc), fill = 'gray93') +
  geom_boxplot(aes(group = site4_gc), alpha = 0.5, outlier.alpha = 0, 
               size = 0.25) +
  xlab('Local GC content') +
  ylab('Average normalized\nexpression (a.u.)') +
  figurefont_theme

table2csv(aov_site4, '../tables/aov_site4.csv', digits = 16)

ggsave('../plots/p_site4_localgc.pdf', p_site4_localgc, width = 2.5, height = 2,
       units = 'in')

#site 3

aov_site3 <- tidy(aov(ave_med_ratio_norm_4 ~ as.factor(site3_gc), 
                      data = filter(sixsite_newback_norm_localgc, 
                                    consensus == 6))) %>%
  mutate(term = gsub('as\\.factor\\(site3\\_gc\\)', 'GC content', term))

p_site3_localgc <- sixsite_newback_norm_localgc %>%
  filter(consensus == 6) %>%
  ggplot(aes(site3_gc, ave_med_ratio_norm_4)) +
  geom_violin(aes(group = site3_gc), fill = 'gray93') +
  geom_boxplot(aes(group = site3_gc), alpha = 0.5, outlier.alpha = 0, 
               size = 0.25) +
  xlab('Local GC content') +
  ylab('Average normalized\nexpression (a.u.)') +
  figurefont_theme 

table2csv(aov_site3, '../tables/aov_site3.csv', digits = 16)

ggsave('../plots/p_site3_localgc.pdf', p_site3_localgc, width = 2.5, height = 2,
       units = 'in')

#site 2

aov_site2 <- tidy(aov(ave_med_ratio_norm_4 ~ as.factor(site2_gc), 
                      data = filter(sixsite_newback_norm_localgc, 
                                    consensus == 6))) %>%
  mutate(term = gsub('as\\.factor\\(site2\\_gc\\)', 'GC content', term))

p_site2_localgc <- sixsite_newback_norm_localgc %>%
  filter(consensus == 6) %>%
  ggplot(aes(site2_gc, ave_med_ratio_norm_4)) +
  geom_violin(aes(group = site2_gc), fill = 'gray93') +
  geom_boxplot(aes(group = site2_gc), alpha = 0.5, outlier.alpha = 0, 
               size = 0.25) +
  xlab('Local GC content') +
  ylab('Average normalized\nexpression (a.u.)') +
  figurefont_theme

table2csv(aov_site2, '../tables/aov_site2.csv', digits = 16)

ggsave('../plots/p_site2_localgc.pdf', p_site2_localgc, width = 2.5, height = 2,
       units = 'in')

#site 1

aov_site1 <- tidy(aov(ave_med_ratio_norm_4 ~ as.factor(site1_gc), 
                      data = filter(sixsite_newback_norm_localgc, 
                                    consensus == 6))) %>%
  mutate(term = gsub('as\\.factor\\(site1\\_gc\\)', 'GC content', term))

p_site1_localgc <- sixsite_newback_norm_localgc %>%
  filter(consensus == 6) %>%
  ggplot(aes(site1_gc, ave_med_ratio_norm_4)) +
  geom_violin(aes(group = site1_gc), fill = 'gray93') +
  geom_boxplot(aes(group = site1_gc), alpha = 0.5, outlier.alpha = 0, 
               size = 0.25) +
  xlab('Local GC content') +
  ylab('Average normalized\nexpression (a.u.)') +
  figurefont_theme

table2csv(aov_site1, '../tables/aov_site1.csv', digits = 16)

ggsave('../plots/p_site1_localgc.pdf', p_site1_localgc, width = 2.5, height = 2,
       units = 'in')

#Figure 6, Supplemental Figure 6, and Supplemental Table 2----------------------

#Make linear regression predicting episomal expression from genomic expression
#using variants containing only consensus CREs. Use this regression to plot 
#correlation line and determine residuals to relationship

s5_cons_log10 <- s5_gen_epi %>%
  var_log10() %>%
  filter(site_combo == 'consensus')

cons_int_epi_lm <- lm(ave_ratio_22_norm ~ ave_med_ratio_norm, 
                      data = s5_cons_log10)

#Add model predictions and residuals to consensus only (R2 value) and to all 
#variants in library

pred_resid <- function(df1, x) {
  df2 <- df1 %>%
    add_predictions(x)
  df3 <- df2 %>%
    add_residuals(x)
  return(df3)
  print('processed pred_resid(df1, df2) in order of (data, model)')
}

s5_gen_epi_cons_lm <- pred_resid(s5_cons_log10, cons_int_epi_lm)

round(cor(s5_gen_epi_cons_lm$ave_ratio_22,
          s5_gen_epi_cons_lm$pred,
          use = "pairwise.complete.obs", 
          method = "pearson")^2, 2)

s5_gen_epi_all_lm <- pred_resid(var_log10(s5_gen_epi), cons_int_epi_lm)

#Figure 6A, plot genomic vs. episomal expression and linear regression as
#reference. Not plotting variants wihtout any CRE (backgrounds) as they do not
#fall into CRE categories, the three points don't seem biased in expression
#around line though.

p_s5_int_trans_site_combo <- s5_gen_epi_all_lm %>%
  filter(site_combo != 'none') %>%
  mutate(site_combo = factor(site_combo, 
                             levels = c('consensus', 'weak', 'mixed'))) %>%
  ggplot(aes(ave_med_ratio_norm, ave_ratio_22_norm)) +
  facet_grid(. ~ site_combo) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_line(aes(ave_med_ratio_norm, pred), color = 'red', size = 0.5) +
  annotation_logticks() +
  scale_y_continuous(breaks = seq(from = 0, to = 2, by =1)) +
  xlab('Average normalized log10 genomic expression (a.u.)') +
  ylab('Average normalized log10\nepisomal expression (a.u.)') +
  panel_border(colour = 'black') +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

ggsave('../plots/p_s5_int_trans_site_combo.pdf', p_s5_int_trans_site_combo,
       width = 4, height = 2.25, units = 'in')

#determine percent below line in mixed

mixed <- s5_gen_epi_all_lm %>%
  filter(site_combo == 'mixed')

m <- count(mixed)
n <- count(filter(mixed, resid < 0))
n/m

#Figure 6B

#Plot residual of each combination of CREs per variant to linear relationship 
#between MPRA expression of variants with consensus CREs only. Linear 
#relationship indicated with a red reference line. Negative residuals indicate 
#higher relative expression of variant in the genomic MPRA and positive
#indicates higher relative expression in episomal MPRA

p_s5_gen_epi_site_combo_resid <- s5_gen_epi_all_lm %>%
  ggplot(aes(as.factor(consensus), resid, fill = as.factor(weak))) +
  geom_boxplot(outlier.size = 1, size = 0.3, 
               outlier.shape = 21, outlier.alpha = 1, 
               position = position_dodge(0.75)) +
  scale_fill_manual(name = 'number of\nweak CREs', 
                    values = cbPalette7_grad_light) +
  xlab("consensus CREs") +
  ylab("residual") +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5, color = 'red') +
  theme(legend.position = 'top', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

ggsave('../plots/p_s5_gen_epi_site_combo_resid.pdf', 
       p_s5_gen_epi_site_combo_resid, width = 3.75, height = 3, units = 'in')


#anova for Supplemental Table 2

aov_s5_gen_epi_all_lm <- aov(resid ~ as.factor(consensus) * as.factor(weak), 
                         data = s5_gen_epi_all_lm)

test <- s5_gen_epi_all_lm %>%
  add_residuals(aov_s5_gen_epi_all_lm)

write_csv(tidy(aov_s5_gen_epi_all_lm),
          '../tables/aov_s5_gen_epi_all_lm.csv')

plot(aov_s5_gen_epi_all_lm, 1)

library(car)

levene_s5_gen_epi_all_lm <- leveneTest(resid ~ as.factor(consensus) * as.factor(weak), 
                                       data = s5_gen_epi_all_lm)

write_csv(tidy(levene_s5_gen_epi_all_lm),
           '../tables/levene_s5_gen_epi_all_lm.csv')

qq_aov_s5_gen_epi_all_lm <- plot(aov_s5_gen_epi_all_lm, 2)

aov_s5_gen_epi_all_lm_resid <- residuals(object =  aov_s5_gen_epi_all_lm)

#Supplemental Figure 6

#Make linear regression predicting episomal expression from genomic expression
#using variants containing only consensus CREs. Use this regression to plot 
#correlation line and determine residuals to relationship

abline_lm <- lm(ave_ratio_22_norm ~ ave_med_ratio_norm, s5_gen_epi)

#Add model predictions and residuals to consensus only (R2 value) and to all 
#variants in library

pred_resid <- function(df1, x) {
  df2 <- df1 %>%
    add_predictions(x)
  df3 <- df2 %>%
    add_residuals(x)
  return(df3)
  print('processed pred_resid(df1, df2) in order of (data, model)')
}

s5_gen_epi_lm <- pred_resid(s5_gen_epi, abline_lm)

round(cor(s5_gen_epi_lm$ave_ratio_22_norm,
          s5_gen_epi_lm$pred,
          use = "pairwise.complete.obs", 
          method = "pearson")^2, 3)

#Supplemental Fig. 6A

#Plot genomic vs. episomal expression and linear regression as reference. Not 
#plotting variants wihtout any CRE (backgrounds) as they do not fall into CRE 
#categories, the three points don't seem biased in expression around line.

p_s5_int_trans_site_combo <- s5_gen_epi_lm %>%
  filter(site_combo != 'none') %>%
  ggplot(aes(ave_med_ratio_norm, ave_ratio_22_norm)) +
  geom_point(alpha = 0.1, size = 0.5)  +
  geom_line(aes(ave_med_ratio_norm, pred), color = 'red') +
  xlab('Average normalized\ngenomic expression (a.u.)') +
  ylab('Average normalized\nepisomal expression (a.u.)') +
  theme(legend.position = 'right',
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

ggsave('../plots/p_s5_gen_epi_abline_linear.png', p_s5_int_trans_site_combo,
       width = 2.3, height = 2.2, units = 'in')

#Supplemental Figure 6B

#Plot residual of each combination of CREs per variant to linear relationship 
#between MPRA expression of variants with consensus CREs only. Linear 
#relationship indicated with a red reference line. Negative residuals indicate 
#higher relative expression of variant in the genomic MPRA and positive
#indicates higher relative expression in episomal MPRA

p_s5_gen_epi_site_combo_resid <- s5_gen_epi_lm %>%
  filter(site_combo != 'none') %>%
  ggplot(aes(as.factor(consensus), resid, fill = as.factor(weak))) +
  geom_boxplot(outlier.size = 0.5, size = 0.3, 
               outlier.shape = 21, outlier.alpha = 1, 
               position = position_dodge(0.75)) +
  scale_fill_manual(name = 'number of\nweak CREs', 
                    values = cbPalette7_grad_light) +
  xlab("consensus CREs") +
  ylab("residual") +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5, color = 'red') +
  theme(legend.position = 'top', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

ggsave('../plots/p_s5_gen_epi_site_combo_resid_abline_linear.pdf', 
       p_s5_gen_epi_site_combo_resid, width = 3.75, height = 3, units = 'in')


aov_s5_gen_epi_lm <- aov(resid ~ as.factor(consensus) * as.factor(weak), 
                         data = s5_gen_epi_lm)

plot(aov_s5_gen_epi_lm, 1)

library(car)

leveneTest(resid ~ as.factor(consensus) * as.factor(weak), 
           data = s5_gen_epi_lm)

#We get no homogeneity of variance

#Reviewer-specific plots barcode effects----------------------------------------

#looking at barcode effects based on % match to perfect variant

SP3_SP5_bc <- read_tsv('../BCMap/uniqueSP2345_bcdetail.txt', 
                        col_names = c(
                          'fluff', 'barcode', 'num_unique', 'total_reads', 
                          'mapped_variant_reads', 'most_common', 'name'), 
                        skip = 1) %>%
  select(-fluff) %>%
  mutate(subpool = ifelse(startsWith(name, 'subpool'), 
                          substr(name, 1, 8), 
                          'control')) %>%
  filter(subpool != 'subpool2' & subpool != 'subpool4')

#join to variants analyzed in episomal MPRA

bc_DNA_percent <- function(df1, df2) {
  bc_count_med <- df1 %>%
    group_by(subpool, name, most_common) %>%
    mutate(barcodes_DNA = n()) %>%
    filter(barcodes_DNA > 7) %>%
    mutate(med_ratio = median(ratio)) %>%
    filter(med_ratio > 0)
  BCmap_join <- left_join(bc_count_med, df2) %>%
    mutate(percent_associated = (mapped_variant_reads/total_reads)*100) %>%
    mutate(other_assoc = total_reads - mapped_variant_reads)
  return(BCmap_join)
} 
  
epi_bc_percent <- bc_DNA_percent(bc_DNA_RNA_22A, SP3_SP5_bc)

bc_assoc_count_0 <- filter(epi_bc_percent, other_assoc == 0)

bc_assoc_count_1 <- filter(epi_bc_percent, other_assoc == 1)

bc_assoc_count_2 <- filter(epi_bc_percent, other_assoc == 2)

bc_assoc_count_3 <- filter(epi_bc_percent, other_assoc == 3)

bc_assoc_count_4 <- filter(epi_bc_percent, other_assoc > 3)

bc_assoc_count <- tibble(reads_to_nonmap = c('0', '1', '2', '3', '>3'),
                         percent = c(nrow(bc_assoc_count_0)/nrow(epi_bc_percent),
                                     nrow(bc_assoc_count_1)/nrow(epi_bc_percent),
                                     nrow(bc_assoc_count_2)/nrow(epi_bc_percent),
                                     nrow(bc_assoc_count_3)/nrow(epi_bc_percent),
                                     nrow(bc_assoc_count_4)/nrow(epi_bc_percent))) %>%
  mutate(percent = percent*100)

p_bc_assoc_count <- bc_assoc_count %>%
  mutate(reads_to_nonmap = factor(reads_to_nonmap, 
                                  levels = c('0', '1', '2', '3', '>3'))) %>%
  ggplot(aes(reads_to_nonmap, percent)) +
  geom_col() +
  ylab('% of barcodes') +
  xlab('Reads mapping to non-variant') +
  figurefont_theme

ggsave('../plots/p_bc_assoc_count.pdf', 
       p_bc_assoc_count, height = 1.5, width = 2, units = 'in')

p_bc_percent_exp <- epi_bc_percent %>%
  ggplot(aes(percent_associated, med_ratio)) +
  geom_point(alpha = 0.1, size = 0.2) +
  scale_y_log10('Variant expression (a.u.)') +
  annotation_logticks(sides = 'l') +
  xlab('% barcode map to variant') +
  figurefont_theme

ggsave('../plots/p_bc_percent_exp.png', 
       p_bc_percent_exp, height = 1.75, width = 2.5, units = 'in')

#reviewer-specific plots, hill plots--------------------------------------------

#subtract expression at 0 conc per each variant to make fitting easier

trans_back_0_norm_conc <- epi_back_norm_pc_spGl4 %>%
  mutate(ave_ratio_2_5_norm = ave_ratio_2_5_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_2_4_norm = ave_ratio_2_4_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_2_3_norm = ave_ratio_2_3_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_2_2_norm = ave_ratio_2_2_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_2_1_norm = ave_ratio_2_1_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_20_norm = ave_ratio_20_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_22_norm = ave_ratio_22_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_0_norm = ave_ratio_0_norm - ave_ratio_0_norm) %>%
  var_conc_exp()

#non-linear fit functions, using the hill equation in the form of reaction 
#velocity, fitting n and "EC50"

library(minpack.lm)

m_m_model_nlslm <- function(df) {
  m_m_nlslm <- nlsLM(
    ave_ratio_norm ~ (max_ave_ratio_norm * conc^n)/(conc_half_max^n + conc^n),
    data = df, 
    start = list(conc_half_max = (2^-3), max_ave_ratio_norm = 2, n = 1))
  return(m_m_nlslm)
}

m_m_nest_coef <- function(df1) {
  add_coef_unnest <- df1 %>%
    mutate(results = map(m_m_fit, tidy)) %>%
    select(-m_m_fit, -data) %>%
    unnest()
}

m_m_nest_pred_resid <- function(df1) {
  pred <- df1 %>%
    mutate(m_m_pred = map2(data, m_m_fit, add_predictions)) %>%
    select(-data, -m_m_fit) %>%
    unnest()
  resid <- df1 %>%
    mutate(m_m_resids = map2(data, m_m_fit, add_residuals)) %>%
    select(-data, -m_m_fit) %>%
    unnest()
  data <- df1 %>%
    select(-m_m_fit) %>%
    unnest()
  data_pred <- left_join(data, pred, 
                         by = c('subpool', 'name', 'most_common', 'background',
                                'conc', 'ave_ratio_norm', 'induction', 
                                'ratio_A_norm', 'ratio_B_norm'))
  data_pred_resid <- left_join(data_pred, resid,
                               by = c('subpool', 'name', 'most_common', 
                                      'background', 'conc', 'ave_ratio_norm', 
                                      'induction', 'ratio_A_norm', 
                                      'ratio_B_norm')) %>%
    ungroup()
  return(data_pred_resid)
}

#filter out variants with higher expression in the two site library (sp3) and
#six site library (sp5) to make more robust model fits

trans_back_0_norm_conc_nest_sp3 <- trans_back_0_norm_conc %>%
  select(-ave_barcode) %>%
  mutate(conc = 2^conc) %>%
  filter(subpool == 'subpool3') %>%
  arrange(desc(ave_ratio_norm)) %>%
  group_by(subpool, name, most_common, background) %>%
  nest() %>%
  slice(1:590)

m_m_nest_fit_nlslm_sp3 <- trans_back_0_norm_conc_nest_sp3 %>%
  mutate(m_m_fit = map(trans_back_0_norm_conc_nest_sp3$data, m_m_model_nlslm))

m_m_coef_nlslm_sp3 <- m_m_nest_coef(m_m_nest_fit_nlslm_sp3)


trans_back_0_norm_conc_nest_sp5 <- trans_back_0_norm_conc %>%
  select(-ave_barcode) %>%
  mutate(conc = 2^conc) %>%
  filter(!grepl('^subpool5_no_site_no_site_no_site_no_site_no_site_no_site', 
                name)) %>%
  filter(subpool == 'subpool5') %>%
  arrange(desc(ave_ratio_norm)) %>%
  group_by(subpool, name, most_common, background) %>%
  nest() %>%
  slice(1:1270)

m_m_nest_fit_nlslm_sp5 <- trans_back_0_norm_conc_nest_sp5 %>%
  mutate(m_m_fit = map(trans_back_0_norm_conc_nest_sp5$data, m_m_model_nlslm))

m_m_coef_nlslm_sp5 <- m_m_nest_coef(m_m_nest_fit_nlslm_sp5)

#Filter out n and EC50 estimates with lower than 0.25 relative standard error.
#Join back to main df to retaion sequence info

hillcoef_sp3 <- m_m_coef_nlslm_sp3 %>%
  filter(term == 'n') %>%
  mutate(n_rse = std.error/estimate) %>%
  filter(n_rse <= 0.25) %>%
  rename(n = estimate) %>%
  select(-term, -std.error, -statistic, -p.value) %>%
  left_join(s3_epi_back_norm_conc,
            by = c('most_common', 'background')) %>%
  select(-subpool) %>%
  ungroup()

hill_EC50_sp3 <- m_m_coef_nlslm_sp3 %>%
  filter(term == 'conc_half_max') %>%
  mutate(EC50_rse = std.error/estimate) %>%
  filter(EC50_rse <= 0.25) %>%
  rename(EC50 = estimate) %>%
  select(-term, -std.error, -statistic, -p.value) %>%
  inner_join(hillcoef_sp3, by = c('name', 'most_common', 'background')) %>%
  ungroup() %>%
  mutate(conc = 2^conc)


hillcoef_sp5 <- m_m_coef_nlslm_sp5 %>%
  filter(term == 'n') %>%
  mutate(n_rse = std.error/estimate) %>%
  filter(n_rse <= 0.25) %>%
  rename(n = estimate) %>%
  select(-term, -std.error, -statistic, -p.value) %>%
  left_join(s5_epi_back_norm_conc,
            by = c('most_common', 'background')) %>%
  select(-subpool) %>%
  ungroup()

hill_EC50_sp5 <- m_m_coef_nlslm_sp5 %>%
  filter(term == 'conc_half_max') %>%
  mutate(EC50_rse = std.error/estimate) %>%
  filter(EC50_rse <= 0.25) %>%
  rename(EC50 = estimate) %>%
  select(-term, -std.error, -statistic, -p.value) %>%
  inner_join(hillcoef_sp5, by = c('name', 'most_common', 'background')) %>%
  ungroup() %>%
  mutate(conc = 2^conc)

#determine model predictions and residuals

m_m_p_r_sp3 <- m_m_nest_pred_resid(m_m_nest_fit_nlslm_sp3) %>%
  rename(ave_ratio_norm_0 = ave_ratio_norm)

m_m_EC50_n_p_r_sp3 <- left_join(hill_EC50_sp3, m_m_p_r_sp3, 
                            by = c('subpool', 'name', 'most_common', 
                                   'background', 'conc', 'induction', 
                                   'ratio_A_norm', 'ratio_B_norm'))


m_m_p_r_sp5 <- m_m_nest_pred_resid(m_m_nest_fit_nlslm_sp5) %>%
  rename(ave_ratio_norm_0 = ave_ratio_norm)

m_m_EC50_n_p_r_sp5 <- left_join(hill_EC50_sp5, m_m_p_r_sp5, 
                                by = c('subpool', 'name', 'most_common', 
                                       'background', 'conc', 'induction', 
                                       'ratio_A_norm', 'ratio_B_norm'))

#plot residuals and relative standare errors

p_resid_dens_sp3 <- m_m_EC50_n_p_r_sp3 %>%
  ggplot(aes(resid)) +
  geom_density(kernel = 'gaussian') +
  xlab('expression residual')

p_resid_distr_sp3 <- ggplot(m_m_EC50_n_p_r_sp3, aes(ave_ratio_norm_0, resid)) +
  geom_point(alpha = 0.1, size = 0.75) +
  xlab('Average normalized\nexpression (a.u.)\n-exp. at 0 µM')

p_EC50_rse_hist_sp3 <- m_m_EC50_n_p_r_sp3 %>%
  filter(conc == 4) %>%
  ggplot(aes(EC50_rse)) +
  geom_density(kernel = 'gaussian') +
  scale_x_continuous(breaks = c(0, 0.1, 0.2))

p_n_rse_hist_sp3 <- m_m_EC50_n_p_r_sp3 %>%
  filter(conc == 4) %>%
  ggplot(aes(n_rse)) +
  geom_density(kernel = 'gaussian') +
  scale_x_continuous(breaks = c(0, 0.1, 0.2))


p_resid_dens_sp5 <- m_m_EC50_n_p_r_sp5 %>%
  ggplot(aes(resid)) +
  geom_density(kernel = 'gaussian') +
  xlab('expression residual')

p_resid_distr_sp5 <- ggplot(m_m_EC50_n_p_r_sp5, aes(ave_ratio_norm_0, resid)) +
  geom_point(alpha = 0.1, size = 0.75) +
  xlab('Average normalized\nexpression (a.u.)\n-exp. at 0 µM')

p_EC50_rse_hist_sp5 <- m_m_EC50_n_p_r_sp5 %>%
  filter(conc == 4) %>%
  ggplot(aes(EC50_rse)) +
  geom_density(kernel = 'gaussian') +
  scale_x_continuous(breaks = c(0, 0.1, 0.2))

p_n_rse_hist_sp5 <- m_m_EC50_n_p_r_sp5 %>%
  filter(conc == 4) %>%
  ggplot(aes(n_rse)) +
  geom_density(kernel = 'gaussian') +
  scale_x_continuous(breaks = c(0, 0.1, 0.2))

#Plot trends in n and EC50 based on library features

p_s3_spacing_n <- m_m_EC50_n_p_r_sp3 %>%
  filter(spacing != 0 & spacing != 70) %>%
  filter(conc == 4) %>%
  ggplot(aes(spacing, n)) +
  facet_grid(. ~ background) +
  geom_jitter(aes(group = spacing), alpha = 0.25, size = 0.5) +
  geom_boxplot(aes(group = spacing), alpha = 0.5, outlier.alpha = 0, 
               size = 0.3) +
  panel_border(colour = 'black') +
  ylab('Hill n') +
  xlab('Spacing') +
  figurefont_theme +
  theme(strip.background = element_rect(colour="black", fill="white"))

p_s3_spacing_ec50 <- m_m_EC50_n_p_r_sp3 %>%
  filter(spacing != 0 & spacing != 70) %>%
  filter(conc == 4) %>%
  ggplot(aes(spacing, EC50)) +
  facet_grid(. ~ background) +
  geom_jitter(aes(group = spacing), alpha = 0.25, size = 0.5) +
  geom_boxplot(aes(group = spacing), alpha = 0.5, outlier.alpha = 0,
               size = 0.3) +
  panel_border(colour = 'black') +
  xlab('Spacing') +
  figurefont_theme +
  theme(strip.background = element_rect(colour="black", fill="white"))


p_s5_cons_n <- m_m_EC50_n_p_r_sp5 %>%
  filter(conc == 4) %>%
  ggplot(aes(as.factor(consensus), n, fill = as.factor(weak))) +
  facet_grid(. ~ background) +
  geom_boxplot(outlier.size = 0.7, size = 0.3, 
               outlier.shape = 21, outlier.alpha = 1, 
               position = position_dodge(0.75)) +
  scale_fill_manual(name = 'number of\nweak CREs', 
                    values = cbPalette7_grad_light) +
  panel_border(colour = 'black') +
  ylab('Hill n') +
  xlab('Number of consensus CREs') +
  figurefont_theme +
  theme(strip.background = element_rect(colour="black", fill="white"))

p_s5_cons_ec50 <- m_m_EC50_n_p_r_sp5 %>%
  filter(conc == 4) %>%
  ggplot(aes(as.factor(consensus), EC50, fill = as.factor(weak))) +
  facet_grid(. ~ background) +
  geom_boxplot(outlier.size = 0.7, size = 0.3, 
               outlier.shape = 21, outlier.alpha = 1, 
               position = position_dodge(0.75)) +
  scale_fill_manual(name = 'number of\nweak CREs', 
                    values = cbPalette7_grad_light) +
  panel_border(colour = 'black') +
  xlab('Number of consensus CREs') +
  figurefont_theme +
  theme(strip.background = element_rect(colour="black", fill="white"))

ggsave('../plots/p_s3_spacing_n.pdf', p_s3_spacing_n, 
       width = 2.75, height = 1.5, units = 'in')

ggsave('../plots/p_s3_spacing_ec50.pdf', p_s3_spacing_ec50, 
       width = 2.75, height = 1.5, units = 'in')

ggsave('../plots/p_s5_cons_n.pdf', p_s5_cons_n, 
       width = 6.5, height = 1.75, units = 'in')

ggsave('../plots/p_s5_cons_ec50.pdf', p_s5_cons_ec50, 
       width = 6.5, height = 1.75, units = 'in')


#Reviewer-specific plots feed forward attempt-----------------------------------

#using caret

library(fastDummies)

test_data <- s5_gen_epi %>%
  filter(MPRA == 'episomal') %>%
  select(site1, site2, site3, site4, site5, site6, background, consensus, weak,
         ave_ratio_22)

test_data_dummy <- dummy_cols(test_data) %>%
  select(-site1, -site2, -site3, -site4, -site5, -site6, -background)

control <- trainControl(method = 'repeatedcv', number = 10, repeats = 3)

model <- train(ave_ratio_22 ~ ., test_data_dummy, method = 'lm')

importance <- varImp(model, scale = FALSE)

#caret doesn't play nicely with categorical variables.....


#Extras-------------------------------------------------------------------------

#Determining total unique sequences analyzed in paper

#Single CRE assay only looked at half of the designed library (with sequence 
#flanks)

s2_flanks <- s2_epi_rep_0_25 %>%
  filter(site == 'consensusflank') %>%
  filter(conc == 0)

most_common <- function(df) {
  df <- df %>%
    select(most_common)
}

count_seq <- function(dfepi, dfgen, df25, dfnew) {
  dfepi <- most_common(dfepi)
  dfgen <- most_common(dfgen)
  df25 <- most_common(df25)
  dfnew <- most_common(dfnew)
  total <- rbind(dfepi, dfgen, df25, dfnew)
  return(total)
  }

test <- count_seq(med_rep_0_22_A_B, gen_rep_1_2, s2_flanks, med_rep_followup)

#Remove 20 from test for controls not analyzed in paper


