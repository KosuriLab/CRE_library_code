#Written for the analysis of a range of inductions in the episomal library
#DNA_BC: DNA
#R0A_BC and R0B_BC: RNA at 0 µM Forsk replicate A and B
#R2#A_BC and R2#B_BC: RNA at 2^# µM Forsk replicate A and B, includes negative 
#and positive #'s

#Tested concentrations: 0, 2^-5, 2^-4, 2^-3, 2^-2, 2^-1, 2^0, 2^2 µM Forsk

#All figures made with this dataset require sections "Load index and bcmap 
#files" to "Median BC expression".

#Some figures require section "Background-normalize expression" and it is
#recommended to run this section before plotting any figure

#Most figures require both genomic and episomal MPRAs, the genomic MPRA data
#processing is performed in 20171129_genlib.R and imported here in section
#"Import genomic MPRA and combine with episomal". Other episomal MPRAs used in
#supplemental figures are perfomred in 20170320_epilib_analysis_0_25 and 
#20170631_epilib_analysis_0_64 and are imported per supplemental figure section.

#In all but figure 1, the section "Separate into sublibraries" is required for
#figure generation.

#Establish workspace------------------------------------------------------------

options(stringsAsFactors = F)

#import necessary libraries

library(tidyverse)
library(lemon)
library(viridis)
library(cowplot)
library(caTools)
library(broom)
library(modelr)

#install packages again below this
library(reshape2)
library(stringr)
library(ggExtra)
library(lazyeval)
library(splines)
library(GGally)
library(updateR)
library(ggsignif)

#General figure customizations

cbPalette7 <- c('#440154FF', '#39568CFF', '#287D8EFF', '#20A387FF', '#73D055FF',
                '#B8DE29FF', '#FDE725FF')

cbPalette7_grad_light <- c('white', '#FDE725FF', '#B8DE29FF', '#55C667FF', 
                           '#1F968BFF', '#39568CFF', '#482677FF')

spacing_5_20_cbpalette <- c('gray20', 'deepskyblue2', 'orangered3', 'sandybrown')

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

#Load index and bcmap files-----------------------------------------------------

bc_DNA <- read_tsv('BCreads_txts/DNA_BC.txt')
bc_R0A <- read_tsv('BCreads_txts/R0A_BC.txt')
bc_R0B <- read_tsv('BCreads_txts/R0B_BC.txt')
bc_R2_5A <- read_tsv('BCreads_txts/R2-5A_BC.txt')
bc_R2_5B <- read_tsv('BCreads_txts/R2-5B_BC.txt')
bc_R2_4A <- read_tsv('BCreads_txts/R2-4A_BC.txt')
bc_R2_4B <- read_tsv('BCreads_txts/R2-4B_BC.txt')
bc_R2_3A <- read_tsv('BCreads_txts/R2-3A_BC.txt')
bc_R2_3B <- read_tsv('BCreads_txts/R2-3B_BC.txt')
bc_R2_2A <- read_tsv('BCreads_txts/R2-2A_BC.txt')
bc_R2_2B <- read_tsv('BCreads_txts/R2-2B_BC.txt')
bc_R2_1A <- read_tsv('BCreads_txts/R2-1A_BC.txt')
bc_R2_1B <- read_tsv('BCreads_txts/R2-1B_BC.txt')
bc_R20A <- read_tsv('BCreads_txts/R20A_BC.txt')
bc_R20B <- read_tsv('BCreads_txts/R20B_BC.txt')
bc_R22A <- read_tsv('BCreads_txts/R22A_BC.txt')
bc_R22B <- read_tsv('BCreads_txts/R22B_BC.txt')

#Load barcode mapping table, sequences (most_common) are rcomp due to sequencing
#format. Pick out controls, SP3, and SP5 in the bcmap that were used in this 
#assay

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

#Join BC reads to BC mapping, keeping the reads only appearing in barcode 
#mapping and replacing na with 0 reads. Determine normalized reads per million

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
#of RNA/DNA normalizedreads per million

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
#determine the median absolute deviation of all barcode ratios. 

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
    summarize(med_ratio = median(ratio))
  mad_ratio <- bc_min_8_df %>%
    group_by(subpool, name, most_common) %>%
    summarize(mad = mad(ratio))
  med_mad <- inner_join(med_ratio, mad_ratio, 
                        by = c('subpool', 'name', 'most_common')) %>%
    mutate(mad_over_med = as.double(mad/med_ratio)) %>%
    mutate(mad_over_med = if_else(
      is.na(mad_over_med),
      as.double(0), 
      mad_over_med))
  bc_med <- inner_join(med_mad, bc_DNA_RNA, 
                       by = c('subpool', 'name', 'most_common')) %>%
    ungroup() %>%
    filter(med_ratio > 0)
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

#Combine biological replicates

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

#Make untidy df with conc as a variable and normalized expressions as single 
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

#Import genomic MPRA and combine with episomal----------------------------------

#Genomic MPRA expression determination and data processing performed in
#"20171129_genlib_analysis/20171129_genlib.R". Processed df with expression
#values is exported from there and imported here to plot together

gen_rep_1_2 <- read_tsv('../20171129_genlib_analysis/rep_1_2.txt')

#Combine genomic and episomal dfs, only comparing to expression at 2^2 µM 
#forskolin in the episomal dataset. Here med_ratio_br# and ave_med_ratio refers 
#to genomic expression, either across biological replicates or averaged.
#Episomal expression is represented by the annotation <name>_22 as representing
#the forskolin concentration the sample was incubated with

gen_epi <- epi_back_norm_pc_spGl4 %>%
  mutate(ave_ratio_22 = (med_ratio_22A + med_ratio_22B)/2) %>%
  select(subpool, name, most_common, barcodes_DNA, med_ratio_22A,
         barcodes_RNA_22A, ratio_22A_norm, med_ratio_22B, barcodes_RNA_22B, 
         ratio_22B_norm, ave_ratio_22, ave_ratio_22_norm) %>%
  inner_join(gen_rep_1_2, by = c('subpool', 'name', 'most_common'))

#Make untidy df where MPRA format is a variable according to average variant
#expression and average #barcodes per variant in the RNA samples, also create 
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

#Figure 1D----------------------------------------------------------------------

#Plot subfigure C, normalized variant expression curves across forskolin
#concentrations. Expression curves for backgrounds and control overlayed

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

#Figure 1E----------------------------------------------------------------------

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
                           name == 'pGL4.29 Promega 1-63 + 1-87'), 
             fill = 'red', shape = 21, size = 1.75) +
  annotation_logticks(scaled = TRUE) +
  xlab("Expression (a.u.) replicate 1") +
  ylab("Expression (a.u.) replicate 2") +
  scale_x_log10(limits = c(0.03, 20), breaks = c(0.1, 1, 10)) + 
  scale_y_log10(limits = c(0.03, 20), breaks = c(0.1, 1, 10)) +
  theme(strip.background = element_rect(colour="black", fill="white"),
        axis.line.y = element_line(), panel.spacing.x=unit(1, "lines")) +
  figurefont_theme

p_fig1_int_med_rep <- gen_epi %>%
  ggplot(aes(med_ratio_br1, med_ratio_br2)) +
  geom_point(alpha = 0.1, size = 1) +
  geom_point(data = filter(gen_epi, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             fill = 'orange', shape = 21, size = 1.75) + 
  geom_point(data = filter(gen_epi, name == 'pGL4.29 Promega 1-63 + 1-87'), 
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

ggsave('../plots/p_fig1_int_med_rep.png', p_fig1_int_med_rep, 
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
                              'bottom',
                              'top')) %>%
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

titration_luc <- read_csv("plate_reader/170726_epi_gen_R.csv") %>%
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

#Replicability beyond 4 µM and Supp. Fig. 2-------------------------------------

#Import expression from MPRA done across 0-64 µM forskolin concentrations,
#concentrations given in powers of 2

epi_rep_0_64 <- read_tsv('../20170631_epilib_analysis_0_64/rep_0_64_med.txt')
epi_rep_0_64_back <- read_tsv(
  '../20170631_epilib_analysis_0_64/rep_0_64_med_back.txt'
  )
epi_rep_0_64_conc <- read_tsv(
  '../20170631_epilib_analysis_0_64/rep_0_64_med_norm.txt'
  ) %>%
  mutate(conc = log2(conc))

#plot background-normalized replicability similar to Supplemental figure 1E

p_var_epi_all_0_64 <- ggplot(epi_rep_0_64_conc, 
                        aes(ratio_A_norm, ratio_B_norm)) +
  facet_rep_wrap(~ conc, nrow = 3) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_point(data = filter(epi_rep_0_64_conc, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             fill = 'orange', shape = 21, size = 1.25) + 
  geom_point(data = filter(epi_rep_0_64_conc, name == 'pGL4.29 Promega 1-63 + 1-87_back_55'), 
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

ggsave('../plots/p_var_epi_all_0_64.png', p_var_epi_all_0_64, 
       width = 4, height = 4.3, units = 'in')

log10_epi_rep_0_64_back <- var_log10(epi_rep_0_64_back)

pearsons_epi_conc_0_64 <- tibble(
  conc = c(0, 1, 2, 4, 8, 16, 25, 32, 64),
  pearsons = c(round(cor(log10_epi_rep_0_64_back$ratio_0A_norm, 
                         log10_epi_rep_0_64_back$ratio_0B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_epi_rep_0_64_back$ratio_1A_norm, 
                         log10_epi_rep_0_64_back$ratio_1B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_epi_rep_0_64_back$ratio_2A_norm, 
                         log10_epi_rep_0_64_back$ratio_2B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_epi_rep_0_64_back$ratio_4A_norm, 
                         log10_epi_rep_0_64_back$ratio_4B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_epi_rep_0_64_back$ratio_8A_norm, 
                         log10_epi_rep_0_64_back$ratio_8B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_epi_rep_0_64_back$ratio_16A_norm, 
                         log10_epi_rep_0_64_back$ratio_16B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_epi_rep_0_64_back$ratio_25A_norm, 
                         log10_epi_rep_0_64_back$ratio_25B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_epi_rep_0_64_back$ratio_32A_norm, 
                         log10_epi_rep_0_64_back$ratio_32B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_epi_rep_0_64_back$ratio_64A_norm, 
                         log10_epi_rep_0_64_back$ratio_64B_norm, 
                         use = "pairwise.complete.obs", method = "pearson"), 3))
)

#Plot titration similar to Figure 1D

p_titr_pc_back_0_64 <- epi_rep_0_64_conc %>%
  ggplot(aes(conc, ave_ratio_norm)) +
  geom_line(aes(group = name), alpha = 0.1) +
  geom_point(data = filter(epi_rep_0_64_conc, 
                           startsWith(name, 
                                      'subpool5_no_site_no_site_no_site_no_site_no_site_no_site')),
             color = 'darkgoldenrod1', shape = 19, stroke = 0.75) +
  geom_line(data = filter(epi_rep_0_64_conc, 
                          startsWith(name, 
                                     'subpool5_no_site_no_site_no_site_no_site_no_site_no_site')),
            color = 'darkgoldenrod1', size = 1) +
  geom_point(data = filter(epi_rep_0_64_conc, 
                           startsWith(name, 
                                      'pGL4.29 Promega 1-63 + 1-87')),
             color = 'firebrick2', shape = 19, stroke = 0.75) +
  geom_line(data = filter(epi_rep_0_64_conc, 
                          startsWith(name, 
                                     'pGL4.29 Promega 1-63 + 1-87')),
            color = 'firebrick2', size = 1) +
  ylab('Average normalized\nexpression (a.u.)') +
  annotation_logticks(sides = 'b', short = unit(0.05, 'cm'), 
                      mid = unit(0.1, 'cm'), long = unit(0.15, 'cm')) +
  scale_x_continuous(breaks = (-2:6), 'log2 forskolin (µM)') +
  figurefont_theme

ggsave('../plots/p_titr_pc_back_0_64.pdf', p_titr_pc_back_0_64, 
       width = 2.8, height = 2, units = 'in')

#Compare average expression between episomal MPRAs between repeated 
#concentrations 0, 1 and 4 µM forskolin

join_epi <- function(df_0_22, df_0_64) {
  df_0_22 <- df_0_22 %>%
    mutate(ave_med_ratio_0 = (med_ratio_0A + med_ratio_0B)/2) %>%
    mutate(ave_med_ratio_1 = (med_ratio_20A + med_ratio_20B)/2) %>%
    mutate(ave_med_ratio_4 = (med_ratio_22A + med_ratio_22B)/2) %>%
    select(subpool, name, most_common, ave_med_ratio_0, ave_med_ratio_1, 
           ave_med_ratio_4)
  df_0_64 <- df_0_64 %>%
    mutate(ave_med_ratio_0 = (med_ratio_0A + med_ratio_0B)/2) %>%
    mutate(ave_med_ratio_1 = (med_ratio_1A + med_ratio_1B)/2) %>%
    mutate(ave_med_ratio_4 = (med_ratio_4A + med_ratio_4B)/2) %>%
    select(subpool, name, most_common, ave_med_ratio_0, ave_med_ratio_1, 
           ave_med_ratio_4)
  join <- inner_join(df_0_22, df_0_64, by = c('subpool', 'name', 'most_common'),
                     suffix = c('_0_22', '_0_64'))
  return(join)
}

epi_epi <- join_epi(med_rep_0_22_A_B, epi_rep_0_64)

p_epi_epi_0 <- ggplot(epi_epi, 
                      aes(ave_med_ratio_0_0_22, ave_med_ratio_0_0_64)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_point(data = filter(epi_epi, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             fill = 'orange', shape = 21, size = 1.25) + 
  geom_point(data = filter(epi_epi, name == 'pGL4.29 Promega 1-63 + 1-87'), 
             fill = 'red', shape = 21, size = 1.25) +
  annotation_logticks(scaled = TRUE) +
  xlab("Average expression (a.u.)\n0-4 µM forskolin MPRA") +
  ylab("Average expression (a.u.)\n0-64 µM forskolin MPRA") +
  scale_x_log10(limits = c(0.05, 20), breaks = c(0.1, 1, 10)) + 
  scale_y_log10(limits = c(0.05, 20), breaks = c(0.1, 1, 10)) +
  background_grid(major = 'xy', minor = 'none') + 
  panel_border(colour = 'black') +
  theme(strip.background = element_rect(colour="black", fill="white"),
        axis.line.y = element_line(), panel.spacing.x=unit(1, "lines")) +
  figurefont_theme

p_epi_epi_1 <- ggplot(epi_epi, 
                      aes(ave_med_ratio_1_0_22, ave_med_ratio_1_0_64)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_point(data = filter(epi_epi, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             fill = 'orange', shape = 21, size = 1.25) + 
  geom_point(data = filter(epi_epi, name == 'pGL4.29 Promega 1-63 + 1-87'), 
             fill = 'red', shape = 21, size = 1.25) +
  annotation_logticks(scaled = TRUE) +
  xlab("Average expression (a.u.)\n0-4 µM forskolin MPRA") +
  ylab("Average expression (a.u.)\n0-64 µM forskolin MPRA") +
  scale_x_log10(limits = c(0.05, 20), breaks = c(0.1, 1, 10)) + 
  scale_y_log10(limits = c(0.05, 20), breaks = c(0.1, 1, 10)) +
  background_grid(major = 'xy', minor = 'none') + 
  panel_border(colour = 'black') +
  theme(strip.background = element_rect(colour="black", fill="white"),
        axis.line.y = element_line(), panel.spacing.x=unit(1, "lines")) +
  figurefont_theme

p_epi_epi_4 <- ggplot(epi_epi, 
                      aes(ave_med_ratio_4_0_22, ave_med_ratio_4_0_64)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_point(data = filter(epi_epi, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             fill = 'orange', shape = 21, size = 1.25) + 
  geom_point(data = filter(epi_epi, name == 'pGL4.29 Promega 1-63 + 1-87'), 
             fill = 'red', shape = 21, size = 1.25) +
  annotation_logticks(scaled = TRUE) +
  xlab("Average expression (a.u.)\n0-4 µM forskolin MPRA") +
  ylab("Average expression (a.u.)\n0-64 µM forskolin MPRA") +
  scale_x_log10(limits = c(0.05, 20), breaks = c(0.1, 1, 10)) + 
  scale_y_log10(limits = c(0.05, 20), breaks = c(0.1, 1, 10)) +
  background_grid(major = 'xy', minor = 'none') + 
  panel_border(colour = 'black') +
  theme(strip.background = element_rect(colour="black", fill="white"),
        axis.line.y = element_line(), panel.spacing.x=unit(1, "lines")) +
  figurefont_theme

ggsave('../plots/p_epi_epi_0.png', p_epi_epi_0, width = 1.75, height = 1.6, 
       units = 'in')

ggsave('../plots/p_epi_epi_1.png', p_epi_epi_1, width = 1.75, height = 1.6, 
       units = 'in')

ggsave('../plots/p_epi_epi_4.png', p_epi_epi_4, width = 1.75, height = 1.6, 
       units = 'in')

log10_epi_epi <- var_log10(epi_epi)

pearsons_epi_epi <- tibble(
  conc = c(0, 1, 4),
  pearsons = c(round(cor(log10_epi_epi$ave_med_ratio_0_0_22, 
                         log10_epi_epi$ave_med_ratio_0_0_64, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_epi_epi$ave_med_ratio_1_0_22, 
                         log10_epi_epi$ave_med_ratio_1_0_64, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(log10_epi_epi$ave_med_ratio_4_0_22, 
                         log10_epi_epi$ave_med_ratio_4_0_64,
                         use = "pairwise.complete.obs", method = "pearson"), 3))
)


#Separate into sublibraries-----------------------------------------------------

#Subpool 3 corresponds to the CRE Spacing and Distance Library. This library
#contains 2 consensus CRE sites with flanks (ATTGACGTCAGC) that vary in 
#CRE Spacing from one another by 0 (no inner flanks), 5, 10, 15, 20, and 
#70 bp (all but 0 appear as spacing - 4 bp due to 4 bp of sequence flanking both
#CREs). Both CREs are then placed at variable distances to the 3' end of the 
#backgrounds. Separation lists the CRE spacing between sites and CRE distance. 
#Distances measured from the end of the background to the CRE proximal to the 
#promoter. Added 2 bp to all distances here to measure to start of CRE without 
#the flanks and then added 64 bp to measure to the minimal promoter. Added 4 to 
#all spacings but 0 to measure difference between start of sites without flanks.

subpool3 <- function(df) {
  df <- df %>%
    filter(subpool == "subpool3") %>%
    ungroup() %>%
    select(-subpool) %>%
    mutate(name = gsub('2BS ', '', name), 
           name = gsub(' bp spacing ', '_', name)) %>%
    separate(name, 
             into = c("subpool", "spacing", "fluff1", "fluff2", "dist", "fluff3", "fluff4"),
             sep = "_", convert = TRUE) %>%
    select(-subpool, -fluff1, -fluff2, -fluff3, -fluff4) %>%
    mutate(dist = as.integer(dist + 2 + 64)) %>%
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

#3bp moving window function used in a lot of the CRE Spacing and Distance
#library figures

moveavg_dist3 <- function(df) {
  df <- df %>%
    mutate(ave_3 = runmean(ave_ratio_norm, 3, alg = 'R', endrule = 'NA'))
}

#Subpool 5 corresponds to the CRE Number and Affinity Library. This library
#contains 6 equally spaced sites with 17 bp CRE Spacing. Per variant, each sites
#is one of: the consensus CRE, the weak CRE, or no CRE. Both the weak and 
#consensus CREs have the same flanking sequence. Here site 1, 2, 3, 4, 5, and 6 
#equate to -191, -166, -141, -116, -91 and -66 site distances to the downstream
#promoter. Separation lists identity of CRE per site, the total CRE sites per
#variants (consensus + weak), and if variants contain only consenus CREs, 
#mixed consensus and weak, or no CREs.

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
  scale_x_continuous("Distance to minimal promoter (bp)",
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
  scale_x_continuous("Distance to minimal promoter (bp)",
                     breaks = c(seq(from = 60, to = 190, by = 10))) +
  scale_y_continuous("Average normalized\nexpression (a.u.)") +
  background_grid(major = 'x', minor = 'none', colour.major = 'grey70') + 
  panel_border(colour = 'black') +
  figurefont_theme  +
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
  range_change <- left_join(med_low, med_high, 
                            by = c('spacing', 'background', 'MPRA')) %>%
    mutate(range_change = median_6796/median_147176) %>%
    filter(spacing != 0, spacing != 70) %>%
    group_by(background, MPRA) %>%
    summarize(med_fold_change = median(range_change))
  return(range_change)
}
  
s3_gen_epi_med_range <- med_range_dist(s3_gen_epi)

s3_gen_epi_med_change_dist <- fold_change_range(s3_gen_epi_med_range) %>%
  filter(background != 41)
  
p_s3_gen_epi_med_change_dist_55 <- s3_gen_epi_med_range %>%
  filter(spacing != 0 & spacing != 70 & background == '55') %>%
  mutate(range = factor(range, levels = c('67-96', '147-176'))) %>%
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
  mutate(range = factor(range, levels = c('67-96', '147-176'))) %>%
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

ggsave('../plots/p_s3_gen_epi_med_change_dist_55.pdf', 
       p_s3_gen_epi_med_change_dist_55, width = 2.7, height = 2.3, unit = 'in')

ggsave('../plots/p_s3_gen_epi_med_change_dist_52.pdf', 
       p_s3_gen_epi_med_change_dist_52, width = 2.7, height = 2.3, unit = 'in')


#Supplemental Figure 3----------------------------------------------------------

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
  scale_y_log10(limits = c(0.5, 4)) +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 66, to = 206, by = 20)) +
  theme(legend.position = 'right',
        strip.background = element_rect(colour="black", fill="white")) +
  panel_border(colour = 'black') + 
  ylab('Average normalized expression (a.u.)') +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'x', colour.major = 'grey90',
                  colour.minor = 'grey95') +
  figurefont_theme

ggsave('../plots/subpool2_dist_0_25.pdf', p_subpool2_dist_0_25, units = 'in',
       width = 5, height = 2.66)

#Supp. Fig, 4, Figure 3 and Supp. Fig. 5----------------------------------------

#Plot average background-normalized expression per MPRA, CRE Spacing, and CRE
#background acros CRE Distances. Because so few variants were retained in this 
#library with background 41 in the genomic MPRA, the join results in too few
#variants to analyze in the episomal MPRA as well. Thus here I show the variants
#retained in the genomic MPRA for completeness but rbind variants retained in 
#the episomal MPRA to show the full distance profile in the supplement

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

#fit 3 bp moving window function

s3_tidy_moveavg3_MPRA <- s3_gen_epi_rbind41 %>%
  select(background, spacing, dist, MPRA, ave_ratio_norm) %>%
  group_by(background, spacing, MPRA) %>%
  arrange(dist, .by_group = TRUE) %>%
  moveavg_dist3()

#plot supplemental figure 4

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
  scale_y_log10() +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'x', colour.major = 'grey90',
                  colour.minor = 'grey95') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 70, to = 190, by = 20)) +
  figurefont_theme +
  theme(legend.position = 'top',
        strip.background = element_rect(colour="black", fill="white"))

ggsave('../plots/space_dist_gen_epi.pdf', p_space_dist_gen_epi, units = 'in',
       width = 6.4, height = 5)

#Plot figure 3

p_subpool3_spa_4_vchr9_5_10 <- s3_tidy_moveavg3_MPRA %>%
  filter(MPRA == 'episomal' & background == 41 & dist < 124 & (spacing == 5 | spacing == 10)) %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c('gray20', 'dodgerblue2'), name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  geom_vline(xintercept = c(78, 88), color = 'gray20', linetype = 2, 
             alpha = 0.5) +
  geom_vline(xintercept = c(83, 92), color = 'dodgerblue2', linetype = 2, 
             alpha = 0.5) +
  scale_y_log10(limits = c(1, 20), breaks = c(1, 10)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 66, to = 126, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

p_subpool3_spa_4_vchr9_5_15 <- s3_tidy_moveavg3_MPRA %>%
  filter(MPRA == 'episomal' & background == 41 & dist < 124 & (spacing == 5 | spacing == 15)) %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c('gray20', 'orangered3'), 
                     name = 'spacing (bp)') +
  scale_fill_manual(values = c('gray20', 'orangered3'), 
                    name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  geom_vline(xintercept = c(78, 88), color = 'gray20', linetype = 2, 
             alpha = 0.5) +
  scale_y_log10(limits = c(1, 20), breaks = c(1, 10)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
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
  geom_line(aes(y = ave_3), size = 0.4) +
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
  scale_y_log10(limits = c(1, 20), breaks = c(1, 10)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 66, to = 126, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

ggsave('../plots/p_subpool3_med_spa_4_vchr9_5_10.pdf', 
       p_subpool3_spa_4_vchr9_5_10, height = 1.5, width = 4.9, units = 'in')

ggsave('../plots/p_subpool3_med_spa_4_vchr9_5_15.pdf', 
       p_subpool3_spa_4_vchr9_5_15, height = 1.5, width = 4.9, units = 'in')

ggsave('../plots/p_subpool3_med_spa_4_vchr9_10_20.pdf', 
       p_subpool3_spa_4_vchr9_10_20, height = 1.5, width = 4.9, units = 'in')

#Supplemental Figure 5

#Overlay of CRE expression profiles following CRE Distance and between the CRE 
#Spacings 5 and 10, 5 and 15, and 10 and 20 bp. Plots similar to figure 3A are
#shown with back 52 and back 55 and across both MPRAs.

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
  scale_y_log10(limits = c(1, 6), breaks = c(1,5)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
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
  scale_y_log10(limits = c(1, 6), breaks = c(1,5)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

p_subpool3_spa_spgl4_trans_10_20 <- s3_tidy_moveavg3_MPRA %>%
  filter(background == 55 & (spacing == 10 | spacing == 20) & MPRA == 'episomal') %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_vline(xintercept = c(99, 109), color = 'dodgerblue2', linetype = 2, 
             alpha = 0.5) +
  scale_color_manual(values = c('dodgerblue2', 'sandybrown'),
                     name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_log10(limits = c(1, 6), breaks = c(1,5)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10),
                     limits = c(67, 191)) +
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
  scale_y_log10(limits = c(0.5, 20), breaks = c(1,10)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
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
  scale_y_log10(limits = c(0.5, 20), breaks = c(1,10)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

p_subpool3_spa_spgl4_int_10_20 <- s3_tidy_moveavg3_MPRA %>%
  filter(background == 55 & (spacing == 10 | spacing == 20) & MPRA == 'genomic') %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c('dodgerblue2', 'sandybrown'),
                     name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_log10(limits = c(0.5, 20), breaks = c(1,10)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10),
                     limits = c(67, 191)) +
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
  scale_y_log10(limits = c(2, 30), breaks = c(10)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
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
  scale_y_log10(limits = c(2, 30), breaks = c(10)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

p_subpool3_spa_vchr5_trans_10_20 <- s3_tidy_moveavg3_MPRA %>%
  filter(background == 52 & (spacing == 10 | spacing == 20) & MPRA == 'episomal') %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c('dodgerblue2', 'sandybrown'),
                     name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_log10(limits = c(2, 30), breaks = c(10)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10),
                     limits = c(67, 191)) +
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
  scale_y_log10(limits = c(1, 40), breaks = c(1, 10)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
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
  scale_y_log10(limits = c(1, 40), breaks = c(1, 10)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10)) +
  theme(legend.position = 'right', axis.ticks.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  figurefont_theme

p_subpool3_spa_vchr5_int_10_20 <- s3_tidy_moveavg3_MPRA %>%
  filter(background == 52 & (spacing == 10 | spacing == 20) & MPRA == 'genomic') %>%
  ggplot(aes(x = dist, y = ave_ratio_norm, color = as.factor(spacing))) +
  geom_line(aes(y = ave_3), size = 0.4) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c('dodgerblue2', 'sandybrown'),
                     name = 'spacing (bp)') +
  ylab('Average normalized expression (a.u.)') + 
  panel_border(colour = 'black') +
  scale_y_log10(limits = c(1, 40), breaks = c(1, 10)) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'x', minor = 'none') +
  scale_x_continuous("Distance to minimal promoter (bp)", 
                     breaks = seq(from = 64, to = 194, by = 10),
                     limits = c(67, 191)) +
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

#Fit linear model to the different locations of CRE (site 1-6), allowing weight
#per weak and consensus CRE per site, and to the different background.

pred_resid <- function(df1, x) {
  df2 <- df1 %>%
    add_predictions(x)
  df3 <- df2 %>%
    add_residuals(x)
  return(df3)
  print('processed pre_res_trans_int(df1, df2) in order of (data, model)')
}

ind_site_ind_back <- function(df) {
  model <- lm(ave_ratio ~ background + site1 + site2 + site3 + site4 + site5 + site6 + 1, 
              data = df)
}

#change order of site types to fit consensus and weak CRE expression relative to
#no site and the weakest-expressing background

subpool5_ncw <- s5_gen_epi %>%
  mutate(site1 = gsub('nosite', 'anosite', site1)) %>%
  mutate(site2 = gsub('nosite', 'anosite', site2)) %>%
  mutate(site3 = gsub('nosite', 'anosite', site3)) %>%
  mutate(site4 = gsub('nosite', 'anosite', site4)) %>%
  mutate(site5 = gsub('nosite', 'anosite', site5)) %>%
  mutate(site6 = gsub('nosite', 'anosite', site6)) %>%
  mutate(background = gsub('v chr9', 'av chr9', background)) %>%
  var_log10()

#Make sum of errors output plottable, y-intercept represents the model without
#weights per variable, so the weight for no CRE sites and background 41

ind_site_ind_back_sumtidy <- function(df) {
  df <- tidy(df)
  sites <- df %>%
    filter(str_detect(term, '^site')) %>%
    mutate(term = gsub('consensus', '_consensus', term)) %>%
    mutate(term = gsub('weak', '_weak', term)) %>%
    separate(term, into = c('variable', 'type'), sep = "_")
  background <- df %>%
    filter(str_detect(term, '^background')) %>%
    mutate(term = gsub('background', 'background_', term)) %>%
    separate(term, into = c('variable', 'type'), sep = '_')
  back41 <- data.frame()
  sum <- rbind(sites, background)
  return(sum)
}

#Fit linear model to episomal data

ind_site_ind_back_epi <- subpool5_ncw %>%
  filter(MPRA == 'episomal') %>%
  ind_site_ind_back()

ind_site_ind_back_sum_epi <- ind_site_ind_back_sumtidy(ind_site_ind_back_epi)

ind_site_ind_back_anova_epi <- tidy(anova(ind_site_ind_back_epi)) %>%
  mutate(term_fctr = factor(term, levels = term)) %>%
  mutate(total_sumsq = sum(sumsq)) %>%
  mutate(per_sumsq = sumsq/total_sumsq)

ind_site_ind_back_p_r_epi <- pred_resid(filter(subpool5_ncw, MPRA == 'episomal'), 
                                        ind_site_ind_back_epi)

lessthan1_2color <- c('red', 'black', 'black', 'black', 'black', 'black', 'black')

p_ind_site_ind_back_epi <- ggplot(ind_site_ind_back_p_r_epi, 
                                  aes(ave_ratio, pred, 
                                      color = as.factor(consensus))) +
  geom_point(alpha = 0.2, size = 1, show.legend = FALSE) +
  scale_color_manual(values = lessthan1_2color) +
  scale_x_continuous(name = 'Measured expression', breaks = c(-1:1),
                     limits = c(-1.5, 1.8)) + 
  scale_y_continuous(name = 'Predicted expression', breaks = c(-1:1),
                     limits = c(-1.5, 1.8)) +
  annotation_logticks(sides = 'bl') +
  annotate("text", x = -0.5, y = 1, 
           label = paste('r =', 
                         round(cor(ind_site_ind_back_p_r_epi$pred,
                                   ind_site_ind_back_p_r_epi$ave_ratio,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2)))

(cor.test(ind_site_ind_back_p_r_epi$pred,
          ind_site_ind_back_p_r_epi$ave_ratio)$estimate)^2

p_ind_site_ind_back_sum_epi <- ind_site_ind_back_sum_epi %>%
  mutate(type = factor(type, 
                       levels = c('v chr9', 's pGl4', 'v chr5', 'consensus', 
                                  'weak'))) %>%
  ggplot(aes(variable, estimate, fill = type)) + 
  geom_bar(stat = 'identity', position = 'dodge', color = 'gray60', 
           size = 0.3) + 
  geom_hline(yintercept = 0, size = 0.25) +
  scale_x_discrete(position = 'bottom') + 
  scale_fill_viridis(discrete = TRUE) + 
  theme(axis.ticks.x = element_blank(), legend.position = 'top',
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank()) +
  ylab('Weight')

p_ind_site_ind_back_anova_epi <- ind_site_ind_back_anova_epi %>%
  ggplot(aes(term_fctr, per_sumsq)) + 
  geom_bar(stat = 'identity') + 
  ylab('Proportion of\nvariance explained') +
  scale_y_continuous(limits = c(0, 0.26), 
                     breaks = seq(from = 0, to = 0.25, by = 0.05)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank(), axis.title.x = element_blank())

ggsave('plots/p_ind_site_ind_back_epi.pdf', p_ind_site_ind_back_epi, 
       width = 2.3, height = 2.3, units = 'in')

ggsave('plots/p_ind_site_ind_back_sum_epi.pdf', p_ind_site_ind_back_sum_epi,
       width = 3, height = 2.5, units = 'in')

ggsave('plots/p_ind_site_ind_back_anova_epi.pdf', p_ind_site_ind_back_anova_epi,
       width = 2.5, height = 2.5)



