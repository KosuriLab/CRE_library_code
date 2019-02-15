#Written for the analysis of a range of inductions in the episomal library
#DNA_BC: DNA
#R0A_BC and R0B_BC: RNA at 0 µM Forsk replicate A and B
#R2#A_BC and R2#B_BC: RNA at 2^# µM Forsk replicate A and B, includes negative 
#and positive #'s

#Tested concentrations: 0, 2^-5, 2^-4, 2^-3, 2^-2, 2^-1, 2^0, 2^2 µM Forsk

#All figures made with this dataset require sections "Load index and bcmap 
#files" to "Median BC expression"

#Most figures require both genomic and episomal MPRAs, the genomic MPRA data
#processing is performed in 20171129_genlib.R and imported here in section
#"Import genomic MPRA and combine with episomal" and is used for figures 1E, 2, 
#4 and supplemental figures 2B, 3, and 4.

#Establish workspace------------------------------------------------------------

#import necessary libraries

library(tidyverse)
library(stringr)
library(viridis)
library(cowplot)
library(ggExtra)
library(modelr)
library(lazyeval)
library(splines)
library(broom)
library(GGally)
library(lemon)
library(devtools)
library(updateR)
library(ggsignif)

#General figure customizations

cbPalette7 <- c('#440154FF', '#39568CFF', '#287D8EFF', '#20A387FF', '#73D055FF',
                '#B8DE29FF', '#FDE725FF')

cbPalette7_grad_light <- c('white', '#FDE725FF', '#B8DE29FF', '#55C667FF', 
                           '#1F968BFF', '#39568CFF', '#482677FF')

spacing_5_20_palette <- c('gray20', 'dodgerblue3', 'indianred2', '#55C667FF')

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

SP3_SP5_map <- read_tsv('BCMap/uniqueSP2345.txt', 
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
                               as.integer(0), 
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

#Join DNA BC reads > 6 reads and join to RNA. Take ratio of RNA/DNA normalized
#reads per million

dna7_join_rna_rep <- function(df1, df2) {
  filter_DNA <- filter(df1, num_reads > 6)
  DNA_RNA_join <- left_join(filter_DNA, df2,
                            by = c("barcode", "name", "subpool", 
                                   "most_common"), 
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

#Figure 1D, episomal------------------------------------------------------------

#Plot replicability with backgrounds in orange and positive control in red

p_fig1_epi_med_rep <- med_rep_0_22_A_B %>%
  ggplot(aes(med_ratio_22A, med_ratio_22B)) +
  geom_point(alpha = 0.1, size = 0.75) +
  geom_point(data = filter(med_rep_0_22_A_B, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             fill = 'orange', shape = 21, size = 1.75) + 
  geom_point(data = filter(med_rep_0_22_A_B, 
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

ggsave('../plots/p_fig1_epi_med_rep.png', p_fig1_epi_med_rep,
       width = 2, height = 2, units = 'in')

log10_med_rep_0_22_A_B <- var_log10(med_rep_0_22_A_B)

pearsons <- tibble(
  sample = c('all', 'subpool3', 'subpool5'),
  pearsons = c(cor(log10_med_rep_0_22_A_B$med_ratio_22A, 
                   log10_med_rep_0_22_A_B$med_ratio_22B, 
                   use = "pairwise.complete.obs", method = "pearson"),
               cor(filter(log10_med_rep_0_22_A_B, 
                          subpool == 'subpool3')$med_ratio_22A,
                   filter(log10_med_rep_0_22_A_B, 
                          subpool == 'subpool3')$med_ratio_22B,
                   use = "pairwise.complete.obs", method = "pearson"),
               cor(filter(log10_med_rep_0_22_A_B, 
                          subpool == 'subpool5')$med_ratio_22A,
                   filter(log10_med_rep_0_22_A_B, 
                          subpool == 'subpool5')$med_ratio_22B,
                   use = "pairwise.complete.obs", method = "pearson")))

#Supplemental Figure 1B---------------------------------------------------------

#Supplemental Figure 1B, subfigures A and B

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

#Supplemental Figure 1B, subfigure C

#Normalize MPRA expression of all variants and 1 control to variant without CRE
#sites (backgrounds). Each variant normalized to the background used in its 
#design

back_norm <- function(df1) {
  gsub_0_22 <- df1 %>%
    ungroup () %>%
    filter(subpool != 'control') %>%
    mutate(
      name = gsub('Smith R. Vista chr9:83712599-83712766', 'v chr9', name),
      name = gsub('Vista Chr5:88673410-88674494', 'v chr5', name),
      name = gsub('scramble pGL4.29 Promega 1-63 \\+ 1-87', 's pGl4', name)
    ) %>%
    mutate(background = name) %>%
    mutate(background = str_sub(background, 
                                nchar(background)-5, 
                                nchar(background)))
  backgrounds <- gsub_0_22 %>%
    filter(startsWith(name, 
                      'subpool5_no_site_no_site_no_site_no_site_no_site_no_site')) %>%
    select(background, med_ratio_0A, med_ratio_0B, med_ratio_2_5A, med_ratio_2_5B, med_ratio_2_4A, 
           med_ratio_2_4B, med_ratio_2_3A, med_ratio_2_3B, med_ratio_2_2A, med_ratio_2_2B, 
           med_ratio_2_1A, med_ratio_2_1B, med_ratio_20A, med_ratio_20B, med_ratio_22A, 
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

#Add back in control that can be normalized to a background for plotting

epi_back_norm_pc_spGl4 <- med_rep_0_22_A_B %>%
  filter(name == 'pGL4.29 Promega 1-63 + 1-87') %>%
  mutate(name = str_c(name, '_scramble pGL4.29 Promega 1-63 + 1-87')) %>%
  mutate(subpool = 'subpool3') %>%
  rbind(med_rep_0_22_A_B) %>%
  back_norm()

#Make untidy df with conc as a variable and expression as a single column. Conc
#is indicated as log2 and 0 µM forskolin concentration is represented as 2^-7
#for plotting

var_conc_exp <- function(df) {
  df_0 <- df %>%
    mutate(ave_barcode_0 = (barcodes_RNA_0A + barcodes_RNA_0B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_0, 
           ave_ratio_0_norm) %>%
    mutate(conc = -7) %>%
    rename(ave_ratio_norm = ave_ratio_0_norm) %>%
    rename(ave_barcode = ave_barcode_0)
  df_2_5 <- df %>%
    mutate(ave_barcode_2_5 = (barcodes_RNA_2_5A + barcodes_RNA_2_5B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_2_5, 
           ave_ratio_2_5_norm) %>%
    mutate(conc = -5) %>%
    rename(ave_ratio_norm = ave_ratio_2_5_norm) %>%
    rename(ave_barcode = ave_barcode_2_5)
  df_2_4 <- df %>%
    mutate(ave_barcode_2_4 = (barcodes_RNA_2_4A + barcodes_RNA_2_4B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_2_4, 
           ave_ratio_2_4_norm) %>%
    mutate(conc = -4) %>%
    rename(ave_ratio_norm = ave_ratio_2_4_norm) %>%
    rename(ave_barcode = ave_barcode_2_4)
  df_2_3 <- df %>%
    mutate(ave_barcode_2_3 = (barcodes_RNA_2_3A + barcodes_RNA_2_3B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_2_3, 
           ave_ratio_2_3_norm) %>%
    mutate(conc = -3) %>%
    rename(ave_ratio_norm = ave_ratio_2_3_norm) %>%
    rename(ave_barcode = ave_barcode_2_3)
  df_2_2 <- df %>%
    mutate(ave_barcode_2_2 = (barcodes_RNA_2_2A + barcodes_RNA_2_2B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_2_2, 
           ave_ratio_2_2_norm) %>%
    mutate(conc = -2) %>%
    rename(ave_ratio_norm = ave_ratio_2_2_norm) %>%
    rename(ave_barcode = ave_barcode_2_2)
  df_2_1 <- df %>%
    mutate(ave_barcode_2_1 = (barcodes_RNA_2_1A + barcodes_RNA_2_1B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_2_1, 
           ave_ratio_2_1_norm) %>%
    mutate(conc = -1) %>%
    rename(ave_ratio_norm = ave_ratio_2_1_norm) %>%
    rename(ave_barcode = ave_barcode_2_1)
  df_20 <- df %>%
    mutate(ave_barcode_20 = (barcodes_RNA_20A + barcodes_RNA_20B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_20, 
           ave_ratio_20_norm) %>%
    mutate(conc = 0) %>%
    rename(ave_ratio_norm = ave_ratio_20_norm) %>%
    rename(ave_barcode = ave_barcode_20)
  df_22 <- df %>%
    mutate(ave_barcode_22 = (barcodes_RNA_22A + barcodes_RNA_22B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_22, 
           ave_ratio_22_norm) %>%
    mutate(conc = 2) %>%
    rename(ave_ratio_norm = ave_ratio_22_norm) %>%
    rename(ave_barcode = ave_barcode_22)
  df_0_22 <- rbind(df_0, df_2_5, df_2_4, df_2_3, df_2_2, df_2_1, df_20, df_22)
  return(df_0_22)
}

#Subtract expression at 0 µM forskolin for expression at each concentration

epi_back_0_norm_conc <- epi_back_norm_pc_spGl4 %>%
  mutate(ave_ratio_2_5_norm = ave_ratio_2_5_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_2_4_norm = ave_ratio_2_4_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_2_3_norm = ave_ratio_2_3_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_2_2_norm = ave_ratio_2_2_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_2_1_norm = ave_ratio_2_1_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_20_norm = ave_ratio_20_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_22_norm = ave_ratio_22_norm - ave_ratio_0_norm) %>%
  mutate(ave_ratio_0_norm = ave_ratio_0_norm - ave_ratio_0_norm) %>%
  var_conc_exp()

#Plot subfigure C, normalized variant expression curves across forskolin
#concentrations. Expression curves for backgrounds and control overlayed

p_titr_pc_back <- epi_back_0_norm_conc %>%
  ggplot(aes(conc, ave_ratio_norm)) +
  geom_line(aes(group = name), alpha = 0.1) +
  geom_point(data = filter(epi_back_0_norm_conc, 
                           startsWith(name, 
                                      'subpool5_no_site_no_site_no_site_no_site_no_site_no_site')),
             color = 'darkgoldenrod1', shape = 19, stroke = 0.75) +
  geom_line(data = filter(epi_back_0_norm_conc, 
                          startsWith(name, 
                                     'subpool5_no_site_no_site_no_site_no_site_no_site_no_site')),
            color = 'darkgoldenrod1', size = 1) +
  geom_point(data = filter(epi_back_0_norm_conc, 
                           startsWith(name, 
                                      'pGL4.29 Promega 1-63 + 1-87')),
             color = 'firebrick2', shape = 19, stroke = 0.75) +
  geom_line(data = filter(epi_back_0_norm_conc, 
                          startsWith(name, 
                                     'pGL4.29 Promega 1-63 + 1-87')),
            color = 'firebrick2', size = 1) +
  ylab('Average normalized\nexpression (a.u.)') +
  annotation_logticks(sides = 'b', short = unit(0.05, 'cm'), 
                      mid = unit(0.1, 'cm'), long = unit(0.15, 'cm')) +
  scale_x_continuous(breaks = (-7:2), 'log2 forskolin (µM)') +
  figurefont_theme

ggsave('../plots/p_titr_pc_back.pdf', p_titr_pc_back, width = 3.6, height = 2,
       units = 'in')

#Import genomic MPRA and combine with episomal----------------------------------

#Genomic MPRA expression determination and data processing performed in
#"20171129_genlib_analysis/20171129_genlib.R". Processed df with expression
#values is exported from there and imported here to plot together

gen_rep_1_2 <- read_tsv('../20171129_genlib_analysis/rep_1_2.txt') %>%
  mutate(ave_med_ratio = (med_ratio_br1 + med_ratio_br2)/2)

#Combine genomic and episomal dfs, only comparing to expression at 2^2 µM 
#forskolin in the episomal dataset. Here med_ratio_br# and ave_med_ratio refers 
#to genomic expression, either across biological replicates or averaged.
#Episomal expression is represented by the annotation <name>_22 as representing
#the forskolin concentration the sample was incubated with

gen_epi <- med_rep_0_22_A_B %>%
  select(subpool, name, most_common, barcodes_DNA, med_ratio_22A, 
         barcodes_RNA_22A, med_ratio_22B, barcodes_RNA_22B) %>%
  mutate(ave_ratio_22 = (med_ratio_22A + med_ratio_22B)/2) %>%
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
         med_ratio_br1, med_ratio_br2, ave_med_ratio, barcodes_RNA_22A, 
         barcodes_RNA_22B, med_ratio_22A, med_ratio_22B, ave_ratio_22) %>%
  mutate(genomic = (barcodes_RNA_br1 + barcodes_RNA_br2)/2) %>%
  mutate(episomal = (barcodes_RNA_22A + barcodes_RNA_22B)/2) %>%
  select(-barcodes_RNA_br1, -barcodes_RNA_br2, -barcodes_RNA_22A, 
         -barcodes_RNA_22B) %>%
  gather(genomic, episomal, key = 'MPRA', value = 'barcodes') %>%
  mutate(genomic = ave_med_ratio) %>%
  mutate(episomal = ave_ratio_22) %>%
  gather(genomic, episomal, key = 'MPRA2', value = 'ave_ratio') %>%
  filter((MPRA == 'genomic' & MPRA2 == 'genomic') | (MPRA == 'episomal' & MPRA2 == 'episomal')) %>%
  select(-MPRA2)

#Figure 1E----------------------------------------------------------------------

#Plot correlation between MPRAs with backgrounds in orange and positive control 
#in red

p_gen_epi_rep <- gen_epi %>%
  ggplot(aes(ave_med_ratio, ave_ratio_22)) +
  geom_point(alpha = 0.1, size = 0.75) +
  geom_point(data = filter(gen_epi, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             fill = 'orange', shape = 21, size = 1.75) +
  geom_point(data = filter(gen_epi, 
                           name == 'pGL4.29 Promega 1-63 + 1-87'), 
             fill = 'red', shape = 21, size = 1.75) +
  xlab("Average genomic expression (a.u.)") +
  ylab("Average episomal expression (a.u.)") +
  scale_x_log10(limits = c(0.01, 20), breaks = c(0.01, 0.1, 1, 10)) + 
  scale_y_log10(limits = c(0.05, 20), breaks = c(0.1, 1, 10)) +
  annotation_logticks(scaled = TRUE) +
  figurefont_theme

ggsave('../plots/p_gen_epi_rep.png', p_gen_epi_rep, width = 2, height = 2, 
       units = 'in')

gen_epi_log10 <- var_log10(gen_epi)

gen_epi_pearsons <- tibble(
  sample = c('all', 'subpool3', 'subpool5'),
  pearsons = c(round(cor(gen_epi_log10$ave_med_ratio, 
                         gen_epi_log10$ave_ratio_22, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(filter(gen_epi_log10, 
                                subpool == 'subpool3')$ave_med_ratio,
                         filter(gen_epi_log10, 
                                subpool == 'subpool3')$ave_ratio_22,
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(filter(gen_epi_log10, 
                                subpool == 'subpool5')$ave_med_ratio,
                         filter(gen_epi_log10, 
                                subpool == 'subpool5')$ave_ratio_22,
                         use = "pairwise.complete.obs", method = "pearson"), 
                     3)))

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
  filter(subpool == 'subpool3') %>%
  subpool3()


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
  filter(subpool == 'subpool5') %>%
  subpool5()

#Figure 2A----------------------------------------------------------------------

#Plot CRE distance vs. variant expression in the CRE Spacing and Distance 
#Library. Distance is binned in 22 bp intervals and plotted per MPRA and
#background, colored according to CRE Spacing.

s3_gen_epi_bin20bp <- s3_gen_epi %>%
  filter(dist <= 176) %>%
  mutate(bin = cut(dist, seq(from = 66, to = 176, by = 22),
                   labels = c('67-88', '89-110', '111-132', '133-154',
                              '155-176')))

p_s3_dist_gen_bin20bp <- s3_gen_epi_bin20bp %>%
  mutate(background = factor(background, levels = c('55', '52', '41'))) %>%
  filter(spacing != 0 & spacing != 70 & background != '41' & MPRA == 'genomic') %>%
  ggplot(aes(bin, ave_ratio)) +
  geom_signif(comparisons = list(c('67-88', '89-110')), y_position = log10(5),
              textsize = 2.4) +
  geom_signif(comparisons = list(c('67-88', '111-132')), y_position = log10(15),
              textsize = 2.4) +
  facet_grid(. ~ background) +
  geom_jitter(aes(color = as.factor(spacing)), 
              position=position_jitter(width=0.3, height=0), alpha = 0.75,
              size = 0.5) +
  geom_boxplot(outlier.shape=NA, size = 0.3, position = position_dodge(1),
               show.legend = FALSE, alpha = 0) +
  scale_color_manual(values = spacing_5_20_palette, name = 'CRE spacing (bp)') +
  theme(legend.position = 'right', axis.ticks.x = element_blank(), 
        strip.background = element_rect(colour="black", fill="white"),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_log10(limits = c(0.015, 30)) +
  annotation_logticks(sides = 'l') +
  panel_border(colour = 'black') +
  ylab('Average expression (a.u.)') +
  xlab('Distance to minimal promoter (bp)') +
  figurefont_theme

p_s3_dist_epi_bin20bp <- s3_gen_epi_bin20bp %>%
  mutate(background = factor(background, levels = c('55', '52', '41'))) %>%
  filter(spacing != 0 & spacing != 70 & background != '41' & MPRA == 'episomal') %>%
  ggplot(aes(bin, ave_ratio)) +
  geom_signif(comparisons = list(c('67-88', '89-110')), y_position = log10(6),
              textsize = 2.4) +
  geom_signif(comparisons = list(c('67-88', '111-132')), y_position = log10(15),
              textsize = 2.4) +
  facet_grid(. ~ background) +
  geom_jitter(aes(color = as.factor(spacing)), 
              position=position_jitter(width=0.3, height=0), alpha = 0.75,
              size = 0.5) +
  geom_boxplot(outlier.shape=NA, size = 0.3, position = position_dodge(1),
               show.legend = FALSE, alpha = 0) +
  scale_color_manual(values = spacing_5_20_palette, name = 'CRE spacing (bp)') +
  theme(legend.position = 'right', axis.ticks.x = element_blank(), 
        strip.background = element_rect(colour="black", fill="white"),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_log10(limits = c(0.1, 30)) +
  annotation_logticks(sides = 'l') +
  panel_border(colour = 'black') +
  ylab('Average expression (a.u.)') +
  xlab('Distance to minimal promoter (bp)') +
  figurefont_theme

ggsave('../plots/p_s3_dist_gen_bin20bp.pdf', p_s3_dist_gen_bin20bp,
       width = 3.375, height = 1.75, unit = 'in')

ggsave('../plots/p_s3_dist_epi_bin20bp.pdf', p_s3_dist_epi_bin20bp,
       width = 3.375, height = 1.75, unit = 'in')

#Figure 2B----------------------------------------------------------------------

#Take median expression of variants with a certain background, MPRA format,
#CRE Spacing per 10 bp CRE Distance increment. Plot change in median expression
#between 67-76 and 167-176 bp CRE Distance increments.

bin_10bp_change_dist <- function(df) {
  df <- df %>%
    filter(dist <= 176) %>%
    mutate(bin = cut(dist, seq(from = 66, to = 177, by = 10),
                     labels = c('67-76', '77-86', '87-96', '97-106', '107-116',
                                '117-126', '127-136', '137-146', '147-156', 
                                '157-166', '167-176'))) %>%
    group_by(background, MPRA, bin, spacing) %>%
    summarize(median_bin = median(ave_ratio)) %>%
    ungroup()
  df_6776 <- df %>%
    filter(bin == '67-76') %>%
    mutate(median_6776 = median_bin) %>%
    select(-bin, -median_bin)
  df_167176 <- df %>%
    filter(bin == '167-176') %>%
    mutate(median_167176 = median_bin) %>%
    select(-bin, -median_bin)
  bin_change <- left_join(df_6776, df_167176, by = c('background', 'MPRA',
                                                     'spacing')) %>%
    mutate(bin_change = median_6776/median_167176)
  return(bin_change)
}
  
s3_gen_epi_bin_10bp_change_dist <- bin_10bp_change_dist(s3_gen_epi)
  
p_s3_gen_epi_bin10bp_med_change <- s3_gen_epi_bin_10bp_change_dist %>%
  mutate(background = factor(background, levels = c('55', '52', '41'))) %>%
  filter(spacing != 0 & spacing != 70 & background != '41') %>%
  ggplot(aes(MPRA, bin_change)) +
  facet_grid(. ~ background) +
  geom_jitter(aes(color = as.factor(spacing)), 
              position=position_jitter(width=0.3, height=0), alpha = 0.75,
              size = 1) +
  geom_boxplot(outlier.shape=NA, size = 0.3, position = position_dodge(1),
               show.legend = FALSE, alpha = 0) +
  scale_color_manual(values = spacing_5_20_palette, name = 'CRE spacing (bp)') +
  theme(legend.position = 'top', axis.ticks.x = element_blank(), 
        strip.background = element_rect(colour="black", fill="white"),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  panel_border(colour = 'black') +
  scale_y_continuous(limits = c(1, 6), breaks = c(1:6)) +
  ylab('Change in median\nexpression from\n167-176 to 67-76 bp') +
  xlab('MPRA') +
  figurefont_theme

ggsave('../plots/p_s3_gen_epi_bin10bp_med_change.pdf', 
       p_s3_gen_epi_bin10bp_med_change, width = 3.25, height = 2.75, unit = 'in')




