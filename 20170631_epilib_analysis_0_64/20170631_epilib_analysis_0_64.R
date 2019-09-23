#Written for the analysis of the episomal library from 0-64 µM forskolin
#D3_BC: DNA
#R#A_BC: RNA at # µM Forsk replicate A
#R#B_BC: RNA at # µM Forsk replicate B

#tested concentrations: 0, 1, 2, 4, 8, 16, 25, 32, 64

library(tidyverse)

#Load index and bcmap files-----------------------------------------------------

bc_DNA <- read_tsv('BCreads_txts/20170631_DNA_BC.txt')
bc_0A <- read_tsv('BCreads_txts/20170631_R0A_BC.txt')
bc_0B <- read_tsv('BCreads_txts/20170631_R0B_BC.txt')
bc_1A <- read_tsv('BCreads_txts/20170631_R1A_BC.txt')
bc_1B <- read_tsv('BCreads_txts/20170631_R1B_BC.txt')
bc_2A <- read_tsv('BCreads_txts/20170631_R2A_BC.txt')
bc_2B <- read_tsv('BCreads_txts/20170631_R2B_BC.txt')
bc_4A <- read_tsv('BCreads_txts/20170631_R4A_BC.txt')
bc_4B <- read_tsv('BCreads_txts/20170631_R4B_BC.txt')
bc_8A <- read_tsv('BCreads_txts/20170631_R8A_BC.txt')
bc_8B <- read_tsv('BCreads_txts/20170631_R8B_BC.txt')
bc_16A <- read_tsv('BCreads_txts/20170631_R16A_BC.txt')
bc_16B <- read_tsv('BCreads_txts/20170631_R16B_BC.txt')
bc_25A <- read_tsv('BCreads_txts/20170631_R25A_BC.txt')
bc_25B <- read_tsv('BCreads_txts/20170631_R25B_BC.txt')
bc_32A <- read_tsv('BCreads_txts/20170631_R32A_BC.txt')
bc_32B <- read_tsv('BCreads_txts/20170631_R32B_BC.txt')
bc_64A <- read_tsv('BCreads_txts/20170631_R64A_BC.txt')
bc_64B <- read_tsv('BCreads_txts/20170631_R64B_BC.txt')

#Load barcode mapping table, sequences (most_common) are rcomp due to sequencing
#format. Pick out controls, SP3, and SP5 in the bcmap that were used in this 
#assay

SP3_SP5_map <- read_tsv('../BCMap/uniqueSP2345.txt', 
                        col_names = c(
                          'fluff', 'barcode', 'name', 'most_common'
                        ), 
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
bc_join_0A <- bc_map_join_bc(SP3_SP5_map, bc_0A)
bc_join_0B <- bc_map_join_bc(SP3_SP5_map, bc_0B)
bc_join_1A <- bc_map_join_bc(SP3_SP5_map, bc_1A)
bc_join_1B <- bc_map_join_bc(SP3_SP5_map, bc_1B)
bc_join_2A <- bc_map_join_bc(SP3_SP5_map, bc_2A)
bc_join_2B <- bc_map_join_bc(SP3_SP5_map, bc_2B)
bc_join_4A <- bc_map_join_bc(SP3_SP5_map, bc_4A)
bc_join_4B <- bc_map_join_bc(SP3_SP5_map, bc_4B)
bc_join_8A <- bc_map_join_bc(SP3_SP5_map, bc_8A)
bc_join_8B <- bc_map_join_bc(SP3_SP5_map, bc_8B)
bc_join_16A <- bc_map_join_bc(SP3_SP5_map, bc_16A)
bc_join_16B <- bc_map_join_bc(SP3_SP5_map, bc_16B)
bc_join_25A <- bc_map_join_bc(SP3_SP5_map, bc_25A)
bc_join_25B <- bc_map_join_bc(SP3_SP5_map, bc_25B)
bc_join_32A <- bc_map_join_bc(SP3_SP5_map, bc_32A)
bc_join_32B <- bc_map_join_bc(SP3_SP5_map, bc_32B)
bc_join_64A <- bc_map_join_bc(SP3_SP5_map, bc_64A)
bc_join_64B <- bc_map_join_bc(SP3_SP5_map, bc_64B)

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

bc_DNA_RNA_0A <- dna7_join_rna_rep(bc_join_DNA, bc_join_0A)
bc_DNA_RNA_0B <- dna7_join_rna_rep(bc_join_DNA, bc_join_0B)
bc_DNA_RNA_1A <- dna7_join_rna_rep(bc_join_DNA, bc_join_1A)
bc_DNA_RNA_1B <- dna7_join_rna_rep(bc_join_DNA, bc_join_1B)
bc_DNA_RNA_2A <- dna7_join_rna_rep(bc_join_DNA, bc_join_2A)
bc_DNA_RNA_2B <- dna7_join_rna_rep(bc_join_DNA, bc_join_2B)
bc_DNA_RNA_4A <- dna7_join_rna_rep(bc_join_DNA, bc_join_4A)
bc_DNA_RNA_4B <- dna7_join_rna_rep(bc_join_DNA, bc_join_4B)
bc_DNA_RNA_8A <- dna7_join_rna_rep(bc_join_DNA, bc_join_8A)
bc_DNA_RNA_8B <- dna7_join_rna_rep(bc_join_DNA, bc_join_8B)
bc_DNA_RNA_16A <- dna7_join_rna_rep(bc_join_DNA, bc_join_16A)
bc_DNA_RNA_16B <- dna7_join_rna_rep(bc_join_DNA, bc_join_16B)
bc_DNA_RNA_25A <- dna7_join_rna_rep(bc_join_DNA, bc_join_25A)
bc_DNA_RNA_25B <- dna7_join_rna_rep(bc_join_DNA, bc_join_25B)
bc_DNA_RNA_32A <- dna7_join_rna_rep(bc_join_DNA, bc_join_32A)
bc_DNA_RNA_32B <- dna7_join_rna_rep(bc_join_DNA, bc_join_32B)
bc_DNA_RNA_64A <- dna7_join_rna_rep(bc_join_DNA, bc_join_64A)
bc_DNA_RNA_64B <- dna7_join_rna_rep(bc_join_DNA, bc_join_64B)

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
med_ratio_R1A <- ratio_bc_med_var(bc_DNA_RNA_1A)
med_ratio_R1B <- ratio_bc_med_var(bc_DNA_RNA_1B)
med_ratio_R2A <- ratio_bc_med_var(bc_DNA_RNA_2A)
med_ratio_R2B <- ratio_bc_med_var(bc_DNA_RNA_2B)
med_ratio_R4A <- ratio_bc_med_var(bc_DNA_RNA_4A)
med_ratio_R4B <- ratio_bc_med_var(bc_DNA_RNA_4B)
med_ratio_R8A <- ratio_bc_med_var(bc_DNA_RNA_8A)
med_ratio_R8B <- ratio_bc_med_var(bc_DNA_RNA_8B)
med_ratio_R16A <- ratio_bc_med_var(bc_DNA_RNA_16A)
med_ratio_R16B <- ratio_bc_med_var(bc_DNA_RNA_16B)
med_ratio_R25A <- ratio_bc_med_var(bc_DNA_RNA_25A)
med_ratio_R25B <- ratio_bc_med_var(bc_DNA_RNA_25B)
med_ratio_R32A <- ratio_bc_med_var(bc_DNA_RNA_32A)
med_ratio_R32B <- ratio_bc_med_var(bc_DNA_RNA_32B)
med_ratio_R64A <- ratio_bc_med_var(bc_DNA_RNA_64A)
med_ratio_R64B <- ratio_bc_med_var(bc_DNA_RNA_64B)

#Combine biological replicates

var_conc_rep_med <- function(df0A, df0B, df1A, df1B, df2A, df2B, df4A, df4B, 
                             df8A, df8B, df16A, df16B, df25A, df25B, df32A, 
                             df32B, df64A, df64B) {
  join_0 <- inner_join(df0A, df0B, 
                       by = c("name", "subpool", "most_common", "barcodes_DNA"),
                       suffix = c("_0A", "_0B"))
  join_1 <- inner_join(df1A, df1B, 
                       by = c("name", "subpool", "most_common", "barcodes_DNA"),
                       suffix = c("_1A", "_1B"))
  join_2 <- inner_join(df2A, df2B, 
                       by = c("name", "subpool", "most_common", "barcodes_DNA"),
                       suffix = c("_2A", "_2B"))
  join_4 <- inner_join(df4A, df4B, 
                       by = c("name", "subpool", "most_common", "barcodes_DNA"),
                       suffix = c("_4A", "_4B"))
  join_8 <- inner_join(df8A, df8B, 
                       by = c("name", "subpool", "most_common", "barcodes_DNA"),
                       suffix = c("_8A", "_8B"))
  join_16 <- inner_join(df16A, df16B, 
                        by = c("name", "subpool", "most_common", 
                               "barcodes_DNA"),
                        suffix = c("_16A", "_16B"))
  join_25 <- inner_join(df25A, df25B, 
                        by = c("name", "subpool", "most_common", 
                               "barcodes_DNA"),
                        suffix = c("_25A", "_25B"))
  join_32 <- inner_join(df32A, df32B, 
                        by = c("name", "subpool", "most_common", 
                               "barcodes_DNA"),
                        suffix = c("_32A", "_32B"))
  join_64 <- inner_join(df64A, df64B, 
                        by = c("name", "subpool", "most_common", 
                               "barcodes_DNA"),
                        suffix = c("_64A", "_64B"))
  join_0_1 <- inner_join(join_0, join_1, 
                         by = c("name", "subpool", "most_common", 
                                "barcodes_DNA"))
  join_0_2 <- inner_join(join_0_1, join_2, 
                         by = c("name", "subpool", "most_common", 
                                "barcodes_DNA"))
  join_0_4 <- inner_join(join_0_2, join_4, 
                         by = c("name", "subpool", "most_common", 
                                "barcodes_DNA"))
  join_0_8 <- inner_join(join_0_4, join_8, 
                         by = c("name", "subpool", "most_common", 
                                "barcodes_DNA"))
  join_0_16 <- inner_join(join_0_8, join_16, 
                          by = c("name", "subpool", "most_common", 
                                 "barcodes_DNA"))
  join_0_25 <- inner_join(join_0_16, join_25, 
                          by = c("name", "subpool", "most_common", 
                                 "barcodes_DNA"))
  join_0_32 <- inner_join(join_0_25, join_32, 
                          by = c("name", "subpool", "most_common", 
                                 "barcodes_DNA"))
  join_0_64 <- inner_join(join_0_32, join_64, 
                          by = c("name", "subpool", "most_common", 
                                 "barcodes_DNA"))
  print('processed dfs in order of samples: 0A, 0B, 1A, 1B, 2A, 2B, 4A, 4B, 8A, 
        8B, 16A, 16B, 25A, 25B, 32A, 32B, 64A, 64B')
  return(join_0_64)
}

med_rep_0_64_A_B <- var_conc_rep_med(med_ratio_R0A, med_ratio_R0B, 
                                     med_ratio_R1A, med_ratio_R1B, 
                                     med_ratio_R2A, med_ratio_R2B, 
                                     med_ratio_R4A, med_ratio_R4B, 
                                     med_ratio_R8A, med_ratio_R8B, 
                                     med_ratio_R16A, med_ratio_R16B,
                                     med_ratio_R25A, med_ratio_R25B, 
                                     med_ratio_R32A, med_ratio_R32B, 
                                     med_ratio_R64A, med_ratio_R64B
)

write.table(med_rep_0_64_A_B,"rep_0_64_med.txt", sep = '\t', row.names = FALSE)

#Background-normalize expression------------------------------------------------

#Normalize MPRA expression of all variants and 1 control to variant without CRE
#sites (backgrounds). Each variant normalized to the background used in its 
#design.

back_norm <- function(df1) {
  gsub_0_64 <- df1 %>%
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
  backgrounds <- gsub_0_64 %>%
    filter(startsWith(name, 
                      'subpool5_no_site_no_site_no_site_no_site_no_site_no_site'
    )) %>%
    select(background, med_ratio_0A, med_ratio_0B, med_ratio_1A, med_ratio_1B, 
           med_ratio_2A, med_ratio_2B, med_ratio_4A, med_ratio_4B, med_ratio_8A,
           med_ratio_8B, med_ratio_16A, med_ratio_16B, med_ratio_25A, 
           med_ratio_25B, med_ratio_32A, med_ratio_32B, med_ratio_64A, 
           med_ratio_64B) %>%
    rename(ratio_0A_back = med_ratio_0A) %>%
    rename(ratio_0B_back = med_ratio_0B) %>%
    rename(ratio_1A_back = med_ratio_1A) %>%
    rename(ratio_1B_back = med_ratio_1B) %>%
    rename(ratio_2A_back = med_ratio_2A) %>%
    rename(ratio_2B_back = med_ratio_2B) %>%
    rename(ratio_4A_back = med_ratio_4A) %>%
    rename(ratio_4B_back = med_ratio_4B) %>%
    rename(ratio_8A_back = med_ratio_8A) %>%
    rename(ratio_8B_back = med_ratio_8B) %>%
    rename(ratio_16A_back = med_ratio_16A) %>%
    rename(ratio_16B_back = med_ratio_16B) %>%
    rename(ratio_25A_back = med_ratio_25A) %>%
    rename(ratio_25B_back = med_ratio_25B) %>%
    rename(ratio_32A_back = med_ratio_32A) %>%
    rename(ratio_32B_back = med_ratio_32B)%>%
    rename(ratio_64A_back = med_ratio_64A) %>%
    rename(ratio_64B_back = med_ratio_64B)
  back_join_norm <- left_join(gsub_0_64, backgrounds, by = 'background') %>%
    mutate(ratio_0A_norm = med_ratio_0A/ratio_0A_back) %>%
    mutate(ratio_0B_norm = med_ratio_0B/ratio_0B_back) %>%
    mutate(ratio_1A_norm = med_ratio_1A/ratio_1A_back) %>%
    mutate(ratio_1B_norm = med_ratio_1B/ratio_1B_back) %>%
    mutate(ratio_2A_norm = med_ratio_2A/ratio_2A_back) %>%
    mutate(ratio_2B_norm = med_ratio_2B/ratio_2B_back) %>%
    mutate(ratio_4A_norm = med_ratio_4A/ratio_4A_back) %>%
    mutate(ratio_4B_norm = med_ratio_4B/ratio_4B_back) %>%
    mutate(ratio_8A_norm = med_ratio_8A/ratio_8A_back) %>%
    mutate(ratio_8B_norm = med_ratio_8B/ratio_8B_back) %>%
    mutate(ratio_16A_norm = med_ratio_16A/ratio_16A_back) %>%
    mutate(ratio_16B_norm = med_ratio_16B/ratio_16B_back) %>%
    mutate(ratio_25A_norm = med_ratio_25A/ratio_25A_back) %>%
    mutate(ratio_25B_norm = med_ratio_25B/ratio_25B_back) %>%
    mutate(ratio_32A_norm = med_ratio_32A/ratio_32A_back) %>%
    mutate(ratio_32B_norm = med_ratio_32B/ratio_32B_back) %>%
    mutate(ratio_64A_norm = med_ratio_64A/ratio_64A_back) %>%
    mutate(ratio_64B_norm = med_ratio_64B/ratio_64B_back) %>%
    mutate(ave_ratio_0_norm = (ratio_0A_norm + ratio_0B_norm)/2) %>%
    mutate(ave_ratio_1_norm = (ratio_1A_norm + ratio_1B_norm)/2) %>%
    mutate(ave_ratio_2_norm = (ratio_2A_norm + ratio_2B_norm)/2) %>%
    mutate(ave_ratio_4_norm = (ratio_4A_norm + ratio_4B_norm)/2) %>%
    mutate(ave_ratio_8_norm = (ratio_8A_norm + ratio_8B_norm)/2) %>%
    mutate(ave_ratio_16_norm = (ratio_16A_norm + ratio_16B_norm)/2) %>%
    mutate(ave_ratio_25_norm = (ratio_25A_norm + ratio_25B_norm)/2) %>%
    mutate(ave_ratio_32_norm = (ratio_32A_norm + ratio_32B_norm)/2) %>%
    mutate(ave_ratio_64_norm = (ratio_64A_norm + ratio_64B_norm)/2) %>%
    mutate(induction = ave_ratio_64_norm/ave_ratio_0_norm)
  return(back_join_norm)
}

#Add back to df CRE control that can be normalized to a background for plotting

epi_back_norm_pc_spGl4 <- med_rep_0_64_A_B %>%
  filter(name == 'pGL4.29 Promega 1-63 + 1-87') %>%
  mutate(name = str_c(name, '_scramble pGL4.29 Promega 1-63 + 1-87')) %>%
  mutate(subpool = 'subpool3') %>%
  rbind(med_rep_0_64_A_B) %>%
  back_norm()

write.table(epi_back_norm_pc_spGl4, "rep_0_64_med_back.txt", 
            sep = '\t', row.names = FALSE)

#Make untidy df with conc as a variable and normalized expressions as single 
#columns. Conc is indicated in log2 and 0 µM forskolin concentration is 
#represented as 2^-1 for plotting on log scales.

var_conc_exp <- function(df) {
  df_0 <- df %>%
    mutate(ave_barcode_0 = (barcodes_RNA_0A + barcodes_RNA_0B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_0, 
           ratio_0A_norm, ratio_0B_norm, ave_ratio_0_norm) %>%
    mutate(conc = 0.25) %>%
    rename(ratio_A_norm = ratio_0A_norm) %>%
    rename(ratio_B_norm = ratio_0B_norm) %>%
    rename(ave_ratio_norm = ave_ratio_0_norm) %>%
    rename(ave_barcode = ave_barcode_0)
  df_1 <- df %>%
    mutate(ave_barcode_1 = (barcodes_RNA_1A + barcodes_RNA_1B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_1, 
           ratio_1A_norm, ratio_1B_norm, ave_ratio_1_norm) %>%
    mutate(conc = 1) %>%
    rename(ratio_A_norm = ratio_1A_norm) %>%
    rename(ratio_B_norm = ratio_1B_norm) %>%
    rename(ave_ratio_norm = ave_ratio_1_norm) %>%
    rename(ave_barcode = ave_barcode_1)
  df_2 <- df %>%
    mutate(ave_barcode_2 = (barcodes_RNA_2A + barcodes_RNA_2B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_2, 
           ratio_2A_norm, ratio_2B_norm, ave_ratio_2_norm) %>%
    mutate(conc = 2) %>%
    rename(ratio_A_norm = ratio_2A_norm) %>%
    rename(ratio_B_norm = ratio_2B_norm) %>%
    rename(ave_ratio_norm = ave_ratio_2_norm) %>%
    rename(ave_barcode = ave_barcode_2)
  df_4 <- df %>%
    mutate(ave_barcode_4 = (barcodes_RNA_4A + barcodes_RNA_4B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_4, 
           ratio_4A_norm, ratio_4B_norm, ave_ratio_4_norm) %>%
    mutate(conc = 4) %>%
    rename(ratio_A_norm = ratio_4A_norm) %>%
    rename(ratio_B_norm = ratio_4B_norm) %>%
    rename(ave_ratio_norm = ave_ratio_4_norm) %>%
    rename(ave_barcode = ave_barcode_4)
  df_8 <- df %>%
    mutate(ave_barcode_8 = (barcodes_RNA_8A + barcodes_RNA_8B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_8, 
           ratio_8A_norm, ratio_8B_norm, ave_ratio_8_norm) %>%
    mutate(conc = 8) %>%
    rename(ratio_A_norm = ratio_8A_norm) %>%
    rename(ratio_B_norm = ratio_8B_norm) %>%
    rename(ave_ratio_norm = ave_ratio_8_norm) %>%
    rename(ave_barcode = ave_barcode_8)
  df_16 <- df %>%
    mutate(ave_barcode_16 = (barcodes_RNA_16A + barcodes_RNA_16B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_16, 
           ratio_16A_norm, ratio_16B_norm, ave_ratio_16_norm) %>%
    mutate(conc = 16) %>%
    rename(ratio_A_norm = ratio_16A_norm) %>%
    rename(ratio_B_norm = ratio_16B_norm) %>%
    rename(ave_ratio_norm = ave_ratio_16_norm) %>%
    rename(ave_barcode = ave_barcode_16)
  df_25 <- df %>%
    mutate(ave_barcode_25 = (barcodes_RNA_25A + barcodes_RNA_25B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_25, 
           ratio_25A_norm, ratio_25B_norm, ave_ratio_25_norm) %>%
    mutate(conc = 25) %>%
    rename(ratio_A_norm = ratio_25A_norm) %>%
    rename(ratio_B_norm = ratio_25B_norm) %>%
    rename(ave_ratio_norm = ave_ratio_25_norm) %>%
    rename(ave_barcode = ave_barcode_25)
  df_32 <- df %>%
    mutate(ave_barcode_32 = (barcodes_RNA_32A + barcodes_RNA_32B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_32, 
           ratio_32A_norm, ratio_32B_norm, ave_ratio_32_norm) %>%
    mutate(conc = 32) %>%
    rename(ratio_A_norm = ratio_32A_norm) %>%
    rename(ratio_B_norm = ratio_32B_norm) %>%
    rename(ave_ratio_norm = ave_ratio_32_norm) %>%
    rename(ave_barcode = ave_barcode_32)
  df_64 <- df %>%
    mutate(ave_barcode_64 = (barcodes_RNA_64A + barcodes_RNA_64B)/2) %>%
    select(subpool, name, most_common, background, induction, ave_barcode_64, 
           ratio_64A_norm, ratio_64B_norm, ave_ratio_64_norm) %>%
    mutate(conc = 64) %>%
    rename(ratio_A_norm = ratio_64A_norm) %>%
    rename(ratio_B_norm = ratio_64B_norm) %>%
    rename(ave_ratio_norm = ave_ratio_64_norm) %>%
    rename(ave_barcode = ave_barcode_64)
  df_0_64 <- rbind(df_0, df_1, df_2, df_4, df_8, df_16, df_25, df_32, df_64)
  return(df_0_64)
}

#Episomal MPRA expression across tested forskolin concentrations

epi_back_norm_conc <- epi_back_norm_pc_spGl4 %>%
  var_conc_exp()

write.table(epi_back_norm_conc, "rep_0_64_med_norm.txt", 
            sep = '\t', row.names = FALSE)



