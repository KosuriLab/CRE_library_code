#Written for the processing of the first episomal MPRA dataset. MPRA performed 
#at 0 and 25 µM forskolin with library designs corresponding to the CRE 
#Distance, CRE Spacing and Distance, CRE Number and Affinity library, and 
#another library based on the CRE Spacing and Distance library but with 8 
#combinations of site affinities across the two CREs and less Distances tested

#bc_JD02: RNA 25 µM forskolin rep. 1
#bc_JD03: RNA 25 µM forskolin rep. 2
#bc_JD04: RNA 0 µM forskolin rep. 1
#bc_JD05: RNA 0 µM forskolin rep. 2
#bc_JD11: DNA pool

library(tidyverse)

#Load index files---------------------------------------------------------------

bc_R25A <- read_tsv('BCreads_txts/bc_JD02.txt')
bc_R25B <- read_tsv('BCreads_txts/bc_JD03.txt')
bc_R0A <- read_tsv('BCreads_txts/bc_JD04.txt')
bc_R0B <- read_tsv('BCreads_txts/bc_JD05.txt')
bc_DNA <- read_tsv('BCreads_txts/bc_JD11.txt')


#Load barcode mapping table, remember sequences are rcomp

barcode_map <- read_tsv('../BCMap/uniqueSP2345.txt', 
                        col_names = c(
                          'fluff', 'barcode', 'name', 'most_common'), 
                        skip = 1) %>%
  select(-fluff) %>%
  mutate(subpool = ifelse(startsWith(name, 'subpool'), 
                          substr(name, 1, 8), 
                          'control')) 


#Join reads---------------------------------------------------------------------

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

bc_join_R25A <- bc_map_join_bc(barcode_map, bc_R25A)
bc_join_R25B <- bc_map_join_bc(barcode_map, bc_R25B)
bc_join_R0A <- bc_map_join_bc(barcode_map, bc_R0A)
bc_join_R0B <- bc_map_join_bc(barcode_map, bc_R0B)
bc_join_DNA <- bc_map_join_bc(barcode_map, bc_DNA)

#Median BC expression-----------------------------------------------------------

#Retain BCs with > 5 reads in DNA sample, left-join RNA to DNA BCs. Take ratio 
#of RNA/DNA normalizedreads per million

bc_dna_join_rna <- function(df1, df2) {
  filter_DNA <- filter(df1, num_reads > 5)
  DNA_RNA_join <- left_join(filter_DNA, df2,
                            by = c("barcode", "name", "subpool", "most_common"), 
                            suffix = c('_DNA', '_RNA')) %>%
    mutate(ratio = normalized_RNA/normalized_DNA)
  print('processed dfs in order of (DNA, RNA) in bc_dna_join_rna(df1, df2)')
  return(DNA_RNA_join)
}

RNA_DNA_bc_R0A <- bc_dna_join_rna(bc_join_DNA, bc_join_R0A)
RNA_DNA_bc_R0B <- bc_dna_join_rna(bc_join_DNA, bc_join_R0B)
RNA_DNA_bc_R25A <- bc_dna_join_rna(bc_join_DNA, bc_join_R25A)
RNA_DNA_bc_R25B <- bc_dna_join_rna(bc_join_DNA, bc_join_R25B)

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

med_RNA_DNA_R0A <- ratio_bc_med_var(RNA_DNA_bc_R0A)
med_RNA_DNA_R0B <- ratio_bc_med_var(RNA_DNA_bc_R0B)
med_RNA_DNA_R25A <- ratio_bc_med_var(RNA_DNA_bc_R25A)
med_RNA_DNA_R25B <- ratio_bc_med_var(RNA_DNA_bc_R25B)


#Combine all samples

med_var_rep <- function(df0A, df0B, df25A, df25B) {
  join_0 <- inner_join(df0A, df0B, 
                       by = c("name", "subpool", "most_common", 'barcodes_DNA'), 
                       suffix = c("_0A", "_0B"))
  join_25 <- inner_join(df25A, df25B, 
                        by = c("name", "subpool", "most_common", 
                               'barcodes_DNA'), 
                        suffix = c("_25A", "_25B"))
  join_0_25 <- inner_join(join_0, join_25, 
                          by = c("name", "subpool", "most_common", 
                                 'barcodes_DNA'))
  print('processed dfs in order: df0A, df0B, df25A, df25B')
  return(join_0_25)
}

med_rep_0_25 <- med_var_rep(med_RNA_DNA_R0A, med_RNA_DNA_R0B, 
                            med_RNA_DNA_R25A, med_RNA_DNA_R25B)

#Background-normalize expression------------------------------------------------

back_norm <- function(df1) {
  gsub_df1 <- df1 %>%
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
  backgrounds <- gsub_df1 %>%
    filter(startsWith(name,
                      'subpool5_no_site_no_site_no_site_no_site_no_site_no_site')
    ) %>%
    select(background, med_ratio_0A, med_ratio_0B, 
           med_ratio_25A, med_ratio_25B) %>%
    rename(med_ratio_0A_back = med_ratio_0A) %>%
    rename(med_ratio_0B_back = med_ratio_0B) %>%
    rename(med_ratio_25A_back = med_ratio_25A) %>%
    rename(med_ratio_25B_back = med_ratio_25B)
  back_join_norm <- left_join(gsub_df1, backgrounds, by = 'background') %>%
    mutate(med_ratio_0A_norm = med_ratio_0A/med_ratio_0A_back) %>%
    mutate(med_ratio_0B_norm = med_ratio_0B/med_ratio_0B_back) %>%
    mutate(med_ratio_25A_norm = med_ratio_25A/med_ratio_25A_back) %>%
    mutate(med_ratio_25B_norm = med_ratio_25B/med_ratio_25B_back) %>%
    mutate(ave_ratio_0_norm = (med_ratio_0A_norm + med_ratio_0B_norm)/2) %>%
    mutate(ave_ratio_25_norm = (med_ratio_25A_norm + med_ratio_25B_norm)/2) %>%
    mutate(induction = ave_ratio_25_norm/ave_ratio_0_norm)
  return(back_join_norm)
}

med_rep_0_25_back_norm <- back_norm(med_rep_0_25)

var_conc_exp_med <- function(df) {
  df_0 <- df %>%
    mutate(ave_barcode_0 = (barcodes_RNA_0A + barcodes_RNA_0B)/2) %>%
    select(subpool, name, most_common, background, ave_barcode_0, 
           ave_ratio_0_norm) %>%
    mutate(conc = 0) %>%
    rename(ave_ratio_norm = ave_ratio_0_norm) %>%
    rename(ave_barcode = ave_barcode_0)
  df_25 <- df %>%
    mutate(ave_barcode_25 = (barcodes_RNA_25A + barcodes_RNA_25B)/2) %>%
    select(subpool, name, most_common, background, ave_barcode_25, 
           ave_ratio_25_norm) %>%
    mutate(conc = 25) %>%
    rename(ave_ratio_norm = ave_ratio_25_norm) %>%
    rename(ave_barcode = ave_barcode_25)
  df_0_25 <- rbind(df_0, df_25)
  return(df_0_25)
}

med_rep_0_25_back_norm_conc <- var_conc_exp_med(med_rep_0_25_back_norm)

write.table(med_rep_0_25_back_norm_conc, "rep_0_25_med_norm.txt", 
            sep = '\t', row.names = FALSE)

