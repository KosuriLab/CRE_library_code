#Written for the analysis of genomic library expression at 8 µM. 2 biological
#replicates, each with 2 DNA technical replicates

#D##_BC: DNA biological replicate # technical replicate #
#R#8_BC: RNA biological replicate # at 8 µM

library(tidyverse)

#Load index and bcmap files-----------------------------------------------------

bc_DNA_1_1 <- read_tsv('BCreads_txts/D11_BC.txt')
bc_DNA_1_2 <- read_tsv('BCreads_txts/D12_BC.txt')
bc_DNA_2_1 <- read_tsv('BCreads_txts/D21_BC.txt')
bc_DNA_2_2 <- read_tsv('BCreads_txts/D22_BC.txt')
bc_RNA_1 <- read_tsv('BCreads_txts/R18_BC.txt')
bc_RNA_2 <- read_tsv('BCreads_txts/R28_BC.txt')

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
#mapping and replacing na with 0 reads.

bc_map_join_bc <- function(df1, df2) {
  df2 <- df2 %>%
    mutate(norm = as.numeric((num_reads * 1000000) / (sum(num_reads))))
  keep_bc <- left_join(df1, df2, by = 'barcode') %>%
    mutate(norm = if_else(is.na(norm), 
                                0, 
                                norm)) %>%
    mutate(num_reads = if_else(is.na(num_reads), 
                               0, 
                               num_reads))
  return(keep_bc)
}

bc_join_DNA_1_1 <- bc_map_join_bc(SP3_SP5_map, bc_DNA_1_1)
bc_join_DNA_1_2 <- bc_map_join_bc(SP3_SP5_map, bc_DNA_1_2)
bc_join_DNA_2_1 <- bc_map_join_bc(SP3_SP5_map, bc_DNA_2_1)
bc_join_DNA_2_2 <- bc_map_join_bc(SP3_SP5_map, bc_DNA_2_2)
bc_join_RNA_1 <- bc_map_join_bc(SP3_SP5_map, bc_RNA_1)
bc_join_RNA_2 <- bc_map_join_bc(SP3_SP5_map, bc_RNA_2)


#Median BC expression-----------------------------------------------------------

#Retain BCs with > 6 reads in both DNA samples, determine average normalized DNA
#BC reads between replicates. Left_join RNA BCs to DNA BCs and determine ratio
#of expression of each BC per biological replicate

ave_dna_join_rna_rep <- function(df1, df2, df3) {
  filter_reads_1 <- filter(df1, num_reads > 6)
  filter_reads_2 <- filter(df2, num_reads > 6)
  DNA_join <- inner_join(filter_reads_1, filter_reads_2, 
                         by = c("barcode", "name", "subpool", "most_common"), 
                         suffix = c("_DNA_tr1", "_DNA_tr2")) %>%
    mutate(ave_norm_DNA = (norm_DNA_tr1 + norm_DNA_tr2)/2)
  DNA_RNA_join <- left_join(DNA_join, df3,
                             by = c("barcode", "name", "subpool", 
                                    "most_common")) %>%
    rename(num_reads_RNA = num_reads) %>%
    rename(norm_RNA = norm) %>%
    mutate(ratio = norm_RNA/ave_norm_DNA)
  print('processed dfs in order of (DNA tr1, DNA tr2, RNA) in 
        bc_dna_biol_rep(df1, df2, df3)')
  return(DNA_RNA_join)
}

bc_ave_DNA_RNA_1 <- ave_dna_join_rna_rep(bc_join_DNA_1_1, bc_join_DNA_1_2, 
                                         bc_join_RNA_1)
bc_ave_DNA_RNA_2 <- ave_dna_join_rna_rep(bc_join_DNA_2_1, bc_join_DNA_2_2, 
                                         bc_join_RNA_2)


#Count barcodes per variant per DNA and RNA, set minimum of 8 BC's per variant 
#in DNA, take median RNA/DNA per variant, then per variant determine 
#the median absolute deviation of all barcode ratios. Then filter out variants 
#with 0 median expression

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

med_ratio_1 <- ratio_bc_med_var(bc_ave_DNA_RNA_1)
med_ratio_2 <- ratio_bc_med_var(bc_ave_DNA_RNA_2)


#combine biological replicates

rep_1_2 <- inner_join(med_ratio_1, med_ratio_2,
             by = c("name", "subpool", "most_common"),
             suffix = c('_br1', '_br2'))

#After combining rename backgrounds to simplified names, make background column. 
#Separate out background values in each dataset and left join to original 
#dataset. Normalize expression of each variant to its background in that 
#biological replicate. Determine average expression and average 
#background-normalized expression across biological replicates.

back_norm <- function(df1) {
  gsub_1_2 <- df1 %>%
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
                                nchar(background)))
  backgrounds <- gsub_1_2 %>%
    filter(startsWith(name, 'subpool5_no_site_no_site_no_site_no_site_no_site_no_site')) %>%
    select(background, med_ratio_br1, med_ratio_br2) %>%
    rename(med_ratio_br1_back = med_ratio_br1) %>%
    rename(med_ratio_br2_back = med_ratio_br2) 
  back_join_norm <- left_join(gsub_1_2, backgrounds, by = 'background') %>%
    mutate(ave_med_ratio = (med_ratio_br1 + med_ratio_br2)/2) %>%
    mutate(med_ratio_br1_norm = med_ratio_br1/med_ratio_br1_back) %>%
    mutate(med_ratio_br2_norm = med_ratio_br2/med_ratio_br2_back) %>%
    mutate(ave_med_ratio_norm = (med_ratio_br1_norm + med_ratio_br2_norm)/2)
}

rep_1_2_back_norm <- rep_1_2 %>%
  filter(name == 'pGL4.29 Promega 1-63 + 1-87') %>%
  mutate(name = str_c(name, '_scramble pGL4.29 Promega 1-63 + 1-87')) %>%
  mutate(subpool = 'subpool3') %>%
  rbind(rep_1_2) %>%
  back_norm()

output_int <- rep_1_2_back_norm %>%
  write.table(
    "rep_1_2.txt", 
    sep = '\t', row.names = FALSE)

