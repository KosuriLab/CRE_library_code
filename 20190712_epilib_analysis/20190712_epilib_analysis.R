#Written for the analysis of the follow-up CRE library in the episomal MPRA


#D3_BC: DNA
#R0A_BC and R0B_BC: RNA at 0 µM Forsk replicate A and B
#R4A_BC and R4B_BC: RNA at 4 µM Forsk replicate A and B


#Establish workspace------------------------------------------------------------

options(stringsAsFactors = F)

#import necessary libraries

library(tidyverse)
library(viridis)
library(cowplot)

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

bc_DNA <- read_tsv('BCreads_txts/D3_BC.txt')
bc_R0A <- read_tsv('BCreads_txts/R0A_BC.txt')
bc_R0B <- read_tsv('BCreads_txts/R0B_BC.txt')
bc_R4A <- read_tsv('BCreads_txts/R4A_BC.txt')
bc_R4B <- read_tsv('BCreads_txts/R4B_BC.txt')

#Load barcode mapping table, sequences (most_common) are rcomp due to sequencing
#format

Follow_up_map <- read_tsv('../BCMap/cre_followup_barcode_statistics.txt') %>%
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

#The twosite_distal library consists of a fixed proximal CRE and moving distal 
#CRE. The offset of the proximal CRE from the 3' end of the background is 
#indicated in 0, 5, 10, and 15 bp increments. The spacing between the CREs is
#also indicated. The distance of the distal CRE as it is placed further from the
#3' end is indicated as the distance from the start of this site to the minimal 
#promoter. Distance stops at 60 bp. Only the 3 original backgrounds were used here.

twosite_distal <- med_rep_0_4_A_B %>%
  filter(grepl('twositedistdistal', name)) %>%
  separate(name, 
           into = c("subpool", "fluff1", "offset", "fluff2", "spacing", 
                    "fluff3", "distaldist", "fluff4", "background"),
           sep = "_", convert = TRUE) %>%
  select(-subpool, -fluff1, -fluff2, -fluff3, -fluff4) %>%
  mutate(distaldist = distaldist + 64)

#The twosite_proximal library consists of a fixed distal CRE and moving proximal 
#CRE. The distal CRE is placed 60 bp from the 3' end of the background - offset.
#The offset of the distal CRE from the 5' end of the background is 
#indicated in 0, 5, 10, and 15 bp increments. The spacing between the CREs is
#also indicated. The distance of the proximal CRE as it is placed closer to the
#3' end is indicated as the distance from the start of this site to the minimal 
#promoter. Only the 3 original backgrounds were used here.

twosite_proximal <- med_rep_0_4_A_B %>%
  filter(grepl('twositedistproximal', name)) %>%
  separate(name, 
           into = c("subpool", "fluff1", "offset", "fluff2", "spacing", 
                    "fluff3", "proximaldist", "fluff4", "background"),
           sep = "_", convert = TRUE) %>%
  select(-subpool, -fluff1, -fluff2, -fluff3, -fluff4) %>%
  mutate(proximaldist = proximaldist + 64)


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
twosite_distal_norm <- back_norm(oldback, twosite_distal)
twosite_proximal_norm <- back_norm(oldback, twosite_proximal)

#Export dfs---------------------------------------------------------------------

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

twosite_distal_norm %>%
  write.table(
    "twosite_distal_norm.txt", 
    sep = '\t', row.names = FALSE)

twosite_proximal_norm %>%
  write.table(
    "twosite_proximal_norm.txt", 
    sep = '\t', row.names = FALSE)

#plots--------------------------------------------------------------------------

p_rep_0 <- med_rep_0_4_A_B %>%
  ggplot(aes(med_ratio_0A, med_ratio_0B)) +
  geom_point(alpha = 0.1, size = 0.75) +
  geom_point(data = filter(med_rep_0_4_A_B, 
                           grepl(
                             'sixsite_nosite_nosite_nosite_nosite_nosite_nosite',
                             name)), 
             fill = 'orange', shape = 21, size = 1.5) + 
  geom_point(data = filter(med_rep_0_4_A_B, 
                           name == 'control_pGL4.29 Promega 1-63 + 1-87'), 
             fill = 'red', shape = 21, size = 1.75) +
  annotation_logticks(scaled = TRUE) +
  xlab("Expression (a.u.) replicate 1") +
  ylab("Expression (a.u.) replicate 2") +
  scale_x_log10(limits = c(0.02, 20), breaks = c(0.1, 1, 10)) + 
  scale_y_log10(limits = c(0.02, 20), breaks = c(0.1, 1, 10)) +
  theme(strip.background = element_rect(colour="black", fill="white"),
        axis.line.y = element_line(), panel.spacing.x=unit(1, "lines")) +
  figurefont_theme

p_rep_4 <- med_rep_0_4_A_B %>%
  ggplot(aes(med_ratio_4A, med_ratio_4B)) +
  geom_point(alpha = 0.1, size = 0.75) +
  geom_point(data = filter(med_rep_0_4_A_B, 
                           grepl(
                             'sixsite_nosite_nosite_nosite_nosite_nosite_nosite',
                             name)), 
             fill = 'orange', shape = 21, size = 1.5) + 
  geom_point(data = filter(med_rep_0_4_A_B, 
                           name == 'control_pGL4.29 Promega 1-63 + 1-87'), 
             fill = 'red', shape = 21, size = 1.75) +
  annotation_logticks(scaled = TRUE) +
  xlab("Expression (a.u.) replicate 1") +
  ylab("Expression (a.u.) replicate 2") +
  scale_x_log10(limits = c(0.02, 20), breaks = c(0.1, 1, 10)) + 
  scale_y_log10(limits = c(0.02, 20), breaks = c(0.1, 1, 10)) +
  theme(strip.background = element_rect(colour="black", fill="white"),
        axis.line.y = element_line(), panel.spacing.x=unit(1, "lines")) +
  figurefont_theme

log10_med_rep_0_4_A_B <- var_log10(med_rep_0_4_A_B)

pearsons_epi <- tibble(
  sample = c('0 µM', '4 µM'),
  pearsons = c(cor(log10_med_rep_0_4_A_B$med_ratio_0A, 
                   log10_med_rep_0_4_A_B$med_ratio_0B, 
                   use = "pairwise.complete.obs", method = "pearson"),
               cor(log10_med_rep_0_4_A_B$med_ratio_4A, 
                   log10_med_rep_0_4_A_B$med_ratio_4B, 
                   use = "pairwise.complete.obs", method = "pearson")))

