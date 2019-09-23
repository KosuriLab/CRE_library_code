#Combines library-BC associations from separate BC mappings into one df,
#removing duplicates

library(tidyverse)

SP24 <- read_tsv("SP2_4/barcode_statistics.txt",
                 col_names = c('barcode', 'num_unique_constructs', 'num_reads',
                               'num_reads_most_common', 'most_common', 
                               'name'))

SP35 <- read_tsv("SP3_5/barcode_statistics.txt", 
                 col_names = c('barcode', 'num_unique_constructs', 'num_reads',
                               'num_reads_most_common', 'most_common', 
                               'name'))

varSP2345 <- rbind(SP24, SP35)

uniqueSP2345 <- varSP2345[!(duplicated(varSP2345$barcode) | duplicated(varSP2345$barcode, 
                                                                       fromLast = TRUE)),]

write.table(uniqueSP2345, "uniqueSP2345_bcdetail.txt", sep = "\t", col.names = TRUE, quote = FALSE)
