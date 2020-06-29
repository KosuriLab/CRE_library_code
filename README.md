# CRE_library_code
Combined and finalized code and data required for figure generation used in the paper:

"Dissection of c-AMP Response Element Architecture by Using Genomic and Episomal Massively Parallel Reporter Assays"

https://doi.org/10.1016/j.cels.2020.05.011

Included in this project is code required to:

-generate CRE variants (cre_lib_generation and cre_followuplib_generation)

-assign barcodes to CRE variants (BCMap)

-analyze CRE MPRAs (20170320_epilib_analysis_0_25, 20170631_epilib_analysis_0_64, 20170921_epi_lib_analysis, 20171129_genlib_analysis, and 20190712_epilib_analysis) of which contains code for analyzing dfs generated from other experiments (plate_reader)

All files too large to be uploaded here are uploaded onto the Gene Expression Omnibus: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137922.

The main MPRA analysis script and the one used for all figure generation is 20170921_epi_lib_analysis. All other MPRA processing scripts need to be run first since their outputs are imported into 20170921_epi_lib_analysis for analysis and figure generation.

For processing all MPRA data, you must generate BCreads_txts (folder containing txt files per index of reverse compliment BCs and their summed reads) per each MPRA from sequencing data using BC mapping. Per each MPRA folder, organize the DNA/RNA barcode reads from GEO into each folder per MPRA as indicated by the Github organization (indicated by MPRA date) and place in a sub-folder of BCreads_txts within each MPRA folder.

