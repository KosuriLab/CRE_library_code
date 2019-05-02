# CRE_library_code
Combined and finalized code and data required for figure generation used in the paper:

"Multiplexed dissection of a model human transcription factor binding site architecture"

This paper is currently not peer-reviewed and may change throughout this process.

Included in this project is code required to:

-generate CRE variants (cre_lib_generation)

-assign barcodes to CRE variants (BCMap)

-analyze CRE MPRAs (20170320_epilib_analysis_0_25, 20170631_epilib_analysis_0_64, 20170921_epi_lib_analysis, and 20171129_genlib_analysis) of which contains code for analyzing dfs generated from other experiments (plate_reader)

The outputs of cre_lib_generation and BCMap are already uploaded here so that MPRA figures can be generated.

The main MPRA analysis script and the one used for all figure generation is 20170921_epi_lib_analysis. All other MPRA processing scripts need to be run first since their outputs are imported into 20170921_epi_lib_analysis for analysis and figure generation.

