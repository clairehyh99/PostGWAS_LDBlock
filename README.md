# PostGWAS_LDBlock
##Pre-requests
1. Do GWAS using GEMMA - otherwise need to adjust some inputs (p_value col no. etc...)
2. Get LD data using plink (plink --blocks...), get .det
##Ready for processing
1. Extract Significant SNPs, loop over all traits
2. Read in LD block files (.det)
3. Format .det, and using BP1 (start physical pos) as BlockID
4. Find LD blocks containing significant SNPs
5. Write in the corresponding traits
