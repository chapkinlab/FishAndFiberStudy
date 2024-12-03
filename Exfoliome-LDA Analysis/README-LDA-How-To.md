# ** Preparing and running LDA code for running 1, 2 and 3 feature classification**

1. RCodes/1_LDA-input-file*.R has code for subsetting data and properly formatting it for LDA.
2. After creating input files, run the LDA code in the LDA-C-code folder to get LDA results.
3. After running LDA code, RCodes/2_LDA-Annotation*.R has the code to annotate the results.
4. RCodes/3_LDA_High_Frequency_Classifiers.R code can use used to find high frequency genes.