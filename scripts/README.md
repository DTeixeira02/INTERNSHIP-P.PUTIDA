# README.md `scripts`

Scripts currently operate with absolute file paths rather than relative due to inexperience. May be updated in the future

IGR-Gene association finder takes IGR_presence_absence.csv from PIGGY, outputs text file listing every cluster and the genes theyre associated with in format [igr, genome, gene]

Operon predictor uses .gbk files from PROKKA and associates genes by their locus tag on the basis that genes are 1. <50 bases apart and 2. Coded on the same strand

Operon IGR graphing script takes outputs of IGR-gene association finder and Operon predictor to assemble a directed graph with networkx for each genome, edges formed between operons and IGRs based on the IGR-gene associations for genes within the operon.
