## Lampe F&F Study
## Filtering, Normalization & DE Analysis.
## Last updated December 2, 2024 by Destiny Mullens

library(dplyr)
library(tibble)

## Import counts
counts <- read.csv("../Working_Data/Counts/LFF-raw-counts-all-combined.csv") %>% 
  column_to_rownames(var = "Ensembl.ID") %>%              ## Move Ensembl ID's to row name
  dplyr::select(.,-Gene.Name) %>%                         ## Remove column with Gene names
  filter(.,rowSums(.)!=0) %>%                             ## Filter out genes with no counts in any sample
  dplyr::select(-contains(c('.S1','.S2','.S4','.S5','NA')))    ## Remove all time points except endpoints (S3 & S6)


## Import metadata
meta <- readRDS("../Working_Data/LFF-Full_Metadata.rds") %>% 
  filter(Sample.Name %in% names(counts)) %>%      ## Only import metadata for samples in counts df
  dplyr::select(.,c(Sample.Name,Participant,Group))           ## Select only the Sample.Name & Group columns

library(edgeR)
## Convert counts to cpm for filtering
counts.cpm <- counts %>% cpm() %>% as.data.frame()

### Set cpm limit for required number of reads to keep a gene (over 20 cpm)
filter.limit <- 20
## Set number of samples in each group that must be at/over filter limit to keep gene
sample.num <- 19

## Count number of times a gene is over >= the filter limit
features.count <- left_join(rownames_to_column(as.data.frame(t(counts.cpm))), meta,by=c("rowname"="Sample.Name")) %>%
  group_by(Group) %>%                                                             ## Group by Group in metadata
  dplyr::summarise(across(everything(), ~ sum(.x >= filter.limit))) %>%           ## Summarize and count number of samples with count >= filter limit
  select_if(.,is.numeric) %>%                                                     ## Select numeric columns
  t() %>% as.data.frame()                                                         ## transpose and make a df

## Select Ensembl.IDs for any gene that that was over the filter limit at least the "sample.num" amount of times
features.to.keep <- rownames(features.count)[features.count$V1 > sample.num] %>%  ## Select row names if more than sample num have gene for group 1
  append(rownames(features.count)[features.count$V2 > sample.num]) %>%            ## Select row name if more than sample num have gene for group 2
  unique()                                                                        ## Select unique Ensembl.ID's only

## Keep genes meeting above criteria for analysis
filtered.counts <- counts[features.to.keep,] %>%                                  ## Filter counts by keep only genes selected above
  na.omit() %>% as.data.frame()                                                   ## Omit any NA's and make a dataframe

## Save raw, filtered counts
write.csv(filtered.counts,"../Working_Data/Counts/Raw-filtered-counts-endpoints-73percent-20cpm.csv")

###################### DE Analysis ######################

## Create DGE list with filtered counts and groups in meta data
y <- DGEList(counts=filtered.counts, group = meta$Group)

## Create another group in DGEList for Participant number
y$samples$participant <- meta$Participant

## Check samples to make sure everything looks okay
y$samples

## Create factored vector for participant number
participant <- factor(meta$Participant)

## Create factored vector for treatment group
treatment <- factor(meta$Group,levels = c("B","A"))

## Create design model for DE analysis 
design <- model.matrix(~participant+treatment) %>% as.data.frame()
rownames(design) <- colnames(y)

## Check samples to make sure everything looks okay
y$samples

## Normalize counts
y <- calcNormFactors(y, method = "TMM")

## Extract normalized counts from DGEList
TMM.norm <- cpm(y,log = FALSE) %>% round() %>% as.data.frame()

## Save normalized counts for other analyses
write.csv(TMM.norm,"../Working_Data/Counts/Norm-Counts-endpoints-73percent-20cpm.csv",row.names = TRUE)

## Estimate dispersion
y <- estimateDisp(y, design, robust = TRUE)

## Fit model
fitq <- glmQLFit(y, design)

## Test for treatment effect
qlf <- glmQLFTest(fitq)
## View top results
topTags(qlf,n=50)

out1q <- topTags(qlf,n=nrow(filtered.counts),sort.by="PValue")$table
out2q <- out1q[order(out1q$FDR),] # Sort by ascending FDR
out2q$FC <- 2^(out1q$logFC)
genes <- counts <- read.csv("../Working_Data/Counts/LFF-raw-counts-all-combined.csv")[,c(1:2)]
out2q <- left_join(rownames_to_column(out2q,var = "Ensembl.ID"),genes)
out2q <- out2q[,c("Ensembl.ID","Gene.Name","FC","logFC","logCPM","F","PValue","FDR")]
write.csv(out2q,"LFF-paired_AvB-endpoints-QLF.csv")