## Lampe Fish & Fiber Study
## This code is for combining all of the libraries for each sample.
## Code last updated December 2nd, 2024 by Destiny Mullens

library(dplyr)
## Combine raw data from separate sequencing runs
counts <- read.csv("../../Batch2/Lampe_FF_Batch2-optimized/summary/Lampe_FF_Batch2-optimized-counts.csv")

counts2 <- read.csv("../../Batch3_resequenced-samples/Lampe_FF_Batch3-optimized/summary/Lampe_FF_Batch3-optimized-counts.csv", sep = "\t")[,c(2:9)] %>% 
  setNames(.,c("Ensembl.ID","LFF.53061.S2.lib2","LFF.53061.S3_S145","LFF.53061.S6.lib1","LFF.53061.S6.lib2","LFF.53087.S3.lib1","LFF.53096.S6_S176","LFF.53126.S1.lib1")) %>% 
  as.data.frame()

counts3 <- read.csv("../../Batch4_resequenced-sample/Resequenced_Samples_Jan2024/summary/Resequenced_Samples_Jan2024-counts.csv", sep = "\t")[,c(2,4)] %>% 
  setNames(.,c("Ensembl.ID","LFF.53002.S4.lib1"))

## Removed re sequenced samples from original counts df
counts <- counts %>% 
  select(!any_of(names(counts2[2:8]))) %>% 
  select(.,!any_of(names(counts3[2])))

## Merge together original sequencing runs with new ones
counts.new <- counts %>% 
  left_join(.,counts2) %>% 
  left_join(.,counts3)

## Order samples by name
counts.new <- counts.new[,order(names(counts.new))]

## Fix names so that all samples have the same formatting with .lib# at the end
## Get sample names
sample.list <- names(counts.new) %>% as.data.frame() %>% setNames(.,c("Names"))
## Separate Sample ID to get library information
sample.list$library <- substring(sample.list$Names, 14, 18)
## Replace those with S## with lib1
sample.list$library[grep("^S",sample.list$library)] <- c('lib1')
## Get the portion of sample name prior to .lib
sample.list$Sample.Name <- substring(sample.list$Names, 1, 12)
## Paste together Sample Name & corrected .lib information
sample.list$New.Names <- paste(sample.list$Sample.Name,sample.list$library,sep = '.')
## Fix these two so they don't have a '.' at the end
sample.list$New.Names[1:2] <- c("Ensembl.ID","Gene.Name")
## Repace names
names(counts.new) <- sample.list$New.Names

write.csv(counts.new,"../Working_Data/Counts/LFF-uncombined-raw-counts-all-batches-.csv",row.names = FALSE)
###########################################################
## Combine libraries for each sample
counts <- read.csv("../Working_Data/Counts/LFF-uncombined-raw-counts-all-batches.csv",sep = ',')
## Get Gene Name & Ensembl.ID info for new data frame for combined libaries
new.counts <- data.frame(counts[,1:2])

## Create list of samples
samples <- unique(substring(names(counts), 1, 12))

for (i in 3:length(samples)) {
  # Get column names that start with the current sample name
  sample_cols <- grep(paste0('^', samples[i]), names(counts))
  if(length(sample_cols)=='1'){
    ## If there is only one library, that library is added to new data frame
    sample.vector <- counts %>% select(.,all_of(sample_cols)) %>% setNames(.,c(paste(samples[i])))
    new.counts <- cbind(new.counts,sample.vector)
  } else {
    ## If there are multiple libaries, they are then added together and added to new data frame
    new.counts[[paste0(samples[i])]] <- rowSums(counts[,sample_cols])
  }
}

## Check to make sure the counts match
summary(rowSums(counts[,3:length(counts)])==rowSums(new.counts[,3:length(new.counts)]))

## Remove samples that are dropped for analysis
samples.to.drop <- c('LFF.53029','LFF.53048','LFF.53152','LFF.53096','LFF.53061')

## Remove samples that will be dropped
new.counts <- new.counts %>% dplyr::select(.,-starts_with(samples.to.drop))

## Save data frame with combined libraries
write.csv(new.counts,"../Working_Data/Counts/LFF-raw-counts-all-combined.csv",row.names = FALSE)
saveRDS(new.counts,"../Working_Data/Counts/LFF-raw-counts-all-combined.rds")
