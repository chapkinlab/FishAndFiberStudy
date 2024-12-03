## This code is to format the input files for to run LDA for Lampe F&F Study endpoint data
## Last edited by Destiny Mullens December 2, 2024
library(dplyr)
library(tibble)
## Import & Cleanup Normalized counts
counts <- read.csv("../Working_Data/Counts/Norm-Counts-endpoints-73percent-20cpm-12022024.csv") %>% remove_rownames %>% column_to_rownames(var = "X")
## Import gene names for later from the raw counts data (norm counts only have Ensembl ID)
gene.names <- readRDS("../Working_Data/Counts/LFF-raw-counts-all-combined-12022024.rds")[,1:2]

## Function to make data noisy...most likely it needs it and it saves a lot of time to do it first
noisy.data <- function(x){
  set.seed(123)
  noise.df <- as.data.frame(matrix(replicate(nrow(x)*ncol(x),rnorm(1,0,0.00001)), ncol = ncol(x), byrow = TRUE))
  data.w.noise <- x + noise.df
  return(data.w.noise)  
}

## Import & Cleanup Metadata
meta <- readRDS("../Working_Data/LFF-Full_Metadata.rds")[,c(1,9)] %>% ## Importing only the columns I know is needed
  setNames(.,c("Sample.Name","Group.TP")) %>% lapply(., as.factor) %>% as.data.frame() %>% 
  subset(., (Sample.Name %in% names(counts))) %>% arrange(Group.TP) 

LDA.key <- 0:(nrow(counts)-1) %>% 
  as.data.frame() %>% 
  cbind(.,row.names(counts)) %>% 
  setNames(.,c("Gene.Number","Ensembl.ID")) %>% 
  left_join(.,gene.names)

write.csv(LDA.key,"Input_Files/LDA-key-LFF-endpoints-12022024.csv",row.names = FALSE)

gene.count <- paste(nrow(counts)) ## Number of genes for header later
counts <- noisy.data(counts) %>% t() %>% as.data.frame() ## Add noise to data
counts <- left_join(meta,rownames_to_column(counts),by=c("Sample.Name"="rowname")) ## Join with metadata to make sure rows are in correct order
                                                                                   ## All rows from group A3 & B3 should be grouped together
labels <- counts[,1:2]
counts <- counts[,3:length(counts)]

########################
## Write Output files ##
########################
## File path to folder to save all input files
fp <- c("Input_Files")
## File name
filename <- paste(fp,"LFF-endpoints-12022024.txt",sep = "/")
## Adds the numbers at the top of input file indicating the number of samples.
sample.count <- paste(length(which(labels$Group.TP=="A3")),",",length(which(labels$Group.TP=="B3")),sep = "")
header <- paste(sample.count,gene.count, sep = "\n")
#Save header & data frame to file
writeLines(header,filename)
write.table(counts,filename,append=TRUE,quote=FALSE,row.names = FALSE,col.names = FALSE, eol = "\n")