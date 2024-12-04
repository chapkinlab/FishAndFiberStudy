## LDA Frequency calculations for Lampe F&F
## Last updated by Destiny Mullens on December 2, 2024
library(ggplot2)
library(dplyr)

## Importing only gene names & bresub
midpoints <- read.csv("Results/LFF-midpoints-3_feature_annotation-FULL.csv")[,5:8]
endpoints <- read.csv("Results/LFF-endpoints-3_feature_annotation-FULL.csv")[,5:8]

## Order by bresub and tak the top 100
midpoints <- midpoints %>% slice_min(., order_by = bresub, n=100)
endpoints <- endpoints %>% slice_min(., order_by = bresub, n=100)

## Run both functions first
probability_feature <- function(num.features, expected.num) {
  prob <- dpois(num.features,expected.num)
}

find.high.freq.num <- function(num.features,observations,criteria){
  expected.num <- 3/num.features*observations
  k <- round(expected.num)
  while (probability_feature(k, expected.num) >= 0.01) {
    k <- k + 1
  }
  print(k)  
}

## Input the number of features, the total number of triplicates & desired probability
## Output will be the minimum number of times a gene can appear to reach the desired probability 
midpoint.freq.num <- find.high.freq.num(1020,100,0.01)
endpoint.freq.num <- find.high.freq.num(1020,100,0.01)

midpoints <- data.frame(table(unlist(midpoints[,2:4]))) %>% setNames(.,c("Gene","Midpoint"))
endpoints <- data.frame(table(unlist(endpoints[,2:4]))) %>% setNames(.,c("Gene","Endpoint"))

## Create list of all gene frequencies 
all.frequencies <- full_join(midpoints,endpoints)

## Add zero value for NA's
all.frequencies[is.na(all.frequencies)] <- 0
## Save gene frequencies
write.csv(all.frequencies,"LDA_Frequencies-top100.csv",row.names = FALSE)

library(reshape)
graph.data <- melt(all.frequencies, id.vars = c("Gene"), measure.vars = c("Midpoint","Endpoint")) %>% setNames(., c("Gene","Dataset","Frequency"))
graph.data$Dataset <- as.factor(graph.data$Dataset)

## Remove genes from the list with lower frequencies
freq_min <- 2
all.reduced <- all.frequencies[which(all.frequencies$Midpoint > freq_min | all.frequencies$Endpoint > freq_min),]
write.csv(all.reduced,"LDA_Frequencies-high-freq-genes-only.csv")

graph.data <- melt(all.reduced, id.vars = c("Gene"), measure.vars = c("Midpoint","Endpoint")) %>% setNames(., c("Gene","Dataset","Frequency"))
graph.data$Gene <- as.character(graph.data$Gene)
## Remove genes that are 2 genes not of interest and replace name of another before plotting
graph.data <- filter(graph.data,graph.data$Gene!="AC009967.1")
graph.data <- filter(graph.data,graph.data$Gene!="AL158163.1")
graph.data$Gene[graph.data$Gene=='AC084757.3'] <- c('FBN1-DT')
graph.data$Dataset <- as.factor(graph.data$Dataset)

library(ggplot2)
## Manuscript color
ggplot(data=graph.data, aes(x=Gene, y=Frequency, fill=Dataset)) + geom_bar(stat="identity", position = "dodge") +
  scale_fill_manual(values=c("#56B4E9","#CC79A7"), labels=c("Middle Time Point", "End Time Point")) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 90,size = 12,hjust=0.95,vjust=.5), panel.grid.major = element_blank(), 
        legend.position = "bottom", 
        legend.text=element_text(size=10), 
        legend.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=10, face="bold", colour = "black"), 
        axis.title.y = element_text(size=10, face="bold", colour = "black"), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "gray90"),
        axis.line = element_line(colour = "black"), 
        legend.background = element_rect(linewidth=0.5, linetype="solid",colour ="black")) + 
  labs(x="Gene Name", y="Gene Frequency")
ggsave("LDA-High-Frequency-Genes-Color-300dpi.tiff",plot = last_plot(),width = 8.5, height = 5, units = c("in"),dpi=300)

## Manuscript BW
ggplot(data=graph.data, aes(x=Gene, y=Frequency, fill=Dataset)) + geom_bar(stat="identity", position = "dodge") +
  scale_fill_manual(values=c("gray70","gray20"), labels=c("Middle Time Point", "End Time Point")) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 90,size = 12,hjust=0.95,vjust=.5), panel.grid.major = element_blank(), 
        legend.position = "bottom", 
        legend.text=element_text(size=10), 
        legend.title = element_text(size=10, face = "bold"),
        axis.title.x = element_text(size=10, face="bold", colour = "black"), 
        axis.title.y = element_text(size=10, face="bold", colour = "black"), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"), 
        legend.background = element_rect(linewidth=0.5, linetype="solid",colour ="black")) +
  labs(x="Gene Name", y="Gene Frequency")
ggsave("LDA-High-Frequency-Genes-BW-300dpi.tiff",plot = last_plot(),width = 8.5, height = 5, units = c("in"),dpi=300)
