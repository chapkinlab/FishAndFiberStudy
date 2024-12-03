## This code is to annotate and perform delta calculations for LDA results for the LFF midpoint data
## Code last updated by Destiny Mullens December 2, 2024

##Import & Format Gene List. I only use the columns for the gene index and gene name and drop the ensembl ID for this
gene.key <- read.csv("Input_Files/LDA-key-LFF-midpoints.csv")[,c(1,3)]

##Import LDA Results as output by the LDA code. It imports the columns with the gene index and bresub only as written.
feature.1 <- read.csv("Results/LFF-midpoints-1_feature.txt",sep = "\t")[,-c(2:4)]
feature.2 <- read.csv("Results/LFF-midpoints-2_feature.txt",sep = "\t")[,-c(2:4)]
feature.3 <- read.csv("Results/LFF-midpoints-3_feature.txt",sep = "\t")[,-c(2:4)]

## 1 feature annotation
library(dplyr)
library(stringr)
feature.1.annotation <- feature.1 %>% mutate_at("gene.index", str_replace, "\\[", "") %>% mutate_at("gene.index", str_replace, "\\]", "") %>% 
  mutate_all(.,as.numeric) %>% setNames(.,c("Gene.Number","bresub")) %>% left_join(.,gene.key)
write.csv(feature.1.annotation, file="Results/LFF-midpoints-1_feature_annotation.csv",row.names=FALSE)

library(tidyr)
## 2 feature annotation
feature.2.annotation <- feature.2 %>% mutate_at("gene.index", str_replace, "\\[", "") %>% mutate_at("gene.index", str_replace, "\\]", "") %>% 
  separate(gene.index,c("Gene.A.Number","Gene.B.Number")) %>% mutate_all(.,as.numeric) %>% 
  left_join(.,feature.1.annotation, by=c("Gene.A.Number"="Gene.Number")) %>% left_join(.,feature.1.annotation, by=c("Gene.B.Number"="Gene.Number")) %>% 
  replace(is.na(.),max(feature.1.annotation$bresub)) %>% mutate(Delta.Gene.A=bresub.y-bresub.x) %>% mutate(Delta.Gene.B=bresub-bresub.x) %>% 
  mutate(Smallest.Delta=pmin(Delta.Gene.A,Delta.Gene.B)) %>% select(.,-c(bresub.y,bresub)) %>% 
  setNames(.,c("Gene.A.Number","Gene.B.Number","bresub","Gene.A.Name","Gene.B.Name","Delta-Gene.A","Delta-Gene.B","Min.Delta")) 
write.csv(feature.2.annotation, file="Results/LFF-midpoints-2_feature_annotation.csv",row.names=FALSE)

## Create Pair-Index on 2 feature annotation helping annotate 3-feature delta calculations
feature.pair.index <- feature.2.annotation %>% mutate(pair.index=str_c(Gene.A.Number,'-',Gene.B.Number)) %>% select(.,c("pair.index","bresub"))

## 3 feature annotation
feature.3.annotation <- feature.3 %>% 
  mutate_at("gene.index", str_replace, "\\[", "") %>% mutate_at("gene.index", str_replace, "\\]", "") %>% 
  separate(gene.index,c("Gene.A.Number","Gene.B.Number","Gene.C.Number")) %>% mutate_all(.,as.numeric) %>%
  left_join(.,feature.1.annotation,by=c("Gene.A.Number"="Gene.Number")) %>% 
  left_join(.,feature.1.annotation,by=c("Gene.B.Number"="Gene.Number")) %>% 
  left_join(.,feature.1.annotation,by=c("Gene.C.Number"="Gene.Number")) %>% 
  replace(is.na(.),max(feature.2.annotation$bresub)) %>% 
  mutate(Delta.Gene.A=bresub.y-bresub.x) %>% 
  mutate(Delta.Gene.B=bresub.x.x-bresub.x) %>% 
  mutate(Delta.Gene.C=bresub.y.y-bresub.x) %>%  
  mutate(pair.index.AB=str_c(Gene.A.Number,'-',Gene.B.Number)) %>% 
  mutate(pair.index.AC=str_c(Gene.A.Number,'-',Gene.C.Number)) %>% 
  mutate(pair.index.BC=str_c(Gene.B.Number,'-',Gene.C.Number)) %>% 
  left_join(.,feature.pair.index,by=c("pair.index.AB"="pair.index")) %>% 
  left_join(.,feature.pair.index,by=c("pair.index.AC"="pair.index")) %>% 
  left_join(.,feature.pair.index,by=c("pair.index.BC"="pair.index")) %>% 
  mutate(Delta.Pair.AB=bresub.x.x.x-bresub.x) %>% 
  mutate(Delta.Pair.AC=bresub.y.y.y-bresub.x) %>% 
  mutate(Delta.Pair.BC=bresub-bresub.x) %>% 
  mutate(Min.Delta=pmin(Delta.Gene.A,Delta.Gene.B,Delta.Gene.C,Delta.Pair.AB,Delta.Pair.AC,Delta.Pair.BC)) %>% 
  select(.,-c(bresub.y,bresub.x.x,bresub.y.y,pair.index.AB,pair.index.AC,pair.index.BC,bresub.x.x.x,bresub.y.y.y,bresub)) %>% 
  setNames(.,c("Gene.A.Number","Gene.B.Number","Gene.C.Number","bresub","Gene.A.Name","Gene.B.Name","Gene.C.Name","Delta-Gene.A","Delta-Gene.B","Delta-Gene.C",
               "Delta-Pair.AB","Delta-Pair.AC","Delta-Pair.BC","Min.Delta"))

write.csv(feature.3.annotation,"Results/LFF-midpoints-3_feature_annotation.csv")