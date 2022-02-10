########################################################################
# Author    : E. Gerrits
# Date      : 10-jul-2020
# Purpose   : Analyse RNAseq data (lexogen) from Wouter
# Dataset   : Human cell lines
# Updated on 10-07-2020

# if packages are not installed yet try
# install.packages("ggplot2") 
# or google how to install

# lines starting with # are to add comments to your script and won't run in the console

library(edgeR)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(ggrepel)
library(ape)
library(Biostrings)
library(ggplot2)
library(phylogram)
library(readxl)

########################################################################
#							            	Import dataset	              					   #
########################################################################
# import countfile
countmtx <- read.csv("~/Documents/Projects/14. DNA damage Wouter (BioRxiv)/2021/countmtx_wouter_aug2021.csv", row.names = 1)
colnames(countmtx)

# split countmtx into datasets from U2OS and HEK293T samples. 1,10,11,12,13,14,15,16 are U2OS; 17,2,3,4,5,6,7,8 are HEK293T
# Emma: ik zou hem niet splitten als je U2OS met HEK wilt vergelijken nu, want dan ze niet neer beide in de dgList straks
split1 <- countmtx[,c(1:8)]
#split1 <- countmtx

# import the target file
target <- as.data.frame(read_excel("targetfile.xlsx")) # add more information to the targetfile in excel
rownames(target) <- target$sampleID

# make sure order of the rows in targetfile are the same as the column in the countmtx
target <- target[colnames(split1),] 

########################################################################
#							            EnsembleID to genesymbol	       					   #
########################################################################

# a txt file is used of a gtf file (used during the alignment).
geneExp <- read.delim("./Homo_sapiens.GRCh38.91_shorter.txt", header=FALSE)
colnames(geneExp) <- c("Ensembl_ID", "Gene_sym")

# countmtx
geneMat <- geneExp[match(rownames(split1), geneExp$Ensembl_ID),]
rownames(split1) <- make.names(geneMat$Gene_sym, unique=TRUE)

########################################################################
#							                   Functions  	            					   #
########################################################################

obtainFilteredDGE <- function(countMat, targetInfo, genes) {
  # Making sure that a DGEList is made.
  dge <- DGEList(counts=countMat, samples=targetInfo, genes = genes)
  # Determining the expressed genes.
  isExpr <- rowSums(cpm(dge)>1) >= 2
  # These expressed genes will be saved within the dge dataset.
  dge <- dge[isExpr, ]
  # Calculate the normalization factors.
  dge <- calcNormFactors(dge)
  # The dge list is returned.
  return (dge)
}

########################################################################
#					             	   Make a filtered object         					   #
########################################################################

# Creating a DGEList object
dgList <- obtainFilteredDGE(split1, target, rownames(split1))
dgList$samples # so if you have added extra columns to the target file you will also see them here
length(dgList$genes$genes) # number of expressed genes in the dataset. 15000-20000 is normal

# dgList object is ook input voor edgeR differentiele expressie analyse
saveRDS(dgList, "dgList.rds")
#dgList <- readRDS("dgList.rds") # inladen

plotMDS(dgList, labels = rownames(dgList$samples)) # similar to PCA plot

# barplot of library sizes (is the colSum of counts, so the total number of counts per sample)
df <- as.data.frame(dgList$samples)
ggplot(df, aes(x = as.factor(sampleID), y = lib.size)) + # basic of ggplot, define x and y
  geom_bar(stat="identity", fill = "darksalmon") + # add a barplot
  ggtitle("Barplot of library sizes") + # add title
  theme_classic() + # select color theme
  ylab("Library size") + # y axis label
  xlab("Sample") + # x axis label
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) # rotate x labels


# barplot of number of expressed genes per sample
nGene <- data.frame(matrix(NA, nrow = length(colnames(split1)), ncol = 1)) # maak leeg dataframe
rownames(nGene) <- target$sample # rownames word sample names
colnames(nGene) <- "genes"
for (i in 1:ncol(split1)){ # per sample
  nGene[i,1] <- length(which(split1[,i] > 0)) # tellen hoeveel genen een expressie hebben van > 0
}

ggplot(nGene, (aes(x = rownames(nGene), y = nGene$genes))) +
  geom_bar(stat = "identity", fill = "green4") +
  xlab("sample") +
  ylab("Number of genes expressed") +
  ggtitle("Number of genes with count > 0") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

########################################################################
#							               NORMALIZATION		            			       #
########################################################################

# counts per millon
# includes the use of using normalized lib sizes
countsPerMillion <- cpm(dgList) # do CPM
head(countsPerMillion)
#write.csv(countsPerMillion, "CPMsU2OS.csv") # if you want to save the CPM file

# log of the counts per million
logCPMcounts <- cpm(dgList,log=TRUE)
colnames(logCPMcounts) <- target$sample

#--------------------------------------------------------------------#

# Check distributions of samples using boxplots
# here the aim is to get all the boxplots on 1 line (which is perfect here)
# so this means the data is very nicely normalized by taking the logCPMs
boxplot(logCPMcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logCPMcounts),col="blue")
title("Boxplots of logCPMs")

########################################################################
#							                   PCA                    		    	     #
########################################################################
# calculate PCs
project.pca <- prcomp(t(logCPMcounts)) # logCPM matrix moet getransposed worden, vandaar de t
summary(project.pca)

# plot
df <- project.pca$x
df <- df[,1:2]
df <- as.data.frame(df)
df <- cbind(df, dgList$samples)
ggplot(df, aes(PC1, PC2)) + # selecteer x en y
  geom_point(size = 3) +  
  scale_shape_manual(values =  c(21, 23)) +
  geom_text_repel(aes(label = sampleID), nudge_x = 5) +
  xlab(paste("PC1 (", (round((project.pca$sdev[1]^2)/(sum(project.pca$sdev^2))*100, digits = 2)), "%)", sep = "")) +
  ylab(paste("PC2 (", (round((project.pca$sdev[2]^2)/(sum(project.pca$sdev^2))*100, digits = 2)), "%)", sep = "")) +
  theme_classic() +
  ggtitle("PCA logCPM")

# PCA plot colored by library size shows that the separation of the samples is not caused by the library size
# So that is good because that indicates that the samples segregate based on biological variation instead of technical variation
# so i hope these 3 groups represent your experimental groups :) 
ggplot(df, aes(PC1, PC2, colour = interaction(df$cell.line, df$treatment))) + # you can change lib.size by any column in your target file
  geom_point(size = 3) + 
  scale_shape_manual(values =  c(21, 23)) +
  geom_text_repel(aes(label = sampleID), nudge_x = 5) +
  xlab(paste("PC1 (", (round((project.pca$sdev[1]^2)/(sum(project.pca$sdev^2))*100, digits = 2)), "%)", sep = "")) +
  ylab(paste("PC2 (", (round((project.pca$sdev[2]^2)/(sum(project.pca$sdev^2))*100, digits = 2)), "%)", sep = "")) +
  theme_classic() +
  ggtitle("PCA logCPM")

########################################################################
#							             Differential expression      		    	     #
########################################################################
# this is with edgeR. 
# Some people also like DESeq2 or limma. Should give more or less the same results.

# Creating a design 
design <- model.matrix(~0 + dgList$samples$treatment, data=dgList$samples) # basically making a binary matrix that contains the information which sample is in which experimental group

# estimate dispersions 
dge.full <- estimateDisp(dgList, design)
plotBCV(dge.full) 
fit <- glmFit(dge.full, design) # fitting of the glm

# differential gene expression function, contains all settings you need
# change p- and logFC thresholds if you like
DiffExpr <- function(contrast, name){
  # do the actual DE analysis
  lrt <- glmLRT(fit, contrast = contrast)
  tags <- topTags(lrt, n=dim(dge.full[[1]])[1], adjust.method="BH", sort.by="logFC") #this is the full table for all genes
  assign(x = paste("tags", name, sep = ""), value = tags, .GlobalEnv)
  
  # extract significant genes (logFC > 1 and p < 0.05)
  toptable <- tags[[1]][which(abs(tags[[1]]$logFC) > 0.5 & tags[[1]]$FDR < 0.05),] # this is the table with only the significant genes
  plotSmear(lrt, de.tags = toptable$genes)
  abline(h=c(-1, 1), col=2)
  title(paste("DEGS for ", name, sep = ""))
  assign(x = paste("topTable", name, sep = ""), value = toptable, .GlobalEnv)
  
  # print info
  return(summary(decideTests(lrt, p = 0.05, lfc = 0.5, adjust.method = "BH")))
}

# setting of the contrast is based on the colnames of your design. In this example there are 3 groups.
# but if you have more than you should also add more values in de contrast vector. 
colnames(design)
contrast <- c(1, -1) #Emma: moet je even kijken welke kolommen je precies moet hebben

# for example group3 vs group1 then contrast would be:
# contrast <- c(-1, 0, 1) 

# now you have defined everything you need so you can use the DiffExpr function
DiffExpr(contrast, "Group1vsGroup2") # alles wat rood is in the plot is een DEG
write.csv(tagsGroup1vsGroup2$table, "DE_U2OSCPTvsU2OSDMSO.csv") # save as csv file

# plot DE genes in a heatmap
n_allgenes <- logCPMcounts[rownames(topTableGroup1vsGroup2), ] # here you can also opt to extract the values from the countspermillion matrix, but logCPM is usually nicer
colnames(n_allgenes) <- target$sample
colfunc <- colorRampPalette(c("steelblue", "slategray2", "white", "tan1", "firebrick3")) # these are nice colors
heatmap.2(n_allgenes, 
          col=colfunc(50), # define number of sections in color scale
          keysize = 1,
          scale="row", # does a row z-score
          trace="none", # is really ugly
          cexRow = 0.7, # font size of the row labels
          cexCol = 0.9, # font size of the col labels
          #Colv = "NA", # do unsupervised, meaning you don't do column clustering
          dendrogram = 'both', # or choose 'row' or 'col or 'none
          key.title = NA,
          main = "Differentially expressed genes")

# volcano plot
# all dots that have a color are significantly differentially expressed between the groups
# grey dots are genes that are not significant DE (p > 0,05 or logFC > 1)
ggplot(tagsGroup1vsGroup2$table, aes(x = logFC, y = -log(FDR))) +
  geom_point(colour = "grey80", size = 0.5) + 
  geom_hline(yintercept = -log(0.05), colour = "grey80") + # these are the sifnigicane thresholds (p < 0.05)
  geom_vline(xintercept = c(-1,1), colour = "grey80") + # logFC thresholds (logFC > 1 or < -1)
  geom_point(data = topTableGroup1vsGroup2, aes(x = logFC, y = -log(FDR)), shape = 21, 
             colour = "black", 
             fill = ifelse(topTableGroup1vsGroup2$logFC > 1, "violetred3", "mediumaquamarine")) +
  geom_text_repel(data = topTableGroup1vsGroup2, label = ifelse(-log(topTableGroup1vsGroup2$FDR) > 50, topTableGroup1vsGroup2$genes, "")) + # set treshold from when you want to show gene names, or put a # for this line if you dont want to show gene names
  theme_classic() +
  xlim(-15,15) +
  ggtitle("Group1 vs Group2")

#####
genenlijst <- read_excel("~/Documents/Wouter/Postdoc/RNAseq/DataSeq/DEcommonIncrease.xlsx")
subset_table <- topTableGroup1vsGroup2[genenlijst,]
# replace all 0 by 0.0000000001
#countsPerMillion[countsPerMillion == 0] <- 0.000000001


