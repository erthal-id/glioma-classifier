# Differential DNAmet - finding probe markers for gliomas fine-tuning Ceccarelli's classification
# dmpFinder() function 
# https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html#identifying-dmrs-and-dmps

# Install ChAMP
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install('ChAMP')

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install('minfi', force = TRUE)

library(minfi)
library(ChAMP)
library(dplyr)
library(data.table)
library(stringr)
library(readr)
library(caret)
install.packages('caret')

# Load DNAmtx
setwd("~/Documents/TCGA_glioma_classifier")
DNAmet_mtx <- read.delim("LGG.GBM.meth.txt")
DNAmet_mtx[1:10, 1:10]
class(DNAmet_mtx)
rownames(DNAmet_mtx) = DNAmet_mtx$Composite.Element.REF
dim(DNAmet_mtx)
DNAmet_mtx = DNAmet_mtx[, 5:936]
head(DNAmet_mtx)
DNAmet_mtx = as.matrix(DNAmet_mtx)
DNAmet_mtx[1:10, 1:10]
range(na.omit(DNAmet_mtx)) #0.00180811 0.99766672

# Load metadata
setwd("~/Documents/TCGA_glioma_classifier")
library(readr)
metadata <- read_csv("mmc2.csv")
metadata = data.frame(metadata)
metadata %>% head()

# Checking column ID
colnames(DNAmet_mtx)[1:10]
length(colnames(DNAmet_mtx)) #932
metadata$Case[1:10]
length(metadata$Case) #1122

# Fixing sample ID in DNAmet_mtx
colnames(DNAmet_mtx) <- str_split_fixed(as.character(colnames(DNAmet_mtx)), "[.][0-9][0-9][A-Z]", 2)[, 1]
colnames(DNAmet_mtx)  <- gsub('[.]', '-', colnames(DNAmet_mtx))
table(duplicated(colnames(DNAmet_mtx)))
# FALSE 
# 932 


# Filtering metadata by DNAmet=TRUE patients
length(intersect(colnames(DNAmet_mtx),
                 metadata$Case)) #932

metadata <- metadata %>% 
  filter(Case %in% colnames(DNAmet_mtx))
dim(metadata) # 932  51
dim(DNAmet_mtx) #25978   932


# Exclude mask and chrs probes
print(load("hm450.anno.Rda")) # hm450.anno 
probes_retain <- subset(hm450.anno, !chrm_A %in% c("chrX","chrY", "chrM") & MASK_general == FALSE)$probeID

DNAmet_mtx_sub = DNAmet_mtx[rownames(DNAmet_mtx) %in% probes_retain, ]
dim(DNAmet_mtx_sub) #23493   932


# Updating labels
metadata$Ceccarelli_six_subtypes <- as.character(metadata$Supervised.DNA.Methylation.Cluster)
metadata[metadata$Ceccarelli_six_subtypes %in% 'PA-like', ]$Ceccarelli_six_subtypes <- 'LGm6-GBM'
metadata$Ceccarelli_six_subtypes %>% table()

# ChAMP - find differentially methylated positions (DMPs)
identical(colnames(DNAmet_mtx_sub), metadata$Case) #FALSE
rownames(metadata) <- metadata$Case
metadata <- metadata[colnames(DNAmet_mtx_sub), ]
identical(colnames(DNAmet_mtx_sub), metadata$Case) #TRUE (okay)
metadata <- metadata[!metadata$Supervised.DNA.Methylation.Cluster %in% NA, ]
DNAmet_mtx_sub <- DNAmet_mtx_sub[, colnames(DNAmet_mtx_sub) %in% metadata$Case]
identical(colnames(DNAmet_mtx_sub), metadata$Case) #TRUE (still true, ok)

# Test 1 - downsampling 
set.seed(42)
Downsample <- createDataPartition(metadata$Ceccarelli_six_subtypes, p=0.3, list=FALSE, times=1) #Downsampling
Downsample <- DNAmet_mtx_sub[Downsample, ]

subtype <- metadata[metadata$Case %in% colnames(Downsample), ]
subtype <- subtype$Ceccarelli_six_subtypes 

dmp <- dmpFinder(dat = Downsample, 
                 pheno = subtype, 
                 type = "categorical",
                 qCutoff = 0.01)
head(dmp)
hist(dmp$qval)

# Test 2 - all saples
subtype <- metadata$Ceccarelli_six_subtypes
dmp <- dmpFinder(dat = DNAmet_mtx_sub, 
                 pheno = subtype, 
                 type = "categorical",
                 qCutoff = 0.01)
head(dmp)
hist(dmp$qval)
dim(dmp)
dmp$probeID <- rownames(dmp)
dmp
dmp[, ]
dmp[!grepl("NA", dmp$probeID), ] %>% dim() #18402     5
dmp <- dmp[!grepl("NA", dmp$probeID), ]
probes_dmpFinder <- dmp$probeID

# probes Diff Methylated by ANOVA 
load("~/Documents/TCGA_glioma_classifier/anova.pvalue.Rda")
w.p.values.adj <- data.frame
w.p.values.adj$probesID <- rownames(w.p.values.adj)
probes_ANOVA <- w.p.values.adj[w.p.values.adj$w.p.values.adj < 0.01, ]$probesID

length(intersect(probes_dmpFinder,
                 probes_ANOVA)) #17242

library(UpSetR)
listInput <- list( probes_dmpFinder = probes_dmpFinder,
                         probes_ANOVA = probes_ANOVA)
upset(fromList(listInput), order.by = "freq", nsets = 2) 



