### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
### Differential Methylation across Ceccareli's glioma subtype ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

#@@@ TCGA glioma samples
#@@@ Aim: tune DNAmet probes (features) to classify those 6 glioma subtypes based on the label establish by Ceccarelli et al 2016
#@@@ probes from 450k and 27k illumina array 
#**IMPORTANT**#
#@@@ If we don't get a good fine tuned model from that we should try to add EPIC array. Maybe only intersections within EPIC and 450k ? 

# Load packages 
library(dplyr)


# Load DNAmtx 
load("/media/hd/maycon/Glioma_classifier/DNAmet_mtx_no_maskprobes_neither_chrprobes.rda")
DNAmet_mtx_sub[1:4, 1:4]
dim(DNAmet_mtx_sub) #23493   932

# Load metadata
library(readr)
metadata <- read_csv("/media/hd/maycon/Glioma_classifier/mmc2.csv")
metadata = data.frame(metadata)
metadata %>% head()
dim(metadata) #1122 51

metadata <- metadata %>% 
  filter(Case %in% colnames(DNAmet_mtx_sub)) #keep only DNAmet samples
dim(metadata) #932  51



# Anova (Tathi's code) -------------------------

# Aggregating 'PA-like' subtype into 'LGm6-GBM' (because now there are 6 subtypes instead of 7)
metadata$Ceccarelli_six_subtypes <- as.character(metadata$Supervised.DNA.Methylation.Cluster)
metadata[metadata$Ceccarelli_six_subtypes %in% 'PA-like', ]$Ceccarelli_six_subtypes <- 'LGm6-GBM'
metadata$Ceccarelli_six_subtypes %>% table()

# Handling NAs 
table(is.na(DNAmet_mtx_sub))
# FALSE     TRUE 
# 20091130  1804346 
DNAmet_mtx_sub <- na.omit(DNAmet_mtx_sub) #removing probes (rows) within at least one beta-value NA
dim(DNAmet_mtx_sub) #18789 932
table(is.na(DNAmet_mtx_sub))
# FALSE 
# 17511348

# Ordering datasets
meta <- metadata
rownames(meta) <- meta$Case
data <- DNAmet_mtx_sub
data <- data[, rownames(meta)]
head(data)

identical(colnames(data), rownames(meta)) #T (samples in the same order. data = your data. meta = meta data)

require(parallel)
values <- as.data.frame(t(data)) #[samples, features]
w.p.values <- unlist(mclapply(values,
                              function(probe) {
                                probe <- data.frame(probe)
                                probe$Ceccarelli_six_subtypes <- meta$Ceccarelli_six_subtypes #column with your groups
                                colnames(probe)[1] <- "value"
                                if(nrow(na.omit(probe)) > 1){ #se você só tiver NA no seu objeto ele não realiza o teste
                                  test <- summary(aov(value ~ Ceccarelli_six_subtypes, data=probe)) #faz o teste anova
                                  if(ncol(test[[1]]) > 3)  #dimensão onde o p-value está armazenado
                                    return(test[[1]][[5]][[1]])
                                  else
                                    return(NA)
                                }
                                else
                                  return(NA)
                                
                              }, mc.cores=8))
w.p.values.adj <- p.adjust(w.p.values, method = "BH")
save(w.p.values, w.p.values.adj, file="/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/anova.pvalue.Rda")
# what about the tukey to correct ANOVA test? ??? Answer: BH method has already been choose

identical(rownames(data),
          names(w.p.values.adj)) #TRUE
data <- as.data.frame(data)
data$p_val_adj <- NA
data$p_val_adj <- as.vector(w.p.values.adj)
data[data$p_val_adj < 0.01, ]

meta$Ceccarelli_six_subtypes %>% table()
# Classic-like            Codel      G-CIMP-high       G-CIMP-low 
# 148                     174         249               25 
# LGm6-GBM                Mesenchymal-like 
# 67                      215

Classic_like <-  meta[meta$Ceccarelli_six_subtypes %in% 'Classic-like', ]$Case
Codel <-  meta[meta$Ceccarelli_six_subtypes %in% 'Codel', ]$Case
G_CIMP_high <- meta[meta$Ceccarelli_six_subtypes %in% 'G-CIMP-high', ]$Case
G_CIMP_low <- meta[meta$Ceccarelli_six_subtypes %in% 'G-CIMP-low', ]$Case
LGm6_GBM <- meta[meta$Ceccarelli_six_subtypes %in% 'LGm6-GBM', ]$Case
Mesenchymal_like <- meta[meta$Ceccarelli_six_subtypes %in% 'Mesenchymal-like', ]$Case

# Classic_like Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Codel)],1,mean,na.rm=T) #group n 
data$DiffMean_Classic_Codel_ClassicOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
data$DiffMean_Classic_G_CIMP_high_ClassicOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
data$DiffMean_Classic_G_CIMP_low_ClassicOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
data$DiffMean_Classic_LGm6_GBM_ClassicOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
data$DiffMean_Classic_Mesenchymal_like_ClassicOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]



#**Alert**#
dim(data) #18789   probes
dim(data_sig) #17701   probes - it seems there is too much significant probes... 
# Obs: it's ok to get the same 17701 significant probes each time we're subseting probes for each glioma subtype because we're just filering the same statistics we had made before within all groups together 


data_sig[data_sig$DiffMean_Classic_Codel_ClassicOrientation > 0.3 &
       data_sig$DiffMean_Classic_G_CIMP_high_ClassicOrientation > 0.3 &
       data_sig$DiffMean_Classic_G_CIMP_low_ClassicOrientation > 0.3 &
       data_sig$DiffMean_Classic_LGm6_GBM_ClassicOrientation > 0.3 &
       data_sig$DiffMean_Classic_Mesenchymal_like_ClassicOrientation > 0.3, ] %>% dim() # 4 probes
library(UpSetR)
listInput_hyper <- list(DiffMean_Classic_Codel_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_Codel_ClassicOrientation > 0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_G_CIMP_high_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_G_CIMP_high_ClassicOrientation > 0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_G_CIMP_low_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_G_CIMP_low_ClassicOrientation > 0.3, ] %>% rownames(),
                 
                  DiffMean_Classic_LGm6_GBM_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_LGm6_GBM_ClassicOrientation > 0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_Mesenchymal_like_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_Mesenchymal_like_ClassicOrientation > 0.3, ] %>% rownames()
                  )

upset(fromList(listInput_hyper), order.by = "freq", nsets = 6) 
# 4/5 overlapped comparisons probeset = 85 probes; 6/6 overlapped groups probeset = 4. Take them all (85+4 = 89 probes) to be the Classic-like subtype probe signature/features 


data_sig[data_sig$DiffMean_Classic_Codel_ClassicOrientation < -0.3 &
       data_sig$DiffMean_Classic_G_CIMP_high_ClassicOrientation < -0.3 &
       data_sig$DiffMean_Classic_G_CIMP_low_ClassicOrientation < -0.3 &
       data_sig$DiffMean_Classic_LGm6_GBM_ClassicOrientation < -0.3 &
       data_sig$DiffMean_Classic_Mesenchymal_like_ClassicOrientation < -0.3, ] %>% dim() # 0 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_Classic_Codel_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_Codel_ClassicOrientation < -0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_G_CIMP_high_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_G_CIMP_high_ClassicOrientation < -0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_G_CIMP_low_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_G_CIMP_low_ClassicOrientation < -0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_LGm6_GBM_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_LGm6_GBM_ClassicOrientation < -0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_Mesenchymal_like_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_Mesenchymal_like_ClassicOrientation < -0.3, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 6) # to many groups overlapping on same probes within DeafMean < -0.3


# Extract the 85+4 hyper methylated probes 
x <- upset(fromList(listInput_hyper), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #413 probes somehow overlapped
x1 <- unlist(listInput_hyper, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #413 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # only probes overlapped on 5/5 comparisons 
x1[ rowSums(x$New_data) == 4] # only probes overlapped on 4/5 comparisons

length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret.

Classic_probeset_89 <- c(x1[ rowSums(x$New_data) == 5], x1[ rowSums(x$New_data) == 4])



# Codel Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(Codel) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Classic_like)],1,mean,na.rm=T) #group n 
data$DiffMean_Codel_Classic_CodelOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Codel) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
data$DiffMean_Codel_G_CIMP_high_CodelOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Codel) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
data$DiffMean_Codel_G_CIMP_low_CodelOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Codel) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
data$DiffMean_Codel_LGm6_GBM_CodelOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Codel) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
data$DiffMean_Codel_Mesenchymal_like_CodelOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]



#**Alert**#
dim(data) #18789   probes
dim(data_sig) #17701   probes - it seems there is too much significant probes... 
# Obs: it's ok to get the same 17701 significant probes each time we're subseting probes for each glioma subtype because we're just filering the same statistics we had made before within all groups together 


# Hyper methylated probes
data_sig[data_sig$DiffMean_Codel_Classic_CodelOrientation > 0.3 &
           data_sig$DiffMean_Codel_G_CIMP_high_CodelOrientation > 0.3 &
           data_sig$DiffMean_Codel_G_CIMP_low_CodelOrientation > 0.3 &
           data_sig$DiffMean_Codel_LGm6_GBM_CodelOrientation > 0.3 &
           data_sig$DiffMean_Codel_Mesenchymal_like_CodelOrientation > 0.3, ] %>% dim() # 14 probes
library(UpSetR)
listInput_hyper <- list(DiffMean_Codel_Classic_CodelOrientation = data_sig[data_sig$DiffMean_Codel_Classic_CodelOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_Codel_G_CIMP_high_CodelOrientation = data_sig[data_sig$DiffMean_Codel_G_CIMP_high_CodelOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_Codel_G_CIMP_low_CodelOrientation = data_sig[data_sig$DiffMean_Codel_G_CIMP_low_CodelOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_Codel_LGm6_GBM_CodelOrientation = data_sig[data_sig$DiffMean_Codel_LGm6_GBM_CodelOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_Codel_Mesenchymal_like_CodelOrientation = data_sig[data_sig$DiffMean_Codel_Mesenchymal_like_CodelOrientation > 0.3, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 
# 4/5 overlapped comparisons probeset = 365 probes; 5/5 overlapped comparisons probeset = 14. Take them all (85+4 = 89 probes) to be the Classic-like subtype probe signature/features 


# Hypo methylated probes
data_sig[data_sig$DiffMean_Codel_Classic_CodelOrientation < -0.3 &
           data_sig$DiffMean_Codel_G_CIMP_high_CodelOrientation < -0.3 &
           data_sig$DiffMean_Codel_G_CIMP_low_CodelOrientation < -0.3 &
           data_sig$DiffMean_Codel_LGm6_GBM_CodelOrientation < -0.3 &
           data_sig$DiffMean_Codel_Mesenchymal_like_CodelOrientation < -0.3, ] %>% dim() #  # 0 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_Codel_Classic_CodelOrientation = data_sig[data_sig$DiffMean_Codel_Classic_CodelOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_Codel_G_CIMP_high_CodelOrientation = data_sig[data_sig$DiffMean_Codel_G_CIMP_high_CodelOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_Codel_G_CIMP_low_CodelOrientation = data_sig[data_sig$DiffMean_Codel_G_CIMP_low_CodelOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_Codel_LGm6_GBM_CodelOrientation = data_sig[data_sig$DiffMean_Codel_LGm6_GBM_CodelOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_Codel_Mesenchymal_like_CodelOrientation = data_sig[data_sig$DiffMean_Codel_Mesenchymal_like_CodelOrientation < -0.3, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5) # to many groups overlapping on same probes within DeafMean < -0.3


# Extract the 365+14 hyper methylated probes 
x <- upset(fromList(listInput_hyper), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #2460 probes somehow overlapped
x1 <- unlist(listInput_hyper, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #2460 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # only probes overlapped on 5/5 comparisons 
x1[ rowSums(x$New_data) == 4] # only probes overlapped on 4/5 comparisons

length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret.

Codel_probeset_379p <- c(x1[ rowSums(x$New_data) == 5], x1[ rowSums(x$New_data) == 4])




# G_CIMP_high Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(G_CIMP_high) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Classic_like)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_high_Classic_GcimpHighOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_high) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Codel)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_high_Codel_GcimpHighOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_high) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_high) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_high) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]



#**Alert**#
dim(data) #18789   probes
dim(data_sig) #17701   probes - it seems there is too much significant probes... 
# Obs: it's ok to get the same 17701 significant probes each time we're subseting probes for each glioma subtype because we're just filering the same statistics we had made before within all groups together 


# Hyper methylated probes
data_sig[data_sig$DiffMean_G_CIMP_high_Classic_GcimpHighOrientation > 0.3 &
           data_sig$DiffMean_G_CIMP_high_Codel_GcimpHighOrientation > 0.3 &
           data_sig$DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation > 0.3 &
           data_sig$DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation > 0.3 &
           data_sig$DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation > 0.3, ] %>% dim() # 0 probes
library(UpSetR)
listInput_hyper <- list(DiffMean_G_CIMP_high_Classic_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Classic_GcimpHighOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_high_Codel_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Codel_GcimpHighOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation > 0.3, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 
# 2 4/5 overlapped comparisons probeset = i) 177 probes and ii) 4 probes. Take them all (177+4 = 181 probes) to be the Classic-like subtype probe signature/features 


# Hypo methylated probes
data_sig[data_sig$DiffMean_G_CIMP_high_Classic_GcimpHighOrientation < -0.3 &
           data_sig$DiffMean_G_CIMP_high_Codel_GcimpHighOrientation < -0.3 &
           data_sig$DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation < -0.3 &
           data_sig$DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation < -0.3 &
           data_sig$DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation < -0.3, ] %>% dim() # 0 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_G_CIMP_high_Classic_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Classic_GcimpHighOrientation < -0.3 , ] %>% rownames(),
                       
                       DiffMean_G_CIMP_high_Codel_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Codel_GcimpHighOrientation < -0.3 , ] %>% rownames(),
                       
                       DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation < -0.3 , ] %>% rownames(),
                       
                       DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation < -0.3 , ] %>% rownames(),
                       
                       DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation < -0.3 , ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5) # to many groups overlapping on same probes within DeafMean < -0.3


# Extract the 365+14 hyper methylated probes 
x <- upset(fromList(listInput_hyper), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #2012 probes somehow overlapped
x1 <- unlist(listInput_hyper, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #2012 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # none 
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 groups

# length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret. - don't need to do this

GcimpHigh_probeset_181p <- c(#x1[ rowSums(x$New_data) == 5], - don't need to do this
                         x1[ rowSums(x$New_data) == 4])




# G_CIMP_low Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(G_CIMP_low) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Classic_like)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_low_Classic_GcimpLowOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_low) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Codel)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_low_Codel_GcimpLowOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_low) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_low) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_low) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]



#**Alert**#
dim(data) #18789   probes
dim(data_sig) #17701   probes - it seems there is too much significant probes... 
# Obs: it's ok to get the same 17701 significant probes each time we're subseting probes for each glioma subtype because we're just filering the same statistics we had made before within all groups together 


# Hyper methylated probes
data_sig[data_sig$DiffMean_G_CIMP_low_Classic_GcimpLowOrientation > 0.2 &
           data_sig$DiffMean_G_CIMP_low_Codel_GcimpLowOrientation > 0.2 &
           data_sig$DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation > 0.2 &
           data_sig$DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation> 0.2 &
           data_sig$DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation> 0.2, ] %>% dim() # 5 probes
library(UpSetR)
listInput_hyper <- list(DiffMean_G_CIMP_low_Classic_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Classic_GcimpLowOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_low_Codel_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Codel_GcimpLowOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation > 0.2, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 
# We had to set the threshold to 0.2 otherwise we couldn't find enough probes (~ 90 probes as in the other comparisons we've done)


# Hypo methylated probes
data_sig[data_sig$DiffMean_G_CIMP_low_Classic_GcimpLowOrientation < -0.2 &
           data_sig$DiffMean_G_CIMP_low_Codel_GcimpLowOrientation < -0.2  &
           data_sig$DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation < -0.2  &
           data_sig$DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation < -0.2  &
           data_sig$DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation < -0.2 , ] %>% dim() # 69 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_G_CIMP_low_Classic_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Classic_GcimpLowOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_G_CIMP_low_Codel_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Codel_GcimpLowOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation < -0.2, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5) # In this one we got some  hypo methylated probes ! So for the GcimpLow probeset we're getting both hyper and hypo methylated probes 


# Extract hyper methylated probes 
x <- upset(fromList(listInput_hyper), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #1798 probes somehow overlapped
x1 <- unlist(listInput_hyper, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #1798 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # all probes present it only 5/5 comparisons
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 comparisons

length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret.

hyper_probes <- c(x1[ rowSums(x$New_data) == 5],
                  x1[ rowSums(x$New_data) == 4])

# Extract hypo methylated probes 
x <- upset(fromList(listInput_hypo), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #2426 probes somehow overlapped
x1 <- unlist(listInput_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #2426 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # all probes present it only 5/5 comparisons
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 comparisons

length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret.

hypo_probes <- c(x1[ rowSums(x$New_data) == 5],
                  x1[ rowSums(x$New_data) == 4])

length(intersect(hyper_probes, hypo_probes)) # zero. Corret.
length(hyper_probes) #64
length(hypo_probes) #115
#  115+64 = 179 probes at total
GcimpLow_probeset_179p <- c(hyper_probes, hypo_probes)





# LGm6_GBM Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(LGm6_GBM) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Classic_like)],1,mean,na.rm=T) #group n 
data$DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(LGm6_GBM) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Codel)],1,mean,na.rm=T) #group n 
data$DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(LGm6_GBM) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
data$DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(LGm6_GBM) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
data$DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(LGm6_GBM) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
data$DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]



#**Alert**#
dim(data) #18789   probes
dim(data_sig) #17701   probes - it seems there is too much significant probes... 
# Obs: it's ok to get the same 17701 significant probes each time we're subseting probes for each glioma subtype because we're just filering the same statistics we had made before within all groups together 


# Hyper methylated probes
data_sig[data_sig$DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation > 0.2 &
           data_sig$DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation > 0.2 &
           data_sig$DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation > 0.2 &
           data_sig$DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation> 0.2 &
           data_sig$DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation > 0.2, ] %>% dim() # 0 probes
library(UpSetR)
listInput_hyper <- list(DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation > 0.2, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 
# We had to set the threshold to 0.2 otherwise we couldn't find enough probes (~ 90 probes as in the other comparisons we've done)


# Hypo methylated probes
data_sig[data_sig$DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation < -0.2 &
           data_sig$DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation < -0.2 &
           data_sig$DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation < -0.2 &
           data_sig$DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation < -0.2 &
           data_sig$DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation < -0.2, ] %>% dim() # 28 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation < -0.2, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5) # In this one (LGm6) we got ONLY  hypo methylated probes !   


# Extract hypo methylated probes 
x <- upset(fromList(listInput_hypo), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #4216 probes somehow overlapped
x1 <- unlist(listInput_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #4216 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # all probes present it only 5/5 comparisons
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 comparisons

length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret.
length(x1[ rowSums(x$New_data) == 5]) #28
length(x1[ rowSums(x$New_data) == 4]) #206
# 28 + 206 = 234

LGm6_probeset_234p <- c(x1[ rowSums(x$New_data) == 5],
                        x1[ rowSums(x$New_data) == 4])





# Mesenchymal_like Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(Mesenchymal_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Classic_like)],1,mean,na.rm=T) #group n 
data$DiffMean_Mesenchymal_Classic_MesenchymaOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Mesenchymal_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Codel)],1,mean,na.rm=T) #group n 
data$DiffMean_Mesenchymal_Codel_MesenchymaOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Mesenchymal_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
data$DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Mesenchymal_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
data$DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Mesenchymal_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
data$DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]


#**Alert**#
dim(data) #18789   probes
dim(data_sig) #17701   probes - it seems there is too much significant probes... 
# Obs: it's ok to get the same 17701 significant probes each time we're subseting probes for each glioma subtype because we're just filering the same statistics we had made before within all groups together 


# Hyper methylated probes
data_sig[data_sig$DiffMean_Mesenchymal_Classic_MesenchymaOrientation > 0.3 &
           data_sig$DiffMean_Mesenchymal_Codel_MesenchymaOrientation > 0.3 &
           data_sig$DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation> 0.3 &
           data_sig$DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation > 0.3 &
           data_sig$DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation, ] %>% dim() # 0 probes
library(UpSetR)
listInput_hyper <- list(DiffMean_Mesenchymal_Classic_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_Classic_MesenchymaOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_Mesenchymal_Codel_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_Codel_MesenchymaOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation > 0.2, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 
# too few probes and one comparison fell out for not having anyprobe in common to the other comparisons

# Hypo methylated probes
data_sig[data_sig$DiffMean_Mesenchymal_Classic_MesenchymaOrientation < -0.2 &
           data_sig$DiffMean_Mesenchymal_Codel_MesenchymaOrientation < -0.2 &
           data_sig$DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation < -0.2 &
           data_sig$DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation < -0.2 &
           data_sig$DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation, ] %>% dim() # 19 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_Mesenchymal_Classic_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_Classic_MesenchymaOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_Mesenchymal_Codel_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_Codel_MesenchymaOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation < -0.2, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5) # In this one (Mesenchymal) we got ONLY  hypo methylated probes but still too few probes (20 probes) ! I found more intresting keep few specific probes than picking a lower threshold (eg 0.1) because it might increase non-specific probes into the model


# Extract hypo methylated probes 
x <- upset(fromList(listInput_hypo), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #3828 probes somehow overlapped
x1 <- unlist(listInput_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #3828 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # none
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 comparisons

# length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret. - not necessary 
length(x1[ rowSums(x$New_data) == 4]) #20
# 20 probes at total

Mesenchymal_probeset_20p <- c(x1[ rowSums(x$New_data) == 4])


DMP_anova_diffmean_subtypes_probeset_list <- list(Classic_probeset_89 = Classic_probeset_89,
                                               Codel_probeset_379p = Codel_probeset_379p,
                                               GcimpHigh_probeset_181p= GcimpHigh_probeset_181p,
                                               GcimpLow_probeset_179p = GcimpLow_probeset_179p,
                                               LGm6_probeset_234p = LGm6_probeset_234p,
                                               Mesenchymal_probeset_20p = Mesenchymal_probeset_20p)

saveRDS(DMP_anova_diffmean_subtypes_probeset_list, file = '/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/DMP_anova_diffmean_subtypes_probeset_list.rds')


# Heatmap visualization -----------
head(metadata) # from the beginning of this code 

my_sample_col <- data.frame(row.names = metadata$Case, Subtypes_six = metadata$Ceccarelli_six_subtypes, Subtypes_seven = metadata$Supervised.DNA.Methylation.Cluster, DNAmet_cluster = metadata$IDH.specific.DNA.Methylation.Cluster) # patient ID in row names; any columns are for sample information

data_sig # DNAmtx 
DNAmtx <- data_sig[, colnames(data_sig) %in% metadata$Case] #keep only beta-values

length(DMP_anova_diffmean_subtypes_probeset_list) 
tune_probset <- c(DMP_anova_diffmean_subtypes_probeset_list[[1]],
  DMP_anova_diffmean_subtypes_probeset_list[[2]],
  DMP_anova_diffmean_subtypes_probeset_list[[3]],
  DMP_anova_diffmean_subtypes_probeset_list[[4]],
  DMP_anova_diffmean_subtypes_probeset_list[[5]],
  DMP_anova_diffmean_subtypes_probeset_list[[6]])
  


library(pheatmap)
p = pheatmap(DNAmtx[rownames(DNAmtx) %in% tune_probset, ], 
             #annotation_row = my_probe_col, 
             annotation_col = my_sample_col,
             show_rownames = FALSE,
             main = 'Glioma subtypes tuned probes'); p

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p, "/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/finetunedprobes_heatmap.pdf")

# Ceccareli pan_glioma_probes
library(openxlsx)
pan_glioma_probes <- read.xlsx("/media/hd/maycon/Glioma_classifier/PanGlioma_MethylationSignatures.xlsx", sheet = 1) 
pan_glioma_probes <- pan_glioma_probes$`1,300.pan-glioma.tumor.specific.probes.(Figure.2A)`


library(pheatmap)
p_cecca = pheatmap(DNAmtx[rownames(DNAmtx) %in% pan_glioma_probes, ], 
             #annotation_row = my_probe_col, 
             annotation_col = my_sample_col,
             show_rownames = FALSE,
             main = 'Glioma subtypes tuned probes'); p_cecca

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p_cecca, "/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/ceccarelli_pangliomaprobes_heatmap.pdf")

# We still need to clean up intersected probes across probesets 
# ATTENTION: maybe the intesect probes are methylated in different orientations. Plots a headmap to evaluate the 6 probesets before start removing any probes
# 1. plot upsetplot of all 6 probesets
# 2. extract only probes in x1[ rowSums(x$New_data) == 1]. Rembemer now it's 6 groups

# Unique probes from tuned probeset
library(UpSetR)
listInput_allprobeset <- list(Classic_probeset_89 = DMP_anova_diffmean_subtypes_probeset_list[[1]],
                       
                       Codel_probeset_379p = DMP_anova_diffmean_subtypes_probeset_list[[2]],
                       
                       GcimpHigh_probeset_181p = DMP_anova_diffmean_subtypes_probeset_list[[3]],
                       
                       GcimpLow_probeset_179p =  DMP_anova_diffmean_subtypes_probeset_list[[4]],
                       
                       LGm6_probeset_234p = DMP_anova_diffmean_subtypes_probeset_list[[5]],
                       
                       Mesenchymal_probeset_20p = DMP_anova_diffmean_subtypes_probeset_list[[6]]
                       
)

upset(fromList(listInput_allprobeset), order.by = "freq", nsets = 6)

# Extract non-intersected probes
x <- upset(fromList(listInput_allprobeset), nsets = 6)
x$New_data[1:5, 1:5]
dim(x$New_data) #901 
x1 <- unlist(listInput_allprobeset, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #901 

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 1] # unique probes for each group (glioma subtype)


tune_probset_unique <- x1[ rowSums(x$New_data) == 1] 


library(pheatmap)
p_uniq = pheatmap(DNAmtx[rownames(DNAmtx) %in% tune_probset_unique, ], 
             #annotation_row = my_probe_col, 
             annotation_col = my_sample_col,
             show_rownames = FALSE,
             main = 'Glioma subtypes tuned probes'); p_uniq

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p_uniq, "/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/finetunedprobes_uniqueprobes_heatmap.pdf")



library(ComplexHeatmap)
library(matlab)
#label 'my_sample_col'
#                    Subtypes_six   Subtypes_seven DNAmet_cluster
# TCGA-CS-4938      G-CIMP-high      G-CIMP-high      IDHmut-K2
# TCGA-CS-4941 Mesenchymal-like Mesenchymal-like       IDHwt-K2
# TCGA-CS-4942      G-CIMP-high      G-CIMP-high      IDHmut-K2
# TCGA-CS-4943      G-CIMP-high      G-CIMP-high      IDHmut-K2
# TCGA-CS-4944      G-CIMP-high      G-CIMP-high      IDHmut-K2
# TCGA-CS-5390            Codel            Codel      IDHmut-K3
# Current label

### Pallete 
library(viridis)
n_colors <- length(table(my_sample_col$Subtypes_six))
Subtypes_six_pal <- viridis(n = n_colors, option = "plasma", direction = -1)
n_colors <- length(table(my_sample_col$Subtypes_seven))
Subtypes_seven_pal <- viridis(n = n_colors, option = "inferno", direction = -1)
n_colors <- length(table(my_sample_col$DNAmet_cluster))
DNAmet_cluster_pal <- viridis(n = n_colors, option = "cividis", direction = -1)

Subtypes_six_lvs <- names(table(my_sample_col$Subtypes_six))
Subtypes_seven_lvs <- names(table(my_sample_col$Subtypes_seven))
DNAmet_cluster_lvs <- names(table(my_sample_col$DNAmet_cluster))

# Name yout columns and it's levels
column_name_1 <- "Subtypes_six"
levels_1 <- Subtypes_six_lvs

column_name_2 <- "Subtypes_seven"
levels_2 <- Subtypes_seven_lvs

column_name_3 <- "DNAmet_cluster"
levels_3 <- DNAmet_cluster_lvs

# Generate color vector
colors_1 <- Subtypes_six_pal
colors_2 <- Subtypes_seven_pal
colors_3 <- DNAmet_cluster_pal
# Assign names to color vector
names(colors_1) <- levels_1
names(colors_2) <- levels_2
names(colors_3) <- levels_3
# Create named list
color_mapping <- list(column_name_1 = colors_1, 
                      column_name_2 = colors_2,
                      column_name_3 = colors_3)


top.anno = HeatmapAnnotation(df = my_sample_col, 
                             col= color_mapping,
                             show_annotation_name = T, annotation_name_gp = gpar(fontsize=7),
                             na_col= "white")


#hm_tun_pset <- Heatmap(as.matrix(DNAmtx[rownames(DNAmtx) %in% tune_probset, ]),
#hm_uniq_tun_pset <- Heatmap(as.matrix(DNAmtx[rownames(DNAmtx) %in% tune_probset_unique, ]),
hm_pan_glioma <- Heatmap(as.matrix(DNAmtx[rownames(DNAmtx) %in% pan_glioma_probes, ]),
                               cluster_columns=T,
                               cluster_rows = T,
                               clustering_method_rows="complete",  
                               show_row_names = F,
                               row_names_gp = gpar(fontsize = 7),
                               show_column_names = F,
                               name = "CpG probes methylation",
                               col= jet.colors(75),
                               row_title = "CpG probes",
                               #row_names_gp = gpar(fontsize = 12),
                               #column_title = "Glioma subtypes tuned probes",
                               #column_title = "Glioma subtypes unique tuned probes",
                               column_title = "Glioma subtypes pan_glioma probes (Ceccarelli)",
                               #split=label.EMT.genes$E_ou_M,
                               #row_names_side="left",
                               
                               top_annotation = top.anno,
                               #row_names_gp = gpar(fontsize = 8),
                               #column_names_gp = gpar(fontsize = 12),
                               heatmap_legend_param = list(
                                 color_bar = 'continuous',
                                 legend_direction = 'vertical',
                                 legend_width = unit(8, 'cm'),
                                 legend_height = unit(5.0, 'cm'),
                                 title_position = 'leftcenter-rot',
                                 title_gp=gpar(fontsize = 12, fontface = 'bold'),
                                 labels_gp=gpar(fontsize = 12, fontface = 'bold')))

#print(hm_tun_pset)
#print(hm_uniq_tun_pset)
print(hm_pan_glioma)

### For left side annotation
# RowAnn <- HeatmapAnnotation(df = label, col=list("Epi.Mes.CellC" = c("Epi" = "blue", "Mes" = "black", "c.Cycle" = "cadetblue1")),
#                             show_annotation_name = T, annotation_name_gp = gpar(fontsize=7),
#                             na_col= "white",which = "row",show_legend =T)


### For left side annotation
# draw(hm_EMT.cCycle.MALTA + RowAnn) 

pdf("/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/panglioma_probeset.pdf",width = 5, height = 10)
### For left side annotation
# draw(hm_EMT.cCycle.MALTA, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()




#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#**Tathi's advice**#
# - use knn to replace NA 
# - be MORE stringest on DiffMean threshold. We don't want to many probes
# - run Random Forest with the probe list we got until know including Ceccarelli's (pan glioma probset)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Prepare the data again ------------
# Load DNAmtx 
load("/media/hd/maycon/Glioma_classifier/DNAmet_mtx_no_maskprobes_neither_chrprobes.rda")
DNAmet_mtx_sub[1:4, 1:4]
dim(DNAmet_mtx_sub) #23493   932

# Load metadata
library(readr)
metadata <- read_csv("/media/hd/maycon/Glioma_classifier/mmc2.csv")
metadata = data.frame(metadata)
metadata %>% head()
dim(metadata) #1122 51

metadata <- metadata %>% 
  filter(Case %in% colnames(DNAmet_mtx_sub)) #keep only DNAmet samples
dim(metadata) #932  51


# Aggregating 'PA-like' subtype into 'LGm6-GBM' (because now there are 6 subtypes instead of 7)
metadata$Ceccarelli_six_subtypes <- as.character(metadata$Supervised.DNA.Methylation.Cluster)
metadata[metadata$Ceccarelli_six_subtypes %in% 'PA-like', ]$Ceccarelli_six_subtypes <- 'LGm6-GBM'
metadata$Ceccarelli_six_subtypes %>% table()

# Handling NAs 
table(is.na(DNAmet_mtx_sub))

library(impute)
DNAmet_mtx_sub <- impute.knn(as.matrix(DNAmet_mtx_sub), k = 10, rowmax = 0.8, colmax = 0.8, maxp = 1500, rng.seed=362436069)[[1]]
# Warning message:
#   In knnimp(x, k, maxmiss = rowmax, maxp = maxp) :
#   1919 rows with more than 80 % entries missing;
# mean imputation used for these rows

table(is.na(DNAmet_mtx_sub)) # all NAs have been replaced


# Ordering datasets
meta <- metadata
rownames(meta) <- meta$Case
data <- DNAmet_mtx_sub
data <- data[, rownames(meta)]
head(data)

identical(colnames(data), rownames(meta)) #T (samples in the same order. data = your data. meta = meta data)

require(parallel)
values <- as.data.frame(t(data)) #[samples, features]
w.p.values <- unlist(mclapply(values,
                              function(probe) {
                                probe <- data.frame(probe)
                                probe$Ceccarelli_six_subtypes <- meta$Ceccarelli_six_subtypes #column with your groups
                                colnames(probe)[1] <- "value"
                                if(nrow(na.omit(probe)) > 1){ #se você só tiver NA no seu objeto ele não realiza o teste
                                  test <- summary(aov(value ~ Ceccarelli_six_subtypes, data=probe)) #faz o teste anova
                                  if(ncol(test[[1]]) > 3)  #dimensão onde o p-value está armazenado
                                    return(test[[1]][[5]][[1]])
                                  else
                                    return(NA)
                                }
                                else
                                  return(NA)
                                
                              }, mc.cores=8))
w.p.values.adj <- p.adjust(w.p.values, method = "BH")
save(w.p.values, w.p.values.adj, file="/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/KNN_NAreplace_anova.pvalue.Rda")
# what about the tukey to correct ANOVA test? ??? Answer: BH method has already been choose

identical(rownames(data),
          names(w.p.values.adj)) #TRUE
data <- as.data.frame(data)
data$p_val_adj <- NA
data$p_val_adj <- as.vector(w.p.values.adj)
data[data$p_val_adj < 0.01, ] %>% dim()

meta$Ceccarelli_six_subtypes %>% table()
# Classic-like            Codel      G-CIMP-high       G-CIMP-low 
# 148                     174         249               25 
# LGm6-GBM                Mesenchymal-like 
# 67                      215

Classic_like <-  meta[meta$Ceccarelli_six_subtypes %in% 'Classic-like', ]$Case
Codel <-  meta[meta$Ceccarelli_six_subtypes %in% 'Codel', ]$Case
G_CIMP_high <- meta[meta$Ceccarelli_six_subtypes %in% 'G-CIMP-high', ]$Case
G_CIMP_low <- meta[meta$Ceccarelli_six_subtypes %in% 'G-CIMP-low', ]$Case
LGm6_GBM <- meta[meta$Ceccarelli_six_subtypes %in% 'LGm6-GBM', ]$Case
Mesenchymal_like <- meta[meta$Ceccarelli_six_subtypes %in% 'Mesenchymal-like', ]$Case


# Classic_like Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Codel)],1,mean,na.rm=T) #group n 
data$DiffMean_Classic_Codel_ClassicOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
data$DiffMean_Classic_G_CIMP_high_ClassicOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
data$DiffMean_Classic_G_CIMP_low_ClassicOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
data$DiffMean_Classic_LGm6_GBM_ClassicOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
data$DiffMean_Classic_Mesenchymal_like_ClassicOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]



#**Alert**#
dim(data) #23493   probes
dim(data_sig) #22180   probes - it seems there is too much significant probes... 
# Obs: it's ok to get the same 22180 significant probes each time we're subseting probes for each glioma subtype because we're just filering the same statistics we had made before within all groups together 


data_sig[data_sig$DiffMean_Classic_Codel_ClassicOrientation > 0.3 &
       data_sig$DiffMean_Classic_G_CIMP_high_ClassicOrientation > 0.3 &
       data_sig$DiffMean_Classic_G_CIMP_low_ClassicOrientation > 0.3 &
       data_sig$DiffMean_Classic_LGm6_GBM_ClassicOrientation > 0.3 &
       data_sig$DiffMean_Classic_Mesenchymal_like_ClassicOrientation > 0.3, ] %>% dim() # 4 probes
library(UpSetR)
listInput_hyper <- list(DiffMean_Classic_Codel_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_Codel_ClassicOrientation > 0.4, ] %>% rownames(),
                  
                  DiffMean_Classic_G_CIMP_high_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_G_CIMP_high_ClassicOrientation > 0.4, ] %>% rownames(),
                  
                  DiffMean_Classic_G_CIMP_low_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_G_CIMP_low_ClassicOrientation > 0.4, ] %>% rownames(),
                 
                  DiffMean_Classic_LGm6_GBM_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_LGm6_GBM_ClassicOrientation > 0.4, ] %>% rownames(),
                  
                  DiffMean_Classic_Mesenchymal_like_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_Mesenchymal_like_ClassicOrientation > 0.4, ] %>% rownames()
                  )

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 


data_sig[data_sig$DiffMean_Classic_Codel_ClassicOrientation < -0.3 &
       data_sig$DiffMean_Classic_G_CIMP_high_ClassicOrientation < -0.3 &
       data_sig$DiffMean_Classic_G_CIMP_low_ClassicOrientation < -0.3 &
       data_sig$DiffMean_Classic_LGm6_GBM_ClassicOrientation < -0.3 &
       data_sig$DiffMean_Classic_Mesenchymal_like_ClassicOrientation < -0.3, ] %>% dim() # 0 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_Classic_Codel_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_Codel_ClassicOrientation < -0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_G_CIMP_high_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_G_CIMP_high_ClassicOrientation < -0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_G_CIMP_low_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_G_CIMP_low_ClassicOrientation < -0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_LGm6_GBM_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_LGm6_GBM_ClassicOrientation < -0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_Mesenchymal_like_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_Mesenchymal_like_ClassicOrientation < -0.3, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5)  


# Extract the hyper methylated probes 
x <- upset(fromList(listInput_hyper), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #167 probes somehow overlapped
x1 <- unlist(listInput_hyper, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #167 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # none
x1[ rowSums(x$New_data) == 4] # only probes overlapped on 4/5 comparisons

Classic_probeset_36p <- c(x1[ rowSums(x$New_data) == 4])
# only hyper probes; diffmean > 0.4; pvalue < 0.05; 4/5 comparisons had this probeset


# Codel Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(Codel) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Classic_like)],1,mean,na.rm=T) #group n 
data$DiffMean_Codel_Classic_CodelOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Codel) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
data$DiffMean_Codel_G_CIMP_high_CodelOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Codel) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
data$DiffMean_Codel_G_CIMP_low_CodelOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Codel) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
data$DiffMean_Codel_LGm6_GBM_CodelOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Codel) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
data$DiffMean_Codel_Mesenchymal_like_CodelOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]



#**Alert**#
dim(data) #23493
dim(data_sig) #22180


# Hyper methylated probes
data_sig[data_sig$DiffMean_Codel_Classic_CodelOrientation > 0.5 &
           data_sig$DiffMean_Codel_G_CIMP_high_CodelOrientation > 0.5 &
           data_sig$DiffMean_Codel_G_CIMP_low_CodelOrientation > 0.5 &
           data_sig$DiffMean_Codel_LGm6_GBM_CodelOrientation > 0.5 &
           data_sig$DiffMean_Codel_Mesenchymal_like_CodelOrientation > 0.5, ] %>% dim()

library(UpSetR)
listInput_hyper <- list(DiffMean_Codel_Classic_CodelOrientation = data_sig[data_sig$DiffMean_Codel_Classic_CodelOrientation > 0.5, ] %>% rownames(),
                        
                        DiffMean_Codel_G_CIMP_high_CodelOrientation = data_sig[data_sig$DiffMean_Codel_G_CIMP_high_CodelOrientation > 0.5, ] %>% rownames(),
                        
                        DiffMean_Codel_G_CIMP_low_CodelOrientation = data_sig[data_sig$DiffMean_Codel_G_CIMP_low_CodelOrientation > 0.5, ] %>% rownames(),
                        
                        DiffMean_Codel_LGm6_GBM_CodelOrientation = data_sig[data_sig$DiffMean_Codel_LGm6_GBM_CodelOrientation > 0.5, ] %>% rownames(),
                        
                        DiffMean_Codel_Mesenchymal_like_CodelOrientation = data_sig[data_sig$DiffMean_Codel_Mesenchymal_like_CodelOrientation > 0.5, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 



# Hypo methylated probes
data_sig[data_sig$DiffMean_Codel_Classic_CodelOrientation < -0.3 &
           data_sig$DiffMean_Codel_G_CIMP_high_CodelOrientation < -0.3 &
           data_sig$DiffMean_Codel_G_CIMP_low_CodelOrientation < -0.3 &
           data_sig$DiffMean_Codel_LGm6_GBM_CodelOrientation < -0.3 &
           data_sig$DiffMean_Codel_Mesenchymal_like_CodelOrientation < -0.3, ] %>% dim() #  # 0 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_Codel_Classic_CodelOrientation = data_sig[data_sig$DiffMean_Codel_Classic_CodelOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_Codel_G_CIMP_high_CodelOrientation = data_sig[data_sig$DiffMean_Codel_G_CIMP_high_CodelOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_Codel_G_CIMP_low_CodelOrientation = data_sig[data_sig$DiffMean_Codel_G_CIMP_low_CodelOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_Codel_LGm6_GBM_CodelOrientation = data_sig[data_sig$DiffMean_Codel_LGm6_GBM_CodelOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_Codel_Mesenchymal_like_CodelOrientation = data_sig[data_sig$DiffMean_Codel_Mesenchymal_like_CodelOrientation < -0.3, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5) # to few comparisons overlapped


# Extract the hyper methylated probes 
x <- upset(fromList(listInput_hyper), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #605 probes somehow overlapped
x1 <- unlist(listInput_hyper, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #605 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # only probes overlapped on 5/5 comparisons 
x1[ rowSums(x$New_data) == 4] # only probes overlapped on 4/5 comparisons

length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret.

Codel_probeset_10p <- c(x1[ rowSums(x$New_data) == 5], x1[ rowSums(x$New_data) == 4])
# only hyper probes; diffmean > 0.5; pvalue < 0.05; 5/5 and 4/5 comparisons had this probeset




# G_CIMP_high Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(G_CIMP_high) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Classic_like)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_high_Classic_GcimpHighOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_high) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Codel)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_high_Codel_GcimpHighOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_high) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_high) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_high) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]



#**Alert**#
dim(data) #23493   
dim(data_sig) #22180  


# Hyper methylated probes
data_sig[data_sig$DiffMean_G_CIMP_high_Classic_GcimpHighOrientation > 0.3 &
           data_sig$DiffMean_G_CIMP_high_Codel_GcimpHighOrientation > 0.3 &
           data_sig$DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation > 0.3 &
           data_sig$DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation > 0.3 &
           data_sig$DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation > 0.3, ] %>% dim() 
library(UpSetR)
listInput_hyper <- list(DiffMean_G_CIMP_high_Classic_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Classic_GcimpHighOrientation > 0.4, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_high_Codel_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Codel_GcimpHighOrientation > 0.4, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation > 0.4, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation > 0.4, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation > 0.4, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 

# Hypo methylated probes
data_sig[data_sig$DiffMean_G_CIMP_high_Classic_GcimpHighOrientation < -0.3 &
           data_sig$DiffMean_G_CIMP_high_Codel_GcimpHighOrientation < -0.3 &
           data_sig$DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation < -0.3 &
           data_sig$DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation < -0.3 &
           data_sig$DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation < -0.3, ] %>% dim() # 0 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_G_CIMP_high_Classic_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Classic_GcimpHighOrientation < -0.4 , ] %>% rownames(),
                       
                       DiffMean_G_CIMP_high_Codel_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Codel_GcimpHighOrientation < -0.4  , ] %>% rownames(),
                       
                       DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation < -0.4  , ] %>% rownames(),
                       
                       DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation < -0.4  , ] %>% rownames(),
                       
                       DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation < -0.4  , ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5) # 


# Extract the hyper methylated probes 
x <- upset(fromList(listInput_hyper), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #1079 
x1 <- unlist(listInput_hyper, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #1079 

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # none 
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 groups

# length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret. - don't need to do this

GcimpHigh_probeset_18p <- c(x1[ rowSums(x$New_data) == 4])
# only hyper probes; diffmean > 0.4; pvalue < 0.05; 4/5 comparisons had this probeset


# G_CIMP_low Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(G_CIMP_low) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Classic_like)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_low_Classic_GcimpLowOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_low) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Codel)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_low_Codel_GcimpLowOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_low) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_low) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_low) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]



#**Alert**#
dim(data) #23493   
dim(data_sig) #22180   


# Hyper methylated probes
data_sig[data_sig$DiffMean_G_CIMP_low_Classic_GcimpLowOrientation > 0.2 &
           data_sig$DiffMean_G_CIMP_low_Codel_GcimpLowOrientation > 0.2 &
           data_sig$DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation > 0.2 &
           data_sig$DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation> 0.2 &
           data_sig$DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation> 0.2, ] %>% dim() 
library(UpSetR)
listInput_hyper <- list(DiffMean_G_CIMP_low_Classic_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Classic_GcimpLowOrientation > 0.4, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_low_Codel_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Codel_GcimpLowOrientation > 0.4, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation > 0.4, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation > 0.4, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation > 0.4, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 
# We had to set the threshold to 0.2 otherwise we couldn't find enough probes (~ 90 probes as in the other comparisons we've done)


# Hypo methylated probes
data_sig[data_sig$DiffMean_G_CIMP_low_Classic_GcimpLowOrientation < -0.2 &
           data_sig$DiffMean_G_CIMP_low_Codel_GcimpLowOrientation < -0.2  &
           data_sig$DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation < -0.2  &
           data_sig$DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation < -0.2  &
           data_sig$DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation < -0.2 , ] %>% dim() # 69 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_G_CIMP_low_Classic_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Classic_GcimpLowOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_G_CIMP_low_Codel_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Codel_GcimpLowOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation < -0.3, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5)

# Extract hyper methylated probes 
x <- upset(fromList(listInput_hyper), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #684 probes somehow overlapped
x1 <- unlist(listInput_hyper, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #684 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # none
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 comparisons

length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret.

hyper_probes <- c(x1[ rowSums(x$New_data) == 5],
                  x1[ rowSums(x$New_data) == 4])

# Extract hypo methylated probes 
x <- upset(fromList(listInput_hypo), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #1321 
x1 <- unlist(listInput_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #1321

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # all probes present it only 5/5 comparisons
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 comparisons

length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret.

hypo_probes <- c(x1[ rowSums(x$New_data) == 5],
                  x1[ rowSums(x$New_data) == 4])

length(intersect(hyper_probes, hypo_probes)) # zero. Corret.
GcimpLow_probeset_22p <- c(hyper_probes, hypo_probes)
# hyper probes; diffmean > 0.4; pvalue < 0.05; 4/5 comparisons had this probeset
# hypo probes; diffmean < -0.3; pvalue < 0.05; 5/5 and 4/5 comparisons had this probeset





# LGm6_GBM Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(LGm6_GBM) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Classic_like)],1,mean,na.rm=T) #group n 
data$DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(LGm6_GBM) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Codel)],1,mean,na.rm=T) #group n 
data$DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(LGm6_GBM) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
data$DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(LGm6_GBM) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
data$DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(LGm6_GBM) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
data$DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]



#**Alert**#
dim(data) #23493   
dim(data_sig) #22180   


# Hyper methylated probes
data_sig[data_sig$DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation > 0.2 &
           data_sig$DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation > 0.2 &
           data_sig$DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation > 0.2 &
           data_sig$DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation> 0.2 &
           data_sig$DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation > 0.2, ] %>% dim() # 0 probes
library(UpSetR)
listInput_hyper <- list(DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation > 0.3, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 


# Hypo methylated probes
data_sig[data_sig$DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation < -0.2 &
           data_sig$DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation < -0.2 &
           data_sig$DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation < -0.2 &
           data_sig$DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation < -0.2 &
           data_sig$DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation < -0.2, ] %>% dim() 

library(UpSetR)
listInput_hypo <- list(DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation < -0.3, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5) 


# Extract hypo methylated probes 
x <- upset(fromList(listInput_hypo), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #2661 probes somehow overlapped
x1 <- unlist(listInput_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #2661 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # all probes present it only 5/5 comparisons
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 comparisons

length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret.

LGm6_probeset_45p <- c(x1[ rowSums(x$New_data) == 5],
                        x1[ rowSums(x$New_data) == 4])
# hypo probes; diffmean < -0.4; pvalue < 0.05; 5/5 and 4/5 comparisons had this probeset



# Mesenchymal_like Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(Mesenchymal_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Classic_like)],1,mean,na.rm=T) #group n 
data$DiffMean_Mesenchymal_Classic_MesenchymaOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Mesenchymal_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Codel)],1,mean,na.rm=T) #group n 
data$DiffMean_Mesenchymal_Codel_MesenchymaOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Mesenchymal_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
data$DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Mesenchymal_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
data$DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Mesenchymal_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
data$DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]


#**Alert**#
dim(data) #23493   
dim(data_sig) #22180   

# Hyper methylated probes
data_sig[data_sig$DiffMean_Mesenchymal_Classic_MesenchymaOrientation > 0.3 &
           data_sig$DiffMean_Mesenchymal_Codel_MesenchymaOrientation > 0.3 &
           data_sig$DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation> 0.3 &
           data_sig$DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation > 0.3 &
           data_sig$DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation, ] %>% dim() # 0 probes
library(UpSetR)
listInput_hyper <- list(DiffMean_Mesenchymal_Classic_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_Classic_MesenchymaOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_Mesenchymal_Codel_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_Codel_MesenchymaOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation > 0.2, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 

# Hypo methylated probes
data_sig[data_sig$DiffMean_Mesenchymal_Classic_MesenchymaOrientation < -0.2 &
           data_sig$DiffMean_Mesenchymal_Codel_MesenchymaOrientation < -0.2 &
           data_sig$DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation < -0.2 &
           data_sig$DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation < -0.2 &
           data_sig$DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation, ] %>% dim() # 19 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_Mesenchymal_Classic_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_Classic_MesenchymaOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_Mesenchymal_Codel_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_Codel_MesenchymaOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation < -0.2, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5) 



# Extract hyper methylated probes 
x <- upset(fromList(listInput_hyper), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #327 probes somehow overlapped
x1 <- unlist(listInput_hyper, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #327 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # none
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 comparisons

hyper_probes <- c(x1[ rowSums(x$New_data) == 4])

# Extract hypo methylated probes 
x <- upset(fromList(listInput_hypo), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #4370 probes somehow overlapped
x1 <- unlist(listInput_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #4370 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # none
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 comparisons

# length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret. - not necessary 

hypo_probes <- c(x1[ rowSums(x$New_data) == 4])

Mesenchymal_probeset_47p <- c(hyper_probes, hypo_probes)
# hyoer probes; diffmean > 0.2; pvalue < 0.05; 4/5 comparisons had this probeset
# hypo probes; diffmean < -0.2; pvalue < 0.05; 4/5 comparisons had this probeset




DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list <- list(Classic_probeset_36p = Classic_probeset_36p,
                                                                         Codel_probeset_10p = Codel_probeset_10p,
                                                                         GcimpHigh_probeset_18p= GcimpHigh_probeset_18p,
                                                                         GcimpLow_probeset_22p = GcimpLow_probeset_22p,
                                                                         LGm6_probeset_45p = LGm6_probeset_45p,
                                                                         Mesenchymal_probeset_47p = Mesenchymal_probeset_47p)

saveRDS(DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list, file = '/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/DMP_KNN_StringestThresholdanova_diffmean_subtypes_probeset_list.rds')



tune_probset <- c(DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[1]],
                  DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[2]],
                  DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[3]],
                  DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[4]],
                  DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[5]],
                  DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[6]])

tune_probset_unique <- unique(tune_probset) # we dont need to use that 'upset strategy' to get only unique probes. We just need to remove duplicates like this.


# Heatmap visualization ---------------

my_sample_col <- data.frame(row.names = metadata$Case, Subtypes_six = metadata$Ceccarelli_six_subtypes, Subtypes_seven = metadata$Supervised.DNA.Methylation.Cluster, DNAmet_cluster = metadata$IDH.specific.DNA.Methylation.Cluster) # patient ID in row names; any columns are for sample information

DNAmtx <- data_sig[, colnames(data_sig) %in% metadata$Case] #keep only beta-values

library(ComplexHeatmap)
library(matlab)
#label 'my_sample_col'
#                    Subtypes_six   Subtypes_seven DNAmet_cluster
# TCGA-CS-4938      G-CIMP-high      G-CIMP-high      IDHmut-K2
# TCGA-CS-4941 Mesenchymal-like Mesenchymal-like       IDHwt-K2
# TCGA-CS-4942      G-CIMP-high      G-CIMP-high      IDHmut-K2
# TCGA-CS-4943      G-CIMP-high      G-CIMP-high      IDHmut-K2
# TCGA-CS-4944      G-CIMP-high      G-CIMP-high      IDHmut-K2
# TCGA-CS-5390            Codel            Codel      IDHmut-K3
# Current label

### Pallete 
library(viridis)
n_colors <- length(table(my_sample_col$Subtypes_six))
Subtypes_six_pal <- viridis(n = n_colors, option = "plasma", direction = -1)
n_colors <- length(table(my_sample_col$Subtypes_seven))
Subtypes_seven_pal <- viridis(n = n_colors, option = "inferno", direction = -1)
n_colors <- length(table(my_sample_col$DNAmet_cluster))
DNAmet_cluster_pal <- viridis(n = n_colors, option = "cividis", direction = -1)

Subtypes_six_lvs <- names(table(my_sample_col$Subtypes_six))
Subtypes_seven_lvs <- names(table(my_sample_col$Subtypes_seven))
DNAmet_cluster_lvs <- names(table(my_sample_col$DNAmet_cluster))

# Name yout columns and it's levels
column_name_1 <- "Subtypes_six"
levels_1 <- Subtypes_six_lvs

column_name_2 <- "Subtypes_seven"
levels_2 <- Subtypes_seven_lvs

column_name_3 <- "DNAmet_cluster"
levels_3 <- DNAmet_cluster_lvs

# Generate color vector
colors_1 <- Subtypes_six_pal
colors_2 <- Subtypes_seven_pal
colors_3 <- DNAmet_cluster_pal
# Assign names to color vector
names(colors_1) <- levels_1
names(colors_2) <- levels_2
names(colors_3) <- levels_3
# Create named list
color_mapping <- list(column_name_1 = colors_1, 
                      column_name_2 = colors_2,
                      column_name_3 = colors_3)


top.anno = HeatmapAnnotation(df = my_sample_col, 
                             col= color_mapping,
                             show_annotation_name = T, annotation_name_gp = gpar(fontsize=7),
                             na_col= "white")



hm_unique_probset <- Heatmap(as.matrix(DNAmtx[rownames(DNAmtx) %in% tune_probset_unique, ]),
                         cluster_columns=T,
                         cluster_rows = T,
                         clustering_method_rows="complete",  
                         show_row_names = F,
                         row_names_gp = gpar(fontsize = 7),
                         show_column_names = F,
                         name = "CpG probes methylation",
                         col= jet.colors(75),
                         row_title = "CpG probes",
                         #row_names_gp = gpar(fontsize = 12),
                         #column_title = "Glioma subtypes tuned probes",
                         #column_title = "Glioma subtypes unique tuned probes",
                         #column_title = "Glioma subtypes unique probset (knn + stringest DiffMean)",
                         column_title = "Glioma subtypes unique probset",
                         #split=label.EMT.genes$E_ou_M,
                         #row_names_side="left",
                         
                         top_annotation = top.anno,
                         #row_names_gp = gpar(fontsize = 8),
                         #column_names_gp = gpar(fontsize = 12),
                         heatmap_legend_param = list(
                           color_bar = 'continuous',
                           legend_direction = 'vertical',
                           legend_width = unit(12, 'cm'),
                           legend_height = unit(10, 'cm'),
                           title_position = 'leftcenter-rot',
                           title_gp=gpar(fontsize = 16, fontface = 'bold'),
                           labels_gp=gpar(fontsize = 16, fontface = 'bold')))

print(hm_unique_probset)

### For left side annotation
# RowAnn <- HeatmapAnnotation(df = label, col=list("Epi.Mes.CellC" = c("Epi" = "blue", "Mes" = "black", "c.Cycle" = "cadetblue1")),
#                             show_annotation_name = T, annotation_name_gp = gpar(fontsize=7),
#                             na_col= "white",which = "row",show_legend =T)


### For left side annotation
# draw(hm_EMT.cCycle.MALTA + RowAnn) 

pdf("/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/knn_stringestDiffmean_unique_probeset.pdf",width = 5, height = 10)
### For left side annotation
# draw(hm_EMT.cCycle.MALTA, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()




# Silhouette coef. for a clustering comparison metric -----------
library(cluster)
library(factoextra)
# Didn't worked out ...


# Cifarp Isabela abstract objective information ---------------
# Probeset intersection 
pan_glioma_probes <- read.xlsx("/media/hd/maycon/Glioma_classifier/PanGlioma_MethylationSignatures.xlsx", sheet = 1) 
pan_glioma_probes <- pan_glioma_probes$`1,300.pan-glioma.tumor.specific.probes.(Figure.2A)`

DMP_KNN_StringestThresholdanova_diffmean_subtypes_probeset_list <- readRDS("/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/DMP_KNN_StringestThresholdanova_diffmean_subtypes_probeset_list.rds")
tune_probset <- c(DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[1]],
                  DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[2]],
                  DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[3]],
                  DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[4]],
                  DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[5]],
                  DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[6]])

tune_probset_unique <- unique(tune_probset)


length(pan_glioma_probes) #1300
length(tune_probset_unique) #143
length(intersect(pan_glioma_probes,
                 tune_probset_unique)) #63 


# Catching how many gliomas IDH mut have been miss clustered next to the IDHwt in our probset *143probes

my_sample_col$index <- 1:dim(my_sample_col)[1]
my_sample_col$Case <- rownames(my_sample_col)
rownames(my_sample_col) <- my_sample_col$index
my_sample_col[column_order(hm_unique_probset), ][700:932, ] #part of the dendogram on heatmap which I've seen gliomas IDHmut next to IDHwt
# 12 G-CIMP-high / 249 (4.8%)
# 12 codel / 174 (6.8 %)








