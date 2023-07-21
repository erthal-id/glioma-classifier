### ### ### ### ### ### ####
### IDHmut 163 probe set ###


# Filter IDHmut gliomas and it's probes previous seleceted by Ceccarelli
# Call this object 'RFtrain'
IDHmut_mtx <- DNAmet_mtx_sub[, colnames(DNAmet_mtx_sub) %in% metadata[metadata$IDH.status %in% 'Mutant', ]$Case] 
DefineEacTC <- read.xlsx("/media/hd/maycon/Glioma_classifier/PanGlioma_MethylationSignatures.xlsx", sheet = 5) 
IDH_mut_subtypes_probes <- DefineEacTC$`163.probes.that.define.each.TCGA.IDH-mutant.glioma.subtype.(149.probes.showed.in.Figure.S3C)`
RFtrain = IDHmut_mtx[rownames(IDHmut_mtx) %in% IDH_mut_subtypes_probes, ]


# Transpose DNAmtx to keep probes in the columns and patients in the rows
# Add patient 'barcodes' and 'DNA cluster label' to the transposed matrix 
# Call this object 'trainingdata'

# Remove NAs if its necessary 
table(is.na(RFtrain)) 
# FALSE  TRUE 
# 71092     8
library(impute)
RFtrain <- impute.knn(as.matrix(RFtrain), k = 10, rowmax = 0.8, colmax = 0.8, maxp = 1500, rng.seed=362436069)[[1]]
table(is.na(RFtrain)) 
# FALSE 
# 71100

trainingdata <- t(RFtrain)
# Note: I'm changing the column to 'Supervised.DNA.Methylation.Cluster'
trainingdata <- merge(trainingdata, metadata[,c("Case", "Supervised.DNA.Methylation.Cluster")], by.x=0,by.y="Case", all.x=T) 
rownames(trainingdata) <- as.character(trainingdata$Row.names)
trainingdata <- trainingdata[,-1] 
trainingdata <- droplevels(trainingdata)

# Remove NAs if its necessary (Checking it again)
table(is.na(trainingdata$Supervised.DNA.Methylation.Cluster)) # two samples without label. Remove it in the meanwhile
trainingdata <- trainingdata[!trainingdata$Supervised.DNA.Methylation.Cluster %in% NA, ]

table(is.na(trainingdata))
# FALSE 
# 71232

save(trainingdata, file="/media/hd/maycon/Glioma_classifier/trainingdata_IDHmut_163pset.Rda")


# In order to test the code, let only 50 samples go through the RF
set.seed(666)
trainingdata_test <- createDataPartition(trainingdata$Supervised.DNA.Methylation.Cluster, p=0.15, list=FALSE, times=1)
trainingdata_test <- trainingdata[trainingdata_test,]

table(trainingdata_test$Supervised.DNA.Methylation.Cluster)
# Codel G-CIMP-high  G-CIMP-low 
# 27    38           4 

# Add samples balancing the other groups to be a fair about n of samples
# Getting the same samples used on the model within 1308 probe set (gcimp_low_keep, gcimp_high_keep, and codel_keep)

trainingdata_test = trainingdata[rownames(trainingdata) %in% 
                                   c(gcimp_low_keep,
                                     gcimp_high_keep,
                                     codel_keep), ]

table(trainingdata_test$Supervised.DNA.Methylation.Cluster)
# Codel G-CIMP-high  G-CIMP-low 
# 15          15          15 

trainingdata_45smp <- trainingdata_test
save(trainingdata_45smp, file="/media/hd/maycon/Glioma_classifier/trainingdata_45smp_IDHmut_163pset.Rda")


### Tathi RF code - Run it on terminal -----------------
start_time <- Sys.time()

library(caret)
library(randomForest)
library(doMC)
library(e1071)
load("/media/hd/maycon/Glioma_classifier/trainingdata_45smp_IDHmut_163pset.Rda")
trainingdata <- trainingdata_45smp
# register cores for doMC
registerDoMC(cores = 4) #tathi ran 10
# set up k-fold cross validation
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)
# you may additionally, if you wish use a different method, for validating your
# model parameters, such as oob (Out of Bag).  oob is faster.

# Set your seed so your work is repeatable
set.seed(42)
# Create a subset of your data to train your model on.  This makes sure you have
# equal representation of the 'papercluster' groups in your training set
inTraining <- createDataPartition(trainingdata$Supervised.DNA.Methylation.Cluster, p=0.8, list=FALSE, times=1)
# Training Set
myTrain <- trainingdata[inTraining, ] #checar se tem membros de todos os grupos #table(myTrain$Supervised.DNA.Methylation.Cluster)
# Testing Set
myTest <- trainingdata[-inTraining, ] #checar se tem membros de todos os grupos #table(myTest$Supervised.DNA.Methylation.Cluster)
# Confirm seed is set
set.seed(210)
# set values for mtry
# mtry is the "Number of variables randomly sampled as candidates at each split"
# traditionally for classification you use the sqrt of the number of variables
# but here we try a range of mtry values to find the best parameters for our model
mtryVals <- floor(c(seq(100, 2000, by=100),
                    sqrt(ncol(trainingdata))))
mtryGrid <- data.frame(.mtry=mtryVals)
# Confirm seed again
set.seed(420)
# Set number of cores
registerDoMC(cores = 10)
# Run Training
stemsig <- train(Supervised.DNA.Methylation.Cluster ~ ., # variable to be trained on
                 data = trainingdata, # Data we are using
                 method = "rf", # Method we are using
                 trControl = fitControl, # How we validate
                 # We created this object above
                 ntree = 5000, # number of trees
                 # is dependent on training data size
                 importance = TRUE, # calculate varible importance
                 # can be omitted to speed up calc
                 tuneGrid = mtryGrid # set mtrys
                 #subset = myTrain # define training set #comment when train with 100% of samples
)
end_time <- Sys.time()
runtime = end_time - start_time
print(runtime)
save(list=ls(),file="/media/hd/maycon/Glioma_classifier/RF_Mut_45smp_100perc_163pset.Rda")
# Start: 
# End:  (41.70936 mins)
#save(list=ls(),file="/media/hd/maycon/Glioma_classifier/RF_Mut_all450smp_100perc.Rda")


# back to RStudio

# to be continued

# Compare RF model on 'train set' vs 'test set' ---------------
load('/media/hd/maycon/Glioma_classifier/RF_Mut_45smp_100perc_163pset.Rda')
myTrain <- trainingdata # just because we used all 45 samples in the model
# stemsig: RF model trained on 45 IDHmut samples (model output)
# myTrain: DNAmatrix, IDHmut probes (163 set), Supervised.DNA.Methylation.Cluster column, the 45 IDHmut samples in the model

table(myTrain$Supervised.DNA.Methylation.Cluster)
# Codel    G-CIMP-high  G-CIMP-low 
# 15/174      15/249     15/25


# Testing against all 387 IDHwt samples left out of the training/testing phase -------
load("/media/hd/maycon/Glioma_classifier/trainingdata_IDHmut_163pset.Rda")
myTest_2 = trainingdata[!rownames(trainingdata) %in% rownames(myTrain), ]

train.pred.2 <- predict(stemsig, myTest_2, type="prob")
myTest_2$RFtype <- predict(stemsig, myTest_2)
myTest_2 <- droplevels(myTest_2)
confusionMatrix(data = as.factor(myTest_2$RFtype), 
                reference = as.factor(myTest_2$Supervised.DNA.Methylation.Cluster))

# accuracy = 0.9653
#               Reference
# Prediction    Codel G-CIMP-high G-CIMP-low
# Codel         155           5          0
# G-CIMP-high     3         224          0
# G-CIMP-low      1           5         10
