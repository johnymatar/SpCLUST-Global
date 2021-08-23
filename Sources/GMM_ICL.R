suppressMessages(library(mclust))
args = commandArgs(trailingOnly=TRUE)
#Read the features matrix from the provided text file
features <- data.matrix(read.table(paste("featuresMat",args[1],".txt",sep=""), sep=" "))
#Remove the last column containing NA
features <- features[,-ncol(features)]
#Set the scope of the potential number of clusters to 1:nbElements
G = 1:nrow(features)
#Set the model names to use: the ones allowing unequal volumes
modelNames = c("VII", "VEI", "VVI", "VEE", "VVE", "VEV", "VVV")
#Get the model and the number of clusters scoring the best ICL
features.ICL <- mclustICL(features, G, modelNames)
features.ICL[is.na(features.ICL)] <- 0
nbClusters <- which(features.ICL == max(features.ICL), arr.ind = TRUE)[1]
G = nbClusters:nbClusters
modelName <- modelNames[which(features.ICL == max(features.ICL), arr.ind = TRUE)[2]]
#Get the clustering labels
features.M <- Mclust(features, G, modelName)
cat(features.M$classification, sep=" ")
#summary(features.M, parameters = TRUE)