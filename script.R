
# KRISTEN'S MODELING FRAMEWORK USING THE ROOT METABOLITES
# Microsoft Open R 3.3.0
library(checkpoint)
checkpoint("2015-12-12", checkpointLocation = getwd())

library(caret)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(pROC)
library(ROCR)
library(pls)
library(kernlab)
library(randomForest)
library(e1071)


# Subset the data
rootmet <- read.delim("rootData.txt", row.names = 1)

# Remove highly correlated preds
corrMetabs <- cor(rootmet)
filterCor <- findCorrelation(corrMetabs, cutoff = .85)
rootmet <- rootmet[,-filterCor]

allNames <- rownames(rootmet)
myVec <- vector()
tmp <- strsplit(split = ".", fixed = T, x = allNames)
for (i in 1:length(allNames)){
      
      myVec <- c(myVec, tmp[[i]][1])
      
      
}
rootmet <- t(rootmet)
colnames(rootmet) <- myVec

# Get the unique hybrid genotypes
hNames <- unique(colnames(rootmet)[grep("x", fixed = T, colnames(rootmet))])

index <- vector()
for(i in 1:length(hNames)){
      par <- strsplit(x = hNames[i], split = "x", fixed = T)
      if(par[[1]][1] %in% colnames(rootmet) & par[[1]][2] %in% colnames(rootmet)) {
            index <- c(index, TRUE)
      }else{
            index <- c(index, FALSE)
      }
}
# Update the list to stick only to the hybrids that have both parents in the data set
hNames <- hNames[index]

#List of Hybrids and Parents + replicates
rootList <- list()
for(h in 1:length(hNames)){
      genotyp <- as.matrix(rootmet[,colnames(rootmet) %in% hNames[h]])
      colnames(genotyp) <- rep("H", ncol(genotyp))
      # Get names of the parents
      parents <- strsplit(x = hNames[h], split = "x", fixed = T)
      PA <- as.matrix(rootmet[,colnames(rootmet) %in% parents[[1]][1]])   
      colnames(PA) <- rep("PA", ncol(PA))
      PB <- as.matrix(rootmet[,colnames(rootmet) %in% parents[[1]][2]])
      colnames(PB) <- rep("PB", ncol(PB))
      trio <- cbind(genotyp, PA, PB) #transpose row-bind
      rootList[[h]] <- trio
}
names(rootList) <- hNames
### PREPARE DESIGN ###
designList <- list()
for(i in 1:length(rootList)) {
      designList[[i]] <- table(colnames(rootList[[i]]))
}

## Remove families with a single replicate in at least one of the groups
#because it does not fulfill basis for testing

idx <- !(unlist(lapply(designList, function(x) {1 %in% x})))
designList <- designList[idx]
rootList <- rootList[idx]
hNames <- hNames[idx]


# Create design
design <- list()
for(i in 1:length(rootList)) {
      tmp2 <- model.matrix(~ 0 + factor(rep(0:2, as.numeric(designList[[i]]))))
      colnames(tmp2) <- c("H", "PA", "PB")
      design[[i]] <- tmp2
}

# Create contrast.matrix
contrast.matrix <- list()
for(i in 1:length(rootList)) {
      tmp3 <- makeContrasts(PA - H, PA - PB, PB - H, levels = design[[i]])
      contrast.matrix[[i]] <- tmp3
}


lmfitList <- list()
for(i in 1:length(rootList)) {
      fit <- lmFit(rootList[[i]], design[[i]])
      fit2 <- contrasts.fit(fit, contrast.matrix[[i]])
      fit2 <- eBayes(fit2)
      results <- decideTests(fit2, adjust.method = "none")
      lmfitList[[i]] <- list(fit = fit2, results = results)
}

# Extract results (are the differences between groups significant according to the moderate t-tests?)
finalMatrix <- matrix(nrow = nrow(rootmet), ncol = length(rootList))
colnames(finalMatrix) <- hNames
rownames(finalMatrix) <- rownames(rootmet)
for (i in 1:length(rootList)){
      tmp <- lmfitList[[i]]$results
      store <- vector()
      for(r in 1:nrow(tmp)) {
            if(all(tmp[r,] == c(0, 0, 0)) | all(tmp[r,] == c(0, 1, 0)) | all(tmp[r,] == c(0, -1, 0)) | all(tmp[r,] == c(1, 1, -1))
               | all(tmp[r,] == c(-1, -1, 1))){
                  store <- c(store, "A")
            }else if(all(tmp[r,] == c(1, 1, 0)) | all(tmp[r,] == c(0, -1, 1))){
                  store <- c(store, "-D")
            }else if(all(tmp[r,] == c(1, 0, 1)) | all(tmp[r,] == c(1, 1, 1)) | all(tmp[r,] == c(1, -1, 1)) | all(tmp[r,] == c(1, 0, 0))
                     | all(tmp[r,] == c(0, 0, 1))){
                  store <- c(store, "-OD")
            }else if(all(tmp[r,] == c(0, 1, -1)) | all(tmp[r,] == c(-1, -1, 0))){
                  store <- c(store, "+D")
            }else if(all(tmp[r,] == c(-1, 0, -1)) | all(tmp[r,] == c(-1, 1, -1)) | all(tmp[r,] == c(-1, -1, -1)) | all(tmp[r,] == c(-1, 0, 0))
                     | all(tmp[r,] == c(0, 0, -1))){
                  store <- c(store, "+OD")
            }
      }
      finalMatrix[,i] <- store
}
## Visualize

hybridClasses <- matrix(nrow = nrow(finalMatrix), ncol = ncol(finalMatrix), 0)
rownames(hybridClasses) <- rownames(finalMatrix)
colnames(hybridClasses) <- colnames(finalMatrix)
hybridClasses[finalMatrix == "-OD"] <- -2
hybridClasses[finalMatrix == "-D"] <- -1
hybridClasses[finalMatrix == "A"] <- 0
hybridClasses[finalMatrix == "+D"] <- 1
hybridClasses[finalMatrix == "+OD"] <- 2
# Plot a large heatmap
pdf("hybridClasses1.pdf", height = 8, width = 10)
pheatmap(hybridClasses, cluster_cols = F, cluster_rows = F,
         legend_breaks = (-2:2)*.8,
         legend_labels = c("-OD", "-D", "A", "+D", "+OD"), show_colnames = F,
         fontsize_row = 5,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(5))
dev.off()
colnames(hybridClasses) <- names(rootList)

#### Classification ####

#create Matrix with all Parents
allParents <- matrix(nrow = nrow(rootmet), ncol = 0)
for(e in 1:length(rootList)) {
      PAs <- rootList[[e]][,colnames(rootList[[e]])=="PA"]
      PBs <- rootList[[e]][,colnames(rootList[[e]])=="PB"]
      Pa <- apply(PAs, 1, median)
      Pb <- apply(PBs, 1, median)
      PAPBmed <- as.data.frame(cbind(Pa,Pb))
      names(PAPBmed) <- paste(c("Pa","Pb"), names(rootList)[e], sep = ".")
      allParents <- cbind(allParents, PAPBmed)
}

## Class balance
# Remove metabolites (outcome) with the same class label appearing more than x %
maxRelFrequency <- .75
filter <- apply(hybridClasses, 1, function(x) {if(any((as.vector(table(x)) / sum(as.vector(table(x)))) >= maxRelFrequency)) {return(FALSE)}else{return(TRUE)}})
hybridClasses <- hybridClasses[filter,]

# Cast out metabolites (outcome) that exhibit the same class(classes) only once (each)
outcomeList <- list()
for(c in 1:nrow(hybridClasses)) {
      outcomeList[[c]] <- factor(hybridClasses[c,]) # one outcome per metabolite
}
names(outcomeList) <- rownames(hybridClasses)
freqTable <- lapply(outcomeList, table)
unbalance <- unlist(lapply(freqTable, function(x)
      {sum(x == 1) == 1 & sum(x < 4) == 1 | sum(x == 1) == 0 & sum(x < 4) == 0}))
hybridClasses <- hybridClasses[unbalance,]

# Plot small heatmap
pdf("Figure2.pdf", height = 3, width = 10)
pheatmap(hybridClasses, cluster_cols = F, cluster_rows = F,
         legend_breaks = (-2:2)*.8,
         legend_labels = c("-OD", "-D", "A", "+D", "+OD"), show_colnames = F,
         fontsize_row = 5,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(5))
dev.off()

# Plot class proportion per metab
proportion <- matrix(NA, nrow = 5, ncol = nrow(hybridClasses))
dimnames(proportion) <- list(c("+OD","+D","A","-D","-OD"),
                             rownames(hybridClasses))
for(i in 1:nrow(hybridClasses)){
      proportion[,i] <- rev(table(factor(hybridClasses[i,], levels = -2:2))/length(hybridClasses[i,]))
           
}
pdf("ClassProportion.pdf", height = 3, width = 8)
pheatmap(proportion, scale = "none", cluster_rows = F, cluster_cols = F,
         color = colorRampPalette(brewer.pal(n = 7,name="Greys"))(100),
         breaks = seq(0,1,by=0.01), fontsize_col = 5)
dev.off()

## Organize the data for the modeling, concatenate Pa/Pb * metabolites
predictorSet <- matrix(nrow = ncol(hybridClasses), ncol = nrow(allParents) * 2)
rownames(predictorSet) <- colnames(hybridClasses)
for(n in 1:ncol(hybridClasses)) {
      tmp <- allParents[,grep(colnames(hybridClasses)[n], fixed = T, colnames(allParents))]
      conc <- c(tmp[,1], tmp[,2])
      predictorSet[n,] <- conc
}
colnames(predictorSet) <- paste(rep(c("Pa", "Pb"), each = nrow(tmp)),
                                rep(rownames(tmp), 2), sep = ".")
rownames(predictorSet) <- colnames(hybridClasses)

# Scale the predictor set
scalePredictorSet <- scale(predictorSet, center = T, scale = T)

#### Function for PLS-RF by Max Kuhn ####
PLSRF <- list(label = "PLS-RF",
                  library = c("pls", "randomForest"),
                  type = "Classification",
                  parameters = data.frame(parameter = c('ncomp', 'mtry'),
                                          class = c("numeric", 'numeric'),
                                          label = c('#Components', 
                                                    '#Randomly Selected Predictors')),
                  grid = function(x, y, len = NULL) {
                        grid <- expand.grid(ncomp = seq(1, min(ncol(x) - 1, len), by = 1),
                                            mtry = 1:len)
                        grid <- subset(grid, mtry <= ncomp)
                  },
                  loop = NULL,
                  fit = function(x, y, wts, param, lev, last, classProbs, ...) { 
                        ## First fit the pls model, generate the training set scores,
                        ## then attach what is needed to the random forest object to 
                        ## be used later
                        pre <- plsda(x, y, ncomp = param$ncomp)
                        scores <- pls:::predict.mvr(pre, x, type = "scores")
                        mod <- randomForest(scores, y, mtry = param$mtry, ...)
                        mod$projection <- pre$projection
                        mod
                  },
                  predict = function(modelFit, newdata, submodels = NULL) {       
                        scores <- as.matrix(newdata)  %*% modelFit$projection
                        predict(modelFit, scores)
                  },
                  prob = NULL,
                  varImp = NULL,
                  predictors = function(x, ...) rownames(x$projection),
                  levels = function(x) x$obsLevels,
                  sort = function(x) x[order(x[,1]),])
#### ####
# Wrapper function for classification
hybridClassification <- function(m, perm = F, seed) {
      X <- scalePredictorSet
      Y <- hybridClasses[m,]
      frqTbl <- table(Y)
      if(sum(frqTbl == 1) == 1){
            class <- names(which(frqTbl == 1))
            idx <- which(Y == class)
            X <- X[-idx,]
            Y <- Y[-idx]
      }
      Y <- factor(Y)
      set.seed(seed) # Set seed
      split <- createDataPartition(Y, p = .75)[[1]]
      training <- X[split,]
      test <- X[-split,]
      trainY <- Y[split]
      testY <- Y[-split]
      if(perm == T){
            set.seed(seed) # Set seed
            testY <- sample(testY)
            set.seed(seed) # Set seed
            trainY <- sample(trainY)
      }
      set.seed(seed) # Set seed
      stratFolds <- createMultiFolds(trainY, k = 3, times = 5)
      ctrl <- trainControl(method = "repeatedcv", index = stratFolds)
      
      # Linear models
      # PLS-DA
      plsdaFit <- train(x = training, y = trainY,
                      method = "pls",
                      metric = "Kappa",
                      trControl = ctrl,
                      tuneGrid = expand.grid(.ncomp = 1:30))
      plsdaVal <- predict(plsdaFit, newdata = test)
      plsdaSummary <- defaultSummary(data.frame("obs" = testY, "pred" = plsdaVal))
      
      # Elastic Net
      glmnetGrid <- expand.grid(.alpha = seq(0,1,by=0.1),
                                .lambda = seq(0.01, .2, length = 40))
      glmFit <- train(x = training, y = trainY,
                        method = "glmnet",
                        metric = "Kappa",
                        trControl = ctrl,
                        tuneGrid = glmnetGrid)
      glmVal <- predict(glmFit, newdata = test)
      glmSummary <- defaultSummary(data.frame("obs" = testY, "pred" = glmVal))
      
      
      # Nonlinear models
      # PLS-RF
      plsrfGrid <- expand.grid(.ncomp = 5:20,
                               .mtry = round(seq(5,20, length.out = 3)))
      PLSRFFit <- train(x = training, y = trainY,
                      method = PLSRF,
                      metric = "Kappa",
                      trControl = ctrl,
                      tuneGrid = plsrfGrid)
      PLSRFVal <- predict(PLSRFFit, newdata = test)
      plsrfSummary <- defaultSummary(data.frame("obs" = testY, "pred" = PLSRFVal))
      
      # SVM
    
      svmFit <- train(x = training, y = trainY,
                      method = "svmRadial",
                      metric = "Kappa",
                      trControl = ctrl,
                      tuneLength = 15)
      svmVal <- predict(svmFit, newdata = test)
      svmSummary <- defaultSummary(data.frame("obs" = testY, "pred" = svmVal))
      
      # SVM weighted
      classes <- levels(trainY)
      allClasses <- c("-2" = 10, "-1" = 5, "0" = 1, "1" = 5, "2" = 10)
      relClasses <- allClasses[classes]
      svmWFit <- train(x = training, y = trainY,
                       method = "svmRadial",
                       metric = "Kappa",
                       trControl = ctrl,
                       tuneLength = 15,
                       class.weights = relClasses)
      svmWVal <- predict(svmWFit, newdata = test)
      svmWSummary <- defaultSummary(data.frame("obs" = testY, "pred" = svmWVal))

      # RF
      rfFit <- train(x = training, y = trainY,
                     method = "rf",
                     metric = "Kappa",
                     trControl = ctrl,
                     tuneGrid = expand.grid(.mtry = c(100, 150, 200, 250)),
                     ntree = 1000)
      rfVal <- predict(rfFit, newdata = test)
      rfSummary <- defaultSummary(data.frame("obs" = testY, "pred" = rfVal))
      
      # RF weighted
      rfWFit <- train(x = training, y = trainY,
                      method = "rf",
                      metric = "Kappa",
                      trControl = ctrl,
                      tuneGrid = expand.grid(.mtry = c(100, 150, 200, 250)),
                     ntree = 1000,
                     class.weights= relClasses)
      rfWVal <- predict(rfWFit, newdata = test)
      rfWSummary <- defaultSummary(data.frame("obs" = testY, "pred" = rfWVal))
      return(rbind(plsdaSummary,glmSummary,plsrfSummary,svmSummary,svmWSummary,
                   rfSummary, rfWSummary))
}


metabolitesResults <- list()
for(i in 1:nrow(hybridClasses)) {
      for(j in 1:5){
            metabolitesResults[[j+((i-1)*5)]] <- hybridClassification(i, seed = j)
      }
}


#save(metabolitesResults, file = "classifResults.RData")
load("classifResults.RData")

# Now the permuted version of the classes
metabolitesResultsPerm <- list()
for(i in 1:nrow(hybridClasses)) {
      for(j in 1:5){
            metabolitesResultsPerm[[j+((i-1)*5)]] <- hybridClassification(i, seed = j, perm = T)
      }
}


#save(metabolitesResultsPerm, file = "classifResultsPerm.RData")
load("classifResultsPerm.RData")

# Plot
numberModels <- 7
observedMean <- array(unlist(metabolitesResults), c(numberModels,2,nrow(hybridClasses)*5))
permMean <- array(unlist(metabolitesResultsPerm), c(numberModels,2,nrow(hybridClasses)*5))

pdf("Figure3.pdf", height = 4, width = 5)
modelNames <- c("PLS-DA", "glmnet", "PLS-RF",
                "SVM","SVM-W", "RF", "RF-W")
boxplot(cbind(t(observedMean[,2,]), t(permMean[,2,])),
        names = c(modelNames, paste("Perm", modelNames)), las = 2,
        col = rep(c("azure4", "firebrick2"), each = numberModels), cex.axis = 0.7, ylab = "Kappa")
abline(h = seq(-.2, .8, by=.2), col = "grey", lty = 2)
dev.off()

# Coefficients of variation
cv <- function(x) {return(sd(x)/mean(x))}
cvList <- apply(t(observedMean[,2,]), 2, cv)

# Kappa per analyte + mean and cv
kappaMat <- matrix(ncol = nrow(hybridClasses), nrow = numberModels)
dimnames(kappaMat) <- list(modelNames,rownames(hybridClasses))
cvMat <- matrix(ncol = nrow(hybridClasses), nrow = numberModels)
dimnames(cvMat) <- list(modelNames,rownames(hybridClasses))
for(i in 1:nrow(hybridClasses)){
      cvMat[,i] <- apply(observedMean[,2,(1:5)+(5*(i-1)), drop = T], 1, cv) 
      kappaMat[,i] <- apply(observedMean[,2,(1:5)+(5*(i-1)), drop = T], 1, mean) 
}
write.table(cvMat, file = "cvMat.txt", sep = "\t")
write.table(kappaMat, file = "kappaMat.txt", sep = "\t")
### Choose model w/ lowest cv
# SVM weighted

SVMClassification <- function(x, seed) {
      X <- scalePredictorSet
      Y <- hybridClasses[x,]
      frqTbl <- table(Y)
      if(sum(frqTbl == 1) == 1){
            class <- names(which(frqTbl == 1))
            idx <- which(Y == class)
            X <- X[-idx,]
            Y <- factor(Y[-idx])
      }
      Y <- factor(Y)
      allClasses <- c("-2" = 10, "-1" = 5, "0" = 1, "1" = 5, "2" = 10)
      relClasses <- allClasses[levels(Y)]
      set.seed(seed) # Set seed
      stratFolds <- createMultiFolds(Y, k = 3, times = 10)
      ctrl <- trainControl(method = "repeatedcv", index = stratFolds)
      # SVM
      svmFit <- train(x = X, y = Y,
                      method = "svmRadial",
                      metric = "Kappa",
                      trControl = ctrl,
                      tuneLength = 15,
                      class.weights = relClasses)
      res <- varImp(svmFit)$importance
      return(apply(res, 1, mean))
}

vip <- list()
for (m in 1:nrow(hybridClasses)){
  vip[[m]] <- SVMClassification(m, seed = m)
  message("Metabolite No. ", m)
}

array <- array(unlist(vip), c(ncol(scalePredictorSet), 1, nrow(hybridClasses)))
averageVIP <- apply(array, 1:2, median)
rownames(averageVIP) <- colnames(scalePredictorSet)

rankMetabs <- averageVIP[order(averageVIP, decreasing = T),]
rankMetabsNames <- names(rankMetabs)

#save(rankMetabsNames, file = "rankMetabsNames.RData")
load("rankMetabsNames.RData")

# Plot boxplot following the ranking
preRank <- t(as.matrix(array[,,]))
colnames(preRank) <- colnames(scalePredictorSet)
preRank <- preRank[,rev(rankMetabsNames)]
colBP <- unlist(lapply(strsplit(rev(rankMetabsNames), split = ".", fixed = T), "[", 1))
colBP[colBP == "Pa"] <- "red"
colBP[colBP == "Pb"] <- "grey"
pdf("Figure4.pdf", height = 4, width = 12)
boxplot(preRank, col = colBP, axes = F, ylab = "Relative importance to the model (%)",
        xlab = "Parental analytes")
axis(1, at = 1:272,tck = -0.02, labels = F)
axis(1, at=c(1, 91, 182, 272), tck = -0.05)
axis(2, at = seq(0,100, by = 20), las = 2)
axis(4, at = seq(0,100, by = 20), labels = F)
dev.off()

# Cor between paternal and maternal ranks
mother <- gsub(rankMetabsNames[grep("Pa", rankMetabsNames, fixed=T)], pattern = "Pa.", replacement = "")
father <- gsub(rankMetabsNames[grep("Pb", rankMetabsNames, fixed=T)], pattern = "Pb.", replacement = "")
cor.test(1:length(mother), match(mother, father), method = "kendall") # No significance

## SVR using the parental metabolites as predictors (from scaledPredictorSet)
# to predict field fresh weight (FW) of the hybrids
FW <- read.delim("biomassData.txt", row.names = 1)
idx <- intersect(rownames(scalePredictorSet), rownames(FW))
relevantBiomass <- FW[idx,"FW"]
relevantTrial <- FW[idx,"trial"]
newPredictorSet <- scalePredictorSet[idx,]
all(rownames(newPredictorSet) == rownames(FW[idx,]))

## SVR
SVRmodel <- function(listOfMetabs, perm = F, seed) {
      set.seed(seed) # Set seed
      split <- createDataPartition(y = relevantBiomass, p = 0.8)[[1]]
      train <- newPredictorSet[split,listOfMetabs]
      test <- newPredictorSet[-split,listOfMetabs]
      trainY <- relevantBiomass[split]
      testY <- relevantBiomass[-split]
      set.seed(seed) # Set seed
      stratFolds <- createMultiFolds(trainY, k = 3, times = 10)
      ctrl <- trainControl(method = "repeatedcv", index = stratFolds)
      if(perm == T) {
            trainY <- sample(trainY)
            testY <- sample(testY)
      }
      SRVFit <- train(x = train, y = trainY,
                        trControl = ctrl,
                        method = "svmRadial",
                        metric = "Rsquared",
                        tuneLength = 15)
      testprediction <- predict(SRVFit, newdata = test)
      result <- defaultSummary(data.frame("obs" = testY, "pred" = testprediction))
      return(result)
}

best5 <- data.frame(ncol = 2, nrow = 0)
best10 <- data.frame(ncol = 2, nrow = 0)
best20<- data.frame(ncol = 2, nrow = 0)
best50 <- data.frame(ncol = 2, nrow = 0)
best100 <- data.frame(ncol = 2, nrow = 0)
worst100 <- data.frame(ncol = 2, nrow = 0)
worst50<- data.frame(ncol = 2, nrow = 0)
worst20<- data.frame(ncol = 2, nrow = 0)
worst10<- data.frame(ncol = 2, nrow = 0)
worst5<- data.frame(ncol = 2, nrow = 0)
all <- data.frame(ncol = 2, nrow = 0)
allPerm <- data.frame(ncol = 2, nrow = 0)
random5 <- data.frame(ncol = 2, nrow = 0)


for(r in 1:25){
      best5 <- rbind(best5, SVRmodel(rankMetabsNames[1:5], seed = r))
      best10 <- rbind(best10, SVRmodel(rankMetabsNames[1:10][sample(1:10, 5)], seed = r))
      best20<- rbind(best20,SVRmodel(rankMetabsNames[1:20][sample(1:20, 5)], seed = r))
      best50 <- rbind(best50,SVRmodel(rankMetabsNames[1:50][sample(1:50, 5)], seed = r))
      best100 <- rbind(best100,SVRmodel(rankMetabsNames[1:100][sample(1:100, 5)], seed = r))
      worst100 <- rbind(worst100 ,SVRmodel(rev(rankMetabsNames)[1:100][sample(1:100, 5)], seed = r))
      worst50<- rbind(worst50,SVRmodel(rev(rankMetabsNames)[1:50][sample(1:50, 5)], seed = r))
      worst20<- rbind(worst20, SVRmodel(rev(rankMetabsNames)[1:20][sample(1:20, 5)], seed = r))
      worst10<- rbind(worst10, SVRmodel(rev(rankMetabsNames)[1:10][sample(1:10, 5)], seed = r))
      worst5 <- rbind(worst5, SVRmodel(rev(rankMetabsNames)[1:5], seed = r))
      all <- rbind(all, SVRmodel(rankMetabsNames, seed = r))
      allPerm <- rbind(allPerm, SVRmodel(rankMetabsNames, perm = T, seed = r))
      set.seed(r)
      random <- rankMetabsNames[sample(1:length(rankMetabsNames), 5)]
      random5 <- rbind(random5, SVRmodel(random, seed = r))
}

pdf("Figure5.pdf", height = 4.5, width = 5)
par(mar = c(10.1, 4.1, 4.1, 2.1))
boxplot(best5[-1,2], best10[-1,2],best20[-1,2],
             best50[-1,2],best100[-1,2],worst100[-1,2],
             worst50[-1,2],worst20[-1,2],worst10[-1,2],
             worst5[-1,2], random5[-1,2], all[-1,2],allPerm[-1,2],
        names = c("Top 5", "5 out of top 10", "5 out of top 20",
                  "5 out of top 50", "5 out of top 100", "5 out of bottom 100",
                  "5 out of bottom 50", "5 out of bottom 20", "5 out of bottom 10",
                  "Bottom 5", "Random 5", "All", "All (perm)"), las = 2,
        ylab = expression(paste(italic(R)^2)), cex.axis = 0.8)
abline(h = (7:0)/10, lty = 2, col = "grey")
dev.off()

### Assess influence of the metabolites from Pa compared to those from Pb
# Pa is DENT and Pb is FLINT
comparePs <- rep("Pa", length(rankMetabsNames))
comparePs[grep("Pb.", fixed = T, rankMetabsNames)] <- "Pb"
comparePs <- factor(comparePs)
countsPa <- vector()
countsPb <- vector()
for(i in 1:length(comparePs)){
      countsPa <- c(countsPa, table(comparePs[1:i])["Pa"]/table(comparePs)["Pa"])
      countsPb <- c(countsPb, table(comparePs[1:i])["Pb"]/table(comparePs)["Pb"])
      
}

pdf("FigureS5.pdf", height = 4, width = 4)
plot(x = 1:length(comparePs), y = countsPa, type = "l", col = "red",
     ylab = "Cumulative frequency", xlab = "Ranking", lwd = 2)
lines(x = 1:length(comparePs), y = countsPb, type = "l", col = "azure4")
legend("topleft", legend = c("Dent", "Flint"), text.col = c("red", "azure4"),
       bty = "n")
abline(a = 0, b = 1/length(comparePs), col= "grey")
axis(side = 4, at = seq(0,1, by = .2), labels = F)
dev.off()

## SVM (classification)
classes <- factor(ifelse(relevantBiomass < mean(relevantBiomass), "Bad","Good"))
SVMmodel <- function(listOfMetabs, perm = F, seed) {
      set.seed(seed) # Set seed
      split <- createDataPartition(y = classes, p = 0.8)[[1]]
      train <- newPredictorSet[split,listOfMetabs]
      test <- newPredictorSet[-split,listOfMetabs]
      trainY <- classes[split]
      testY <- classes[-split]
      set.seed(seed) # Set seed
      stratFolds <- createMultiFolds(trainY, k = 3, times = 10)
      ctrl <- trainControl(method = "repeatedcv", index = stratFolds)
      if(perm == T) {
            trainY <- sample(trainY)
            testY <- sample(testY)
      }
      SVMFit <- train(x = train, y = trainY,
                      trControl = ctrl,
                      method = "svmRadial",
                      metric = "Accuracy",
                      tuneLength = 15)
      testprediction <- predict(SVMFit, newdata = test)
      result <- defaultSummary(data.frame("obs" = testY, "pred" = testprediction))
      return(result)
}

best5Clas <- data.frame(ncol = 2, nrow = 0)
best10Clas <- data.frame(ncol = 2, nrow = 0)
best20Clas<- data.frame(ncol = 2, nrow = 0)
best50Clas <- data.frame(ncol = 2, nrow = 0)
best100Clas <- data.frame(ncol = 2, nrow = 0)
worst100Clas <- data.frame(ncol = 2, nrow = 0)
worst50Clas<- data.frame(ncol = 2, nrow = 0)
worst20Clas<- data.frame(ncol = 2, nrow = 0)
worst10Clas<- data.frame(ncol = 2, nrow = 0)
worst5Clas<- data.frame(ncol = 2, nrow = 0)
allClas <- data.frame(ncol = 2, nrow = 0)
allPermClas <- data.frame(ncol = 2, nrow = 0)
random5Clas <- data.frame(ncol = 2, nrow = 0)


for(r in 1:25){
      best5Clas <- rbind(best5Clas, SVMmodel(rankMetabsNames[1:5], seed = r))
      best10Clas <- rbind(best10Clas, SVMmodel(rankMetabsNames[1:10][sample(1:10, 5)], seed = r))
      best20Clas<- rbind(best20Clas,SVMmodel(rankMetabsNames[1:20][sample(1:20, 5)], seed = r))
      best50Clas <- rbind(best50Clas,SVMmodel(rankMetabsNames[1:50][sample(1:50, 5)], seed = r))
      best100Clas <- rbind(best100Clas,SVMmodel(rankMetabsNames[1:100][sample(1:100, 5)], seed = r))
      worst100Clas <- rbind(worst100Clas ,SVMmodel(rev(rankMetabsNames)[1:100][sample(1:100, 5)], seed = r))
      worst50Clas<- rbind(worst50Clas,SVMmodel(rev(rankMetabsNames)[1:50][sample(1:50, 5)], seed = r))
      worst20Clas<- rbind(worst20Clas, SVMmodel(rev(rankMetabsNames)[1:20][sample(1:20, 5)], seed = r))
      worst10Clas<- rbind(worst10Clas, SVMmodel(rev(rankMetabsNames)[1:10][sample(1:10, 5)], seed = r))
      worst5Clas <- rbind(worst5Clas, SVMmodel(rev(rankMetabsNames)[1:5], seed = r))
      allClas <- rbind(allClas, SVMmodel(rankMetabsNames, seed = r))
      allPermClas <- rbind(allPermClas, SVMmodel(rankMetabsNames, perm = T, seed = r))
      set.seed(r)
      random <- rankMetabsNames[sample(1:length(rankMetabsNames), 5)]
      random5Clas <- rbind(random5Clas, SVMmodel(random, seed = r))
}

pdf("FigureS2.pdf", height = 4.5, width = 5)
par(mar = c(10.1, 4.1, 4.1, 2.1))
boxplot(best5Clas[-1,1], best10Clas[-1,1],best20Clas[-1,1],
        best50Clas[-1,1],best100Clas[-1,1],worst100Clas[-1,1],
        worst50Clas[-1,1],worst20Clas[-1,1],worst10Clas[-1,1],
        worst5Clas[-1,1], random5Clas[-1,1], allClas[-1,1],allPermClas[-1,1],
        names = c("Top 5", "5 out of top 10", "5 out of top 20",
                  "5 out of top 50", "5 out of top 100", "5 out of bottom 100",
                  "5 out of bottom 50", "5 out of bottom 20", "5 out of bottom 10",
                  "Bottom 5", "Random 5", "All", "All (perm)"), las = 2,
        ylab = "Accuracy", cex.axis = 0.8)
abline(h = seq(4, 9)/10, lty = 2, col = "grey")
dev.off()

## SVM (classification, reporting AUC)
classes <- factor(ifelse(relevantBiomass < mean(relevantBiomass), "Bad","Good"))
SVMmodel <- function(listOfMetabs, perm = F, seed) {
      set.seed(seed) # Set seed
      split <- createDataPartition(y = classes, p = 0.8)[[1]]
      train <- newPredictorSet[split,listOfMetabs]
      test <- newPredictorSet[-split,listOfMetabs]
      trainY <- classes[split]
      testY <- classes[-split]
      set.seed(seed) # Set seed
      stratFolds <- createMultiFolds(trainY, k = 3, times = 10)
      ctrl <- trainControl(method = "repeatedcv", index = stratFolds, classProbs = T,
                           summaryFunction = twoClassSummary)
      if(perm == T) {
            trainY <- sample(trainY)
            testY <- sample(testY)
      }
      SVMFit <- train(x = train, y = trainY,
                      trControl = ctrl,
                      method = "svmRadial",
                      metric = "ROC",
                      tuneLength = 15)
      testprediction <- predict(SVMFit, newdata = test, type = "prob")$Good
      result <- roc(response = testY, predictor = testprediction)
      return(as.numeric(result$auc))
}

best5Clas <- vector()
best10Clas <- vector()
best20Clas<- vector()
best50Clas <- vector()
best100Clas <- vector()
worst100Clas <- vector()
worst50Clas<- vector()
worst20Clas<- vector()
worst10Clas<- vector()
worst5Clas<- vector()
allClas <- vector()
allPermClas <- vector()
random5Clas <- vector()


for(r in 1:25){
      best5Clas <- c(best5Clas, SVMmodel(rankMetabsNames[1:5], seed = r))
      best10Clas <- c(best10Clas, SVMmodel(rankMetabsNames[1:10][sample(1:10, 5)], seed = r))
      best20Clas<- c(best20Clas,SVMmodel(rankMetabsNames[1:20][sample(1:20, 5)], seed = r))
      best50Clas <- c(best50Clas,SVMmodel(rankMetabsNames[1:50][sample(1:50, 5)], seed = r))
      best100Clas <- c(best100Clas,SVMmodel(rankMetabsNames[1:100][sample(1:100, 5)], seed = r))
      worst100Clas <- c(worst100Clas ,SVMmodel(rev(rankMetabsNames)[1:100][sample(1:100, 5)], seed = r))
      worst50Clas<- c(worst50Clas,SVMmodel(rev(rankMetabsNames)[1:50][sample(1:50, 5)], seed = r))
      worst20Clas<- c(worst20Clas, SVMmodel(rev(rankMetabsNames)[1:20][sample(1:20, 5)], seed = r))
      worst10Clas<- c(worst10Clas, SVMmodel(rev(rankMetabsNames)[1:10][sample(1:10, 5)], seed = r))
      worst5Clas <- c(worst5Clas, SVMmodel(rev(rankMetabsNames)[1:5], seed = r))
      allClas <- c(allClas, SVMmodel(rankMetabsNames, seed = r))
      allPermClas <- c(allPermClas, SVMmodel(rankMetabsNames, perm = T, seed = r))
      set.seed(r)
      random <- rankMetabsNames[sample(1:length(rankMetabsNames), 5)]
      random5Clas <- c(random5Clas, SVMmodel(random, seed = r))
}

pdf("FigureS4.pdf", height = 4.5, width = 5)
par(mar = c(10.1, 4.1, 4.1, 2.1))
boxplot(best5Clas, best10Clas,best20Clas,
        best50Clas,best100Clas,worst100Clas,
        worst50Clas,worst20Clas,worst10Clas,
        worst5Clas,random5Clas,allClas,allPermClas,
        names = c("Top 5", "5 out of top 10", "5 out of top 20",
                  "5 out of top 50", "5 out of top 100", "5 out of bottom 100",
                  "5 out of bottom 50", "5 out of bottom 20", "5 out of bottom 10",
                  "Bottom 5", "Random 5", "All", "All (perm)"), las = 2,
        ylab = "Area under the ROC curve (AUC)", cex.axis = 0.8)
abline(h = seq(-2, 9, by = 1)/10, lty = 2, col = "grey")
dev.off()

## Precision-Recall curves
classes <- factor(ifelse(relevantBiomass < mean(relevantBiomass), "Bad","Good"))
cols <- rainbow(11,alpha = 0.7)
SVMmodelROC <- function(listOfMetabs, perm = F, seed, col = "black", add = T, cutoff = NULL) {
      set.seed(seed) # Set seed
      split <- createDataPartition(y = classes, p = 0.8)[[1]]
      train <- newPredictorSet[split,listOfMetabs]
      test <- newPredictorSet[-split,listOfMetabs]
      trainY <- classes[split]
      testY <- classes[-split]
      set.seed(seed) # Set seed
      stratFolds <- createMultiFolds(trainY, k = 3, times = 10)
      ctrl <- trainControl(method = "repeatedcv", index = stratFolds, classProbs = T)
      if(perm == T) {
            trainY <- sample(trainY)
            testY <- sample(testY)
      }
      SVMFit <- train(x = train, y = trainY,
                      trControl = ctrl,
                      method = "svmRadial",
                      metric = "Accuracy",
                      tuneLength = 15)
      probs <- predict(SVMFit, newdata = test, type = "prob")$Good
      pred <- prediction(predictions = probs, labels = testY, label.ordering = c("Bad", "Good"))
      perf <- performance(pred, measure = "prec", x.measure = "rec")
      if(add == T){
            if(is.null(cutoff)){
                  plot(perf, add = T, col = col, lwd = 2) 
            }else{
                  plot(perf, add = T, col = col, print.cutoffs.at = cutoff, lwd = 2)
            }
            
      }else{
            if(is.null(cutoff)){
                  plot(perf, add = F, col = col, xlim = c(0,1), ylim = c(0,1), lwd = 2)   
            }else{
                  plot(perf, add = F, col = col, xlim = c(0,1), ylim = c(0,1),
                       print.cutoffs.at = cutoff, lwd = 2)   
            }
            
      }
      
}
pdf("Figure6.pdf", 5, 5)
SVMmodelROC(rankMetabsNames[1:5], seed = 99, col = cols[1], add = F, cutoff = 0.5)
SVMmodelROC(rankMetabsNames[1:10][sample(1:10, 5)], seed = 99, col = cols[2], add = T)
SVMmodelROC(rankMetabsNames[1:20][sample(1:20, 5)], seed = 99, col = cols[3], add = T)
SVMmodelROC(rankMetabsNames[1:50][sample(1:50, 5)], seed = 99, col = cols[4], add = T)
SVMmodelROC(rankMetabsNames[1:100][sample(1:100, 5)], seed = 99, col = cols[5], add = T)
SVMmodelROC(rev(rankMetabsNames)[1:100][sample(1:100, 5)], seed = 99, col = cols[6], add = T)
SVMmodelROC(rev(rankMetabsNames)[1:50][sample(1:50, 5)], seed = 99, col = cols[7], add = T)
SVMmodelROC(rev(rankMetabsNames)[1:20][sample(1:20, 5)], seed = 99, col = cols[8], add = T)
SVMmodelROC(rev(rankMetabsNames)[1:10][sample(1:10, 5)], seed = 99, col = cols[9], add = T)
SVMmodelROC(rev(rankMetabsNames)[1:5], seed = 99,  col = cols[10], add = T)
set.seed(99)
random <- rankMetabsNames[sample(1:length(rankMetabsNames), 5)]
SVMmodelROC(random, seed = 99, col = cols[11], add = T)
SVMmodelROC(rankMetabsNames, seed = 99, col = rgb(0,0,0, alpha = 0.7), add = T)
SVMmodelROC(rankMetabsNames, perm = T, seed = 99, col = rgb(.5,.5,.5,alpha = 0.7), add = T)

legend("topright", xpd = T, legend = c("Top 5", "5 out of top 10", "5 out of top 20",
                  "5 out of top 50", "5 out of top 100", "5 out of bottom 100",
                  "5 out of bottom 50", "5 out of bottom 20", "5 out of bottom 10",
                  "Bottom 5", "Random 5", "All", "All (perm)"), lty = 1, col = c(rainbow(11,alpha = 0.7), rgb(0,0,0, alpha = 0.3), rgb(.5,.5,.5,alpha = 0.3)), cex = 0.5)
dev.off()
# Cross-predictability across trials
SVMmodelCP <- function(listOfMetabs, perm = F, train = "2012", seed) {
      set.seed(seed) # Set seed
      split <- which(relevantTrial == train)
      train <- newPredictorSet[split,listOfMetabs]
      test <- newPredictorSet[-split,listOfMetabs]
      trainY <- classes[split]
      testY <- classes[-split]
      set.seed(seed) # Set seed
      stratFolds <- createMultiFolds(trainY, k = 3, times = 10)
      ctrl <- trainControl(method = "repeatedcv", index = stratFolds)
      if(perm == T) {
            trainY <- sample(trainY)
            testY <- sample(testY)
      }
      SVMFit <- train(x = train, y = trainY,
                      trControl = ctrl,
                      method = "svmRadial",
                      metric = "Accuracy",
                      tuneLength = 15)
      
      testprediction <- predict(SVMFit, newdata = test)
      result <- defaultSummary(data.frame("obs" = testY, "pred" = testprediction))
      return(result)
}


CP2012 <- SVMmodelCP(listOfMetabs = rankMetabsNames[1:5], perm = F,train = "2012",
                            seed = 1)

CP2010 <- SVMmodelCP(listOfMetabs = rankMetabsNames[1:5], perm = F,train = "2010",
                     seed = 1)

# Is HP correlated with the frequency of OD or D?
freqFW <- FW[colnames(hybridClasses),"FW"]
freqOD <- apply(hybridClasses, 2, function(x) sum(abs(x) == 2)/length(x))
freqDOD <- apply(hybridClasses, 2, function(x) sum(x != 0)/length(x))
trial <- FW[colnames(hybridClasses), "trial"]
noFW <- which(is.na(freqFW))
pdf("FigureSuppNote.pdf", 4, 4)
plot(freqFW[-noFW], freqDOD[-noFW], pch = 16,
     col = ifelse(trial[-noFW] == 2010, rgb(0,0,0,.3), rgb(1,0,0,.3)),
     xlim = c(350, 850), ylim = c(0,.8),
     xlab = expression(paste("FW [q ha"^"-1","]")),
     ylab = "1 - % additivity")
dev.off()
cor(freqFW[-noFW], freqDOD[-noFW])

# Is the mean FW robust enough to serve as split for defining classes?
hist(relevantBiomass, xlab = expression(paste("FW [dt ha"^"-1","]")), main = NULL)
abline(v=mean(relevantBiomass))
abline(v=median(relevantBiomass), col = "red")
legend("topright", legend = c("mean","median"), pch = NULL, text.col = c("black","red"),
       bty = "n")

# Check if binning of FW w/ different values affects accuracy
FWcutoff <- seq(min(relevantBiomass), max(relevantBiomass), length.out = 11)
FWcutoff <- FWcutoff[-c(1,11)]
binning <- data.frame(nrow=1,ncol=2)
for(i in 1:length(FWcutoff)){
      classes <- factor(ifelse(relevantBiomass < FWcutoff[i], "Bad","Good"))
      for(k in 1:20){
            binning <- rbind(binning, SVMmodel(listOfMetabs = rankMetabsNames[1:5], seed = k))
      }
}
binning <- binning[-1,]
pdf("FigureS3.pdf", 7, 7)
par(mfrow = c(2,1))
boxplot(binning[1:20,1],
        binning[21:40,1],
        binning[41:60,1],
        binning[61:80,1],
        binning[81:100,1],
        binning[101:120,1],
        binning[121:140,1],
        binning[141:160,1],
        binning[161:180,1], ylab = "Accuracy", names = round(FWcutoff), xlab = expression(paste("FW binning cutoff [dt ha"^"-1","]")))

boxplot(binning[1:20,2],
        binning[21:40,2],
        binning[41:60,2],
        binning[61:80,2],
        binning[81:100,2],
        binning[101:120,2],
        binning[121:140,2],
        binning[141:160,2],
        binning[161:180,2], ylab = "Kappa", names = round(FWcutoff), xlab = expression(paste("FW binning cutoff [dt ha"^"-1","]")))
dev.off()
