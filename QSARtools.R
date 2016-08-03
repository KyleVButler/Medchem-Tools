#make a list of matrices for various fps
#fps that work - maccs circular standard estate(simple) kr(slow) shortestpath signature(doesn't work)
fpmatlist <-  as.list(c("circular", "standard", "extended"))

fpreturn <- function(fplist, cmpdlist){
  fptest <- get.fingerprint(cmpdlist[[1]], type = fplist, fp.mode = 'bit')
  fpmatrix <- matrix(data = 0, nrow = length(cmpdlist), ncol = attr(fptest, "nbit"))
  for(i in 1:length(cmpdlist)){
    fp <- get.fingerprint(cmpdlist[[i]], type = fplist, fp.mode = 'bit')
    fpmatrix[i, attr(fp, "bits")] <- 1
  }
  fpmatrix <- fpmatrix[, -nearZeroVar(fpmatrix)]

  return(fpmatrix)
}

#generates a list of fingerprint matrices from input IAtomContainers
generatefpmatlist <- function(containerinput){
  fpmatlist <-  as.list(c("circular", "standard", "extended"))
  fpmatrixlist <- lapply(fpmatlist, fpreturn, cmpdlist = containerinput)
  names(fpmatrixlist) <- c("circular", "standard", "extended")
  return(fpmatrixlist)
}

#finds best model from existing data, given vector of data and vector of smiles for 
#corresponding compounds
getmodelfromfp <- function(ydata, xdata){
  containerstrain <- parse.smiles(xdata)
  fpmatrixlist <- generatefpmatlist(containerstrain)
  trainIndex <- createDataPartition(ydata, p=.75, list=FALSE)
  y.tr <- ydata[trainIndex]
  y.te <- ydata[-trainIndex]
  trainvector <- c("gbm", "svmPoly", "pls", "ranger")
  resultsmatrix <- matrix(data = NA, nrow = length(trainvector), ncol = length(fpmatrixlist))
  rownames(resultsmatrix) <- trainvector
  colnames(resultsmatrix) <- as.vector(names(fpmatrixlist))
  ctrl <- trainControl(method = 'cv', number = 10, repeats = 5, 
                       summaryFunction = defaultSummary)
  #ctrl <- trainControl(method = "boot",
  #                     summaryFunction = defaultSummary) 
  for(n in 1:length(fpmatrixlist)){
    x1.tr <- fpmatrixlist[[n]][trainIndex, ]
    x1.te <- fpmatrixlist[[n]][-trainIndex, ]
    for (i in 1:length(trainvector)){
      svm.fit1 <- train(x1.tr, y.tr, method = trainvector[i], trControl = ctrl, tuneLength = 10,
                        metric = 'RMSE')
      pred.activity <- predict(svm.fit1, newdata = x1.te)
      resultsmatrix[i, n] <- mean((predict(svm.fit1, newdata = x1.te)-y.te)^2)
    }
  }
  
  print(resultsmatrix)
  besttrainidx <- which(resultsmatrix == min(resultsmatrix), arr.ind = TRUE)[1]
  bestfpidx <- which(resultsmatrix == min(resultsmatrix), arr.ind = TRUE)[2]
  x1.tr <- fpmatrixlist[[bestfpidx]][trainIndex, ]
  x1.te <- fpmatrixlist[[bestfpidx]][-trainIndex, ]
  bestfit <- train(x1.tr, y.tr, method = trainvector[besttrainidx], trControl = ctrl, tuneLength = 10,
                   metric = 'RMSE')
  attr(bestfit, "BestTrain") <- trainvector[besttrainidx]
  attr(bestfit, "BestFP") <- names(fpmatrixlist)[bestfpidx]
  print(bestfit)
  dev.new()
  #plot(y.tr, predict(bestfit, x1.tr), xlim = range(y.tr), ylim = range(y.tr))
  #abline(a = 0, b = 1)
  #dev.new()
  #plot(y.te, predict(bestfit, x1.te), xlim = range(y.te), ylim = range(y.te))
  plot(y.tr, predict(bestfit, x1.tr), xlim = range(y.tr), ylim = range(y.tr), 
       main = paste("Best Fingerprint ", attr(bestfit, "BestFP"), "\nBest Training Method ", 
                    attr(bestfit, "BestTrain")), xlab = "Data", ylab = "Predicted")
  abline(a = 0, b = 1)
  points(y.te, predict(bestfit, x1.te), col = "blue", pch = 16)
  return(bestfit)
}



#now to get the model for predictions
predictfromsmiles <- function(data.x, NewSmilesIN, bestmodel){
  containers <- parse.smiles(data.x)
  fptest <- get.fingerprint(containers[[1]], type = attr(bestmodel, "BestFP"), fp.mode = 'bit')
  bitl <- attr(fptest, "nbit")
  fpmatrix <- matrix(data = 0, nrow = length(containers), ncol = bitl)
  for(i in 1:length(containers)){
    fp <- get.fingerprint(containers[[i]], type = attr(bestmodel, "BestFP"), fp.mode = 'bit')
    fpmatrix[i, attr(fp, "bits")] <- 1
  }
  matrixloss <- nearZeroVar(fpmatrix)
  fpmatrix <- fpmatrix[, -matrixloss]
  
  containers.predict <- parse.smiles(as.character(NewSmilesIN$SMILES))
  fptest <- get.fingerprint(containers.predict[[1]], type = attr(bestmodel, "BestFP"), fp.mode = 'bit')
  bitl <- attr(fptest, "nbit")
  fpmatrix.predict <- matrix(data = 0, nrow = length(containers.predict), ncol = bitl)
  for(i in 1:length(containers.predict)){
    fp <- get.fingerprint(containers.predict[[i]], type = attr(bestmodel, "BestFP"), fp.mode = 'bit')
    fpmatrix.predict[i, attr(fp, "bits")] <- 1
  }
  fpmatrix.predict <- fpmatrix.predict[, -matrixloss]  
  
  predictvector <- predict(bestmodel, fpmatrix.predict)
  return(predictvector)
}

rfcircularmodel <- function(ydata, xdata){
  containerstrain <- parse.smiles(xdata)
  fpmatrixlist <- fpreturn("circular", containerstrain)
  ctrl <- trainControl(method = 'cv', number = 10, repeats = 5, 
                       summaryFunction = defaultSummary)
  fit1 <- train(fpmatrixlist, ydata, method = "ranger", trControl = ctrl, tuneLength = 10,
                    metric = 'RMSE')
  attr(fit1, "BestTrain") <- "ranger"
  attr(fit1, "BestFP") <- "circular"
  dev.new()
  plot(ydata, predict(fit1, fpmatrixlist), xlim = range(ydata), ylim = range(ydata))
  abline(a = 0, b = 1)
  return(fit1)
}
  

getclassmodel <- function(ydata, xdata, classmethod){
  containerstrain <- parse.smiles(xdata)
  x <- circularfp(containerstrain)
  y <- ydata
  ctrl <- trainControl(method = 'cv', number = 10, repeats = 5, 
                       classProbs = TRUE, summaryFunction = twoClassSummary)
  #ctrl <- trainControl(method = 'boot', 
  #                     classProbs = TRUE, summaryFunction = twoClassSummary)
  svm.fit1 <- train(x, y, method = classmethod, trControl = ctrl, tuneLength = 10,
                    metric = 'ROC')
  print(svm.fit1)
  return(svm.fit1)
}

gbmclassmodel <- function(ydata, xdata, classmethod){
  containerstrain <- parse.smiles(xdata)
  x <- circularfp(containerstrain)
  y <- ydata
  ctrl <- trainControl(method = 'cv', number = 10, repeats = 5, 
                       classProbs = TRUE, summaryFunction = twoClassSummary)
  #grid <- expand.grid(.n.trees = seq(50, 500, by = 50), .interaction.depth = 5:10, .shrinkage = c(0.001, 0.1),
  #                    .n.minobsinnode = c(10))
  svm.fit1 <- train(x, y, method = "gbm", trControl = ctrl, tuneLength = 10,
                    metric = 'ROC')
  print(svm.fit1)
  return(svm.fit1)
}

circularfp <- function(cmpdlist){
  fptest <- get.fingerprint(cmpdlist[[1]], type = "circular", fp.mode = 'bit')
  fpmatrix <- matrix(data = 0, nrow = length(cmpdlist), ncol = attr(fptest, "nbit"))
  for(i in 1:length(cmpdlist)){
    fp <- get.fingerprint(cmpdlist[[i]], type = "circular", fp.mode = 'bit')
    fpmatrix[i, attr(fp, "bits")] <- 1
  }
  fpmatrix <- fpmatrix[, -nearZeroVar(fpmatrix)]
  
  return(fpmatrix)
}

classpredictfromsmiles <- function(data.x, NewSmilesIN, bestmodel){
  containers <- parse.smiles(data.x)
  fptest <- get.fingerprint(containers[[1]], type = "circular", fp.mode = 'bit')
  bitl <- attr(fptest, "nbit")
  fpmatrix <- matrix(data = 0, nrow = length(containers), ncol = bitl)
  for(i in 1:length(containers)){
    fp <- get.fingerprint(containers[[i]], type = "circular", fp.mode = 'bit')
    fpmatrix[i, attr(fp, "bits")] <- 1
  }
  matrixloss <- nearZeroVar(fpmatrix)
  fpmatrix <- fpmatrix[, -matrixloss]
  
  containers.predict <- parse.smiles(as.character(NewSmilesIN$SMILES))
  fptest <- get.fingerprint(containers.predict[[1]], type = "circular", fp.mode = 'bit')
  bitl <- attr(fptest, "nbit")
  fpmatrix.predict <- matrix(data = 0, nrow = length(containers.predict), ncol = bitl)
  for(i in 1:length(containers.predict)){
    fp <- get.fingerprint(containers.predict[[i]], type = "circular", fp.mode = 'bit')
    fpmatrix.predict[i, attr(fp, "bits")] <- 1
  }
  fpmatrix.predict <- fpmatrix.predict[, -matrixloss]
  
  predictprobs <- predict(bestmodel, newdata = fpmatrix.predict, type = "prob")
  predictvec <- predict(bestmodel, newdata = fpmatrix.predict)
  predictions <- data.frame(predictprobs, predictvec)
  return(predictions)
}


gbmcircularregmodel <- function(ydata, xdata){
  containerstrain <- parse.smiles(xdata)
  fpmatrixlist <- fpreturn("circular", containerstrain)
  ctrl <- trainControl(method = 'cv', number = 10, repeats = 4, 
                       summaryFunction = defaultSummary)
  grid <- expand.grid(n.trees = seq(50, 500, 50), interaction.depth = 5:10, shrinkage = 0.1, n.minobsinnode = 10)
  fit1 <- train(fpmatrixlist, ydata, method = "gbm", trControl = ctrl, tuneGrid = grid,
                metric = 'RMSE')
  dev.new()
  plot(ydata, predict(fit1, fpmatrixlist), xlim = range(ydata), ylim = range(ydata))
  abline(a = 0, b = 1)
  return(fit1)
}

#data.x is the smiles string used to build bestmodel, it must be added so that the same nearZeroVar are removed from NewSmilesIN

gbmcircularregpredict <- function(data.x, NewSmilesIN, bestmodel){
  containers <- parse.smiles(data.x)
  fptest <- get.fingerprint(containers[[1]], type = "circular", fp.mode = 'bit')
  bitl <- attr(fptest, "nbit")
  fpmatrix <- matrix(data = 0, nrow = length(containers), ncol = bitl)
  for(i in 1:length(containers)){
    fp <- get.fingerprint(containers[[i]], type = "circular", fp.mode = 'bit')
    fpmatrix[i, attr(fp, "bits")] <- 1
  }
  matrixloss <- nearZeroVar(fpmatrix)
  fpmatrix <- fpmatrix[, -matrixloss]
  
  containers.predict <- parse.smiles(as.character(NewSmilesIN))
  fptest <- get.fingerprint(containers.predict[[1]], type = "circular", fp.mode = 'bit')
  bitl <- attr(fptest, "nbit")
  fpmatrix.predict <- matrix(data = 0, nrow = length(containers.predict), ncol = bitl)
  for(i in 1:length(containers.predict)){
    fp <- get.fingerprint(containers.predict[[i]], type = "circular", fp.mode = 'bit')
    fpmatrix.predict[i, attr(fp, "bits")] <- 1
  }
  fpmatrix.predict <- fpmatrix.predict[, -matrixloss]  
  
  predictvector <- predict(bestmodel, fpmatrix.predict)
  return(predictvector)
}