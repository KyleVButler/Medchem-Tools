#This will plot multiple IC50 curves, but will need modification to work with variable response (Emax values)

setwd("/Users/kylebutler/Desktop/R Files")
testdata <- read.csv("MultiDRCTestData.csv", header = TRUE)
head(testdata)
CompoundList <- levels(testdata$Compound)
#Normalize each compounds data to maximum 100
for(i in 1:length(CompoundList)){
  weight <- (100/(max(testdata[testdata$Compound == CompoundList[i], 3])))
  testdata[testdata$Compound == CompoundList[i], 3] <- testdata[testdata$Compound == CompoundList[i], 3]*weight
}
#Change Dose from micromolar to molar
testdata[, 2] <- testdata[ ,2]*0.000001

#Create a matrix to see the results
results <- matrix(data = NA, nrow = length(CompoundList), ncol = 2)
row.names(results) <- CompoundList
colnames(results) <- c("Slope", "EC50")

#create a list that will be used for curves and ribbons
fitslist <- list()
length(fitslist) <- length(CompoundList)

#Create a data set for generating the fitted line - may need to change the log values below to fit the curves properly
for(i in 1:length(CompoundList)){
  fits <- expand.grid(conc=exp(seq(log(1.00e-03), log(1.00e-09), length=200))) 
  LLfit <- drm(data = testdata[testdata$Compound == CompoundList[i], 2:3],Response~Dose,fct=LL.4())
  results[i, 1] <- LLfit$coefficients[1]
  results[i, 2] <- LLfit$coefficients[4]
  pm2 <- predict(LLfit, newdata=fits, interval="confidence") 
  fits$p <- pm2[,1]
  fits$pmin <- pm2[,2]
  fits$pmax <- pm2[,3]
  fits$Compound <- c(CompoundList, rep(NA, nrow(fits)-length(CompoundList)))
  fitslist[[i]] <- fits
}
results


#This will plot all but the curves are the same color
dev.new()
gg1 <- ggplot(data = testdata, aes(x = Dose, y = Response, colour = Compound)) + stat_summary(fun.y = "mean", 
  fun.ymin = min, fun.ymax = max, size=1)
for(i in 1:length(CompoundList)){
  gg1 <- gg1 + geom_line(data=fitslist[[i]], aes(x=conc, y=p, colour = Compound), cex = 1.5)
}
gg1 <- gg1 + scale_x_log10() + theme(text = element_text(size=28)) + labs(x = "Concentration (M)", y = "% Activity") + 
  theme(legend.position = c(0.8, 0.9))
gg1


#If you want curves to be different colors you need to manually duplicate the geom_line for as many elements are in fitslist
dev.new()
gg1 <- ggplot(data = testdata, aes(x = Dose, y = Response, colour = Compound)) + stat_summary(fun.y = "mean", 
  fun.ymin = min, fun.ymax = max, size=1)

gg1 <- gg1 + geom_line(data=fitslist[[1]], aes(x=conc, y=p, colour = CompoundList[1]), cex = 1.5)
gg1 <- gg1 + geom_line(data=fitslist[[2]], aes(x=conc, y=p, colour = CompoundList[2]), cex = 1.5)
gg1 <- gg1 + geom_line(data=fitslist[[3]], aes(x=conc, y=p, colour = CompoundList[3]), cex = 1.5)

gg1 <- gg1 + scale_x_log10() + theme(text = element_text(size=28)) + labs(x = "Concentration (M)", y = "% Activity") + 
  theme(legend.position = c(0.8, 0.9))
gg1

#For Ribbons you also need to manually insert them again
dev.new()
gg1 <- ggplot(data = testdata, aes(x = Dose, y = Response, colour = Compound)) + stat_summary(fun.y = "mean", 
  fun.ymin = min, fun.ymax = max, size=1)

gg1 <- gg1 + geom_line(data=fitslist[[1]], aes(x=conc, y=p, colour = CompoundList[1]), cex = 1.5)
gg1 <- gg1 + geom_line(data=fitslist[[2]], aes(x=conc, y=p, colour = CompoundList[2]), cex = 1.5)
gg1 <- gg1 + geom_line(data=fitslist[[3]], aes(x=conc, y=p, colour = CompoundList[3]), cex = 1.5)
gg1 <- gg1 + geom_ribbon(data=fitslist[[1]], aes(x=conc, y=p, ymin=pmin, ymax=pmax, colour = CompoundList[1]), alpha=0.15)
gg1 <- gg1 + geom_ribbon(data=fitslist[[2]], aes(x=conc, y=p, ymin=pmin, ymax=pmax, colour = CompoundList[2]), alpha=0.15)
gg1 <- gg1 + geom_ribbon(data=fitslist[[3]], aes(x=conc, y=p, ymin=pmin, ymax=pmax, colour = CompoundList[3]), alpha=0.15)

gg1 <- gg1 + scale_x_log10() + theme(text = element_text(size=28)) + labs(x = "Concentration (M)", y = "% Activity") + 
  theme(legend.position = c(0.8, 0.9))
gg1
