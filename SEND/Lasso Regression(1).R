#install.packages("glmnet")
#install.packages("Rcpp")
#install.packages("ggplot2")
#install.packages("psych")
#install.packages("dplyr")
#install.packages("caret")
#install.packages("mlbench")
library(caret)
library(glmnet)
library(psych)
library(mlbench)
library(ggplot2)
library(dplyr)
#library(doSNOW)
#library(Rcpp)
#library(mclust)

China=read.csv(file="C_China_SexNation",header=T)
India=read.csv(file="C_India_SexNation",header=T)
China=China[-1]
India=India[-1]
China[,c('Sex','nation')]=lapply(China[,c('Sex','nation')],factor)
India[,c('Sex','nation')]=lapply(India[,c('Sex','nation')],factor)
CombinedData=rbind(China,India)
CombinedData[,"V91"]=factor(CombinedData[,"V91"])
CombinedData[,"V146"]=factor(CombinedData[,"V146"])
str(CombinedData)
summary(CombinedData)



#Partition data into training and testing data set
set.seed(1234)
indices=createDataPartition(CombinedData$nation,times=1,p=0.7,list=FALSE)
CombinedData.train=CombinedData[indices,]
CombinedData.test=CombinedData[-indices,]

#Control parameter
train.control=trainControl(method = "repeatedcv",
                           number = 10, 
                           repeats = 1,
                           verboseIter = T)

#fit the model

set.seed(1234)

lasso=train(V91~.,
            CombinedData.train,
            method="glmnet",
            tuneGrid=expand.grid(alpha = 1,
                                 lambda=seq(0.0001,0.01,length=40)),
            trControl=train.control)
plot(lasso)
plot(lasso$finalModel,xvar="lambda",label = T)
plot(lasso$finalModel,xvar="dev",label = T)
plot(varImp(lasso,scale=T))

#CombinedData[,"nation"]=as.numeric(CombinedData[,"nation"])
fit=glmnet(x=as.matrix(CombinedData.train[,-1]),y=CombinedData.train[,1],family="binomial",alpha=1,lambda=lasso$bestTune$lambda)
fit$beta

pred=predict(lasso,CombinedData.train)
sum(CombinedData.train[,1]==pred)/length(CombinedData.train[,1])

pred2=predict(lasso,CombinedData.test)
sum(CombinedData.test[,1]==pred2)/length(CombinedData.test[,1])
