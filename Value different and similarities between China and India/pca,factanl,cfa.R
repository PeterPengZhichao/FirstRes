#install.packages("GPArotation")
install.packages("lavaan")
install.packages("semPlot")
install.packages("BDgraph")
install.packages("Rmpi")
install.packages("OpenMx",dependencies=TRUE)

library(GPArotation)
library(psych)
library(lavaan)
library(BDgraph)
library(Rmpi)
library(OpenMx)
library(semPlot)

#read in data
China=read.csv(file="C_China_SexNation",header=T)
India=read.csv(file="C_India_SexNation",header=T)
China=China[-1]
India=India[-1]
China[,c('Sex','nation')]=lapply(China[,c('Sex','nation')],factor)
India[,c('Sex','nation')]=lapply(India[,c('Sex','nation')],factor)
CombinedData=rbind(China,India)
str(CombinedData)
summary(CombinedData)

drData=CombinedData[,-c(1,5,29,30,31,32)]
pca=princomp(drData)
summary(pca)
print(pca$loadings, digits=3, sort=TRUE)
plot(pca)
screeplot(pca,type = "l")

#Exploratory factor analysis
scree(drData)
fa.parallel(drData)#suggest 8, look into it and group into 5
fa1=factanal(drData,8,rotation = "varimax",scores = "regression")
fa2=factanal(drData,5,rotation = "varimax",scores = "regression")
print(fa1, digits=3, cutoff=.3, sort=TRUE)




#confirmtary factor analysis
five.model='
  EthicWithInterest=~V318+V319+V315+V301+V314
  MeaningOfLife=~ V313+V310+V312+V316+V309+V303+V311
  EthicWOInterest=~ V298+V296+V299+V302+V297+V300+V306
  SexRelated=~ V305+V304+V307+V308
  Happiness=~ V96+V116+V132'

five.fit=cfa(five.model,data=drData)
coef(five.fit)
summary(five.fit,standardized=TRUE,rsquare=TRUE)
semPaths(five.fit, whatLabels = "std", layout = "tree")

#recode
CombinedData$V91[CombinedData$V91==1]=0
CombinedData$V91[CombinedData$V91==2]=1
CombinedData$V146[CombinedData$V146==1]=0
CombinedData$V146[CombinedData$V146==2]=1

#compute latent variable should update each time
CombinedData$EthicWithInterest=with(CombinedData,0.774*V318+0.708*V319+0.590*V315+0.608*V301+0.605*V314)
CombinedData$MeaningOfLife=with(CombinedData,0.560*V313+0.680*V310+0.565*V312+0.507*V316+0.546*V309+0.497*V303+0.528*V311)
CombinedData$EthicWOInterest=with(CombinedData,0.716*V298+0.436*V296+0.761*V299+0.404*V302+0.669*V297+0.676*V300+0.588*V306)
CombinedData$SexRelated=with(CombinedData,0.686*V305+0.713*V304+0.485*V307+0.495*V308)
CombinedData$Happiness=with(CombinedData,0.68*V96+0.645*V116+0.696*V132)

model1=lm(Happiness~SexRelated+EthicWithInterest+EthicWOInterest+MeaningOfLife+CombinedData$V91+CombinedData$V146+CombinedData$Sex+CombinedData$nation+CombinedData$V355, data=CombinedData)
summary(model1)

model2=glm(CombinedData$V146~MeaningOfLife+Happiness+SexRelated+EthicWOInterest+EthicWithInterest+CombinedData$V91+CombinedData$Sex+CombinedData$nation+CombinedData$V355, data=CombinedData,family = binomial(link = "logit"))
summary(model2)




