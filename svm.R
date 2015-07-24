#Argument1: training data file [checked.bed]
#Argument2: random regions variance data [regions_20_02.csv]
#Argument3: file name to save results in

library(e1071)
args <- commandArgs(trailingOnly = TRUE)
if(file.exists(args[1])){
  data0<-read.table(args[1])

#data0<-read.table('checked.bed')
#print(data2)
if(file.exists(args[2])){
  data1<-read.table(args[2], sep=",", header = FALSE)}

#data1<-read.table('result_20_02.csv', sep=",", header = FALSE)
#print(data3)
data2 <- data0[,3:5]
data3 <- data1[,3:5]
index<-1:nrow(data2)
testindex<-sample(index,trunc(length(index)/3))
testset<-data2[testindex,]
trainset<-data2[-testindex,]
#print(trainset)

model1<-tune.svm(V3 ~ .,data=trainset,gamma=10^(-6:1),cost=10^(1:3))
#print(summary(model1))
bestGamma<-model1$best.parameters[[1]]
bestC<-model1$best.parameters[[2]]

svm.model<-svm(V3 ~ .,data=trainset,method='C-classification', kernel='radial',
               cost=bestC,gamma=bestGamma,verbose=TRUE,cross=10)

print(summary(svm.model))
svm.pred<-predict(svm.model,testset[,2:3],decision.values = TRUE)
#print('yo')

svm.pred1<-predict(svm.model,data3[,2:3],decision.values = TRUE)
#print('yo1')

#if(file.exists('regions1.csv')){
#  file.remove('regions1.csv')}

if(file.exists(args[3])){
  file.remove(args[3])}

for(i in 1:length(svm.pred)){
if(svm.pred[i]<0)
  svm.pred[i]=-1
else
  svm.pred[i]=1
}

for(i in 1:length(svm.pred1)){
  if(svm.pred1[i]<0)
    svm.pred1[i]=-1
  else
    svm.pred1[i]=1
  if((svm.pred1[i]==1) && (data1[i,4]>1000) && (data1[i,5]>100))
  {
    df<-data.frame(data1[i,])
   #print(df)
    write.table(df, file = args[3], append = TRUE, quote = TRUE, sep = " , ",
                eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                col.names = FALSE, qmethod = c("escape", "double"),
                fileEncoding = "")
  }
}

print(table(pred=svm.pred,true=testset[,1]))
}
