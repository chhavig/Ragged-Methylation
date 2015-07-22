#Argument1: file with hyperbola parameters [best1.txt]
#Argument2: methylation variance data of all random data points (has status 1, for all given by Anastasiya, and -1 for all others) [results_20_02.csv]
#Argument3: methylation variance data of training data points (has status 1 for ragged, and -1 for non ragged) [checked.bed]
#Argument4: file where it'll save all regions identified as ragged [regions.csv]

args <- commandArgs(trailingOnly = TRUE)
#conn <- file('best1.txt',open='r')
conn <- file(args[1],open='r')
linn=readLines(conn)
words <- strsplit(linn[1], "\t")[[1]]
print(words)
#data1<-read.csv('result_20_02.csv', header = FALSE)
data1<-read.csv(args[2], header = FALSE)
data2<-NULL
#data2<-read.table('checked.bed')
data2<-read.table(args[3])
#print(data1)
#print(nrow(data1))
xlow=min((data1[,'V4']))
xhigh=max((data1[,'V4']))
ylow=min((data1[,'V5']))
#ylow=-100
yhigh=max((data1[,'V5']))
#yhigh=100
xr=c(xlow,xhigh)
yr=c(ylow,yhigh)
#print(xr)
#print(yr)
#if(file.exists('regions.csv')){
#  file.remove('regions.csv')
#}
if(file.exists(args[4])){
  file.remove(args[4])
}
if(data1[1,'V3']==-1)
{
  plot((data1[1,'V4']),(data1[1,'V5']),main="Ragged methylation smoothed", 
       xlab="Smoothed Variance within replicates",ylab="Variance across replicates", 
       pch=18, col='green',xlim=xr,ylim=yr)
}
if(data1[1,'V3']==1)
{
  plot((data1[1,'V4']),(data1[1,'V5']),main="Ragged methylation smoothed", 
       xlab="Smoothed Variance within replicates", ylab="Variance across replicates", 
       pch=20, col='red',xlim=xr,ylim=yr)
}
#print(nrow(data1))
for(i in 2:nrow(data1))
{
  if(data1[i,'V3']==-1)
  {
    #print('here1')
    points((data1[i,'V4']),(data1[i,'V5']),pch=18,col='green')
  }
}

if(!is.null(data2)){

for(i in 1:nrow(data2))
{
  if(data2[i,'V3']==1)
  {
    #print('here2')
    points((data2[i,'V4']),(data2[i,'V5']),pch=20,col='red')
  }
  if(data2[i,'V3']==-1)
  {
  #print('here2')
    points((data2[i,'V4']),(data2[i,'V5']),pch=20,col='blue')
  }
}
#points(data2[,'V4'],data2[,'V5'],pch=20,col='red')

}

f=function(i)(as.double(words[4])/(as.double(i)-as.double(words[2])))+as.double(words[3])
 
for(i in 1:nrow(data1))
{
  if(data1[i,'V3']==1)
  {
    #print('here2')
    points((data1[i,'V4']),(data1[i,'V5']),pch=15,col='red')
  }
  
  if(data1[i,'V5']>f(data1[i,'V4']))
  {
    df<-data.frame(data1[i,])
    #print(df)
    write.table(df, file = args[4], append = TRUE, quote = TRUE, sep = " , ",
                eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                col.names = FALSE, qmethod = c("escape", "double"),
                fileEncoding = "")
  }
}

curve(f,0.0,xhigh+100,add=TRUE,n=5000000)
identify(data1[,'V4'],data1[,'V5'], labels=data1[,'V1'])
close(conn)