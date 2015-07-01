args <- commandArgs(trailingOnly = TRUE)
conn=file(args[1],open='r')
linn=readLines(conn)
words <- strsplit(linn[1], "\t")[[1]]
print(words)
data1<-read.csv('output.csv', header = FALSE)
data2<-read.csv('random (3).csv',header=FALSE)
xr<-c(0,0.16)
yr<-c(0,1.0)
plot(data2[,'V8'],data2[,'V6'],main="Ragged methylation smoothed", 
     xlab="Variance",ylab="Mean", pch=18, col=95, xlim=xr, ylim=yr)
#identify(data2[,'V8'],data2[,'V6'], labels=data2[,'V3'])
points(data1[,'V8'],data1[,'V6'],pch=20,col='red')
print(as.double(as.double(words[4])/(-1*as.double(words[3]))))
print(as.double(as.double(words[2])/(-1*as.double(words[3]))))
abline(a=as.double(as.double(words[4])/(-1*as.double(words[3]))),
       b=as.double(as.double(words[2])/(-1*as.double(words[3]))))
plot(data2[,'V7'],data2[,'V5'],main="Ragged methylation raw", 
     xlab="Variance",ylab="Mean", pch=18, col='green',xlim=xr,ylim=yr)
points(data1[,'V7'],data1[,'V5'],pch=20,col='blue')