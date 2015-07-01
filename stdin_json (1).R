library(jsonlite)

f <- file("stdin", "r")

#stop <- FALSE

df <- stream_in(f)
#print(df)
v <- c()

d=list()

for (s in colnames(df)){
	#s = colnames(df)[i]
	d[[s]] <- unlist(df[s][[1]][[1]])
}
print(d)

#while(!stop)

x<-names(d)
#print(x)
comp=list()
sum_var=0
for(i in names(d)){
  if(i=='REGION_INFO'){
    region<-c(strsplit(x," "))
  }
  else if (i=='STATUS'){
    status<-as.numeric(d[[i]])
  }
  else{
    DF<-NULL
    #print(length(d[[i]]))
    for(j in 1:length(d[[i]])){
      #print(strsplit (d[[i]][j], " ")[[1]])
      rbind(DF,c(strsplit (d[[i]][j], " ")[[1]]))->DF
    }
    #print (DF)
    colnames(DF) <- c('chr', 'start', 'end', 'score', 'M', 'U')
    comp[[i]]<-DF
  }
}

#smoothing

sum_var=0
if(FALSE){
for(i in 1:names(comp)){
  data<-comp[[i]]
  #print(data)
  M <- matrix (data [,'M'], ncol=1 )
  Cov <- matrix (data [,'M']+ data [,'U'], ncol=1 )
  
  if ((var((data[,'M'])/(data[,'M']+data[,'U'])))==0)
  {
    df1<-data.frame(region,i,var((data[,'M'])/(data[,'M']+data[,'U'])),var((data[,'M'])/(data[,'M']+data[,'U'])))
  }
  
  if ((var((data[,'M'])/(data[,'M']+data[,'U'])))!=0) {
    
    BStmp <- BSseq (chr = c(data[,'chr']), pos = data [,'start'], M=M, Cov=Cov)
    BStemp<-BSmooth(BStmp, ns=10, h=300, maxGap=3000, keep.se=TRUE, verbose=FALSE)
    
    cover=data[,'M']+data[,'U']
    print(summary(getMeth(BStemp)))
    #print(var(getMeth(BStemp)))
    
    if (((var(getMeth(BStemp)))/(var((data[,'M'])/(data[,'M']+data[,'U']))))>=1.45)
    {
      df1<-data.frame(region,i,var((data[,'M'])/(data[,'M']+data[,'U'])),var((data[,'M'])/(data[,'M']+data[,'U'])))
    }
    
    else
    {
      sum_var=sum_var+var(df1[2])
      df1<-data.frame(region,i,var((data[,'M'])/(data[,'M']+data[,'U'])),var(getMeth(BStemp)))
    }
    
  }
  write.table(df1, file = 'random.csv', append = TRUE, quote = TRUE, sep = " , ",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = FALSE, qmethod = c("escape", "double"),
              fileEncoding = "")
}
}

#for(i in 1:length(comp)){
#  df1<-comp[[i]]
#  sum_var=sum_var+var(df1[,'score'])
#}

df2<-comp[[1]]
print(df2)
num<-c()

for(j in 1:nrow(df2)){
  search=paste0(df2[j,'chr'],df2[j,'start'],' ')
  #print(search)
  num1<-c()
  #print(length(comp))
  for(k in 1:length(comp)){
    df3<-comp[[k]]
    #print(df3)
    present=FALSE
    #print(search == df3[,'chr'])
    
    for(l in 1:nrow(df3)){
      chck=paste0(df3[l,'chr'],df3[l,'start'],' ')
      if(chck==search){
        #print('Yep')
        present=TRUE
        n<-as.numeric(df3[l,'score'])
        #print(n)
        break;
      }
    }
    
    if(present==FALSE){
      break;
    }
    else{
      num1<-c(num1,n)
      #print(num1)
    }
  }
  
  #print('Out here')
  
  if(length(num1)==length(comp)){
    num<-c(num,num1)
    #print(num)
  }
  #print(num)
}

print(num)

print(sum_var)
print(sum_var/length(comp))

get=matrix(c(num), nrow=3)
print(get)

va<-c()

for(go in 1:ncol(get)){
  va<-c(va,var(get[,go]))
}

print(mean(va))

#print(comp)

