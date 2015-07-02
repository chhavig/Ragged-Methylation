library(jsonlite)
library(bsseq)

f <- file("stdin", "r")

stop <- FALSE
while (!stop){
  line <- readLines(f,n=1)
  if (length(line)>0){
    d <- fromJSON(line)
#df <- stream_in(f)
#print(df)

#d=list()
    #print(d)
    
    #while(!stop)
    
    x<-names(d)
    #print(x)
    comp=list()
    sum_var=0
    region<-0
    status<-0
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
        
        colnames(DF) <- c('chr', 'start', 'end', 'score', 'M', 'U')
        #print(DF)
        #DF[,'start']<-as.numeric(DF[,'start'])
        #print(DF[1,'start']+DF[2,'start'])
        #print(DF)
        #DF$end<-as.numeric(DF$end)
        #DF$score<-as.numeric(DF$score)
        #DF$M<-as.numeric(DF$M)
        #DF$U<-as.numeric(DF$U)
        #print (DF)
        comp[[i]]<-DF
      }
    }
    
    #smoothing
    
    sum_var=0
    #if(FALSE){
      for(i in 1:length(comp)){
        data<-comp[[i]]
        #print(data)
        M <- matrix (as.numeric(data [,'M']), ncol=1 )
        #print(M)
        Cov <- matrix (as.numeric(data [,'M'])+as.numeric(data [,'U']), ncol=1 )
        #print(Cov)
        
        cover=as.numeric(data [,'M'])+ as.numeric(data [,'U'])
        sc=((as.numeric(data [,'M']))/(cover))*100
        
        if ((var(sc))==0)
        {
          df1<-data.frame(region,i,var(sc),var(sc))
          sum_var=sum_var+var(sc)
        }
        
        if ((var(sc))!=0) {
          
          BStmp <- BSseq (chr = c(data[,'chr']), pos = as.numeric(data [,'start']), M=M, Cov=Cov)
          BStemp<-BSmooth(BStmp, ns=10, h=300, maxGap=3000, keep.se=TRUE, verbose=FALSE)
          
          
          print(summary(getMeth(BStemp)))
          print(var(getMeth(BStemp)*100))
          
          if (((var(getMeth(BStemp)*100))/(var(sc)))>=1.45)
          {
            df1<-data.frame(region,i,var(sc),var(sc))
            sum_var=sum_var+var(sc)
          }
          
          else
          {
            #sum_var=sum_var+var(df1[2])
            df1<-data.frame(region,i,var(sc),var(getMeth(BStemp)*100))
            sum_var=sum_var+var(getMeth(BStemp)*100)
          }
          
        }
        #write.table(df1, file = 'random.csv', append = TRUE, quote = TRUE, sep = " , ",
        #            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
        #            col.names = FALSE, qmethod = c("escape", "double"),
        #            fileEncoding = "")
      }
    #}
    
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
            n<-(df3[l,'score'])
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
    
    #print(sum_var)
    print('Region:')
    print((region))
    print('Result 1, mean of variance different replicates:')
    print(mean(sum_var))
    
    get=matrix(c(num), nrow=3)
    #print(get)
    
    va<-c()
    
    for(go in 1:ncol(get)){
      va<-c(va,var(get[,go]))
    }
    
    print('Result 2, mean of variance of different CpGs across replicates:')
    print(mean(va))
    
    df1<-data.frame(region,mean(sum_var),mean(va))
    
    write.table(df1, file = 'results.csv', append = TRUE, quote = TRUE, sep = " , ",
                eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                col.names = FALSE, qmethod = c("escape", "double"),
                fileEncoding = "")
    
  } else {
    stop <- TRUE
    print("Finished")
  }
}

#for (s in colnames(df)){
	#s = colnames(df)[i]
#	d[[s]] <- unlist(df[s][[1]][[1]])
#}

#print(comp)

