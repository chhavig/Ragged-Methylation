#argument 1: file name to save the variance data

library(jsonlite)
library(bsseq)

args <- commandArgs(trailingOnly = TRUE)

f <- file("stdin", "r")

if(file.exists(args[1])){
  file.remove(args[1])
}

lost=0
cnt=0
stop <- FALSE
while (!stop){
  line <- readLines(f,n=1)
  if (length(line)>0){
    d <- fromJSON(line)
cnt=cnt+1
print(cnt)
    
    x<-names(d)
    #print(x)
    comp=list()
    sum_var=0
    region<-0
    status<-0
    for(i in names(d)){
      if(i=='REGION_INFO'){
        region<-d[[i]]
      }
      
      else if (i=='STATUS'){
        status<-as.numeric(d[[i]])
      }
      
      else{
        
        #saving all replicates methylation details in comp
        DF<-NULL
        #print(length(d[[i]]))
        for(j in 1:length(d[[i]])){
          #print(strsplit (d[[i]][j], " ")[[1]])
          rbind(DF,c(strsplit (d[[i]][j], " ")[[1]]))->DF
        }
        
        colnames(DF) <- c('chr', 'start', 'end', 'score', 'M', 'U')
        comp[[i]]<-DF
      }
    }
    
    print('Region:')
    print(region)
    words <- strsplit(region, " ")[[1]]
    
    #print(length(comp))
    sum_var=0
    #loop 1, smoothing methylation data and calculating variance within replicates
    if(length(comp)>0){
      for(i in 1:length(comp)){
        data<-comp[[i]]
        M <- matrix (as.numeric(data [,'M']), ncol=1 )
        Cov <- matrix (as.numeric(data [,'M'])+as.numeric(data [,'U']), ncol=1 )
        
        cover=as.numeric(data [,'M'])+ as.numeric(data [,'U'])
        sc=((as.numeric(data [,'M']))/(cover))*100
        if(!is.na(var(sc))){
          
        if ((var(sc))==0)
        {
          df1<-data.frame(region,i,var(sc),var(sc))
          sum_var=sum_var+var(sc)
        }
          
        if ((var(sc))!=0) {
          
        BStmp <- BSseq (chr = c(data[,'chr']), pos = as.numeric(data [,'start']), M=M, Cov=Cov)
        BStemp<-BSmooth(BStmp, ns=10, h=300, maxGap=3000, keep.se=TRUE, verbose=FALSE)
        
        #calculating regions lost when var is NA
        if(is.na(var(getMeth(BStemp))))
          lost=lost+1
        
        if(!is.na(var(getMeth(BStemp)))){
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
            data[,'score']<-(getMeth(BStemp)*100)
            comp[[i]]<-data
          }
        }
        }
        }
      }

    df2<-comp[[1]]
    #print(df2)
    num<-c()
    
    #loop2, calculating variance across replicates
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
            present=TRUE
            n<-(df3[l,'score'])
            break;
          }
        }
        
        if(present==FALSE){
          break;
        }
        else{
          num1<-c(num1,n)
        }
      }
      
      #print('Out here')
      
      if(length(num1)>1){
        va1<-var(num1) #variance across replicates for a CpG
        num<-c(num,va1)
      }
    }
    
    #print(num)
    print('Result 1, mean of variance different replicates:')
    print(mean(sum_var))
    
    if(!is.null(num)){
      print('Result 2, mean of variance of different CpGs across replicates:')
      print(mean(num))
      dff<-data.frame(words[[1]],region,status,mean(sum_var),mean(num))
    }
    
    else{
      print('Cant calculate Result 2')
      dff<-data.frame(words[[1]],region,status,mean(sum_var),0)
    }
    write.table(dff, file = args[1], append = TRUE, quote = TRUE, sep = " , ",
                eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                col.names = FALSE, qmethod = c("escape", "double"),
                fileEncoding = "")
    }
    
  } else {
    stop <- TRUE
    print("Finished")
  }
}

print("Regions lost:")
print(lost)