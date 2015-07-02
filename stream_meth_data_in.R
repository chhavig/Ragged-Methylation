library(jsonlite)

f <- file("stdin", "r")

stop <- FALSE

while (!stop){
	line <- readLines(f,n=1)

	if (length(line)>0){
		cat("Line found. Reading line.\n")
		dataFrameIn <- fromJSON(line)

		methDataByRep <- list()

		for (i in names(dataFrameIn)){
			if(i == 'REGION_INFO'){
				cat("Found region information\n")
				region <- c(strsplit(dataFrameIn[[i]]," ")[[1]])
			}
			else if (i == 'STATUS'){
				cat("Found status\n")
				status <- as.numeric(dataFrameIn[[i]])
			}
			else {
				cat("Found methylation data for replicate", i, "\n")
				methDataFrame <- NULL
				for (cpg in dataFrameIn[[i]]){
					methDataFrame <- rbind(methDataFrame,c(strsplit (cpg, " ")[[1]]))
				}

			colnames(methDataFrame) <- c('chr', 'start', 'end', 'score', 'M', 'U')
			newDF <- NULL
			newDF <- cbind(newDF, methDataFrame[,'chr'])
			newDF <- cbind(newDF, as.numeric(methDataFrame[,'start']))
			colnames(newDF) <- c('chr','start')
			print(newDF)

			#print(methDataFrame)
			#print(methDataFrame[,'start'])
			#print(typeof(methDataFrame[,'start']))
			#tmp <- as.numeric(methDataFrame[,'start'])
			#print(typeof(tmp))
			#methDataFrame[,'start'] <- tmp
			#print(typeof(methDataFrame[,'start']))
			#print(methDataFrame[,'start'])
			#methDataFrame$end <- as.numeric(methDataFrame$end)
			#methDataFrame$score <- as.numeric(methDataFrame$score)
			#methDataFrame$M <- as.numeric(methDataFrame$M)
			#methDataFrame$U <- as.numeric(methDataFrame$U)
			
			#print(methDataFrame)
			methDataByRep[[i]] <- methDataFrame
			}
		}
	}
	else {
		stop <- TRUE
	}

	print(length(methDataByRep))
}

