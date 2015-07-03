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
				region <- c(strsplit(dataFrameIn[[i]]," ")[[1]])
				cat("Found region information",dataFrameIn[[i]],'\n')
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
			####################################################
			#At the moment methDataFrame is a matrix
			methDataFrame <- data.frame(methDataFrame)
			#Now it's a data frame
			for (n in colnames(methDataFrame)){
				if (n != 'chr'){
					methDataFrame[,n] <- as.numeric(methDataFrame[,n])
				} else {
					methDataFrame[,n] <- as.character(methDataFrame[,n])
				}
			}
			#Now all the columns are the right type
			#####################################################
			methDataByRep[[i]] <- methDataFrame
			}
		}

		#print(methDataByRep)
	}
	else {
		stop <- TRUE
	}

	
}

