library(jsonlite)

f <- file("stdin", 'r')

stop <- FALSE
while (!stop){
	line <- readLines(f,n=1)
	if (length(line)>0){
		l <- fromJSON(line)
		print(l)
		#Analysis goes here 
		

	} else {
		stop <- TRUE
		print("Finished")
	}
}
