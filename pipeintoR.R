f <- file("stdin", "r")
#open(f)
x <- c()
y <- c()
print("Here first")
stop = FALSE
while (!stop){
	line = readLines(f, n=1)
	if (length(line) > 0){
		print(line)
		l <- strsplit(line, ',')[[1]]
		x <- c(x, as.numeric(l[1]))
		y <- c(y, as.numeric(l[2]))
	} else {
		stop = TRUE
		print("Finished")
	}
}
print(x)
print(y)
plot(x,y)
