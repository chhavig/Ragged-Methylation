library(jsonlite)

f <- file("stdin", "r")

df <- stream_in(f)
v <- c()

for (s in colnames(df)){
	#s = colnames(df)[i]
	v <- cbind(v, unlist(df[s][[1]][[1]]))
}
print(v)

v <- data.frame(v)

colnames(v) = colnames(df)

print(v)
