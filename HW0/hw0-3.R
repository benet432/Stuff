#a
snip="Hello, my name is Bob. I am a statistician. I like statistics very much."
split=unlist(strsplit(snip, ""))
j=1:length(split)
for (i in 1:9) {j[i]=paste("0", j[i], sep="")}
for(i in 1:length(split)) {write(split[i], paste('out_',j[i],'.txt', sep=""))}

#b
filelist= list.files(pattern = ".*.txt") 
allfile= lapply(filelist,function(x){try(read.table(x,colClasses='character'))}) 
empty = sapply(allfile, function(x){inherits(x, "try-error")})
char= character(72)
char[empty]=" "
char[!empty]=allfile[!empty] 
write(paste(do.call(rbind, char) $V1, collapse = ''), "newfile.txt")




