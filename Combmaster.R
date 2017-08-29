combmaster <- function(dir){
library(plyr)
file_list <- list.files(path=dir,
                                recursive=T,
                                pattern="Master*"
                                ,full.names=T)
master.df <<- ldply(file_list, read.table, header=TRUE, sep="\t")
}