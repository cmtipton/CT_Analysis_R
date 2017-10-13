
CTseq <- function(directory){

folders <- list.files(directory)

for (i in 1:length(folders)){
	location <- paste0(directory, '/', folders[i])
	seqanalysis(location)
	}
}