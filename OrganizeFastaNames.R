
library(plyr)
library(stringr)

sample_list <- list.files(path="/Users/cmtipto/Desktop/Analysis-Done/",
                                recursive=T,
                                pattern="*.txz"
                                ,full.names=T)


a <- str_split_fixed(sample_list, "_", 5)
b <- str_split_fixed(a[,1], "/", 8)
c <- a[,c(2:4)]
d <- str_split_fixed(a[,5], "\\.", 2)

b1 <- as.data.frame(b)
c1 <- as.data.frame(c)
d1 <- as.data.frame(d)

#colnames(b1)[8] <- "Run"
#colnames(c1) <- c("Subject","Tissue", "Population")
#colnames(d1)[1] <- "Lane"

df1 <- cbind(b1[,8],c1)
df <- cbind(df1, d1[,1])
colnames(df) <- c("Run","Subject","Tissue","Population", "Lane")
#df <- merge(df1,d1)[,c("Run","Subject","Tissue","Population","Lane")]
