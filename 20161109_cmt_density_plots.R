library(plyr)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(RSQLite)


#function to load SQLite db file

	loadsql <- function(filename){



	sqlite.driver <- dbDriver("SQLite")
	db <- dbConnect(sqlite.driver,
                dbname = filename)

	mytable <<- dbReadTable(db, 'sequence')
	}

loadsql('/Users/cmtipto/Desktop/Working/Alessia/674/674_2017_02_17_16_10_24.db')


df <- cbind(substr(df[,1], 1, 1), df)
colnames(df)[1] <- 'Isotype'

df <- filter(df, Isotype != 'U')

df$Mutation_Rate <- df$Vquest_8_V.REGION.Nb.of.mutations/df$Vquest_8_V.REGION.Nb.of.nucleotides * 100


mu <- ddply(df, "Isotype", summarise, grp.median=median(Mutation_Rate))


colourCount = length(unique(df$Isotype))
getPalette = colorRampPalette(brewer.pal(4, "Set1"))(colourCount)


g <- ggplot(df, aes(x = Mutation_Rate, fill = Isotype)) + geom_density(alpha = .3, adjust = 5) + geom_vline(data=mu, aes(xintercept=grp.median, color=Isotype), linetype="dashed")
g <- g + theme_bw()
g <- g + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title = element_text(size = 16), axis.title.y = element_text(vjust = 1.5)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(name = "Mutation Rate", limits = c(0, 35)) + scale_y_continuous(name = "Sequence Density")
g + facet_grid(Isotype ~ .) + theme(strip.text.y = element_text(size = 16, angle = 0)) + scale_fill_manual(values = getPalette) + scale_color_manual(values = getPalette) + theme(legend.position="none")
