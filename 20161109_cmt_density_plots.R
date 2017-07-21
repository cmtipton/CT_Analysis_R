library(plyr)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(data.table)


#674
dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/674/674_NP_NAIVE_withIdenticals.txt')
df.N.NP.674 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/674/674_PBL_NAIVE_withIdenticals.txt')
df.N.PBL.674 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/674/674_NP_POP3_A_withIdenticals.txt')
df.POP3.A.674 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/674/674_NP_POP3_G_withIdenticals.txt')
df.POP3.G.674 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/674/674_NP_POP3_E_withIdenticals.txt')
df.POP3.E.674 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/674/674_NP_SM_withIdenticals.txt')
df.SM.NP.674 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/674/674_PBL_SM_withIdenticals.txt')
df.SM.PBL.674 <- data.frame(dt)

df.N.NP.674$Population <- rep('NP Naive')
df.N.PBL.674$Population <- rep('PBL Naive')
df.POP3.A.674$Population <- rep('NP IgA ASC')
df.POP3.G.674$Population <- rep('NP IgG ASC')
df.POP3.E.674$Population <- rep('NP IgE ASC')
df.SM.NP.674$Population <- rep('NP Memory')
df.SM.PBL.674$Population <- rep('PBL Memory')

#1097
dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1097/1097_NP_NAIVE_withIdenticals.txt')
df.N.NP.1097 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1097/1097_PBL_NAIVE_withIdenticals.txt')
df.N.PBL.1097 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1097/1097_NP_POP3_A_withIdenticals.txt')
df.POP3.A.1097 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1097/1097_NP_POP3_G_withIdenticals.txt')
df.POP3.G.1097 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1097/1097_NP_POP3_E_withIdenticals.txt')
df.POP3.E.1097 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1097/1097_PBL_SM_withIdenticals.txt')
df.SM.PBL.1097 <- data.frame(dt)

df.N.NP.1097$Population <- rep('NP Naive')
df.N.PBL.1097$Population <- rep('PBL Naive')
df.POP3.A.1097$Population <- rep('NP IgA ASC')
df.POP3.G.1097$Population <- rep('NP IgG ASC')
df.POP3.E.1097$Population <- rep('NP IgE ASC')
df.SM.PBL.1097$Population <- rep('PBL Memory')




df.674 <- rbind(df.N.NP.674, df.N.PBL.674, df.POP3.A.674, df.POP3.G.674, df.POP3.E.674, df.SM.NP.674, df.SM.PBL.674)
df.1097 <- rbind(df.N.NP.1097, df.N.PBL.1097, df.POP3.A.1097, df.POP3.G.1097, df.POP3.E.1097, df.SM.PBL.1097)

df.674$Subject <- rep('674')
df.1097$Subject <- rep('1097')

df <- rbind(df.674, df.1097)

df$Mutation_Rate <- df$Vquest_8_V.REGION.Nb.of.mutations/df$Vquest_8_V.REGION.Nb.of.nucleotides *
    100

mu <- ddply(df, c('Population', 'Subject'), summarise, grp.median=median(Mutation_Rate))


colourCount = length(unique(df$Subject))
getPalette = colorRampPalette(brewer.pal(4, "Set1"))(colourCount)


g <- ggplot(df, aes(x = Mutation_Rate, fill = Subject)) + geom_density(alpha = 0.3,
    adjust = 5) + geom_vline(data = mu, aes(xintercept = grp.median, color = Subject),
    linetype = "dashed")
g <- g + theme_bw()
g <- g + scale_x_continuous(name = "Mutation Rate", limits = c(0, 30)) + scale_y_continuous(name = "Sequence Density")
g + facet_grid(Population ~ .) + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 20), axis.title.y = element_text(vjust = 1.5),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(size = 16,
        angle = 0)) + scale_fill_manual(values = getPalette) + scale_color_manual(values = getPalette) +
    theme(legend.position = "none")
