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


#1417
dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1417/1417_np_naive_withIdenticals.txt')
df.N.NP.1417 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1417/1417_pbl_naive_withIdenticals.txt')
df.N.PBL.1417 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1417/1417_np_pop3_A_withIdenticals.txt')
df.POP3.A.1417 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1417/1417_np_pop3_G_withIdenticals.txt')
df.POP3.G.1417 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1417/1417_np_pop3_E_withIdenticals.txt')
df.POP3.E.1417 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1417/1417_np_sm_withIdenticals.txt')
df.SM.NP.1417 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1417/1417_pbl_sm_withIdenticals.txt')
df.SM.PBL.1417 <- data.frame(dt)

df.N.NP.1417$Population <- rep('NP Naive')
df.N.PBL.1417$Population <- rep('PBL Naive')
df.POP3.A.1417$Population <- rep('NP IgA ASC')
df.POP3.G.1417$Population <- rep('NP IgG ASC')
df.POP3.E.1417$Population <- rep('NP IgE ASC')
df.SM.NP.1417$Population <- rep('NP Memory')
df.SM.PBL.1417$Population <- rep('PBL Memory')


#1325
dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1325/1325_NP_NAIVE_withIdenticals.txt')
df.N.NP.1325 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1325/1325_PBL_NAIVE_withIdenticals.txt')
df.N.PBL.1325 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1325/1325_NP_POP3_A_withIdenticals.txt')
df.POP3.A.1325 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1325/1325_NP_POP3_G_withIdenticals.txt')
df.POP3.G.1325 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1325/1325_NP_POP3_E_withIdenticals.txt')
df.POP3.E.1325 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1325/1325_PBL_SM_withIdenticals.txt')
df.SM.PBL.1325 <- data.frame(dt)

df.N.NP.1325$Population <- rep('NP Naive')
df.N.PBL.1325$Population <- rep('PBL Naive')
df.POP3.A.1325$Population <- rep('NP IgA ASC')
df.POP3.G.1325$Population <- rep('NP IgG ASC')
df.POP3.E.1325$Population <- rep('NP IgE ASC')
df.SM.PBL.1325$Population <- rep('PBL Memory')


#1789
dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1789/1789_NP_NAIVE_withIdenticals.txt')
df.N.NP.1789 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1789/1789_PBL_NAIVE_withIdenticals.txt')
df.N.PBL.1789 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1789/1789_NP_POP3_A_withIdenticals.txt')
df.POP3.A.1789 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1789/1789_NP_POP3_G_withIdenticals.txt')
df.POP3.G.1789 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1789/1789_NP_POP3_E_withIdenticals.txt')
df.POP3.E.1789 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1789/1789_NP_SM_withIdenticals.txt')
df.SM.NP.1789 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1789/1789_PBL_SM_withIdenticals.txt')
df.SM.PBL.1789 <- data.frame(dt)

df.N.NP.1789$Population <- rep('NP Naive')
df.N.PBL.1789$Population <- rep('PBL Naive')
df.POP3.A.1789$Population <- rep('NP IgA ASC')
df.POP3.G.1789$Population <- rep('NP IgG ASC')
df.POP3.E.1789$Population <- rep('NP IgE ASC')
df.SM.NP.1789$Population <- rep('NP Memory')
df.SM.PBL.1789$Population <- rep('PBL Memory')


#1792
dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1792/1792_NP_NAIVE_withIdenticals.txt')
df.N.NP.1792 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1792/1792_PBL_NAIVE_withIdenticals.txt')
df.N.PBL.1792 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1792/1792_NP_POP3_A_withIdenticals.txt')
df.POP3.A.1792 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1792/1792_NP_POP3_G_withIdenticals.txt')
df.POP3.G.1792 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1792/1792_NP_POP3_E_withIdenticals.txt')
df.POP3.E.1792 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1792/1792_NP_SM_withIdenticals.txt')
df.SM.NP.1792 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1792/1792_PBL_SM_withIdenticals.txt')
df.SM.PBL.1792 <- data.frame(dt)

df.N.NP.1792$Population <- rep('NP Naive')
df.N.PBL.1792$Population <- rep('PBL Naive')
df.POP3.A.1792$Population <- rep('NP IgA ASC')
df.POP3.G.1792$Population <- rep('NP IgG ASC')
df.POP3.E.1792$Population <- rep('NP IgE ASC')
df.SM.NP.1792$Population <- rep('NP Memory')
df.SM.PBL.1792$Population <- rep('PBL Memory')


#1809
dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1809/1809_NP_NAIVE_withIdenticals.txt')
df.N.NP.1809 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1809/1809_PBL_NAIVE_withIdenticals.txt')
df.N.PBL.1809 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1809/1809_NP_POP3_A_withIdenticals.txt')
df.POP3.A.1809 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1809/1809_NP_POP3_G_withIdenticals.txt')
df.POP3.G.1809 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1809/1809_NP_POP3_E_withIdenticals.txt')
df.POP3.E.1809 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1809/1809_NP_SM_withIdenticals.txt')
df.SM.NP.1809 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1809/1809_PBL_SM_withIdenticals.txt')
df.SM.PBL.1809 <- data.frame(dt)

df.N.NP.1809$Population <- rep('NP Naive')
df.N.PBL.1809$Population <- rep('PBL Naive')
df.POP3.A.1809$Population <- rep('NP IgA ASC')
df.POP3.G.1809$Population <- rep('NP IgG ASC')
df.POP3.E.1809$Population <- rep('NP IgE ASC')
df.SM.NP.1809$Population <- rep('NP Memory')
df.SM.PBL.1809$Population <- rep('PBL Memory')


#1095
dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1095/1095_np_naive_withIdenticals.txt')
df.N.NP.1095 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1095/1095_PBL_NAIVE_withIdenticals.txt')
df.N.PBL.1095 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1095/1095_np_pop3_A_withIdenticals.txt')
df.POP3.A.1095 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1095/1095_np_pop3_G_withIdenticals.txt')
df.POP3.G.1095 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1095/1095_np_pop3_E_withIdenticals.txt')
df.POP3.E.1095 <- data.frame(dt)

dt <- fread('/Users/cmtipto/OneDrive/Working/Alessia/1095/1095_PBL_SM_withIdenticals.txt')
df.SM.PBL.1095 <- data.frame(dt)

df.N.NP.1095$Population <- rep('NP Naive')
df.N.PBL.1095$Population <- rep('PBL Naive')
df.POP3.A.1095$Population <- rep('NP IgA ASC')
df.POP3.G.1095$Population <- rep('NP IgG ASC')
df.POP3.E.1095$Population <- rep('NP IgE ASC')
df.SM.PBL.1095$Population <- rep('PBL Memory')

df.1809 <- rbind(df.N.NP.1809, df.N.PBL.1809, df.POP3.A.1809, df.POP3.G.1809, df.POP3.E.1809, df.SM.NP.1809, df.SM.PBL.1809)
df.1792 <- rbind(df.N.NP.1792, df.N.PBL.1792, df.POP3.A.1792, df.POP3.G.1792, df.POP3.E.1792, df.SM.NP.1792, df.SM.PBL.1792)
df.1789 <- rbind(df.N.NP.1789, df.N.PBL.1789, df.POP3.A.1789, df.POP3.G.1789, df.POP3.E.1789, df.SM.NP.1789, df.SM.PBL.1789)
df.1325 <- rbind(df.N.NP.1325, df.N.PBL.1325, df.POP3.A.1325, df.POP3.G.1325, df.POP3.E.1325, df.SM.PBL.1325)
df.1417 <- rbind(df.N.NP.1417, df.N.PBL.1417, df.POP3.A.1417, df.POP3.G.1417, df.POP3.E.1417, df.SM.NP.1417, df.SM.PBL.1417)
df.674 <- rbind(df.N.NP.674, df.N.PBL.674, df.POP3.A.674, df.POP3.G.674, df.POP3.E.674, df.SM.NP.674, df.SM.PBL.674)
df.1095 <- rbind(df.N.NP.1095, df.N.PBL.1095, df.POP3.A.1095, df.POP3.G.1095, df.POP3.E.1095, df.SM.PBL.1095)

df.674$Subject <- rep('674')
df.1095$Subject <- rep('1095')
df.1417$Subject <- rep('1417')
df.1325$Subject <- rep('1325')
df.1789$Subject <- rep('1789')
df.1792$Subject <- rep('1792')
df.1809$Subject <- rep('1809')




df <- rbind(df.674, df.1095, df.1417, df.1325, df.1789, df.1792, df.1809)

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
