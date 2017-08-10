# location input requires file location of the grouped lineage file from IgSeq

seqanalysis <- function(location) {
    library(RColorBrewer)
    library(ggplot2)
    library(dplyr)
    library(plyr)
    library(data.table)
    # library(reshape)
    library(reshape2)
    library(spaa)
    library(RSQLite)
    library(stringr)

    # function to load SQLite db file

    loadsql <- function(filename) {



        sqlite.driver <- dbDriver("SQLite")
        db <- dbConnect(sqlite.driver, dbname = filename)

        df <<- dbReadTable(db, "sequence")
    }


    # establish input and output file locations

    location <- "~/Desktop/Test"

    location2 <- location
    location <- paste0(location, "/IgSeq")
    # file.dir <- paste0(location, '/lineageOutput/withIdenticals/') file.x <<-
    # paste0(file.dir, 'GROUPED_LINEAGE_COUNTS.txt')
    file.out <- paste0(location, "/CT_Sequencing_Analysis/")
    master.out <- paste0(file.out, "MasterData_", sub(".*/", "", location2), ".txt")
    dir.create(file.path(location, "/CT_Sequencing_Analysis"), showWarnings = FALSE)
    subject <- sub(".*/", "", location2)
    dblocation <- list.files(location, pattern = "\\.db$")
    dbfile <- paste0(location, "/", dblocation[1])
    print("Loading database file...")
    loadsql(dbfile)

    # Add columns of info
    config.location <- list.files(location2, pattern = "config")
    config.file <- paste0(location2, "/", config.location[1])
    config.info <- read.csv(config.file, sep = "\t", header = TRUE)
    time.course <- config.info$Time_Course[1]


    df$Population <- sub(".*_", "", df$population)
    df$Tissue <- str_match(df$population, "_(.*?)_")[, 2]

    if (length(levels(factor(df$Tissue))) != 1) {
        df$Population <- paste0(df$Tissue, "_", df$Population)
    }
    if (time.course == "Y") {
        df$Population <- sub(".*x", "", df$population)
    }
    df$Pop <- sub(".*_", "", df$population)
    df$Subject <- sub("_.*", "", df$population)
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) <
        tol
    df$Singleton <- !is.wholenumber(df$lineageID)



    # get rid of warnings because of data.table call

    # oldw <- getOption('warn') options(warn = -1)



    # establish function to change column names of grouped lineage file

    # changelin <- function(loc){ data <-read.table(loc, sep='\t', header=TRUE)
    # print(paste0('Reading: ', loc)) colnames(data)[1:4] <- c('Lineage', 'IGHV',
    # 'Pops', 'Seqs') for (i in 5:ncol(data)){ colnames(data)[i] <-
    # substr(colnames(data)[i], 21, nchar(colnames(data)[i])) } lineages.union <<-
    # data }

    # call function to change column names in the specified file changelin(file.x)



    # determine number of populations
    poplist <- levels(as.factor(df$Population))
    numpops <- length(poplist)
    print(paste0("Analyzing ", numpops, " populations of cells:"))
    poplist


    # poplist <- list() for (p in 1:numpops){ file.merged <- paste0(location,
    # '/mergeOutput') file.merged.list <- list.dirs(file.merged) j <- sub('.*/', '',
    # file.merged.list[p+1]) #jx <- sub('Run[0-9][0-9]_', '', j) jy <-
    # sub('_S[0-9]+', '', j) poplist[p] <- jy }





    # make data frame that will hold population names and info

    # master.data <- data.frame() seqlin <- data.frame() for (j in 1:numpops){
    # master.data[j,1] <- poplist[j] } colnames(master.data) <- c('Population')
    # master.data


    # Make Grouped Lineage Data Frame


    non.single <- df[df$Singleton == FALSE, ]
    lineage.data <- non.single[, c("Population", "lineageID", "Vgene")]
    table.x <- table(lineage.data$lineageID, lineage.data$Population)
    lineages.union <- as.data.frame.matrix(table.x)


    melt.x <- melt(lineages.union)
    filtered.melt.x <- melt.x[melt.x$value != 0, ]
    lineage.table <- table(filtered.melt.x$variable)
    melt.lineage.table <- as.data.frame(melt(lineage.table))
    colnames(melt.lineage.table)[2] <- "Lineages"

    vgene.table <- table(df$Population, df$Vgene)
    vgene.perc <- prop.table(vgene.table) * 100


    seqs <- as.data.frame(table(df$Population))
    master.data <- cbind(seqs, melt.lineage.table$Lineages)

    colnames(master.data) <- c("Population", "Sequences", "Lineages")


    if (length(levels(factor(df$Tissue))) != 1) {
        master.data$Tissue <- sub("_.*", "", master.data$Population)
    } else {
        master.data$Tissue <- "Unknown"
    }

    if (time.course == "Y") {
        master.data$Tissue <- str_match(master.data$Population, "_(.*?)_")[, 2]
    }

    print(master.data)

    # establish function to calculate clonality based on the Shannon diversity index
    # (1 - Pielou's Evenness)

    clon <- function(x, col = 1) {
        library(vegan)
        div <- diversity(x[, col], index = "shannon")
        H <- div/log(nrow(x))
        cl <<- 1 - H
        return(round(cl, 2))
    }

    # Add clonality values for each population into a new column


    clonality.df <- data.frame()

    for (k in 1:numpops) {
        clonality.df[k, 1] <- clon(lineages.union, col = k)
    }

    master.data <- cbind(master.data, clonality.df)
    colnames(master.data)[5] <- "Clonality"

    df$Mutation_Rate <- df$V_Mutations/df$V_Nucleotides * 100
    master.data$Mutation_Rate_Mean <- ddply(df, .(Population), summarize, mean = mean(Mutation_Rate))[,
        2]
    master.data$Mutation_Rate_Median <- ddply(df, .(Population), summarize, median = median(Mutation_Rate))[,
        2]




    iso.info <- data.frame()

    for (n in 1:numpops) {
        Ig <- df[df$isotype != "U", ]
        IgG <- df[df$isotype == "G", ]
        IgA <- df[df$isotype == "A", ]
        IgM <- df[df$isotype == "M", ]
        IgE <- df[df$isotype == "E", ]
        IgU <- df[df$isotype == "U", ]
        Ig <- Ig[Ig$Population == poplist[n], ]
        IgG <- IgG[IgG$Population == poplist[n], ]
        IgA <- IgA[IgA$Population == poplist[n], ]
        IgM <- IgM[IgM$Population == poplist[n], ]
        IgE <- IgE[IgE$Population == poplist[n], ]
        IgU <- IgU[IgU$Population == poplist[n], ]
        iso.info[n, 1] <- nrow(IgM)
        iso.info[n, 2] <- nrow(IgG)
        iso.info[n, 3] <- nrow(IgA)
        iso.info[n, 4] <- nrow(IgE)
        iso.info[n, 5] <- nrow(IgU)
        iso.info[n, 6] <- round((nrow(IgM)/nrow(Ig) * 100), 2)
        iso.info[n, 7] <- round((nrow(IgG)/nrow(Ig) * 100), 2)
        iso.info[n, 8] <- round((nrow(IgA)/nrow(Ig) * 100), 2)
        iso.info[n, 9] <- round((nrow(IgE)/nrow(Ig) * 100), 2)
        iso.info[n, 10] <- round((nrow(IgU)/nrow(df) * 100), 2)
        iso.info[n, 11] <- round(mean(IgM$V_Mutations), 2)
        iso.info[n, 12] <- round(mean(IgG$V_Mutations), 2)
        iso.info[n, 13] <- round(mean(IgA$V_Mutations), 2)
        iso.info[n, 14] <- round(mean(IgE$V_Mutations), 2)
        iso.info[n, 15] <- round(mean(IgU$V_Mutations), 2)
        iso.info[n, 16] <- round(median(IgM$V_Mutations), 2)
        iso.info[n, 17] <- round(median(IgG$V_Mutations), 2)
        iso.info[n, 18] <- round(median(IgA$V_Mutations), 2)
        iso.info[n, 19] <- round(median(IgE$V_Mutations), 2)
        iso.info[n, 20] <- round(median(IgU$V_Mutations), 2)
        iso.info[n, 21] <- round(mean(IgM$Mutation_Rate), 2)
        iso.info[n, 22] <- round(mean(IgG$Mutation_Rate), 2)
        iso.info[n, 23] <- round(mean(IgA$Mutation_Rate), 2)
        iso.info[n, 24] <- round(mean(IgE$Mutation_Rate), 2)
        iso.info[n, 25] <- round(mean(IgU$Mutation_Rate), 2)
        iso.info[n, 26] <- round(median(IgM$Mutation_Rate), 2)
        iso.info[n, 27] <- round(median(IgG$Mutation_Rate), 2)
        iso.info[n, 28] <- round(median(IgA$Mutation_Rate), 2)
        iso.info[n, 29] <- round(median(IgE$Mutation_Rate), 2)
        iso.info[n, 30] <- round(median(IgU$Mutation_Rate), 2)
    }

    colnames(iso.info) <- c("IgM Sequences", "IgG Sequences", "IgA Sequences", "IgE Sequences",
        "Unknown Isotype Sequences", "Percent IgM", "Percent IgG", "Percent IgA",
        "Percent IgE", "Percent Unknown Isotype", "IgM Mutations (Mean)", "IgG Mutations (Mean)",
        "IgA Mutations (Mean)", "IgE Mutations (Mean)", "Unknown Isotype Mutations (Mean)",
        "IgM Mutations (Median)", "IgG Mutations (Median)", "IgA Mutations (Median)",
        "IgE Mutations (Median)", "Unknown Isotype Mutations (Median)", "IgM Mutation Rate (Mean)",
        "IgG Mutation Rate (Mean)", "IgA Mutation Rate (Mean)", "IgE Mutation Rate (Mean)",
        "Unknown Isotype Mutation Rate (Mean)", "IgM Mutation Rate (Median)", "IgG Mutation Rate (Median)",
        "IgA Mutation Rate (Median)", "IgE Mutation Rate (Median)", "Unknown Isotype Mutation Rate (Median)")

    master.data <- cbind(master.data, iso.info)


    master.data$Disease <- config.info$Disease[1]
    master.data$Vaccine <- config.info$Vaccine[1]
    master.data$Vaccine_Time <- config.info$Vaccine_Time[1]


    vgene.table <- table(df$Population, df$Vgene)
    vgene.perc <- prop.table(vgene.table) * 100

    master.data$IGHV4_34 <- vgene.perc[, "IGHV4-34"]
    V434 <- df[df$Vgene == "IGHV4-34", ]

    substrRight <- function(x, n) {
        substr(x, nchar(x) - n + 1, nchar(x))
    }


    V434$AVY <- substrRight(V434$AA_FR1, 3)
    AVY.table <- table(V434$Population, V434$AVY)
    master.data$AVY_perc <- AVY.table[, "AVY"]/rowSums(AVY.table) * 100

















    # make a data frame with percent of matching clones in each

    percof.df <- data.frame()

    for (m in 1:numpops) {
        df.x <- filter(lineages.union, lineages.union[[m + 4]] > 0)
        lin.x <- nrow(df.x)
        for (n in 1:numpops) {
            df.x_y <- filter(df.x, df.x[[n + 4]] > 0)
            lin.x_y <- nrow(df.x_y)
            perc.x_y <- lin.x_y/lin.x * 100
            percof.df[m, n] <- round(perc.x_y, 2)
        }
        colnames(percof.df)[m] <- colnames(lineages.union)[m + 4]
    }

    # rename columns
    for (l in 1:numpops) {
        colnames(percof.df)[l] <- paste0("Percent_Matching_with_", colnames(percof.df)[l])
    }


    master.data <- cbind(master.data, percof.df)


    # Calculate the Morisita index for each pair of samples into a new matrix

    pops.df <- lineages.union[, 5:ncol(lineages.union)]
    mor <- niche.overlap(pops.df, method = c("morisita"))
    mor.mat <- as.matrix(round(mor, 3))
    for (i in 1:numpops) {
        mor.mat[i, i] <- 1
    }

    # rename columns

    mor.df <- as.data.frame(mor.mat)
    for (l in 1:(numpops)) {
        colnames(mor.mat)[l] <- paste0("MI_", colnames(mor.mat)[l])
    }



    # Add Morisita index values to master.data data frame

    master.data <- cbind(master.data, mor.mat)




    # Calculate the Percent of Each Isotype

    dir.create(file.path(location, "CT_Sequencing_Analysis/IsotypeInfo/"), showWarnings = FALSE)

    iso.out <- paste0(file.out, "IsotypeInfo/")

    dir.create(file.path(iso.out, "Isotype_Pie_Charts"), showWarnings = FALSE)
    out.isopie <- paste0(iso.out, "/Isotype_Pie_Charts/")


    iso.info <<- data.frame()

    for (n in 1:numpops) {
        merge.dir <- paste0(location, "/mergeOutput/")
        merge.list <- list.files(merge.dir)
        merge.file <- paste0(merge.dir, merge.list[n], "/", paste0(merge.list[n],
            "_withIdenticals.txt"))
        print(paste0("Reading: ", merge.list[n]))
        data.m <- read.table(merge.file, sep = "\t", header = TRUE)
        data.m <- cbind(substr(data.m[, 1], 1, 1), data.m)
        # data.m <- cbind(substr(rownames(data.m), 1, 1), data.m)
        data.m$Mutation.Rate <- data.m$Vquest_8_V.REGION.Nb.of.mutations/data.m$Vquest_8_V.REGION.Nb.of.nucleotides *
            100
        colnames(data.m)[1] <- "Isotype"
        Ig <- filter(data.m, Isotype != "U")
        IgG <- filter(data.m, Isotype == "G")
        IgA <- filter(data.m, Isotype == "A")
        IgM <- filter(data.m, Isotype == "M")
        IgE <- filter(data.m, Isotype == "E")
        IgU <- filter(data.m, Isotype == "U")
        iso.info[n, 1] <- nrow(IgM)
        iso.info[n, 2] <- nrow(IgG)
        iso.info[n, 3] <- nrow(IgA)
        iso.info[n, 4] <- nrow(IgE)
        iso.info[n, 5] <- nrow(IgU)
        iso.info[n, 6] <- round((nrow(IgM)/nrow(Ig) * 100), 2)
        iso.info[n, 7] <- round((nrow(IgG)/nrow(Ig) * 100), 2)
        iso.info[n, 8] <- round((nrow(IgA)/nrow(Ig) * 100), 2)
        iso.info[n, 9] <- round((nrow(IgE)/nrow(Ig) * 100), 2)
        iso.info[n, 10] <- round((nrow(IgU)/nrow(data.m) * 100), 2)
        iso.info[n, 11] <- round(mean(IgM$Vquest_8_V.REGION.Nb.of.mutations), 2)
        iso.info[n, 12] <- round(mean(IgG$Vquest_8_V.REGION.Nb.of.mutations), 2)
        iso.info[n, 13] <- round(mean(IgA$Vquest_8_V.REGION.Nb.of.mutations), 2)
        iso.info[n, 14] <- round(mean(IgE$Vquest_8_V.REGION.Nb.of.mutations), 2)
        iso.info[n, 15] <- round(mean(IgU$Vquest_8_V.REGION.Nb.of.mutations), 2)
        iso.info[n, 16] <- round(median(IgM$Vquest_8_V.REGION.Nb.of.mutations), 2)
        iso.info[n, 17] <- round(median(IgG$Vquest_8_V.REGION.Nb.of.mutations), 2)
        iso.info[n, 18] <- round(median(IgA$Vquest_8_V.REGION.Nb.of.mutations), 2)
        iso.info[n, 19] <- round(median(IgE$Vquest_8_V.REGION.Nb.of.mutations), 2)
        iso.info[n, 20] <- round(median(IgU$Vquest_8_V.REGION.Nb.of.mutations), 2)
        iso.info[n, 21] <- round(mean(IgM$Mutation.Rate), 2)
        iso.info[n, 22] <- round(mean(IgG$Mutation.Rate), 2)
        iso.info[n, 23] <- round(mean(IgA$Mutation.Rate), 2)
        iso.info[n, 24] <- round(mean(IgE$Mutation.Rate), 2)
        iso.info[n, 25] <- round(mean(IgU$Mutation.Rate), 2)
        iso.info[n, 26] <- round(median(IgM$Mutation.Rate), 2)
        iso.info[n, 27] <- round(median(IgG$Mutation.Rate), 2)
        iso.info[n, 28] <- round(median(IgA$Mutation.Rate), 2)
        iso.info[n, 29] <- round(median(IgE$Mutation.Rate), 2)
        iso.info[n, 30] <- round(median(IgU$Mutation.Rate), 2)


        a <- filter(Ig, Isotype == "A")
        g <- filter(Ig, Isotype == "G")
        e <- filter(Ig, Isotype == "E")
        m <- filter(Ig, Isotype == "M")

        # slices <- c(nrow(m), nrow(g), nrow(a), nrow(e)) lbs <- c('IgM', 'IgG', 'IgA',
        # 'IgE')

        # pdf(paste0(file.out.isopie, merge.list[n],'_IsotypePie.pdf'), width = 15,
        # height = 9) p <- pie(slices, labels = lbs, col=c('#b0923b', '#8960b3',
        # '#56ae6c', '#ba495b'), cex = 1.5)


        slices <- c(nrow(m), nrow(g), nrow(a))
        lbs <- c("IgM", "IgG", "IgA")

        pdf(paste0(out.isopie, merge.list[n], "_IsotypePie.pdf"), width = 7, height = 7)
        p <- pie(slices, labels = lbs, col = c("#bc5d41", "#84a955", "#965da7"),
            cex = 1.5)


        print(p)
        dev.off()

    }

    colnames(iso.info) <- c("IgM Sequences", "IgG Sequences", "IgA Sequences", "IgE Sequences",
        "Unknown Isotype Sequences", "Percent IgM", "Percent IgG", "Percent IgA",
        "Percent IgE", "Percent Unknown Isotype", "IgM Mutations (Mean)", "IgG Mutations (Mean)",
        "IgA Mutations (Mean)", "IgE Mutations (Mean)", "Unknown Isotype Mutations (Mean)",
        "IgM Mutations (Median)", "IgG Mutations (Median)", "IgA Mutations (Median)",
        "IgE Mutations (Median)", "Unknown Isotype Mutations (Median)", "IgM Mutation Rate (Mean)",
        "IgG Mutation Rate (Mean)", "IgA Mutation Rate (Mean)", "IgE Mutation Rate (Mean)",
        "Unknown Isotype Mutation Rate (Mean)", "IgM Mutation Rate (Median)", "IgG Mutation Rate (Median)",
        "IgA Mutation Rate (Median)", "IgE Mutation Rate (Median)", "Unknown Isotype Mutation Rate (Median)")

    master.data <- cbind(master.data, iso.info)



    # Make Isotype Bar Chart

    df.iso <- master.data[, c("Population", "Percent IgM", "Percent IgG", "Percent IgA")]
    colnames(df.iso) <- c("Population", "IgM", "IgG", "IgA")
    iso.melt <- melt(df.iso)

    colnames(iso.melt) <- c("Population", "Isotype", "Percent")


    # dir.create(file.path(file.out, 'Isotype_Stacked_Bar'), showWarnings = FALSE)
    # file.out.isobar <- paste0(file.out, '/Isotype_Stacked_Bar/')

    pdf(paste0(iso.out, subject, "_Isotype_StackedBar.pdf"))

    p <- ggplot()
    p <- p + geom_bar(data = iso.melt, aes_string(x = "Population", y = "Percent",
        fill = "Isotype"), stat = "identity", position = "stack")
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + ylab("Percent of Sequences") + xlab("")
    p <- p + theme(axis.text.x = element_text(size = 14, vjust = 0.5, angle = 90,
        hjust = 1), axis.text.y = element_text(size = 12), axis.title = element_text(size = 16),
        axis.title.y = element_text(vjust = 1.5))
    p <- p + scale_fill_manual(values = c("#bc5d41", "#84a955", "#965da7"))

    print(p)
    dev.off()



    # Make Morisita heat Map


    pdf(paste0(file.out, subject, "_Morisita_Heatmap.pdf"))
    st <- numpops + 5
    fi <- numpops + (st - 1)
    mor.x <- master.data[, c(1, st:fi)]
    dm <- melt(mor.x, id.vars = "Population")

    p <- ggplot(dm, aes(y = Population, x = variable))
    p <- p + geom_tile(aes(fill = value), colour = "black") + xlab("") + ylab("") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9,
            colour = "black"), axis.text.y = element_text(size = 9, colour = "black")) +
        scale_fill_gradient(low = "white", high = "steelblue") + geom_text(aes(label = round(value,
        2)))

    print(p)
    dev.off()

    # Stacked Barplot

    stacked <- function(data, minimum = 0, maximum = 105, pull = 1, border = "grey30") {

        # Input data should be a file exported from the combined rearrangements view of
        # Analyzer 3.0.  Minimum and maximum alter the scale of the y axis to focus on a
        # certain range of percentages.  Pull number gives column that the data is sorted
        # by before aa are pulled from.



        library(lattice)
        library(tidyr)

        # Sort based on pull value and identify the most abundant aa sequences
        temp <- data
        order.temp <- temp[order(temp[, (pull)], decreasing = TRUE), ]
        aa1 <- order.temp[1, 1]
        aa2 <- order.temp[2, 1]
        aa3 <- order.temp[3, 1]
        aa4 <- order.temp[4, 1]
        aa5 <- order.temp[5, 1]
        aa6 <- order.temp[6, 1]
        aa7 <- order.temp[7, 1]
        aa8 <- order.temp[8, 1]
        aa9 <- order.temp[9, 1]
        aa10 <- order.temp[10, 1]

        # Isolate data needed by removing the sum and present in columns
        df <- data[, c(1, 5:ncol(data))]

        df$Lineage <- factor(df$Lineage)

        # Change the data frame to long form
        df_long <- gather(df, sample, sequences, 2:ncol(df))


        # Identify and change to percent of total templates
        sampleNames <- unique(df_long$sample)

        for (i in 1:length(sampleNames)) {
            df_long$sequences[df_long$sample == sampleNames[i]] <- df_long$sequences[df_long$sample ==
                sampleNames[i]]/sum(df_long$sequences[df_long$sample == sampleNames[i]]) *
                100
        }

        # Identify aa's to be colored
        g <- order(df_long$sequences)
        if (pull == 1) {
            colVec <- ifelse(df_long$Lineage == aa1, "red", ifelse(df_long$Lineage ==
                aa2, "green", ifelse(df_long$Lineage == aa3, "purple", ifelse(df_long$Lineage ==
                aa4, "blue", ifelse(df_long$Lineage == aa5, "yellow", ifelse(df_long$Lineage ==
                aa6, "coral", ifelse(df_long$Lineage == aa7, "pink", ifelse(df_long$Lineage ==
                aa8, "navy", ifelse(df_long$Lineage == aa9, "darkgreen", ifelse(df_long$Lineage ==
                aa10, "maroon", "white"))))))))))[g]
        } else {

            colVec <- ifelse(df_long$Lineage == aa1, "red", ifelse(df_long$Lineage ==
                aa2, "green", ifelse(df_long$Lineage == aa3, "purple", ifelse(df_long$Lineage ==
                aa4, "blue", ifelse(df_long$Lineage == aa5, "yellow", ifelse(df_long$Lineage ==
                aa6, "coral", ifelse(df_long$Lineage == aa7, "pink", ifelse(df_long$Lineage ==
                aa8, "navy", ifelse(df_long$Lineage == aa9, "darkgreen", ifelse(df_long$Lineage ==
                aa10, "maroon", "white"))))))))))[g]
        }

        # Order the columns by name df_long <- df_long[order(df_long$sample),]

        # or specific order df_long$sample <- ordered(df_long$sample, levels =
        # c('PBMC_1', 'PBMC_2', 'PBMC_3', 'PBMC_4', 'R11LN_Fresh', 'S7LN_Fresh',
        # 'Post_Normal_FFPE', 'Post_Tumor_FFPE'))

        # Plot the stacked barplots
        barchart(sequences ~ sample, groups = factor(1:nrow(df_long), levels = g),
            stack = TRUE, data = df_long, border = border, col = colVec, ylim = c(minimum,
                maximum), ylab = "Cumulative Percent of Templates", scales = list(x = list(rot = 90)))
        # can use border = for border color change, box.width or box.ratio for width of
        # bars
    }

    dir.create(file.path(location, "/CT_Sequencing_Analysis/CloneStack/"), showWarnings = FALSE)

    stacked.out <- paste0(file.out, "CloneStack/")

    for (n in 5:(4 + numpops)) {
        c <- (n - 4)
        pdf(paste0(stacked.out, subject, "_pull-Pop-", c, "_Clone_StackedBar.pdf"))

        p <- stacked(lineages.union, pull = (n))

        print(p)
        dev.off()
    }

    pdf(paste0(stacked.out, subject, "_noPull_Clone_StackedBar.pdf"))

    p <- stacked(lineages.union)

    print(p)
    dev.off()

    # return master.data
    write.table(master.data, file = master.out, quote = FALSE, sep = "\t", row.names = FALSE)


    pdf(paste0(file.out, subject, "_Clonality.pdf"))

    p = ggplot()
    p = p + geom_bar(data = master.data, aes_string(x = "Population", y = "Clonality"),
        stat = "identity")
    p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p = p + scale_x_discrete(name = "") + scale_y_continuous(name = "Clonality")
    p = p + theme(axis.text.x = element_text(angle = 90, size = 14, vjust = 0.5),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 16),
        axis.title.y = element_text(vjust = 1.5))

    print(p)
    dev.off()


    dir.create(file.path(location, "/CT_Sequencing_Analysis/MutationInfo/"), showWarnings = FALSE)

    mut.out <- paste0(file.out, "MutationInfo/")

    master.data <- read.table(master.out, sep = "\t", header = TRUE)

    pdf(paste0(mut.out, subject, "_median_IgM_Mutation_rate.pdf"))

    p = ggplot()
    p = p + geom_bar(data = master.data, aes_string(x = "Population", y = "IgM.Mutation.Rate..Median."),
        stat = "identity")
    p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p = p + scale_x_discrete(name = "") + scale_y_continuous(name = "IgM Mutation Rate (Mean)",
        limits = c(0, 20))
    p = p + theme(axis.text.x = element_text(angle = 90, size = 14, vjust = 0.5),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 16),
        axis.title.y = element_text(vjust = 1.5))

    print(p)
    dev.off()

    pdf(paste0(mut.out, subject, "_median_IgG_Mutation_rate.pdf"))

    p = ggplot()
    p = p + geom_bar(data = master.data, aes_string(x = "Population", y = "IgG.Mutation.Rate..Median."),
        stat = "identity")
    p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p = p + scale_x_discrete(name = "") + scale_y_continuous(name = "IgG Mutation Rate (Mean)",
        limits = c(0, 20))
    p = p + theme(axis.text.x = element_text(angle = 90, size = 14, vjust = 0.5),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 16),
        axis.title.y = element_text(vjust = 1.5))

    print(p)
    dev.off()

    pdf(paste0(mut.out, subject, "_median_IgA_Mutation_rate.pdf"))

    p = ggplot()
    p = p + geom_bar(data = master.data, aes_string(x = "Population", y = "IgA.Mutation.Rate..Median."),
        stat = "identity")
    p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p = p + scale_x_discrete(name = "") + scale_y_continuous(name = "IgA Mutation Rate (Mean)",
        limits = c(0, 20))
    p = p + theme(axis.text.x = element_text(angle = 90, size = 14, vjust = 0.5),
        axis.text.y = element_text(size = 12), axis.title = element_text(size = 16),
        axis.title.y = element_text(vjust = 1.5))

    print(p)
    dev.off()

    if (is.na(master.data$IgE.Mutation.Rate..Median.[1]) == FALSE) {

        pdf(paste0(file.out, subject, "_median_IgE_Mutation_rate.pdf"))

        p = ggplot()
        p = p + geom_bar(data = master.data, aes_string(x = "Population", y = "IgE.Mutation.Rate..Median."),
            stat = "identity")
        p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        p = p + scale_x_discrete(name = "") + scale_y_continuous(name = "IgE Mutation Rate (Mean)",
            limits = c(0, 20))
        p = p + theme(axis.text.x = element_text(angle = 90, size = 14, vjust = 0.5),
            axis.text.y = element_text(size = 12), axis.title = element_text(size = 16),
            axis.title.y = element_text(vjust = 1.5))

        print(p)
        dev.off()

    }



    dir.create(file.path(location, "/CT_Sequencing_Analysis/ComparePops/"), showWarnings = FALSE)

    compare.out <- paste0(file.out, "ComparePops/")


    df.x <- lineages.union[, 5:(numpops + 4)]
    df.y <- prop.table(as.matrix(df.x), 2)
    df.z <- as.data.frame(df.y)

    zero <- 1 * 10^-6
    line0 <- 1 * 10^-6
    ticks <- c(zero, 1 * 10^-5, 1 * 10^-4, 1 * 10^-3, 1 * 10^-2, 1 * 10^-1, 1, 10,
        100)
    axisRange <- c(5 * 10^-7, 100)

    for (s in 1:numpops) {
        for (t in 1:numpops) {
            if (s != t) {
                df2 <- df.z[, c(s, t)]
                df.colnames <- colnames(df2)
                df2 <- filter(df2, df.colnames[1] > 0 | df.colnames[2] > 0)
                xx <- df2[, 1]
                yy <- df2[, 2]
                df2[xx == 0, df.colnames[1]] <- zero
                df2[yy == 0, df.colnames[2]] <- zero

                bar = colnames(df2)
                names(df2) <- c("col1", "col2")

                pdf(paste0(compare.out, subject, "_Compare-Pop", s, "_", t, ".pdf"))
                p = ggplot() + coord_fixed()
                p = p + geom_point(data = df2, aes_string(x = "col1", y = "col2"),
                  shape = 21, colour = "black")

                p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
                p = p + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12),
                  axis.title.y = element_text(vjust = 1.05), axis.title.x = element_text(vjust = 0.5))
                p = p + geom_hline(yintercept = c(line0), linetype = "dotted", size = 0.75) +
                  geom_vline(xintercept = c(line0), linetype = "dotted", size = 0.75)
                p = p + scale_x_log10(breaks = ticks, labels = c(0, ticks[2:length(ticks)]),
                  limits = axisRange)
                p = p + scale_y_log10(breaks = ticks, labels = c(0, ticks[2:length(ticks)]),
                  limits = axisRange)
                p = p + annotation_logticks() + labs(x = bar[1], y = bar[2])
                p
                print(p)
                dev.off()
            }
        }
    }


    dir.create(file.path(location, "/CT_Sequencing_Analysis/ComparePops/Zoom"), showWarnings = FALSE)

    compare.out <- paste0(file.out, "ComparePops/Zoom/")



    zero <- 1 * 10^-6
    line0 <- 1 * 10^-6
    ticks <- c(zero, 1 * 10^-5, 1 * 10^-4, 1 * 10^-3, 1 * 10^-2, 1 * 10^-1, 1, 10)
    axisRange <- c(5 * 10^-7, 10)

    for (s in 1:numpops) {
        for (t in 1:numpops) {
            if (s != t) {
                df2 <- df.z[, c(s, t)]
                df.colnames <- colnames(df2)
                df2 <- filter(df2, df.colnames[1] > 0 | df.colnames[2] > 0)
                xx <- df2[, 1]
                yy <- df2[, 2]
                df2[xx == 0, df.colnames[1]] <- zero
                df2[yy == 0, df.colnames[2]] <- zero

                bar = colnames(df2)
                names(df2) <- c("col1", "col2")

                pdf(paste0(compare.out, subject, "_Compare-Pop", s, "_", t, ".pdf"))
                p = ggplot() + coord_fixed()
                p = p + geom_point(data = df2, aes_string(x = "col1", y = "col2"),
                  shape = 21, colour = "black")

                p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
                p = p + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12),
                  axis.title.y = element_text(vjust = 1.05), axis.title.x = element_text(vjust = 0.5))
                p = p + geom_hline(yintercept = c(line0), linetype = "dotted", size = 0.75) +
                  geom_vline(xintercept = c(line0), linetype = "dotted", size = 0.75)
                p = p + scale_x_log10(breaks = ticks, labels = c(0, ticks[2:length(ticks)]),
                  limits = axisRange)
                p = p + scale_y_log10(breaks = ticks, labels = c(0, ticks[2:length(ticks)]),
                  limits = axisRange)
                p = p + annotation_logticks() + labs(x = bar[1], y = bar[2])
                p
                print(p)
                dev.off()
            }
        }
    }

    dir.create(file.path(location, "/CT_Sequencing_Analysis/ComparePops/Zoom/superZoom"),
        showWarnings = FALSE)

    compare.out <- paste0(file.out, "ComparePops/Zoom/superZoom/")



    zero <- 1 * 10^-6
    line0 <- 1 * 10^-6
    ticks <- c(zero, 1 * 10^-5, 1 * 10^-4, 1 * 10^-3, 1 * 10^-2, 1 * 10^-1, 1)
    axisRange <- c(5 * 10^-7, 1)

    for (s in 1:numpops) {
        for (t in 1:numpops) {
            if (s != t) {
                df2 <- df.z[, c(s, t)]
                df.colnames <- colnames(df2)
                df2 <- filter(df2, df.colnames[1] > 0 | df.colnames[2] > 0)
                xx <- df2[, 1]
                yy <- df2[, 2]
                df2[xx == 0, df.colnames[1]] <- zero
                df2[yy == 0, df.colnames[2]] <- zero

                bar = colnames(df2)
                names(df2) <- c("col1", "col2")

                pdf(paste0(compare.out, subject, "_Compare-Pop", s, "_", t, ".pdf"))
                p = ggplot() + coord_fixed()
                p = p + geom_point(data = df2, aes_string(x = "col1", y = "col2"),
                  shape = 21, colour = "black")

                p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
                p = p + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12),
                  axis.title.y = element_text(vjust = 1.05), axis.title.x = element_text(vjust = 0.5))
                p = p + geom_hline(yintercept = c(line0), linetype = "dotted", size = 0.75) +
                  geom_vline(xintercept = c(line0), linetype = "dotted", size = 0.75)
                p = p + scale_x_log10(breaks = ticks, labels = c(0, ticks[2:length(ticks)]),
                  limits = axisRange)
                p = p + scale_y_log10(breaks = ticks, labels = c(0, ticks[2:length(ticks)]),
                  limits = axisRange)
                p = p + annotation_logticks() + labs(x = bar[1], y = bar[2])
                p
                print(p)
                dev.off()
            }
        }
    }

    print("Analysis Done")



    # g <- tableGrob(master.data[c(1, 7:ncol(master.data))], rows = NULL) g <-
    # gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)), t = 1, b =
    # nrow(g), l = 1, r = ncol(g)) plot.new() grid.arrange(g)

    options(warn = oldw)

}
