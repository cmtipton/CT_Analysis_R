3.72# location input requires file location of the grouped lineage file from IgSeq

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

        full.data <<- dbReadTable(db, "sequence")
    }

    # establish input and output file locations
    location <- '/Users/cmtipto/Desktop/Test/Test'
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
    loading <- paste0("Loading database file for subject ", subject, "...")
    print(loading)
    loadsql(dbfile)

    # Add columns of info
    #config.location <- list.files(location2, pattern = "config")
    #config.file <- paste0(location2, "/", config.location[1])
    #config.info <- read.csv(config.file, sep = "\t", header = TRUE)
    #time.course <- config.info$Time_Course[1]
    under <- str_count(full.data$population[1],"\\_")
    dash <- str_count(full.data$population[1],"\\-")
    nolane <- str_count(full.data$population[1],"\\_S[0-9]")

    if (under == 4){
        sample.info <- str_split_fixed(full.data$population, "_", 5)
        sample.info <- as.data.frame(sample.info)
        colnames(sample.info) <- c('Run', 'Subject', 'Tissue', 'Pop', 'Lane')
    } else if (dash == 4) {
        sample.info <- str_split_fixed(full.data$population, "-", 5)
        sample.info <- as.data.frame(sample.info)
        colnames(sample.info) <- c('Run', 'Subject', 'Tissue', 'Pop', 'Lane')
    } else if (nolane == 1 & under == 3) {
        sample.info <- str_split_fixed(full.data$population, "_", 4)
        sample.info <- as.data.frame(sample.info)
        colnames(sample.info) <- c('Run', 'Subject', 'Tissue', 'Pop')
    } else if (nolane == 0 & dash == 3) {
        print('Error: No lane detected')
        sample.info <- str_split_fixed(full.data$population, "_", 4)
        sample.info <- as.data.frame(sample.info)
        colnames(sample.info) <- c('Run', 'Subject', 'Tissue', 'Pop')
    } else if (nolane == 0 & under == 3) {
        print('Error: No lane detected')
        sample.info <- str_split_fixed(full.data$population, "_", 4)
        sample.info <- as.data.frame(sample.info)
        colnames(sample.info) <- c('Run', 'Subject', 'Tissue', 'Pop')
    } else {
        print('Error with input nomenclature')
        sample.info <- str_split_fixed(full.data$population, "_", 4)
        sample.info <- as.data.frame(sample.info)
        colnames(sample.info) <- c('Run', 'Subject', 'Tissue', 'Pop')
    }

    #full.data$Population <- sub(".*_", "", full.data$population)
    #full.data$Tissue <- str_match(full.data$population, "_(.*?)_")[, 2]
    full.data.backup <- full.data
    full.data <- cbind(sample.info, full.data)
    full.data$Population <- paste0(full.data$Tissue, "_", full.data$Pop)

    #if (time.course == "Y") {
    #    full.data$Population <- sub(".*x", "", full.data$population)
    #}
    #full.data$Pop <- sub(".*_", "", full.data$Population)
    #full.data$Subject <- sub("_.*", "", full.data$population)
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) <
        tol
    full.data$Singleton <- !is.wholenumber(full.data$lineageID)



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
    poplist <- levels(as.factor(full.data$population))
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


    non.single <- full.data[full.data$Singleton == FALSE, ]

    Sample.df <- table(full.data$Run, full.data$Lane, full.data$Pop, full.data$Tissue)
    Sample.df <- melt(Sample.df)
    Sample.df <- as.data.frame(Sample.df)
    Sample.df <- Sample.df[Sample.df$value != 0,]
    colnames(Sample.df) <- c('Run', 'Lane', 'Pop', 'Tissue', 'Sequences')

    lineage.data <- non.single[, c("population", "lineageID", "Vgene")]
    table.x <- table(lineage.data$lineageID, lineage.data$population)
    lineages.union <- as.data.frame.matrix(table.x)


    melt.x <- melt(lineages.union)
    filtered.melt.x <- melt.x[melt.x$value != 0, ]
    lineage.table <- table(filtered.melt.x$variable)
    melt.lineage.table <- as.data.frame(melt(lineage.table))
    colnames(melt.lineage.table)[2] <- "Lineages"

    vtable <- table(full.data$population, full.data$Vgene)
    vperc <- prop.table(vtable) * 100

    jtable <- table(full.data$population, full.data$JGene)
    jperc <- prop.table(jtable) * 100

    seqs <- as.data.frame(table(non.single$population))
    master.data <- cbind(seqs, melt.lineage.table$Lineages)

    colnames(master.data) <- c("Population", "Sequences (Non-Single)", "Lineages")
    master.data <- cbind(Sample.df, master.data)

    #if (length(levels(factor(full.data$Tissue))) != 1) {
    #    master.data$Tissue <- sub("_.*", "", master.data$Population)
    # else {
    #    master.data$Tissue <- "Unknown"
    #}

    #if (time.course == "Y") {
    #    master.data$Tissue <- str_match(master.data$Population, "_(.*?)_")[, 2]
    #}

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

    colnames(clonality.df)[1] <- "Clonality"

    master.data <- cbind(master.data, clonality.df)


    full.data$Mutation_Rate <- full.data$V_Mutations/full.data$V_Nucleotides * 100
    master.data$Mutation_Rate_Mean <- ddply(full.data, .(population), summarize, mean = mean(Mutation_Rate))[,
        2]
    master.data$Mutation_Rate_Median <- ddply(full.data, .(population), summarize, median = median(Mutation_Rate))[,
        2]




    iso.info <- data.frame()

    for (n in 1:numpops) {
        Ig <- full.data[full.data$isotype != "U", ]
        IgG <- full.data[full.data$isotype == "G", ]
        IgA <- full.data[full.data$isotype == "A", ]
        IgM <- full.data[full.data$isotype == "M", ]
        IgE <- full.data[full.data$isotype == "E", ]
        IgU <- full.data[full.data$isotype == "U", ]
        Ig <- Ig[Ig$population == poplist[n], ]
        IgG <- IgG[IgG$population == poplist[n], ]
        IgA <- IgA[IgA$population == poplist[n], ]
        IgM <- IgM[IgM$population == poplist[n], ]
        IgE <- IgE[IgE$population == poplist[n], ]
        IgU <- IgU[IgU$population == poplist[n], ]
        iso.info[n, 1] <- nrow(IgM)
        iso.info[n, 2] <- nrow(IgG)
        iso.info[n, 3] <- nrow(IgA)
        iso.info[n, 4] <- nrow(IgE)
        iso.info[n, 5] <- nrow(IgU)
        iso.info[n, 6] <- round((nrow(IgM)/nrow(Ig) * 100), 2)
        iso.info[n, 7] <- round((nrow(IgG)/nrow(Ig) * 100), 2)
        iso.info[n, 8] <- round((nrow(IgA)/nrow(Ig) * 100), 2)
        iso.info[n, 9] <- round((nrow(IgE)/nrow(Ig) * 100), 2)
        iso.info[n, 10] <- round((nrow(IgU)/nrow(full.data) * 100), 2)
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


    non.single$Mutation_Rate <- non.single$V_Mutations/non.single$V_Nucleotides * 100
    master.data$Mutation_Rate_Mean_NS <- ddply(non.single, .(population), summarize, mean = mean(Mutation_Rate))[,
        2]
    master.data$Mutation_Rate_Median_NS <- ddply(non.single, .(population), summarize, median = median(Mutation_Rate))[,
        2]


    nonsingle.iso.info <- data.frame()

    for (n in 1:numpops) {
        Ig <- non.single[non.single$isotype != "U", ]
        IgG <- non.single[non.single$isotype == "G", ]
        IgA <- non.single[non.single$isotype == "A", ]
        IgM <- non.single[non.single$isotype == "M", ]
        IgE <- non.single[non.single$isotype == "E", ]
        IgU <- non.single[non.single$isotype == "U", ]
        Ig <- Ig[Ig$population == poplist[n], ]
        IgG <- IgG[IgG$population == poplist[n], ]
        IgA <- IgA[IgA$population == poplist[n], ]
        IgM <- IgM[IgM$population == poplist[n], ]
        IgE <- IgE[IgE$population == poplist[n], ]
        IgU <- IgU[IgU$population == poplist[n], ]
        nonsingle.iso.info[n, 1] <- nrow(IgM)
        nonsingle.iso.info[n, 2] <- nrow(IgG)
        nonsingle.iso.info[n, 3] <- nrow(IgA)
        nonsingle.iso.info[n, 4] <- nrow(IgE)
        nonsingle.iso.info[n, 5] <- nrow(IgU)
        nonsingle.iso.info[n, 6] <- round((nrow(IgM)/nrow(Ig) * 100), 2)
        nonsingle.iso.info[n, 7] <- round((nrow(IgG)/nrow(Ig) * 100), 2)
        nonsingle.iso.info[n, 8] <- round((nrow(IgA)/nrow(Ig) * 100), 2)
        nonsingle.iso.info[n, 9] <- round((nrow(IgE)/nrow(Ig) * 100), 2)
        nonsingle.iso.info[n, 10] <- round((nrow(IgU)/nrow(non.single) * 100), 2)
        nonsingle.iso.info[n, 11] <- round(mean(IgM$V_Mutations), 2)
        nonsingle.iso.info[n, 12] <- round(mean(IgG$V_Mutations), 2)
        nonsingle.iso.info[n, 13] <- round(mean(IgA$V_Mutations), 2)
        nonsingle.iso.info[n, 14] <- round(mean(IgE$V_Mutations), 2)
        nonsingle.iso.info[n, 15] <- round(mean(IgU$V_Mutations), 2)
        nonsingle.iso.info[n, 16] <- round(median(IgM$V_Mutations), 2)
        nonsingle.iso.info[n, 17] <- round(median(IgG$V_Mutations), 2)
        nonsingle.iso.info[n, 18] <- round(median(IgA$V_Mutations), 2)
        nonsingle.iso.info[n, 19] <- round(median(IgE$V_Mutations), 2)
        nonsingle.iso.info[n, 20] <- round(median(IgU$V_Mutations), 2)
        nonsingle.iso.info[n, 21] <- round(mean(IgM$Mutation_Rate), 2)
        nonsingle.iso.info[n, 22] <- round(mean(IgG$Mutation_Rate), 2)
        nonsingle.iso.info[n, 23] <- round(mean(IgA$Mutation_Rate), 2)
        nonsingle.iso.info[n, 24] <- round(mean(IgE$Mutation_Rate), 2)
        nonsingle.iso.info[n, 25] <- round(mean(IgU$Mutation_Rate), 2)
        nonsingle.iso.info[n, 26] <- round(median(IgM$Mutation_Rate), 2)
        nonsingle.iso.info[n, 27] <- round(median(IgG$Mutation_Rate), 2)
        nonsingle.iso.info[n, 28] <- round(median(IgA$Mutation_Rate), 2)
        nonsingle.iso.info[n, 29] <- round(median(IgE$Mutation_Rate), 2)
        nonsingle.iso.info[n, 30] <- round(median(IgU$Mutation_Rate), 2)
    }

    colnames(nonsingle.iso.info) <- c("IgM Sequences_NS", "IgG Sequences_NS", "IgA Sequences_NS", "IgE Sequences_NS",
        "Unknown Isotype Sequences_NS", "Percent IgM_NS", "Percent IgG_NS", "Percent IgA_NS",
        "Percent IgE_NS", "Percent Unknown Isotype_NS", "IgM Mutations (Mean)_NS", "IgG Mutations (Mean)_NS",
        "IgA Mutations (Mean)_NS", "IgE Mutations (Mean)_NS", "Unknown Isotype Mutations (Mean)_NS",
        "IgM Mutations (Median)_NS", "IgG Mutations (Median)_NS", "IgA Mutations (Median)_NS",
        "IgE Mutations (Median)_NS", "Unknown Isotype Mutations (Median)_NS", "IgM Mutation Rate (Mean)_NS",
        "IgG Mutation Rate (Mean)_NS", "IgA Mutation Rate (Mean)_NS", "IgE Mutation Rate (Mean)_NS",
        "Unknown Isotype Mutation Rate (Mean)_NS", "IgM Mutation Rate (Median)_NS", "IgG Mutation Rate (Median)_NS",
        "IgA Mutation Rate (Median)_NS", "IgE Mutation Rate (Median)_NS", "Unknown Isotype Mutation Rate (Median)_NS")

    master.data <- cbind(master.data, nonsingle.iso.info)


#     ## Summarizes data.
#     ## Gives count, mean, standard deviation, standard error of the mean, and confidence
#     ## interval (default 95%).
#     ##   data: a data frame.
#     ##   measurevar: the name of a column that contains the variable to be summariezed
#     ##   groupvars: a vector containing names of columns that contain grouping variables
#     ##   na.rm: a boolean that indicates whether to ignore NA's
#     ##   conf.interval: the percent range of the confidence interval (default is 95%)
#     summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95) {
#         library(doBy)
#
#         # New version of length which can handle NA's: if na.rm==T, don't count them
#         length2 <- function (x, na.rm=FALSE) {
#             if (na.rm) sum(!is.na(x))
#             else       length(x)
#         }
#
#         # Collapse the data
#         formula <- as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "))
#         datac <- summaryBy(formula, data=data, FUN=c(length2,mean,median), na.rm=na.rm)
#
#         # Rename columns
#         names(datac)[ names(datac) == paste(measurevar, ".mean",    sep="") ] <- "mean"
#         names(datac)[ names(datac) == paste(measurevar, ".median",      sep="") ] <- "median"
#         names(datac)[ names(datac) == paste(measurevar, ".length2", sep="") ] <- "N"
#
#         return(datac)
#     }
#
# tgc <- summarySE(non.single, measurevar="Mutation_Rate", groupvars=c("Population", "LineageID","isotype"))


full.data$V_SilentMutation_Rate <- full.data$V_SilentMutations / full.data$V_Nucleotides * 100
full.data$V_NonsilentMutation_Rate <- full.data$V_NonsilentMutations / full.data$V_Nucleotides * 100
full.data$FR1_Mutation_Rate <- full.data$FR1_Mutations / full.data$FR1_Nucleotides * 100
full.data$FR1_SilentMutation_Rate <- full.data$FR1_SilentMutations / full.data$FR1_Nucleotides * 100
full.data$FR1_NonsilentMutation_Rate <- full.data$FR1_NonsilentMutations / full.data$FR1_Nucleotides * 100
full.data$CDR1_Mutation_Rate <- full.data$CDR1_Mutations / full.data$CDR1_Nucleotides * 100
full.data$CDR1_SilentMutation_Rate <- full.data$CDR1_SilentMutations / full.data$CDR1_Nucleotides * 100
full.data$CDR1_NonsilentMutation_Rate <- full.data$CDR1_NonsilentMutations / full.data$CDR1_Nucleotides * 100
full.data$FR2_Mutation_Rate <- full.data$FR2_Mutations / full.data$FR2_Nucleotides * 100
full.data$FR2_SilentMutation_Rate <- full.data$FR2_SilentMutations / full.data$FR2_Nucleotides * 100
full.data$FR2_NonsilentMutation_Rate <- full.data$FR2_NonsilentMutations / full.data$FR2_Nucleotides * 100
full.data$CDR2_Mutation_Rate <- full.data$CDR2_Mutations / full.data$CDR2_Nucleotides * 100
full.data$CDR2_SilentMutation_Rate <- full.data$CDR2_SilentMutations / full.data$CDR2_Nucleotides * 100
full.data$CDR2_NonsilentMutation_Rate <- full.data$CDR2_NonsilentMutations / full.data$CDR2_Nucleotides * 100
full.data$FR3_Mutation_Rate <- full.data$FR3_Mutations / full.data$FR3_Nucleotides * 100
full.data$FR3_SilentMutation_Rate <- full.data$FR3_SilentMutations / full.data$FR3_Nucleotides * 100
full.data$FR3_NonsilentMutation_Rate <- full.data$FR3_NonsilentMutations / full.data$FR3_Nucleotides * 100

extended.mut.info <- data.frame()

for (n in 1:numpops) {
    df.em <- full.data[full.data$population == poplist[n], ]
    extended.mut.info[n, 1] <- round(mean(df.em$V_SilentMutation_Rate, na.rm = TRUE), 2)
    extended.mut.info[n, 2] <- round(mean(df.em$V_NonsilentMutation_Rate, na.rm = TRUE), 2)
    extended.mut.info[n, 3] <- round(mean(df.em$FR1_Mutation_Rate, na.rm = TRUE), 2)
    extended.mut.info[n, 4] <- round(mean(df.em$FR1_SilentMutation_Rate, na.rm = TRUE), 2)
    extended.mut.info[n, 5] <- round(mean(df.em$FR1_NonsilentMutation_Rate, na.rm = TRUE), 2)
    extended.mut.info[n, 6] <- round(mean(df.em$CDR1_Mutation_Rate, na.rm = TRUE), 2)
    extended.mut.info[n, 7] <- round(mean(df.em$CDR1_SilentMutation_Rate, na.rm = TRUE), 2)
    extended.mut.info[n, 8] <- round(mean(df.em$CDR1_NonsilentMutation_Rate, na.rm = TRUE), 2)
    extended.mut.info[n, 9] <- round(mean(df.em$FR2_Mutation_Rate, na.rm = TRUE), 2)
    extended.mut.info[n, 10] <- round(mean(df.em$FR2_SilentMutation_Rate, na.rm = TRUE), 2)
    extended.mut.info[n, 11] <- round(mean(df.em$FR2_NonsilentMutation_Rate, na.rm = TRUE), 2)
    extended.mut.info[n, 12] <- round(mean(df.em$CDR2_Mutation_Rate, na.rm = TRUE), 2)
    extended.mut.info[n, 13] <- round(mean(df.em$CDR2_SilentMutation_Rate, na.rm = TRUE), 2)
    extended.mut.info[n, 14] <- round(mean(df.em$CDR2_NonsilentMutation_Rate, na.rm = TRUE), 2)
    extended.mut.info[n, 15] <- round(mean(df.em$FR3_Mutation_Rate, na.rm = TRUE), 2)
    extended.mut.info[n, 16] <- round(mean(df.em$FR3_SilentMutation_Rate, na.rm = TRUE), 2)
    extended.mut.info[n, 17] <- round(mean(df.em$FR3_NonsilentMutation_Rate, na.rm = TRUE), 2)
}

    colnames(extended.mut.info) <- c("V_SilentMutation_Rate", "V_NonsilentMutation_Rate", "FR1_Mutation_Rate", "FR1_SilentMutation_Rate", "FR1_NonsilentMutation_Rate", "CDR1_Mutation_Rate", "CDR1_SilentMutation_Rate", "CDR1_NonsilentMutation_Rate", "FR2_Mutation_Rate", "FR2_SilentMutation_Rate", "FR2_NonsilentMutation_Rate", "CDR2_Mutation_Rate", "CDR2_SilentMutation_Rate", "CDR2_NonsilentMutation_Rate", "FR3_Mutation_Rate", "FR3_SilentMutation_Rate", "FR3_NonsilentMutation_Rate")

    master.data <- cbind(master.data, extended.mut.info)
    #master.data$Disease <- config.info$Disease[1]
    #master.data$Vaccine <- config.info$Vaccine[1]
    #master.data$Vaccine_Time <- config.info$Vaccine_Time[1]


    vgene.table <- table(full.data$population, full.data$Vgene)
    vgene.perc <- prop.table(vgene.table) * 100

    master.data$IGHV4_34 <- vgene.perc[, "IGHV4-34"]
    V434 <- full.data[full.data$Vgene == "IGHV4-34", ]

    substrRight <- function(x, n) {
        substr(x, nchar(x) - n + 1, nchar(x))
    }


    V434$AVY <- substrRight(V434$AA_FR1, 3)
    AVY.table <- table(V434$population, V434$AVY)
    master.data$AVY_perc <- AVY.table[, "AVY"]/rowSums(AVY.table) * 100


#Calculate the percent of each VH gene and JH gene

    gene.info <- data.frame()

    for (n in 1:numpops) {
        df.popgenes <- full.data[full.data$population == poplist[n], ]
        IGHV1.2 <- df.popgenes[df.popgenes$Vgene == 'IGHV1-2', ]
        IGHV1.3 <- df.popgenes[df.popgenes$Vgene == 'IGHV1-3', ]
        IGHV1.8 <- df.popgenes[df.popgenes$Vgene == 'IGHV1-8', ]
        IGHV1.18 <- df.popgenes[df.popgenes$Vgene == 'IGHV1-18', ]
        IGHV1.24 <- df.popgenes[df.popgenes$Vgene == 'IGHV1-24', ]
        IGHV1.45 <- df.popgenes[df.popgenes$Vgene == 'IGHV1-45', ]
        IGHV1.46 <- df.popgenes[df.popgenes$Vgene == 'IGHV1-46', ]
        IGHV1.58 <- df.popgenes[df.popgenes$Vgene == 'IGHV1-58', ]
        IGHV1.69 <- df.popgenes[df.popgenes$Vgene == 'IGHV1-69', ]
        IGHV1.f <- df.popgenes[df.popgenes$Vgene == 'IGHV1-f', ]
        IGHV2.5 <- df.popgenes[df.popgenes$Vgene == 'IGHV2-5', ]
        IGHV2.26 <- df.popgenes[df.popgenes$Vgene == 'IGHV2-26', ]
        IGHV2.70 <- df.popgenes[df.popgenes$Vgene == 'IGHV2-70', ]
        IGHV3.7 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-7', ]
        IGHV3.9 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-9', ]
        IGHV3.11 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-11', ]
        IGHV3.13 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-13', ]
        IGHV3.15 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-15', ]
        IGHV3.20 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-20', ]
        IGHV3.21 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-21', ]
        IGHV3.22 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-22', ]
        IGHV3.23 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-23', ]
        IGHV3.30 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-30', ]
        IGHV3.30.3 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-30-3', ]
        IGHV3.33 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-33', ]
        IGHV3.43 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-43', ]
        IGHV3.43D <- df.popgenes[df.popgenes$Vgene == 'IGHV3-43D', ]
        IGHV3.48 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-48', ]
        IGHV3.49 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-49', ]
        IGHV3.53 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-53', ]
        IGHV7.4.1 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-64', ]
        IGHV3.66 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-66', ]
        IGHV3.72 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-72', ]
        IGHV3.73 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-73', ]
        IGHV3.74 <- df.popgenes[df.popgenes$Vgene == 'IGHV3-74', ]
        IGHV3.d <- df.popgenes[df.popgenes$Vgene == 'IGHV3-d', ]
        IGHV3.h <- df.popgenes[df.popgenes$Vgene == 'IGHV3-h', ]
        IGHV4.4 <- df.popgenes[df.popgenes$Vgene == 'IGHV4-4', ]
        IGHV4.28 <- df.popgenes[df.popgenes$Vgene == 'IGHV4-28', ]
        IGHV4.30.2 <- df.popgenes[df.popgenes$Vgene == 'IGHV4-30-2', ]
        IGHV4.30.4 <- df.popgenes[df.popgenes$Vgene == 'IGHV4-30-4', ]
        IGHV4.31 <- df.popgenes[df.popgenes$Vgene == 'IGHV4-31', ]
        IGHV4.34 <- df.popgenes[df.popgenes$Vgene == 'IGHV4-34', ]
        IGHV4.39 <- df.popgenes[df.popgenes$Vgene == 'IGHV4-39', ]
        IGHV4.55 <- df.popgenes[df.popgenes$Vgene == 'IGHV4-55', ]
        IGHV4.59 <- df.popgenes[df.popgenes$Vgene == 'IGHV4-59', ]
        IGHV4.61 <- df.popgenes[df.popgenes$Vgene == 'IGHV4-61', ]
        IGHV4.b <- df.popgenes[df.popgenes$Vgene == 'IGHV4-b', ]
        IGHV5.51 <- df.popgenes[df.popgenes$Vgene == 'IGHV5-51', ]
        IGHV6.1 <- df.popgenes[df.popgenes$Vgene == 'IGHV6-1', ]
        IGHJ1 <- df.popgenes[df.popgenes$JGene == 'IGHJ1', ]
        IGHJ2 <- df.popgenes[df.popgenes$JGene == 'IGHJ2', ]
        IGHJ3 <- df.popgenes[df.popgenes$JGene == 'IGHJ3', ]
        IGHJ4 <- df.popgenes[df.popgenes$JGene == 'IGHJ4', ]
        IGHJ5 <- df.popgenes[df.popgenes$JGene == 'IGHJ5', ]
        IGHJ6 <- df.popgenes[df.popgenes$JGene == 'IGHJ6', ]
        gene.info[n, 1] <- round(nrow(IGHV1.2)/nrow(df.popgenes)*100, 2)
        gene.info[n, 2] <- round(nrow(IGHV1.3)/nrow(df.popgenes)*100, 2)
        gene.info[n, 3] <- round(nrow(IGHV1.8)/nrow(df.popgenes)*100, 2)
        gene.info[n, 4] <- round(nrow(IGHV1.18)/nrow(df.popgenes)*100, 2)
        gene.info[n, 5] <- round(nrow(IGHV1.24)/nrow(df.popgenes)*100, 2)
        gene.info[n, 6] <- round(nrow(IGHV1.45)/nrow(df.popgenes)*100, 2)
        gene.info[n, 7] <- round(nrow(IGHV1.46)/nrow(df.popgenes)*100, 2)
        gene.info[n, 8] <- round(nrow(IGHV1.58)/nrow(df.popgenes)*100, 2)
        gene.info[n, 9] <- round(nrow(IGHV1.69)/nrow(df.popgenes)*100, 2)
        gene.info[n, 10] <- round(nrow(IGHV1.f)/nrow(df.popgenes)*100, 2)
        gene.info[n, 11] <- round(nrow(IGHV2.5)/nrow(df.popgenes)*100, 2)
        gene.info[n, 12] <- round(nrow(IGHV2.26)/nrow(df.popgenes)*100, 2)
        gene.info[n, 13] <- round(nrow(IGHV2.70)/nrow(df.popgenes)*100, 2)
        gene.info[n, 14] <- round(nrow(IGHV3.7)/nrow(df.popgenes)*100, 2)
        gene.info[n, 15] <- round(nrow(IGHV3.9)/nrow(df.popgenes)*100, 2)
        gene.info[n, 16] <- round(nrow(IGHV3.11)/nrow(df.popgenes)*100, 2)
        gene.info[n, 17] <- round(nrow(IGHV3.13)/nrow(df.popgenes)*100, 2)
        gene.info[n, 18] <- round(nrow(IGHV3.15)/nrow(df.popgenes)*100, 2)
        gene.info[n, 19] <- round(nrow(IGHV3.20)/nrow(df.popgenes)*100, 2)
        gene.info[n, 20] <- round(nrow(IGHV3.21)/nrow(df.popgenes)*100, 2)
        gene.info[n, 21] <- round(nrow(IGHV3.22)/nrow(df.popgenes)*100, 2)
        gene.info[n, 22] <- round(nrow(IGHV3.23)/nrow(df.popgenes)*100, 2)
        gene.info[n, 23] <- round(nrow(IGHV3.30)/nrow(df.popgenes)*100, 2)
        gene.info[n, 24] <- round(nrow(IGHV3.30.3)/nrow(df.popgenes)*100, 2)
        gene.info[n, 25] <- round(nrow(IGHV3.33)/nrow(df.popgenes)*100, 2)
        gene.info[n, 26] <- round(nrow(IGHV3.43)/nrow(df.popgenes)*100, 2)
        gene.info[n, 27] <- round(nrow(IGHV3.43D)/nrow(df.popgenes)*100, 2)
        gene.info[n, 28] <- round(nrow(IGHV3.48)/nrow(df.popgenes)*100, 2)
        gene.info[n, 29] <- round(nrow(IGHV3.49)/nrow(df.popgenes)*100, 2)
        gene.info[n, 30] <- round(nrow(IGHV3.53)/nrow(df.popgenes)*100, 2)
        gene.info[n, 31] <- round(nrow(IGHV3.64)/nrow(df.popgenes)*100, 2)
        gene.info[n, 32] <- round(nrow(IGHV3.66)/nrow(df.popgenes)*100, 2)
        gene.info[n, 33] <- round(nrow(IGHV3.72)/nrow(df.popgenes)*100, 2)
        gene.info[n, 34] <- round(nrow(IGHV3.73)/nrow(df.popgenes)*100, 2)
        gene.info[n, 35] <- round(nrow(IGHV3.74)/nrow(df.popgenes)*100, 2)
        gene.info[n, 36] <- round(nrow(IGHV3.d)/nrow(df.popgenes)*100, 2)
        gene.info[n, 37] <- round(nrow(IGHV3.h)/nrow(df.popgenes)*100, 2)
        gene.info[n, 38] <- round(nrow(IGHV4.4)/nrow(df.popgenes)*100, 2)
        gene.info[n, 39] <- round(nrow(IGHV4.28)/nrow(df.popgenes)*100, 2)
        gene.info[n, 40] <- round(nrow(IGHV4.30.2)/nrow(df.popgenes)*100, 2)
        gene.info[n, 41] <- round(nrow(IGHV4.30.4)/nrow(df.popgenes)*100, 2)
        gene.info[n, 42] <- round(nrow(IGHV4.31)/nrow(df.popgenes)*100, 2)
        gene.info[n, 43] <- round(nrow(IGHV4.34)/nrow(df.popgenes)*100, 2)
        gene.info[n, 44] <- round(nrow(IGHV4.39)/nrow(df.popgenes)*100, 2)
        gene.info[n, 45] <- round(nrow(IGHV4.55)/nrow(df.popgenes)*100, 2)
        gene.info[n, 46] <- round(nrow(IGHV4.59)/nrow(df.popgenes)*100, 2)
        gene.info[n, 47] <- round(nrow(IGHV4.61)/nrow(df.popgenes)*100, 2)
        gene.info[n, 48] <- round(nrow(IGHV4.b)/nrow(df.popgenes)*100, 2)
        gene.info[n, 49] <- round(nrow(IGHV5.51)/nrow(df.popgenes)*100, 2)
        gene.info[n, 50] <- round(nrow(IGHV6.1)/nrow(df.popgenes)*100, 2)
        gene.info[n, 51] <- round(nrow(IGHJ1)/nrow(df.popgenes)*100, 2)
        gene.info[n, 52] <- round(nrow(IGHJ2)/nrow(df.popgenes)*100, 2)
        gene.info[n, 53] <- round(nrow(IGHJ3)/nrow(df.popgenes)*100, 2)
        gene.info[n, 54] <- round(nrow(IGHJ4)/nrow(df.popgenes)*100, 2)
        gene.info[n, 55] <- round(nrow(IGHJ5)/nrow(df.popgenes)*100, 2)
        gene.info[n, 56] <- round(nrow(IGHJ6)/nrow(df.popgenes)*100, 2)
    }



colnames(gene.info) <- c("IGHV1-2", "IGHV1-3", "IGHV1-8", "IGHV1-18", "IGHV1-24",
    "IGHV1-45", "IGHV1-46", "IGHV1-58", "IGHV1-69", "IGHV1-f", "IGHV2-5", "IGHV2-26",
    "IGHV2-70", "IGHV3-7", "IGHV3-9", "IGHV3-11", "IGHV3-13", "IGHV3-15", "IGHV3-20",
    "IGHV3-21", "IGHV3-22", "IGHV3-23", "IGHV3-30", "IGHV3-30-3", "IGHV3-33", "IGHV3-43",
    "IGHV3-43D", "IGHV3-48", "IGHV3-49", "IGHV3-53", "IGHV3-64", "IGHV3-66", "IGHV3-72",
    "IGHV3-73", "IGHV3-74", "IGHV3-d", "IGHV3-h", "IGHV4-4", "IGHV4-28", "IGHV4-30-2",
    "IGHV4-30-4", "IGHV4-31", "IGHV4-34", "IGHV4-39", "IGHV4-55", "IGHV4-59", "IGHV4-61",
    "IGHV4-b", "IGHV5-51", "IGHV6-1", "IGHJ1", "IGHJ2", "IGHJ3", "IGHJ4", "IGHJ5",
    "IGHJ6")



master.data <- cbind(master.data, gene.info)


IGHV1.2 <- full.data[full.data$Vgene == "IGHV1-2", ]
IGHV1.2$Last3FR1 <- substrRight(IGHV1.2$AA_FR1, 3)
IGHV1.2.table <- table(IGHV1.2$population, IGHV1.2$Last3FR1)
master.data$IGHV1.2_Last3FR1 <- IGHV1.2.table[, "KAS"]/rowSums(IGHV1.2.table) * 100

IGHV1.3 <- full.data[full.data$Vgene == "IGHV1-3", ]
IGHV1.3$Last3FR1 <- substrRight(IGHV1.3$AA_FR1, 3)
IGHV1.3.table <- table(IGHV1.3$population, IGHV1.3$Last3FR1)
master.data$IGHV1.3_Last3FR1 <- IGHV1.3.table[, "KAS"]/rowSums(IGHV1.3.table) * 100

IGHV1.8 <- full.data[full.data$Vgene == "IGHV1-8", ]
IGHV1.8$Last3FR1 <- substrRight(IGHV1.8$AA_FR1, 3)
IGHV1.8.table <- table(IGHV1.8$population, IGHV1.8$Last3FR1)
master.data$IGHV1.8_Last3FR1 <- IGHV1.8.table[, "KAS"]/rowSums(IGHV1.8.table) * 100

IGHV1.18 <- full.data[full.data$Vgene == "IGHV1-18", ]
IGHV1.18$Last3FR1 <- substrRight(IGHV1.18$AA_FR1, 3)
IGHV1.18.table <- table(IGHV1.18$population, IGHV1.18$Last3FR1)
master.data$IGHV1.18_Last3FR1 <- IGHV1.18.table[, "KAS"]/rowSums(IGHV1.18.table) * 100

IGHV1.24 <- full.data[full.data$Vgene == "IGHV1-24", ]
IGHV1.24$Last3FR1 <- substrRight(IGHV1.24$AA_FR1, 3)
IGHV1.24.table <- table(IGHV1.24$population, IGHV1.24$Last3FR1)
master.data$IGHV1.24_Last3FR1 <- IGHV1.24.table[, "KVS"]/rowSums(IGHV1.24.table) * 100

IGHV1.45 <- full.data[full.data$Vgene == "IGHV1-45", ]
IGHV1.45$Last3FR1 <- substrRight(IGHV1.45$AA_FR1, 3)
IGHV1.45.table <- table(IGHV1.45$population, IGHV1.45$Last3FR1)
master.data$IGHV1.45_Last3FR1 <- IGHV1.45.table[, "KAS"]/rowSums(IGHV1.45.table) * 100

IGHV1.46 <- full.data[full.data$Vgene == "IGHV1-46", ]
IGHV1.46$Last3FR1 <- substrRight(IGHV1.46$AA_FR1, 3)
IGHV1.46.table <- table(IGHV1.46$population, IGHV1.46$Last3FR1)
master.data$IGHV1.46_Last3FR1 <- IGHV1.46.table[, "KAS"]/rowSums(IGHV1.46.table) * 100

IGHV1.58 <- full.data[full.data$Vgene == "IGHV1-58", ]
IGHV1.58$Last3FR1 <- substrRight(IGHV1.58$AA_FR1, 3)
IGHV1.58.table <- table(IGHV1.58$population, IGHV1.58$Last3FR1)
master.data$IGHV1.58_Last3FR1 <- IGHV1.58.table[, "KAS"]/rowSums(IGHV1.58.table) * 100

IGHV1.69 <- full.data[full.data$Vgene == "IGHV1-69", ]
IGHV1.69$Last3FR1 <- substrRight(IGHV1.69$AA_FR1, 3)
IGHV1.69.table <- table(IGHV1.69$population, IGHV1.69$Last3FR1)
master.data$IGHV1.69_Last3FR1 <- IGHV1.69.table[, "KAS"]/rowSums(IGHV1.69.table) * 100

IGHV2.5 <- full.data[full.data$Vgene == "IGHV2-5", ]
IGHV2.5$Last3FR1 <- substrRight(IGHV2.5$AA_FR1, 3)
IGHV2.5.table <- table(IGHV2.5$population, IGHV2.5$Last3FR1)
master.data$IGHV2.5_Last3FR1 <- IGHV2.5.table[, "TFS"]/rowSums(IGHV2.5.table) * 100

IGHV2.26 <- full.data[full.data$Vgene == "IGHV2-26", ]
IGHV2.26$Last3FR1 <- substrRight(IGHV2.26$AA_FR1, 3)
IGHV2.26.table <- table(IGHV2.26$population, IGHV2.26$Last3FR1)
master.data$IGHV2.26_Last3FR1 <- IGHV2.26.table[, "TVS"]/rowSums(IGHV2.26.table) * 100

IGHV2.70 <- full.data[full.data$Vgene == "IGHV2-70", ]
IGHV2.70$Last3FR1 <- substrRight(IGHV2.70$AA_FR1, 3)
IGHV2.70.table <- table(IGHV2.70$population, IGHV2.70$Last3FR1)
master.data$IGHV2.70_Last3FR1 <- IGHV2.70.table[, "TFS"]/rowSums(IGHV2.70.table) * 100

IGHV3.7 <- full.data[full.data$Vgene == "IGHV3-7", ]
IGHV3.7$Last3FR1 <- substrRight(IGHV3.7$AA_FR1, 3)
IGHV3.7.table <- table(IGHV3.7$population, IGHV3.7$Last3FR1)
master.data$IGHV3.7_Last3FR1 <- IGHV3.7.table[, "AAS"]/rowSums(IGHV3.7.table) * 100

IGHV3.9 <- full.data[full.data$Vgene == "IGHV3-9", ]
IGHV3.9$Last3FR1 <- substrRight(IGHV3.9$AA_FR1, 3)
IGHV3.9.table <- table(IGHV3.9$population, IGHV3.9$Last3FR1)
master.data$IGHV3.9_Last3FR1 <- IGHV3.9.table[, "AAS"]/rowSums(IGHV3.9.table) * 100

IGHV3.11 <- full.data[full.data$Vgene == "IGHV3-11", ]
IGHV3.11$Last3FR1 <- substrRight(IGHV3.11$AA_FR1, 3)
IGHV3.11.table <- table(IGHV3.11$population, IGHV3.11$Last3FR1)
master.data$IGHV3.11_Last3FR1 <- IGHV3.11.table[, "AAS"]/rowSums(IGHV3.11.table) * 100

IGHV3.13 <- full.data[full.data$Vgene == "IGHV3-13", ]
IGHV3.13$Last3FR1 <- substrRight(IGHV3.13$AA_FR1, 3)
IGHV3.13.table <- table(IGHV3.13$population, IGHV3.13$Last3FR1)
master.data$IGHV3.13_Last3FR1 <- IGHV3.13.table[, "AAS"]/rowSums(IGHV3.13.table) * 100

IGHV3.15 <- full.data[full.data$Vgene == "IGHV3-15", ]
IGHV3.15$Last3FR1 <- substrRight(IGHV3.15$AA_FR1, 3)
IGHV3.15.table <- table(IGHV3.15$population, IGHV3.15$Last3FR1)
master.data$IGHV3.15_Last3FR1 <- IGHV3.15.table[, "AAS"]/rowSums(IGHV3.15.table) * 100

IGHV3.20 <- full.data[full.data$Vgene == "IGHV3-20", ]
IGHV3.20$Last3FR1 <- substrRight(IGHV3.20$AA_FR1, 3)
IGHV3.20.table <- table(IGHV3.20$population, IGHV3.20$Last3FR1)
master.data$IGHV3.20_Last3FR1 <- IGHV3.20.table[, "AAS"]/rowSums(IGHV3.20.table) * 100

IGHV3.21 <- full.data[full.data$Vgene == "IGHV3-21", ]
IGHV3.21$Last3FR1 <- substrRight(IGHV3.21$AA_FR1, 3)
IGHV3.21.table <- table(IGHV3.21$population, IGHV3.21$Last3FR1)
master.data$IGHV3.21_Last3FR1 <- IGHV3.21.table[, "AAS"]/rowSums(IGHV3.21.table) * 100

IGHV3.23 <- full.data[full.data$Vgene == "IGHV3-23", ]
IGHV3.23$Last3FR1 <- substrRight(IGHV3.23$AA_FR1, 3)
IGHV3.23.table <- table(IGHV3.23$population, IGHV3.23$Last3FR1)
master.data$IGHV3.23_Last3FR1 <- IGHV3.23.table[, "AAS"]/rowSums(IGHV3.23.table) * 100

IGHV3.23D <- full.data[full.data$Vgene == "IGHV3-23D", ]
IGHV3.23D$Last3FR1 <- substrRight(IGHV3.23D$AA_FR1, 3)
IGHV3.23D.table <- table(IGHV3.23D$population, IGHV3.23D$Last3FR1)
master.data$IGHV3.23D_Last3FR1 <- IGHV3.23D.table[, "AAS"]/rowSums(IGHV3.23D.table) * 100

IGHV3.30 <- full.data[full.data$Vgene == "IGHV3-30", ]
IGHV3.30$Last3FR1 <- substrRight(IGHV3.30$AA_FR1, 3)
IGHV3.30.table <- table(IGHV3.30$population, IGHV3.30$Last3FR1)
master.data$IGHV3.30_Last3FR1 <- IGHV3.30.table[, "AAS"]/rowSums(IGHV3.30.table) * 100

IGHV3.30.3 <- full.data[full.data$Vgene == "IGHV3-30-3", ]
IGHV3.30.3$Last3FR1 <- substrRight(IGHV3.30.3$AA_FR1, 3)
IGHV3.30.3.table <- table(IGHV3.30.3$population, IGHV3.30.3$Last3FR1)
master.data$IGHV3.30.3_Last3FR1 <- IGHV3.30.3.table[, "AAS"]/rowSums(IGHV3.30.3.table) * 100

IGHV3.33 <- full.data[full.data$Vgene == "IGHV3-33", ]
IGHV3.33$Last3FR1 <- substrRight(IGHV3.33$AA_FR1, 3)
IGHV3.33.table <- table(IGHV3.33$population, IGHV3.33$Last3FR1)
master.data$IGHV3.33_Last3FR1 <- IGHV3.33.table[, "AAS"]/rowSums(IGHV3.33.table) * 100

IGHV3.43 <- full.data[full.data$Vgene == "IGHV3-43", ]
IGHV3.43$Last3FR1 <- substrRight(IGHV3.43$AA_FR1, 3)
IGHV3.43.table <- table(IGHV3.43$population, IGHV3.43$Last3FR1)
master.data$IGHV3.43_Last3FR1 <- IGHV3.43.table[, "AAS"]/rowSums(IGHV3.43.table) * 100

IGHV3.43D <- full.data[full.data$Vgene == "IGHV3-43D", ]
IGHV3.43D$Last3FR1 <- substrRight(IGHV3.43D$AA_FR1, 3)
IGHV3.43D.table <- table(IGHV3.43D$population, IGHV3.43D$Last3FR1)
master.data$IGHV3.43D_Last3FR1 <- IGHV3.43D.table[, "AAS"]/rowSums(IGHV3.43D.table) * 100

IGHV3.48 <- full.data[full.data$Vgene == "IGHV3-48", ]
IGHV3.48$Last3FR1 <- substrRight(IGHV3.48$AA_FR1, 3)
IGHV3.48.table <- table(IGHV3.48$population, IGHV3.48$Last3FR1)
master.data$IGHV3.48_Last3FR1 <- IGHV3.48.table[, "AAS"]/rowSums(IGHV3.48.table) * 100

IGHV3.49 <- full.data[full.data$Vgene == "IGHV3-49", ]
IGHV3.49$Last3FR1 <- substrRight(IGHV3.49$AA_FR1, 3)
IGHV3.49.table <- table(IGHV3.49$population, IGHV3.49$Last3FR1)
master.data$IGHV3.49_Last3FR1 <- IGHV3.49.table[, "TAS"]/rowSums(IGHV3.49.table) * 100

IGHV3.53 <- full.data[full.data$Vgene == "IGHV3-53", ]
IGHV3.53$Last3FR1 <- substrRight(IGHV3.53$AA_FR1, 3)
IGHV3.53.table <- table(IGHV3.53$population, IGHV3.53$Last3FR1)
master.data$IGHV3.53_Last3FR1 <- IGHV3.53.table[, "AAS"]/rowSums(IGHV3.53.table) * 100

IGHV3.64 <- full.data[full.data$Vgene == "IGHV3-64", ]
IGHV3.64$Last3FR1 <- substrRight(IGHV3.64$AA_FR1, 3)
IGHV3.64.table <- table(IGHV3.64$population, IGHV3.64$Last3FR1)
master.data$IGHV3.64_Last3FR1 <- IGHV3.64.table[, "AAS"]/rowSums(IGHV3.64.table) * 100

IGHV3.64D <- full.data[full.data$Vgene == "IGHV3-64D", ]
IGHV3.64D$Last3FR1 <- substrRight(IGHV3.64D$AA_FR1, 3)
IGHV3.64D.table <- table(IGHV3.64D$population, IGHV3.64D$Last3FR1)
master.data$IGHV3.64D_Last3FR1 <- IGHV3.64D.table[, "SAS"]/rowSums(IGHV3.64D.table) * 100

IGHV3.66 <- full.data[full.data$Vgene == "IGHV3-66", ]
IGHV3.66$Last3FR1 <- substrRight(IGHV3.66$AA_FR1, 3)
IGHV3.66.table <- table(IGHV3.66$population, IGHV3.66$Last3FR1)
master.data$IGHV3.66_Last3FR1 <- IGHV3.66.table[, "AAS"]/rowSums(IGHV3.66.table) * 100

IGHV3.72 <- full.data[full.data$Vgene == "IGHV3-72", ]
IGHV3.72$Last3FR1 <- substrRight(IGHV3.72$AA_FR1, 3)
IGHV3.72.table <- table(IGHV3.72$population, IGHV3.72$Last3FR1)
master.data$IGHV3.72_Last3FR1 <- IGHV3.72.table[, "AAS"]/rowSums(IGHV3.72.table) * 100

IGHV3.73 <- full.data[full.data$Vgene == "IGHV3-73", ]
IGHV3.73$Last3FR1 <- substrRight(IGHV3.73$AA_FR1, 3)
IGHV3.73.table <- table(IGHV3.73$population, IGHV3.73$Last3FR1)
master.data$IGHV3.73_Last3FR1 <- IGHV3.73.table[, "AAS"]/rowSums(IGHV3.73.table) * 100

IGHV3.74 <- full.data[full.data$Vgene == "IGHV3-74", ]
IGHV3.74$Last3FR1 <- substrRight(IGHV3.74$AA_FR1, 3)
IGHV3.74.table <- table(IGHV3.74$population, IGHV3.74$Last3FR1)
master.data$IGHV3.74_Last3FR1 <- IGHV3.74.table[, "AAS"]/rowSums(IGHV3.74.table) * 100

IGHV4.4 <- full.data[full.data$Vgene == "IGHV4-4", ]
IGHV4.4$Last3FR1 <- substrRight(IGHV4.4$AA_FR1, 3)
IGHV4.4.table <- table(IGHV4.4$population, IGHV4.4$Last3FR1)
master.data$IGHV4.4_Last3FR1 <- IGHV4.4.table[, "AVS"]/rowSums(IGHV4.4.table) * 100

IGHV4.28 <- full.data[full.data$Vgene == "IGHV4-28", ]
IGHV4.28$Last3FR1 <- substrRight(IGHV4.28$AA_FR1, 3)
IGHV4.28.table <- table(IGHV4.28$population, IGHV4.28$Last3FR1)
master.data$IGHV4.28_Last3FR1 <- IGHV4.28.table[, "AVS"]/rowSums(IGHV4.28.table) * 100

IGHV4.30.2 <- full.data[full.data$Vgene == "IGHV4-30-2", ]
IGHV4.30.2$Last3FR1 <- substrRight(IGHV4.30.2$AA_FR1, 3)
IGHV4.30.2.table <- table(IGHV4.30.2$population, IGHV4.30.2$Last3FR1)
master.data$IGHV4.30.2_Last3FR1 <- IGHV4.30.2.table[, "AVS"]/rowSums(IGHV4.30.2.table) * 100

IGHV4.30.4 <- full.data[full.data$Vgene == "IGHV4-30-4", ]
IGHV4.30.4$Last3FR1 <- substrRight(IGHV4.30.4$AA_FR1, 3)
IGHV4.30.4.table <- table(IGHV4.30.4$population, IGHV4.30.4$Last3FR1)
master.data$IGHV4.30.4_Last3FR1 <- IGHV4.30.4.table[, "TVS"]/rowSums(IGHV4.30.4.table) * 100

IGHV4.31 <- full.data[full.data$Vgene == "IGHV4-31", ]
IGHV4.31$Last3FR1 <- substrRight(IGHV4.31$AA_FR1, 3)
IGHV4.31.table <- table(IGHV4.31$population, IGHV4.31$Last3FR1)
master.data$IGHV4.31_Last3FR1 <- IGHV4.31.table[, "TVS"]/rowSums(IGHV4.31.table) * 100

IGHV4.34 <- full.data[full.data$Vgene == "IGHV4-34", ]
IGHV4.34$Last3FR1 <- substrRight(IGHV4.34$AA_FR1, 3)
IGHV4.34.table <- table(IGHV4.34$population, IGHV4.34$Last3FR1)
master.data$IGHV4.34_Last3FR1 <- IGHV4.34.table[, "AVY"]/rowSums(IGHV4.34.table) * 100

IGHV4.38 <- full.data[full.data$Vgene == "IGHV4-38", ]
IGHV4.38$Last3FR1 <- substrRight(IGHV4.38$AA_FR1, 3)
IGHV4.38.table <- table(IGHV4.38$population, IGHV4.38$Last3FR1)
master.data$IGHV4.38_Last3FR1 <- IGHV4.38.table[, "TVS"]/rowSums(IGHV4.38.table) * 100

IGHV4.39 <- full.data[full.data$Vgene == "IGHV4-39", ]
IGHV4.39$Last3FR1 <- substrRight(IGHV4.39$AA_FR1, 3)
IGHV4.39.table <- table(IGHV4.39$population, IGHV4.39$Last3FR1)
master.data$IGHV4.39_Last3FR1 <- IGHV4.39.table[, "TVS"]/rowSums(IGHV4.39.table) * 100

IGHV4.59 <- full.data[full.data$Vgene == "IGHV4-59", ]
IGHV4.59$Last3FR1 <- substrRight(IGHV4.59$AA_FR1, 3)
IGHV4.59.table <- table(IGHV4.59$population, IGHV4.59$Last3FR1)
master.data$IGHV4.59_Last3FR1 <- IGHV4.59.table[, "TVS"]/rowSums(IGHV4.59.table) * 100

IGHV4.61 <- full.data[full.data$Vgene == "IGHV4-61", ]
IGHV4.61$Last3FR1 <- substrRight(IGHV4.61$AA_FR1, 3)
IGHV4.61.table <- table(IGHV4.61$population, IGHV4.61$Last3FR1)
master.data$IGHV4.61_Last3FR1 <- IGHV4.61.table[, "TVS"]/rowSums(IGHV4.61.table) * 100

IGHV5.10 <- full.data[full.data$Vgene == "IGHV5-10", ]
IGHV5.10$Last3FR1 <- substrRight(IGHV5.10$AA_FR1, 3)
IGHV5.10.table <- table(IGHV5.10$population, IGHV5.10$Last3FR1)
master.data$IGHV5.10_Last3FR1 <- IGHV5.10.table[, "KGS"]/rowSums(IGHV5.10.table) * 100

IGHV5.51 <- full.data[full.data$Vgene == "IGHV5-51", ]
IGHV5.51$Last3FR1 <- substrRight(IGHV5.51$AA_FR1, 3)
IGHV5.51.table <- table(IGHV5.51$population, IGHV5.51$Last3FR1)
master.data$IGHV5.51_Last3FR1 <- IGHV5.51.table[, "KGS"]/rowSums(IGHV5.51.table) * 100

IGHV6.1 <- full.data[full.data$Vgene == "IGHV6-1", ]
IGHV6.1$Last3FR1 <- substrRight(IGHV6.1$AA_FR1, 3)
IGHV6.1.table <- table(IGHV6.1$population, IGHV6.1$Last3FR1)
master.data$IGHV6.1_Last3FR1 <- IGHV6.1.table[, "AIS"]/rowSums(IGHV6.1.table) * 100

IGHV7.4.1 <- full.data[full.data$Vgene == "IGHV7-4-1", ]
IGHV7.4.1$Last3FR1 <- substrRight(IGHV7.4.1$AA_FR1, 3)
IGHV7.4.1.table <- table(IGHV7.4.1$population, IGHV7.4.1$Last3FR1)
master.data$IGHV7.4.1_Last3FR1 <- IGHV7.4.1.table[, "KAS"]/rowSums(IGHV7.4.1.table) * 100


write.table(master.data, file = master.out, quote = FALSE, sep = "\t", row.names = FALSE)

}







##############################






    # make a data frame with percent of matching clones in each

    percof.df <- data.frame()

    for (m in 1:numpops) {
        df.x <- filter(lineages.union, lineages.union[[m + 4]] > 0)
        lin.x <- nrow(df.x)
        for (n in 1:numpops) {
            df.x_y <- filter(df.x, df.x[[n + 4]] > 0)
            lin.x_y <- nrow(df.x_y)
            perc.x_y <- lin.x_y/lin.x * 100
            percof.full.data[m, n] <- round(perc.x_y, 2)
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

        full.data$Lineage <- factor(full.data$Lineage)

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
