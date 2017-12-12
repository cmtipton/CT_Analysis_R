
# location input requires file location of the grouped lineage file from IgSeq
seqanalysis <- function(location, remove.identicals = 'n') {
    #Load dependencies
    library(RColorBrewer)
    library(ggplot2)
    library(dplyr)
    library(plyr)
    library(data.table)
    library(reshape2)
    library(spaa)
    library(RSQLite)
    library(stringr)
    library(doBy)
    library(vegan)
    library(ggthemes)


    # function to load SQLite db file
    loadsql <- function(filename) {
        sqlite.driver <- dbDriver("SQLite")
        db <- dbConnect(sqlite.driver, dbname = filename)
        full.data <<- dbReadTable(db, "sequence")
    }

    #establish input and output file locations

    #location <- '/Users/cmtipto/Desktop/Test/999'
    location2 <- location
    file.out <- paste0(location, "/CT_Sequencing_Analysis/")
    master.out <- paste0(file.out, "MasterData_", sub(".*/", "", location2), ".txt")
    dir.create(file.path(location, "/CT_Sequencing_Analysis"), showWarnings = FALSE)
    subject <- sub(".*/", "", location2)
    dblocation <- list.files(location, pattern = "\\.db$")
    dbfile <- paste0(location, "/", dblocation[1])
    loading <- paste0("Loading database file for subject ", subject, "...")
    print(loading)
    loadsql(dbfile)


    # Label singletons with boolean column
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) <
        tol
    full.data$Singleton <- !is.wholenumber(full.data$lineageID)

    # To get rid of identical sequences
    if (remove.identicals == 'y'){
        full.data$FR2on <- paste0(full.data$NT_CDR1, full.data$NT_FR2, full.data$NT_CDR2, full.data$NT_FR3, full.data$NT_CDR3, full.data$NT_J)
        full.data <- full.data[!duplicated(full.data$FR2on), ]
        }


    # Separate and generate population/sample columns
    full.data$Sample <- full.data$population

    if (str_count(full.data$population[1], "\\_") == 4) {
        sample.info <- str_split_fixed(full.data$population, "_", 5)
        full.data$population <- paste0(sample.info[, 3], "_", sample.info[, 4])
    }

    full.data$General_Population <- full.data$population


    # determine number of populations
    poplist <- levels(as.factor(full.data$population))
    numpops <- length(poplist)
    samplelist <- levels(as.factor(full.data$Sample))
    numsamples <- length(samplelist)
    if (numpops != numsamples) {
        full.data$population <- full.data$Sample
    }

    print(paste0("Analyzing ", numsamples, " populations of cells:"))
    poplist

    samplelist <- levels(as.factor(full.data$Sample))


    non.single <- full.data[full.data$Singleton == FALSE, ]

    lineage.data <- non.single[, c("population", "lineageID", "Vgene")]
    table.pop <- table(lineage.data$lineageID, lineage.data$population)
    table.vgene <- table(lineage.data$lineageID, lineage.data$Vgene)
    pop.df <- cbind(rownames(table.pop), as.data.frame.matrix(table.pop, row.names = NULL))
    colnames(pop.df)[1] <- "LineageID"
    v.df <- cbind(rownames(table.vgene), as.data.frame.matrix(table.vgene, row.names = NULL))
    colnames(v.df)[1] <- "LineageID"
    table.vgene.melt <- melt(v.df, key = LineageID)
    v.melt <- table.vgene.melt[table.vgene.melt$value != 0, ]
    lineages.union <- cbind(v.melt$variable, pop.df)
    colnames(lineages.union)[1] <- "Vgene"

    melt.x <- melt(lineages.union)
    filtered.melt.x <- melt.x[melt.x$value != 0, ]
    lineage.table <- table(filtered.melt.x$variable)
    melt.lineage.table <- as.data.frame(melt(lineage.table))
    colnames(melt.lineage.table) <- c("Population", "Lineages")

    vtable <- table(full.data$population, full.data$Vgene)
    vperc <- prop.table(vtable) * 100
    jtable <- table(full.data$population, full.data$JGene)
    jperc <- prop.table(jtable) * 100

    seqs <- as.data.frame(table(non.single$population))
    seqs <- cbind(Sample = samplelist, seqs)
    master.data <- cbind(seqs, melt.lineage.table$Lineages)
    colnames(master.data) <- c("Sample", "Population", "Sequences (Non-Single)",
        "Lineages")

    print(master.data)


    # establish function to calculate clonality based on the Shannon diversity index
    # (1 - Pielou's Evenness)
    clon <- function(x, col = 1) {
        div <- diversity(x[, col], index = "shannon")
        H <- div/log(nrow(x))
        cl <<- 1 - H
        return(round(cl, 2))
    }

    # Add clonality values for each population into a new column
    lin.df <- lineages.union[, 3:ncol(lineages.union)]
    clonality.df <- data.frame()
    for (k in 1:numsamples) {
        clonality.df[k, 1] <- clon(lin.df, col = k)
    }
    colnames(clonality.df)[1] <- "Clonality"

    master.data <- cbind(master.data, clonality.df)


    full.data$Mutation_Rate <- full.data$V_Mutations/full.data$V_Nucleotides * 100
    master.data$Mutation_Rate_Mean <- ddply(full.data, .(Sample), summarize, mean = mean(Mutation_Rate))[,
        2]
    master.data$Mutation_Rate_Median <- ddply(full.data, .(Sample), summarize, median = median(Mutation_Rate))[,
        2]


    iso.info <- data.frame()
    for (n in 1:numsamples) {
        Ig <- full.data[full.data$isotype != "U", ]
        IgG <- full.data[full.data$isotype == "G", ]
        IgA <- full.data[full.data$isotype == "A", ]
        IgM <- full.data[full.data$isotype == "M", ]
        IgE <- full.data[full.data$isotype == "E", ]
        IgU <- full.data[full.data$isotype == "U", ]
        Ig <- Ig[Ig$Sample == samplelist[n], ]
        IgG <- IgG[IgG$Sample == samplelist[n], ]
        IgA <- IgA[IgA$Sample == samplelist[n], ]
        IgM <- IgM[IgM$Sample == samplelist[n], ]
        IgE <- IgE[IgE$Sample == samplelist[n], ]
        IgU <- IgU[IgU$Sample == samplelist[n], ]
        switched <- Ig[Ig$isotype == "G" | Ig$isotype == "A" | Ig$isotype == "E",
            ]
        Ig.less3 <- Ig[Ig$Mutation_Rate < 3, ]
        IgM.less3 <- IgM[IgM$Mutation_Rate < 3, ]
        IgG.less3 <- IgG[IgG$Mutation_Rate < 3, ]
        IgA.less3 <- IgA[IgA$Mutation_Rate < 3, ]
        IgE.less3 <- IgE[IgE$Mutation_Rate < 3, ]
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
        iso.info[n, 31] <- round(mean(switched$V_Mutations), 2)
        iso.info[n, 32] <- round(median(switched$V_Mutations), 2)
        iso.info[n, 33] <- round(mean(switched$Mutation_Rate), 2)
        iso.info[n, 34] <- round(median(switched$Mutation_Rate), 2)
        iso.info[n, 35] <- round((nrow(Ig.less3)/nrow(Ig) * 100), 2)
        iso.info[n, 36] <- round((nrow(IgM.less3)/nrow(Ig) * 100), 2)
        iso.info[n, 37] <- round((nrow(IgG.less3)/nrow(Ig) * 100), 2)
        iso.info[n, 38] <- round((nrow(IgA.less3)/nrow(Ig) * 100), 2)
        iso.info[n, 39] <- round((nrow(IgE.less3)/nrow(Ig) * 100), 2)
    }

    colnames(iso.info) <- c("IgM Sequences", "IgG Sequences", "IgA Sequences", "IgE Sequences",
        "Unknown Isotype Sequences", "Percent IgM", "Percent IgG", "Percent IgA",
        "Percent IgE", "Percent Unknown Isotype", "IgM Mutations (Mean)", "IgG Mutations (Mean)",
        "IgA Mutations (Mean)", "IgE Mutations (Mean)", "Unknown Isotype Mutations (Mean)",
        "IgM Mutations (Median)", "IgG Mutations (Median)", "IgA Mutations (Median)",
        "IgE Mutations (Median)", "Unknown Isotype Mutations (Median)", "IgM Mutation Rate (Mean)",
        "IgG Mutation Rate (Mean)", "IgA Mutation Rate (Mean)", "IgE Mutation Rate (Mean)",
        "Unknown Isotype Mutation Rate (Mean)", "IgM Mutation Rate (Median)", "IgG Mutation Rate (Median)",
        "IgA Mutation Rate (Median)", "IgE Mutation Rate (Median)", "Unknown Isotype Mutation Rate (Median)",
        "Switched Mutations (Mean)", "Switched Mutations (Median)", "Switched Mutation Rate (Mean)",
        "Switched Mutation Rate (Median)", "Ig Less than 3", "IgM Less than 3", "IgG Less than 3",
        "IgA Less than 3", "IgE Less than 3")

    master.data <- cbind(master.data, iso.info)


    non.single$Mutation_Rate <- non.single$V_Mutations/non.single$V_Nucleotides *
        100
    master.data$Mutation_Rate_Mean_NS <- ddply(non.single, .(Sample), summarize,
        mean = mean(Mutation_Rate))[, 2]
    master.data$Mutation_Rate_Median_NS <- ddply(non.single, .(Sample), summarize,
        median = median(Mutation_Rate))[, 2]


    nonsingle.iso.info <- data.frame()
    for (n in 1:numsamples) {
        Ig <- non.single[non.single$isotype != "U", ]
        IgG <- non.single[non.single$isotype == "G", ]
        IgA <- non.single[non.single$isotype == "A", ]
        IgM <- non.single[non.single$isotype == "M", ]
        IgE <- non.single[non.single$isotype == "E", ]
        IgU <- non.single[non.single$isotype == "U", ]
        Ig <- Ig[Ig$Sample == samplelist[n], ]
        IgG <- IgG[IgG$Sample == samplelist[n], ]
        IgA <- IgA[IgA$Sample == samplelist[n], ]
        IgM <- IgM[IgM$Sample == samplelist[n], ]
        IgE <- IgE[IgE$Sample == samplelist[n], ]
        IgU <- IgU[IgU$Sample == samplelist[n], ]
        switched <- Ig[Ig$isotype == "G" | Ig$isotype == "A" | Ig$isotype == "E",
            ]
        Ig.less3 <- Ig[Ig$Mutation_Rate < 3, ]
        IgM.less3 <- IgM[IgM$Mutation_Rate < 3, ]
        IgG.less3 <- IgG[IgG$Mutation_Rate < 3, ]
        IgA.less3 <- IgA[IgA$Mutation_Rate < 3, ]
        IgE.less3 <- IgE[IgE$Mutation_Rate < 3, ]
        nonsingle.iso.info[n, 1] <- nrow(IgM)
        nonsingle.iso.info[n, 2] <- nrow(IgG)
        nonsingle.iso.info[n, 3] <- nrow(IgA)
        nonsingle.iso.info[n, 4] <- nrow(IgE)
        nonsingle.iso.info[n, 5] <- nrow(IgU)
        nonsingle.iso.info[n, 6] <- round((nrow(IgM)/nrow(Ig) * 100), 2)
        nonsingle.iso.info[n, 7] <- round((nrow(IgG)/nrow(Ig) * 100), 2)
        nonsingle.iso.info[n, 8] <- round((nrow(IgA)/nrow(Ig) * 100), 2)
        nonsingle.iso.info[n, 9] <- round((nrow(IgE)/nrow(Ig) * 100), 2)
        nonsingle.iso.info[n, 10] <- round((nrow(IgU)/nrow(full.data) * 100), 2)
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
        nonsingle.iso.info[n, 31] <- round(mean(switched$V_Mutations), 2)
        nonsingle.iso.info[n, 32] <- round(median(switched$V_Mutations), 2)
        nonsingle.iso.info[n, 33] <- round(mean(switched$Mutation_Rate), 2)
        nonsingle.iso.info[n, 34] <- round(median(switched$Mutation_Rate), 2)
        nonsingle.iso.info[n, 35] <- round((nrow(Ig.less3)/nrow(Ig) * 100), 2)
        nonsingle.iso.info[n, 36] <- round((nrow(IgM.less3)/nrow(Ig) * 100), 2)
        nonsingle.iso.info[n, 37] <- round((nrow(IgG.less3)/nrow(Ig) * 100), 2)
        nonsingle.iso.info[n, 38] <- round((nrow(IgA.less3)/nrow(Ig) * 100), 2)
        nonsingle.iso.info[n, 39] <- round((nrow(IgE.less3)/nrow(Ig) * 100), 2)
    }

    colnames(nonsingle.iso.info) <- c("IgM Sequences_NS", "IgG Sequences_NS", "IgA Sequences_NS",
        "IgE Sequences_NS", "Unknown Isotype Sequences_NS", "Percent IgM_NS", "Percent IgG_NS",
        "Percent IgA_NS", "Percent IgE_NS", "Percent Unknown Isotype_NS", "IgM Mutations (Mean)_NS",
        "IgG Mutations (Mean)_NS", "IgA Mutations (Mean)_NS", "IgE Mutations (Mean)_NS",
        "Unknown Isotype Mutations (Mean)_NS", "IgM Mutations (Median)_NS", "IgG Mutations (Median)_NS",
        "IgA Mutations (Median)_NS", "IgE Mutations (Median)_NS", "Unknown Isotype Mutations (Median)_NS",
        "IgM Mutation Rate (Mean)_NS", "IgG Mutation Rate (Mean)_NS", "IgA Mutation Rate (Mean)_NS",
        "IgE Mutation Rate (Mean)_NS", "Unknown Isotype Mutation Rate (Mean)_NS",
        "IgM Mutation Rate (Median)_NS", "IgG Mutation Rate (Median)_NS", "IgA Mutation Rate (Median)_NS",
        "IgE Mutation Rate (Median)_NS", "Unknown Isotype Mutation Rate (Median)_NS",
        "Switched Mutations (Mean)_NS", "Switched Mutations (Median)_NS", "Switched Mutation Rate (Mean)_NS",
        "Switched Mutation Rate (Median)_NS", "Ig Less than 3_NS", "IgM Less than 3_NS",
        "IgG Less than 3_NS", "IgA Less than 3_NS", "IgE Less than 3_NS")

    master.data <- cbind(master.data, nonsingle.iso.info)



    full.data$V_SilentMutation_Rate <- full.data$V_SilentMutations/full.data$V_Nucleotides *
        100
    full.data$V_NonsilentMutation_Rate <- full.data$V_NonsilentMutations/full.data$V_Nucleotides *
        100
    full.data$FR1_Mutation_Rate <- full.data$FR1_Mutations/full.data$FR1_Nucleotides *
        100
    full.data$FR1_SilentMutation_Rate <- full.data$FR1_SilentMutations/full.data$FR1_Nucleotides *
        100
    full.data$FR1_NonsilentMutation_Rate <- full.data$FR1_NonsilentMutations/full.data$FR1_Nucleotides *
        100
    full.data$CDR1_Mutation_Rate <- full.data$CDR1_Mutations/full.data$CDR1_Nucleotides *
        100
    full.data$CDR1_SilentMutation_Rate <- full.data$CDR1_SilentMutations/full.data$CDR1_Nucleotides *
        100
    full.data$CDR1_NonsilentMutation_Rate <- full.data$CDR1_NonsilentMutations/full.data$CDR1_Nucleotides *
        100
    full.data$FR2_Mutation_Rate <- full.data$FR2_Mutations/full.data$FR2_Nucleotides *
        100
    full.data$FR2_SilentMutation_Rate <- full.data$FR2_SilentMutations/full.data$FR2_Nucleotides *
        100
    full.data$FR2_NonsilentMutation_Rate <- full.data$FR2_NonsilentMutations/full.data$FR2_Nucleotides *
        100
    full.data$CDR2_Mutation_Rate <- full.data$CDR2_Mutations/full.data$CDR2_Nucleotides *
        100
    full.data$CDR2_SilentMutation_Rate <- full.data$CDR2_SilentMutations/full.data$CDR2_Nucleotides *
        100
    full.data$CDR2_NonsilentMutation_Rate <- full.data$CDR2_NonsilentMutations/full.data$CDR2_Nucleotides *
        100
    full.data$FR3_Mutation_Rate <- full.data$FR3_Mutations/full.data$FR3_Nucleotides *
        100
    full.data$FR3_SilentMutation_Rate <- full.data$FR3_SilentMutations/full.data$FR3_Nucleotides *
        100
    full.data$FR3_NonsilentMutation_Rate <- full.data$FR3_NonsilentMutations/full.data$FR3_Nucleotides *
        100

    extended.mut.info <- data.frame()

    for (n in 1:numsamples) {
        df.em <- full.data[full.data$Sample == samplelist[n], ]
        extended.mut.info[n, 1] <- round(mean(df.em$V_SilentMutation_Rate, na.rm = TRUE),
            2)
        extended.mut.info[n, 2] <- round(mean(df.em$V_NonsilentMutation_Rate, na.rm = TRUE),
            2)
        extended.mut.info[n, 3] <- round(mean(df.em$FR1_Mutation_Rate, na.rm = TRUE),
            2)
        extended.mut.info[n, 4] <- round(mean(df.em$FR1_SilentMutation_Rate, na.rm = TRUE),
            2)
        extended.mut.info[n, 5] <- round(mean(df.em$FR1_NonsilentMutation_Rate, na.rm = TRUE),
            2)
        extended.mut.info[n, 6] <- round(mean(df.em$CDR1_Mutation_Rate, na.rm = TRUE),
            2)
        extended.mut.info[n, 7] <- round(mean(df.em$CDR1_SilentMutation_Rate, na.rm = TRUE),
            2)
        extended.mut.info[n, 8] <- round(mean(df.em$CDR1_NonsilentMutation_Rate,
            na.rm = TRUE), 2)
        extended.mut.info[n, 9] <- round(mean(df.em$FR2_Mutation_Rate, na.rm = TRUE),
            2)
        extended.mut.info[n, 10] <- round(mean(df.em$FR2_SilentMutation_Rate, na.rm = TRUE),
            2)
        extended.mut.info[n, 11] <- round(mean(df.em$FR2_NonsilentMutation_Rate,
            na.rm = TRUE), 2)
        extended.mut.info[n, 12] <- round(mean(df.em$CDR2_Mutation_Rate, na.rm = TRUE),
            2)
        extended.mut.info[n, 13] <- round(mean(df.em$CDR2_SilentMutation_Rate, na.rm = TRUE),
            2)
        extended.mut.info[n, 14] <- round(mean(df.em$CDR2_NonsilentMutation_Rate,
            na.rm = TRUE), 2)
        extended.mut.info[n, 15] <- round(mean(df.em$FR3_Mutation_Rate, na.rm = TRUE),
            2)
        extended.mut.info[n, 16] <- round(mean(df.em$FR3_SilentMutation_Rate, na.rm = TRUE),
            2)
        extended.mut.info[n, 17] <- round(mean(df.em$FR3_NonsilentMutation_Rate,
            na.rm = TRUE), 2)
    }

    colnames(extended.mut.info) <- c("V_SilentMutation_Rate", "V_NonsilentMutation_Rate",
        "FR1_Mutation_Rate", "FR1_SilentMutation_Rate", "FR1_NonsilentMutation_Rate",
        "CDR1_Mutation_Rate", "CDR1_SilentMutation_Rate", "CDR1_NonsilentMutation_Rate",
        "FR2_Mutation_Rate", "FR2_SilentMutation_Rate", "FR2_NonsilentMutation_Rate",
        "CDR2_Mutation_Rate", "CDR2_SilentMutation_Rate", "CDR2_NonsilentMutation_Rate",
        "FR3_Mutation_Rate", "FR3_SilentMutation_Rate", "FR3_NonsilentMutation_Rate")

    master.data <- cbind(master.data, extended.mut.info)


    vgene.table <- table(full.data$Sample, full.data$Vgene)
    vgene.perc <- prop.table(vgene.table, 1) * 100

    master.data$IGHV4_34 <- vgene.perc[, "IGHV4-34"]
    V434 <- full.data[full.data$Vgene == "IGHV4-34", ]

    substrRight <- function(x, n) {
        substr(x, nchar(x) - n + 1, nchar(x))
    }


    V434$AVY <- substrRight(V434$AA_FR1, 3)
    AVY.table <- table(V434$Sample, V434$AVY)
    if (numsamples == nrow(AVY.table)) {
        master.data$AVY_perc <- AVY.table[, "AVY"]/rowSums(AVY.table) * 100
    } else {
        master.data$AVY_perc <- "N/A"
    }

    # Calculate the percent of each VH gene and JH gene

    gene.info <- data.frame()

    for (n in 1:numsamples) {
        df.popgenes <- full.data[full.data$Sample == samplelist[n], ]
        IGHV1.2 <- df.popgenes[df.popgenes$Vgene == "IGHV1-2", ]
        IGHV1.3 <- df.popgenes[df.popgenes$Vgene == "IGHV1-3", ]
        IGHV1.8 <- df.popgenes[df.popgenes$Vgene == "IGHV1-8", ]
        IGHV1.18 <- df.popgenes[df.popgenes$Vgene == "IGHV1-18", ]
        IGHV1.24 <- df.popgenes[df.popgenes$Vgene == "IGHV1-24", ]
        IGHV1.45 <- df.popgenes[df.popgenes$Vgene == "IGHV1-45", ]
        IGHV1.46 <- df.popgenes[df.popgenes$Vgene == "IGHV1-46", ]
        IGHV1.58 <- df.popgenes[df.popgenes$Vgene == "IGHV1-58", ]
        IGHV1.69 <- df.popgenes[df.popgenes$Vgene == "IGHV1-69", ]
        IGHV1.f <- df.popgenes[df.popgenes$Vgene == "IGHV1-f", ]
        IGHV2.5 <- df.popgenes[df.popgenes$Vgene == "IGHV2-5", ]
        IGHV2.26 <- df.popgenes[df.popgenes$Vgene == "IGHV2-26", ]
        IGHV2.70 <- df.popgenes[df.popgenes$Vgene == "IGHV2-70", ]
        IGHV3.7 <- df.popgenes[df.popgenes$Vgene == "IGHV3-7", ]
        IGHV3.9 <- df.popgenes[df.popgenes$Vgene == "IGHV3-9", ]
        IGHV3.11 <- df.popgenes[df.popgenes$Vgene == "IGHV3-11", ]
        IGHV3.13 <- df.popgenes[df.popgenes$Vgene == "IGHV3-13", ]
        IGHV3.15 <- df.popgenes[df.popgenes$Vgene == "IGHV3-15", ]
        IGHV3.20 <- df.popgenes[df.popgenes$Vgene == "IGHV3-20", ]
        IGHV3.21 <- df.popgenes[df.popgenes$Vgene == "IGHV3-21", ]
        IGHV3.22 <- df.popgenes[df.popgenes$Vgene == "IGHV3-22", ]
        IGHV3.23 <- df.popgenes[df.popgenes$Vgene == "IGHV3-23", ]
        IGHV3.30 <- df.popgenes[df.popgenes$Vgene == "IGHV3-30", ]
        IGHV3.30.3 <- df.popgenes[df.popgenes$Vgene == "IGHV3-30-3", ]
        IGHV3.33 <- df.popgenes[df.popgenes$Vgene == "IGHV3-33", ]
        IGHV3.43 <- df.popgenes[df.popgenes$Vgene == "IGHV3-43", ]
        IGHV3.43D <- df.popgenes[df.popgenes$Vgene == "IGHV3-43D", ]
        IGHV3.48 <- df.popgenes[df.popgenes$Vgene == "IGHV3-48", ]
        IGHV3.49 <- df.popgenes[df.popgenes$Vgene == "IGHV3-49", ]
        IGHV3.53 <- df.popgenes[df.popgenes$Vgene == "IGHV3-53", ]
        IGHV7.4.1 <- df.popgenes[df.popgenes$Vgene == "IGHV3-64", ]
        IGHV3.66 <- df.popgenes[df.popgenes$Vgene == "IGHV3-66", ]
        IGHV3.72 <- df.popgenes[df.popgenes$Vgene == "IGHV3-72", ]
        IGHV3.73 <- df.popgenes[df.popgenes$Vgene == "IGHV3-73", ]
        IGHV3.64 <- df.popgenes[df.popgenes$Vgene == "IGHV3-64", ]
        IGHV3.74 <- df.popgenes[df.popgenes$Vgene == "IGHV3-74", ]
        IGHV3.d <- df.popgenes[df.popgenes$Vgene == "IGHV3-d", ]
        IGHV3.h <- df.popgenes[df.popgenes$Vgene == "IGHV3-h", ]
        IGHV4.4 <- df.popgenes[df.popgenes$Vgene == "IGHV4-4", ]
        IGHV4.28 <- df.popgenes[df.popgenes$Vgene == "IGHV4-28", ]
        IGHV4.30.2 <- df.popgenes[df.popgenes$Vgene == "IGHV4-30-2", ]
        IGHV4.30.4 <- df.popgenes[df.popgenes$Vgene == "IGHV4-30-4", ]
        IGHV4.31 <- df.popgenes[df.popgenes$Vgene == "IGHV4-31", ]
        IGHV4.34 <- df.popgenes[df.popgenes$Vgene == "IGHV4-34", ]
        IGHV4.39 <- df.popgenes[df.popgenes$Vgene == "IGHV4-39", ]
        IGHV4.55 <- df.popgenes[df.popgenes$Vgene == "IGHV4-55", ]
        IGHV4.59 <- df.popgenes[df.popgenes$Vgene == "IGHV4-59", ]
        IGHV4.61 <- df.popgenes[df.popgenes$Vgene == "IGHV4-61", ]
        IGHV4.b <- df.popgenes[df.popgenes$Vgene == "IGHV4-b", ]
        IGHV5.51 <- df.popgenes[df.popgenes$Vgene == "IGHV5-51", ]
        IGHV6.1 <- df.popgenes[df.popgenes$Vgene == "IGHV6-1", ]
        IGHJ1 <- df.popgenes[df.popgenes$JGene == "IGHJ1", ]
        IGHJ2 <- df.popgenes[df.popgenes$JGene == "IGHJ2", ]
        IGHJ3 <- df.popgenes[df.popgenes$JGene == "IGHJ3", ]
        IGHJ4 <- df.popgenes[df.popgenes$JGene == "IGHJ4", ]
        IGHJ5 <- df.popgenes[df.popgenes$JGene == "IGHJ5", ]
        IGHJ6 <- df.popgenes[df.popgenes$JGene == "IGHJ6", ]
        gene.info[n, 1] <- round(nrow(IGHV1.2)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 2] <- round(nrow(IGHV1.3)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 3] <- round(nrow(IGHV1.8)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 4] <- round(nrow(IGHV1.18)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 5] <- round(nrow(IGHV1.24)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 6] <- round(nrow(IGHV1.45)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 7] <- round(nrow(IGHV1.46)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 8] <- round(nrow(IGHV1.58)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 9] <- round(nrow(IGHV1.69)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 10] <- round(nrow(IGHV1.f)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 11] <- round(nrow(IGHV2.5)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 12] <- round(nrow(IGHV2.26)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 13] <- round(nrow(IGHV2.70)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 14] <- round(nrow(IGHV3.7)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 15] <- round(nrow(IGHV3.9)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 16] <- round(nrow(IGHV3.11)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 17] <- round(nrow(IGHV3.13)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 18] <- round(nrow(IGHV3.15)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 19] <- round(nrow(IGHV3.20)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 20] <- round(nrow(IGHV3.21)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 21] <- round(nrow(IGHV3.22)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 22] <- round(nrow(IGHV3.23)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 23] <- round(nrow(IGHV3.30)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 24] <- round(nrow(IGHV3.30.3)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 25] <- round(nrow(IGHV3.33)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 26] <- round(nrow(IGHV3.43)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 27] <- round(nrow(IGHV3.43D)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 28] <- round(nrow(IGHV3.48)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 29] <- round(nrow(IGHV3.49)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 30] <- round(nrow(IGHV3.53)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 31] <- round(nrow(IGHV3.64)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 32] <- round(nrow(IGHV3.66)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 33] <- round(nrow(IGHV3.72)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 34] <- round(nrow(IGHV3.73)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 35] <- round(nrow(IGHV3.74)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 36] <- round(nrow(IGHV3.d)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 37] <- round(nrow(IGHV3.h)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 38] <- round(nrow(IGHV4.4)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 39] <- round(nrow(IGHV4.28)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 40] <- round(nrow(IGHV4.30.2)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 41] <- round(nrow(IGHV4.30.4)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 42] <- round(nrow(IGHV4.31)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 43] <- round(nrow(IGHV4.34)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 44] <- round(nrow(IGHV4.39)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 45] <- round(nrow(IGHV4.55)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 46] <- round(nrow(IGHV4.59)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 47] <- round(nrow(IGHV4.61)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 48] <- round(nrow(IGHV4.b)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 49] <- round(nrow(IGHV5.51)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 50] <- round(nrow(IGHV6.1)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 51] <- round(nrow(IGHJ1)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 52] <- round(nrow(IGHJ2)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 53] <- round(nrow(IGHJ3)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 54] <- round(nrow(IGHJ4)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 55] <- round(nrow(IGHJ5)/nrow(df.popgenes) * 100, 2)
        gene.info[n, 56] <- round(nrow(IGHJ6)/nrow(df.popgenes) * 100, 2)
    }



    colnames(gene.info) <- c("IGHV1-2", "IGHV1-3", "IGHV1-8", "IGHV1-18", "IGHV1-24",
        "IGHV1-45", "IGHV1-46", "IGHV1-58", "IGHV1-69", "IGHV1-f", "IGHV2-5", "IGHV2-26",
        "IGHV2-70", "IGHV3-7", "IGHV3-9", "IGHV3-11", "IGHV3-13", "IGHV3-15", "IGHV3-20",
        "IGHV3-21", "IGHV3-22", "IGHV3-23", "IGHV3-30", "IGHV3-30-3", "IGHV3-33",
        "IGHV3-43", "IGHV3-43D", "IGHV3-48", "IGHV3-49", "IGHV3-53", "IGHV3-64",
        "IGHV3-66", "IGHV3-72", "IGHV3-73", "IGHV3-74", "IGHV3-d", "IGHV3-h", "IGHV4-4",
        "IGHV4-28", "IGHV4-30-2", "IGHV4-30-4", "IGHV4-31", "IGHV4-34", "IGHV4-39",
        "IGHV4-55", "IGHV4-59", "IGHV4-61", "IGHV4-b", "IGHV5-51", "IGHV6-1", "IGHJ1",
        "IGHJ2", "IGHJ3", "IGHJ4", "IGHJ5", "IGHJ6")



    master.data <- cbind(master.data, gene.info)



    # make a data frame with percent of matching clones in each

    percof.df <- data.frame()

    for (m in 1:numsamples) {
        df.x <- filter(lineages.union, lineages.union[[m + 2]] > 0)
        lin.x <- nrow(df.x)
        for (n in 1:numsamples) {
            df.x_y <- filter(df.x, df.x[[n + 2]] > 0)
            lin.x_y <- nrow(df.x_y)
            perc.x_y <- lin.x_y/lin.x * 100
            percof.df[m, n] <- round(perc.x_y, 2)
        }
        colnames(percof.df)[m] <- colnames(lineages.union)[m + 2]
    }

    # rename columns
    for (l in 1:numsamples) {
        colnames(percof.df)[l] <- paste0("Percent_Matching_with_", colnames(percof.df)[l])
    }


    master.data <- cbind(master.data, percof.df)




    # Do the same for individual isotypes


    lineage.data <- non.single[, c("population", "lineageID", "Vgene")]
    table.pop <- table(lineage.data$lineageID, lineage.data$population)
    table.vgene <- table(lineage.data$lineageID, lineage.data$Vgene)
    pop.df <- cbind(rownames(table.pop), as.data.frame.matrix(table.pop, row.names = NULL))
    colnames(pop.df)[1] <- "LineageID"
    v.df <- cbind(rownames(table.vgene), as.data.frame.matrix(table.vgene, row.names = NULL))
    colnames(v.df)[1] <- "LineageID"
    table.vgene.melt <- melt(v.df, key = LineageID)
    v.melt <- table.vgene.melt[table.vgene.melt$value != 0, ]
    lineages.union <- cbind(v.melt$variable, pop.df)

iso.df <- data.frame()
iso.df2 <- data.frame()
isovar <- c('M', 'G', 'A', 'E')

z <- 1
        for (n in 1:numsamples){
            for (m in 1:numsamples){

                for (k in 1:4){
                        for (j in 1:4){
                            non.single.x <- non.single[non.single$isotype == isovar[k] & non.single$population == poplist[n],]
                            non.single.y <- non.single[non.single$isotype == isovar[j] & non.single$population == poplist[m],]
                            lineage.data.x <- non.single.x[, c("population", "lineageID", "Vgene")]
                            lineage.data.y <- non.single.y[, c("population", "lineageID", "Vgene")]

                            perc.x_y == 'NA'
                            pop1 <- paste0(poplist[n], ".Ig", isovar[k])
                            pop2 <- paste0(poplist[m], ".Ig", isovar[j])
                            if (nrow(lineage.data.x) > 0 & nrow(lineage.data.y) > 0 & pop1 != pop2){
                                lineage.data.x$Pop <- pop1
                                lineage.data.y$Pop <- pop2
                                lineage.data.z <- rbind(lineage.data.x, lineage.data.y)
                                table.pop.x <- table(lineage.data.z$lineageID, lineage.data.z$Pop)
                                table.vgene.x <- table(lineage.data.z$lineageID, lineage.data.z$Vgene)
                                pop.df.x <- cbind(rownames(table.pop.x), as.data.frame.matrix(table.pop.x, row.names = NULL))
                                colnames(pop.df.x)[1] <- "LineageID"
                                v.df.x <- cbind(rownames(table.vgene.x), as.data.frame.matrix(table.vgene.x, row.names = NULL))
                                colnames(v.df.x)[1] <- "LineageID"
                                table.vgene.melt.x <- melt(v.df.x, key = LineageID)
                                v.melt.x <- table.vgene.melt.x[table.vgene.melt.x$value != 0, ]
                                lineages.union.x <- cbind(v.melt.x$variable, pop.df.x)

                                        df.only4 <- filter(lineages.union.x, lineages.union.x[[4]] > 0)
                                        lin.only4 <- nrow(df.only4)
                                        df.x_y <- filter(df.only4, df.only4[[3]] > 0)
                                        lin.x_y <- nrow(df.x_y)
                                        perc.x_y <- lin.x_y/lin.only4 * 100
                            } else if (pop1 == pop2){
                                    perc.x_y <- 100
                            } else {
                                    perc.x_y <- 'NA'
                            }

                            iso.df2[z, 1] <- paste0(poplist[n], ".Ig", isovar[k])
                            iso.df2[z, 2] <- paste0(poplist[m], ".Ig", isovar[j])
                            iso.df2[z, 3] <- perc.x_y
                            z <- z + 1

                    }
                }
            }
        }
tab.iso <- dcast(iso.df2, V1 ~ V2)
rownames(tab.iso) <- tab.iso[,1]
iso.df <- tab.iso[,-1]
conn.out <- paste0(file.out, "Percent_Conn_", sub(".*/", "", location2), ".txt")
write.table(iso.df, conn.out, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)



    percof.df <- data.frame()

    for (m in 1:numsamples) {
        df.x <- filter(lineages.union, lineages.union[[m + 2]] > 0)
        lin.x <- nrow(df.x)
        for (n in 1:numsamples) {
            df.x_y <- filter(df.x, df.x[[n + 2]] > 0)
            lin.x_y <- nrow(df.x_y)
            perc.x_y <- lin.x_y/lin.x * 100
            percof.df[m, n] <- round(perc.x_y, 2)
        }
        colnames(percof.df)[m] <- colnames(lineages.union)[m + 2]
    }

    # rename columns
    for (l in 1:numsamples) {
        colnames(percof.df)[l] <- paste0("Percent_Matching_with_", colnames(percof.df)[l])
    }


    master.data <- cbind(master.data, percof.df)


    # Calculate the Morisita index for each pair of samples into a new matrix

    pops.df <- lineages.union[, 3:ncol(lineages.union)]
    mor <- niche.overlap(pops.df, method = c("morisita"))
    mor.mat <- as.matrix(round(mor, 3))
    for (i in 1:numsamples) {
        mor.mat[i, i] <- 1
    }

    # rename columns

    mor.df <- as.data.frame(mor.mat)
    for (l in 1:(numsamples)) {
        colnames(mor.mat)[l] <- paste0("MI_", colnames(mor.mat)[l])
    }



    # Add Morisita index values to master.data data frame

    master.data <- cbind(master.data, mor.mat)


    # Add Mextended clone mutation info
    clone.mut <- data.frame()
    for (i in 1:numsamples) {
        df1 <- non.single[non.single$Sample == samplelist[i], ]
        df.lineagedata <- summaryBy(V_Mutations ~ lineageID, data = df1, FUN = list(mean,
            max, min, median, sd))
        clone.mut[i, 1] <- samplelist[i]
        clone.mut[i, 2] <- round(mean(df.lineagedata$V_Mutations.mean), 2)
        clone.mut[i, 3] <- median(df.lineagedata$V_Mutations.median)
        clone.mut[i, 4] <- round(mean(df.lineagedata$V_Mutations.min), 2)
    }
    colnames(clone.mut) <- c("Sample", "Mean Clonal Mutation", "Median Clonal Mutation",
        "Mean Trunk Length")

    master.data <- cbind(master.data, clone.mut[, c(2:4)])

    clone.mut.iso <- data.frame()
    for (i in 1:numsamples) {
        Ig <- non.single[non.single$isotype != "U", ]
        IgG <- non.single[non.single$isotype == "G", ]
        IgA <- non.single[non.single$isotype == "A", ]
        IgM <- non.single[non.single$isotype == "M", ]
        IgE <- non.single[non.single$isotype == "E", ]

        Ig <- Ig[Ig$Sample == samplelist[i], ]
        IgG <- IgG[IgG$Sample == samplelist[i], ]
        IgA <- IgA[IgA$Sample == samplelist[i], ]
        IgM <- IgM[IgM$Sample == samplelist[i], ]
        IgE <- IgE[IgE$Sample == samplelist[i], ]

        clone.mut.iso[i, 1] <- samplelist[i]

        if (nrow(Ig) > 0) {
            df.lineagedata.Ig <- summaryBy(V_Mutations ~ lineageID, data = Ig, FUN = list(mean,
                max, min, median, sd))
            clone.mut.iso[i, 2] <- round(mean(df.lineagedata.Ig$V_Mutations.mean),
                2)
            clone.mut.iso[i, 3] <- median(df.lineagedata.Ig$V_Mutations.median)
            clone.mut.iso[i, 4] <- round(mean(df.lineagedata.Ig$V_Mutations.min),
                2)
        } else {
            clone.mut.iso[i, 2] <- "N/A"
            clone.mut.iso[i, 3] <- "N/A"
            clone.mut.iso[i, 4] <- "N/A"
        }
        if (nrow(IgG) > 0) {
            df.lineagedata.IgG <- summaryBy(V_Mutations ~ lineageID, data = IgG,
                FUN = list(mean, max, min, median, sd))
            clone.mut.iso[i, 5] <- round(mean(df.lineagedata.IgG$V_Mutations.mean),
                2)
            clone.mut.iso[i, 6] <- median(df.lineagedata.IgG$V_Mutations.median)
            clone.mut.iso[i, 7] <- round(mean(df.lineagedata.IgG$V_Mutations.min),
                2)
        } else {
            clone.mut.iso[i, 5] <- "N/A"
            clone.mut.iso[i, 6] <- "N/A"
            clone.mut.iso[i, 7] <- "N/A"
        }
        if (nrow(IgA) > 0) {
            df.lineagedata.IgA <- summaryBy(V_Mutations ~ lineageID, data = IgA,
                FUN = list(mean, max, min, median, sd))
            clone.mut.iso[i, 8] <- round(mean(df.lineagedata.IgA$V_Mutations.mean),
                2)
            clone.mut.iso[i, 9] <- median(df.lineagedata.IgA$V_Mutations.median)
            clone.mut.iso[i, 10] <- round(mean(df.lineagedata.IgA$V_Mutations.min),
                2)
        } else {
            clone.mut.iso[i, 8] <- "N/A"
            clone.mut.iso[i, 9] <- "N/A"
            clone.mut.iso[i, 10] <- "N/A"
        }
        if (nrow(IgM) > 0) {
            df.lineagedata.IgM <- summaryBy(V_Mutations ~ lineageID, data = IgM,
                FUN = list(mean, max, min, median, sd))
            clone.mut.iso[i, 11] <- round(mean(df.lineagedata.IgM$V_Mutations.mean),
                2)
            clone.mut.iso[i, 12] <- median(df.lineagedata.IgM$V_Mutations.median)
            clone.mut.iso[i, 13] <- round(mean(df.lineagedata.IgM$V_Mutations.min),
                2)
        } else {
            clone.mut.iso[i, 11] <- "N/A"
            clone.mut.iso[i, 12] <- "N/A"
            clone.mut.iso[i, 13] <- "N/A"
        }
        if (nrow(IgE) > 0) {
            df.lineagedata.IgE <- summaryBy(V_Mutations ~ lineageID, data = IgE,
                FUN = list(mean, max, min, median, sd))
            clone.mut.iso[i, 14] <- round(mean(df.lineagedata.IgE$V_Mutations.mean),
                2)
            clone.mut.iso[i, 15] <- median(df.lineagedata.IgE$V_Mutations.median)
            clone.mut.iso[i, 16] <- round(mean(df.lineagedata.IgE$V_Mutations.min),
                2)
        } else {
            clone.mut.iso[i, 14] <- "N/A"
            clone.mut.iso[i, 15] <- "N/A"
            clone.mut.iso[i, 16] <- "N/A"
        }
    }
    colnames(clone.mut.iso) <- c("Sample", "Ig Mean Clonal Mutation", "Ig Median Clonal Mutation",
        "Ig Mean Trunk Length", "IgG Mean Clonal Mutation", "IgG Median Clonal Mutation",
        "IgG Mean Trunk Length", "IgA Mean Clonal Mutation", "IgA Median Clonal Mutation",
        "IgA Mean Trunk Length", "IgM Mean Clonal Mutation", "IgM Median Clonal Mutation",
        "IgM Mean Trunk Length", "IgE Mean Clonal Mutation", "IgE Median Clonal Mutation",
        "IgE Mean Trunk Length")

    master.data <- cbind(master.data, clone.mut.iso[, c(2:ncol(clone.mut.iso))])

    write.table(master.data, file = master.out, quote = FALSE, sep = "\t", row.names = FALSE)
}
