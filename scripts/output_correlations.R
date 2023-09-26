# Load data ####
source("load_data.R")

# Load functions ####
source("functions_correlations.R")

# Print some results ####

wrapperCor(data[data$Sampling.date == "2020-01-12", ], virus = "SARS-CoV-2", animal = "Homo.sapiens", bootstrap = TRUE)
wrapperCor(dataJB[dataJB$Sampling.date == "2020-01-12", ], virus = "SARS-CoV-2", animal = "Homo.sapiens", bootstrap = TRUE)

wrapperCor(data[data$Sampling.date == "2020-01-12", ], virus = "SARS-CoV-2", animal = "Nyctereutes.procyonoides", bootstrap = TRUE)
wrapperCor(dataJB[dataJB$Sampling.date == "2020-01-12", ], virus = "SARS-CoV-2", animal = "Nyctereutes.procyonoides", bootstrap = TRUE)

wrapperCor(data, virus = "Influenza A virus H3N2", animal = "Homo.sapiens", bootstrap = TRUE)

wrapperCor(data[data$Sampling.date == "2020-01-01", ], virus = "SARS-CoV-2", animal = "Homo.sapiens", bootstrap = TRUE)
wrapperCor(data[data$Sampling.date == "2020-01-01", ], virus = "SARS-CoV-2", animal = "Homo.sapiens", bootstrap = TRUE, remove0.x = TRUE)
wrapperCor(data[data$Sampling.date == "2020-01-01" & data$SARS2.reads > 0, ], virus = "SARS-CoV-2", animal = "Homo.sapiens", bootstrap = TRUE)

# Prepare output for LaTeX ####
latex <- ""

# JB, Jan 01, SC2 vs HS
tmp <- wrapperCor(dataJB[dataJB$Sampling.date == "2020-01-01", ], virus = "SARS-CoV-2", animal = "Homo.sapiens", bootstrap = TRUE, printOutput = TRUE)
latex <- rbind(latex, paste("\\def \\JBFirstSCHS {", tmp$txt, "}", sep = ""))

# JB, Jan 01, SC2 vs HS, without F13 and F54
tmp <- wrapperCor(dataJB[dataJB$Sampling.date == "2020-01-01" & !is.element(dataJB$Lab.code, c("F13", "F54")), ], virus = "SARS-CoV-2", animal = "Homo.sapiens", bootstrap = TRUE, printOutput = TRUE)
latex <- rbind(latex, paste("\\def \\JBFirstSCHStwo {", tmp$txt, "}", sep = ""))

# JB, Jan 12, SC2 vs HS
tmp <- wrapperCor(dataJB[dataJB$Sampling.date == "2020-01-12", ], virus = "SARS-CoV-2", animal = "Homo.sapiens", bootstrap = TRUE, printOutput = TRUE)
latex <- rbind(latex, paste("\\def \\JBTwelfthSCHS {", tmp$txt, "}", sep = ""))

# JB, Jan 12, SC2 vs RD
tmp <- wrapperCor(dataJB[dataJB$Sampling.date == "2020-01-12", ], virus = "SARS-CoV-2", animal = "Nyctereutes.procyonoides", bootstrap = TRUE)
latex <- rbind(latex, paste("\\def \\JBTwelfthSCRD {", tmp$txt, "}", sep = ""))

# JB, Jan 01, SC2 vs largemouth bass
tmp <- wrapperCor(dataJB[dataJB$Sampling.date == "2020-01-01", ], virus = "SARS-CoV-2", animal = "Micropterus.salmoides", bootstrap = TRUE)
latex <- rbind(latex, paste("\\def \\JBFirstSCLB {", tmp$txt, "}", sep = ""))

# JB, Jan 01, SC2 vs largemouth bass, without F13 and F54
tmp <- wrapperCor(dataJB[dataJB$Sampling.date == "2020-01-01" & !is.element(dataJB$Lab.code, c("F13", "F54")), ], virus = "SARS-CoV-2", animal = "Micropterus.salmoides", bootstrap = TRUE)
latex <- rbind(latex, paste("\\def \\JBFirstSCLBtwo {", tmp$txt, "}", sep = ""))


# ACC, all dates, H3N2 vs HS
tmp <- wrapperCor(data, virus = "Influenza A virus H3N2", animal = "Homo.sapiens", bootstrap = TRUE)
latex <- rbind(latex, paste("\\def \\ACCAllHNHS {", tmp$txt, "}", sep = ""))

# ACC, Jan 12, SC2 vs HS
tmp <- wrapperCor(data[data$Sampling.date == "2020-01-12", ], virus = "SARS-CoV-2", animal = "Homo.sapiens", bootstrap = TRUE)
latex <- rbind(latex, paste("\\def \\ACCTwelfthSCHS {", tmp$txt, "}", sep = ""))

# ACC, Jan 12, SC2 vs RD
tmp <- wrapperCor(data[data$Sampling.date == "2020-01-12", ], virus = "SARS-CoV-2", animal = "Nyctereutes.procyonoides", bootstrap = TRUE)
latex <- rbind(latex, paste("\\def \\ACCTwelfthSCRD {", tmp$txt, "}", sep = ""))

# ACC, Jan 01, SC2 vs HS
tmp <- wrapperCor(data[data$Sampling.date == "2020-01-01", ], virus = "SARS-CoV-2", animal = "Homo.sapiens", bootstrap = TRUE)
latex <- rbind(latex, paste("\\def \\ACCFirstSCHS {", tmp$txt, "}", sep = ""))

# ACC, Jan 01, SC2 vs spotted bass
tmp <- wrapperCor(data[data$Sampling.date == "2020-01-01", ], virus = "SARS-CoV-2", animal = "Micropterus.punctulatus", bootstrap = TRUE)
latex <- rbind(latex, paste("\\def \\ACCFirstSCSB {", tmp$txt, "}", sep = ""))

write(latex, "../paper/commands_outputR.tex")


# View some Theil-Sen results ####
# -> Super wide CIs!

wrapperCor(dataJB[dataJB$Sampling.date == "2020-01-01", ], virus = "SARS-CoV-2", animal = "Micropterus.salmoides", bootstrap = TRUE, printOutput = TRUE, method = "theil-sen")
wrapperCor(dataJB[dataJB$Sampling.date == "2020-01-01" & !is.element(dataJB$Lab.code, c("F13", "F54")), ], virus = "SARS-CoV-2", animal = "Micropterus.salmoides", bootstrap = TRUE, printOutput = TRUE, method = "theil-sen")

wrapperCor(dataJB[dataJB$Sampling.date == "2020-01-01", ], virus = "SARS-CoV-2", animal = "Nyctereutes.procyonoides", bootstrap = TRUE, printOutput = TRUE, method = "theil-sen")
wrapperCor(dataJB[dataJB$Sampling.date == "2020-01-01" & !is.element(dataJB$Lab.code, c("F13", "F54")), ], virus = "SARS-CoV-2", animal = "Nyctereutes.procyonoides", bootstrap = TRUE, printOutput = TRUE, method = "theil-sen")

wrapperCor(dataJB[dataJB$Sampling.date == "2020-01-01", ], virus = "SARS-CoV-2", animal = "Homo.sapiens", bootstrap = TRUE, printOutput = TRUE, method = "theil-sen")
wrapperCor(dataJB[dataJB$Sampling.date == "2020-01-01" & !is.element(dataJB$Lab.code, c("F13", "F54")), ], virus = "SARS-CoV-2", animal = "Homo.sapiens", bootstrap = TRUE, printOutput = TRUE, method = "theil-sen")

