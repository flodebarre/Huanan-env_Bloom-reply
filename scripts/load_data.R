# Script to load and format the data

# Initializations ####

# Install package for Theil-Sen correlation
# (Not used in the main text)
pck <- library("deming", logical.return = TRUE) 
if(!pck){
  install.packages("deming")
}
require("deming")

# Load data ####
## |- Load ACC data ####

data <- read.csv("../data/ACC/mtSC2_by_sample.csv")
viral <- read.csv("../data/ACC/otherviruses.csv")
dataSpecies <- read.csv("../data/ACC/metadata-animals.csv")
metadata <- read.csv("../data/Liu-etal/Liu_etal_metadata.csv")

# Backup all data
data.all <- data
# Subselect data from HSM (exclude warehouses and sewage nearby) and with sequences
data <- data[is.element(data$Sample_location, c("East", "West")) & !is.na(data$Homo.sapiens), ]

## |- Load JB data ####

JB1 <- read.csv("../data/JB/mito_composition_by_sample.csv")
JB2 <- read.csv("../data/JB/sars2_aligned_by_sample.csv")
# Homogenize names
names(JB1)[names(JB1) == "sample"] <- "Lab.code"
names(JB2)[names(JB2) == "sample"] <- "Lab.code"
names(JB1)[names(JB1) == "Sample.name"] <- "Sample.ID"
names(JB2)[names(JB2) == "Sample.name"] <- "Sample.ID"

# Subselect samples from inside the market, and remove amplicon sequencing
table(JB1$Isolation.source)
idsOK <- sort(unique(JB1$Sample.ID[!is.element(JB1$Isolation.source, c("Warehouses related to west wine of HSM", "Sewerage well in the surrounding areas")) & JB1$description != "SARS-CoV-2 Amplicon based SARS-CoV-2 whole genome sequencing for cell supernatants"]))

# Group the data
JBmt <- JB1[is.element(JB1$Sample.ID, idsOK), c("Sample.ID", "Lab.code", "species", "preprocessed_reads", "aligned_reads", "Collection.date")]
names(JBmt)[names(JBmt) == "Collection.date"] <- "Sampling.date"

JBSC2 <- JB2[is.element(JB2$Sample.ID, idsOK), c("Sample.ID", "Lab.code", "preprocessed_reads", "SARS2_aligned_reads", "Collection.date")]
names(JBSC2)[names(JBSC2) == "Collection.date"] <- "Sampling.date"
names(JBSC2)[names(JBSC2) == "SARS2_aligned_reads"] <- "aligned_reads"
JBSC2$species <- "SARS-CoV-2"

cls <- c("Sample.ID", "Lab.code", "species", "preprocessed_reads", "aligned_reads", "Sampling.date")
dataJB <- rbind(JBmt[,  cls], JBSC2[, cls])

# There is an issue with A20: JB included both the metagenomic results and amplicon-based sequencing results
# We need to remove the latter, and first, to identify them
tmp <- dataJB[dataJB$species == "SARS-CoV-2" & dataJB$Lab.code == "A20", ]
# Get the one with max aligned reads
pprA20 <- tmp[tmp$aligned_reads == max(tmp$aligned_reads), "preprocessed_reads"]
# Remove A20 amplicon from the table
dataJB <- dataJB[!(dataJB$Lab.code == "A20" & dataJB$preprocessed_reads == pprA20), ]

rm(pprA20, tmp)

## |- Animal species ####

# Function to capitalize first letter
firstup <- function(x) {
  # https://stackoverflow.com/questions/18509527/first-letter-to-upper-case
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Dictionary to switch between common and latin names
dicoSp <- firstup(dataSpecies$Common_name)
names(dicoSp) <- gsub(" ", "\\.", dataSpecies$Species)

# Sets of species
chordates <- gsub(" ", "\\.", dataSpecies$Species[is.element(dataSpecies$Group, c("Bird", "Fish", "Mammal", "Reptile"))])
mammals <- gsub(" ", "\\.", dataSpecies$Species[is.element(dataSpecies$Group, c("Mammal"))])
birds <- gsub(" ", "\\.", dataSpecies$Species[is.element(dataSpecies$Group, c("Bird"))])

chordatesJB <- gsub(" ", "\\.", unique(dataJB$species))
chordatesJB <- chordatesJB[chordatesJB != "SARS-CoV-2"]

# Remove green monkey (in culture)
chordates <- chordates[chordates != "Chlorocebus.sabaeus"]
mammals <- mammals[mammals != "Chlorocebus.sabaeus"]

# Other ####
# Define color palette
colS <- "black" # Significant
colNS <- gray(0.6) # Non significant
