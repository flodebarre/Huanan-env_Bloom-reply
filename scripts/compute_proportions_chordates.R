# Compute proportions of given species among chordate reads, 
# in JB's data, and order them. 

# Load data ####

source("load_data.R")

# Computations ####

# Compute total reads chordates
tmp <- dataJB[dataJB$species != "SARS-CoV-2", ]
agg <- aggregate(x = as.numeric(tmp$aligned_reads), by = list(Sample.ID = tmp$Sample.ID), FUN = sum)
names(agg)[names(agg) == "x"] <- "tot_reads_chordates"
dataJB <- merge(dataJB, agg, all.x = TRUE)

# Compute proportion of a given species' reads among Chordates
computeProportionSpecies <- function(species = "Nyctereutes procyonoides", proportion = 0.2){

  # Data for the species
  dataJBspecies <- dataJB[dataJB$species == species, ]
  nrow(dataJBspecies)
  sort(table(dataJBspecies$Sample.ID))
  
  # Compute proportions
  dataJBspecies$prop.chordates <- dataJBspecies$aligned_reads / dataJBspecies$tot_reads_chordates
  dataJBspecies$prop.total <- dataJBspecies$aligned_reads / dataJBspecies$preprocessed_reads
  nrow(dataJBspecies)
  dataJBspecies$prop.chordates
  
  # Add information about SC2
  tmp <- dataJB[dataJB$species == "SARS-CoV-2", ]
  names(tmp)[names(tmp) == "aligned_reads"] <- "SARS2.reads"
  dataJBspecies <- merge(dataJBspecies, tmp[, c("Sample.ID", "SARS2.reads")], all.x = TRUE)
  nrow(dataJBspecies)
  dataJBspecies$prop.chordates
  
  # Add metadata
  dataJBspecies <- merge(dataJBspecies, data[, c("Sample.ID", "Sampling.date", "address")], all.x = TRUE)
  
  
  # Order the results
  list(prop.chordates = dataJBspecies[order(dataJBspecies$SARS2.reads >0, dataJBspecies$prop.chordates, decreasing = TRUE), c("Sample.ID", "Lab.code", "address", "Sampling.date", "prop.chordates", "SARS2.reads")], 
       prop.all = dataJBspecies[order(dataJBspecies$SARS2.reads >0, dataJBspecies$prop.total, decreasing = TRUE), c("Sample.ID", "Lab.code", "address", "Sampling.date", "prop.total", "SARS2.reads")], 
       table.chordates = table(speciesProp = dataJBspecies$prop.chordates > proportion, SC2 = dataJBspecies$SARS2.reads > 0), 
       table.all = table(speciesProp = dataJBspecies$prop.total > proportion, SC2 = dataJBspecies$SARS2.reads > 0))
}

# Compute for Raccoon dogs
RD <- computeProportionSpecies("Nyctereutes procyonoides", 0.2)
RD
write.csv(RD$prop.chordates, file = "../results/prop-chordates_raccoon-dog.csv", row.names = FALSE)

# Compute for humans
HS <- computeProportionSpecies("Homo sapiens", 0.2)
HS
write.csv(HS$prop.chordates, file = "../results/prop-chordates_human.csv", row.names = FALSE)
