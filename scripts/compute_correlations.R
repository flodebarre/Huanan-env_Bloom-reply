# Compute correlations for combinations of viruses, hosts, date ranges
# (Not used in the paper but kept for legacy)

# Load data ####

source("load_data.R")

# Load functions ####

source("functions_correlations.R")

# Compute all correlations #### 

# Combinations of parameters
viruses <- c("SARS-CoV-2", "Influenza A virus H3N2", "Bamboo rat circovirus", "Bamboo rat coronavirus", "Raccoon dog amdovirus", "Civet astrovirus", "Civet kobuvirus", "Myocastor coypus polyomavirus 1")
viruses <- unique(viral$Virus_shorter)
parms <- expand.grid(virus = viruses, animal = chordates, 
                     stringsAsFactors = FALSE)

# /!\ This takes a few minutes
res1.all <- lapply(as.list(seq_len(nrow(parms))), function(i) ww(i, data, parms))
res1.Jan01 <- lapply(as.list(seq_len(nrow(parms))), function(i) ww(i, data[data$Sampling.date == "2020-01-01", ], parms))
res1.Jan12 <- lapply(as.list(seq_len(nrow(parms))), function(i) ww(i, data[data$Sampling.date == "2020-01-12", ], parms))
res1.Jan0112 <- lapply(as.list(seq_len(nrow(parms))), function(i) ww(i, data[data$Sampling.date <= "2020-01-12", ], parms))

parmsJB <- expand.grid(virus = "SARS-CoV-2", animal = chordatesJB, stringsAsFactors = FALSE)
JB1.all <- lapply(as.list(seq_len(nrow(parmsJB))), function(i) ww(i, dataJB, parmsJB))
JB1.Jan01 <- lapply(as.list(seq_len(nrow(parmsJB))), function(i) ww(i, dataJB[dataJB$Sampling.date == "2020-01-01", ], parmsJB))
JB1.Jan12 <- lapply(as.list(seq_len(nrow(parmsJB))), function(i) ww(i, dataJB[dataJB$Sampling.date == "2020-01-12", ], parmsJB))
JB1.Jan0112 <- lapply(as.list(seq_len(nrow(parmsJB))), function(i) ww(i, dataJB[dataJB$Sampling.date <= "2020-01-12", ], parmsJB))


# Format results ####

formatRes <- function(res, parms){
  resdf <- parms
  resdf$cors <- unlist(lapply(res, function (x) x$estimate))
  resdf$corsCI1 <- unlist(lapply(res, function (x) x$conf.int[1]))
  resdf$corsCI2 <- unlist(lapply(res, function (x) x$conf.int[2]))
  resdf$pval <- unlist(lapply(res, function (x) x$p.value))
  resdf$commonName <- dicoSp[resdf$animal]
  #resdf$bootCI1 <- unlist(lapply(res, function (x) x$`bootCI.2.5%`))
  #resdf$bootCI2 <- unlist(lapply(res, function (x) x$`bootCI.97.5%`))
  resdf
}

resdf1.all <- formatRes(res1.all, parms)
resdf1.Jan01 <- formatRes(res1.Jan01, parms)
resdf1.Jan12 <- formatRes(res1.Jan12, parms)
resdf1.Jan0112 <- formatRes(res1.Jan0112, parms)

JB1.all <- formatRes(JB1.all, parmsJB)
JB1.Jan01 <- formatRes(JB1.Jan01, parmsJB)
JB1.Jan12 <- formatRes(JB1.Jan12, parmsJB)
JB1.Jan0112 <- formatRes(JB1.Jan0112, parmsJB)

# Export results ####

write.csv(resdf1.all, file = "../results/corACC1_all.csv", row.names = FALSE)
write.csv(resdf1.Jan01, file = "../results/corACC1_Jan01.csv", row.names = FALSE)
write.csv(resdf1.Jan12, file = "../results/corACC1_Jan12.csv", row.names = FALSE)
write.csv(resdf1.Jan0112, file = "../results/corACC1_Jan0112.csv", row.names = FALSE)

write.csv(JB1.all, file = "../results/corJB1_all.csv", row.names = FALSE)
write.csv(JB1.Jan01, file = "../results/corJB1_Jan01.csv", row.names = FALSE)
write.csv(JB1.Jan12, file = "../results/corJB1_Jan12.csv", row.names = FALSE)
write.csv(JB1.Jan0112, file = "../results/corJB1_Jan0112.csv", row.names = FALSE)

