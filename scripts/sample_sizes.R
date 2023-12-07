# Compute sample sizes

# Load data 
source("load_data.R")
metadata$Sampling.date <- as.Date(metadata$Sampling.date)

# Samples around E-10-29-2 ####

# Subset of data from inside the market
tmp <- metadata[is.element(metadata$Sample_location, c("West", "East")), ]

# Stats
aggregate(tmp$SARS.CoV.2.qPCR.result, by = list(tmp$Sampling.date), FUN = function(x) mean(x == "Positive"))
aggregate(tmp$SARS.CoV.2.qPCR.result, by = list(tmp$Sampling.date), FUN = function(x) sum(x == "Positive"))
aggregate(tmp$SARS.CoV.2.qPCR.result, by = list(tmp$Sampling.date), FUN = length)

tmp[tmp$Sampling.date == "2020-02-20" & tmp$SARS.CoV.2.qPCR.result == "Positive", ]
tmp[tmp$Sampling.date == "2020-02-15" & tmp$SARS.CoV.2.qPCR.result == "Positive", ]
tmp[tmp$Sampling.date == "2020-02-20" & tmp$SARS.CoV.2.qPCR.result == "Positive", ]
sum(tmp$Sampling.date > "2020-02-20")
