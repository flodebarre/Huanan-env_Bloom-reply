# Load results ####
ACC1.all <- read.csv(file = "../results/corACC1_all.csv")
ACC1.Jan01 <- read.csv(file = "../results/corACC1_Jan01.csv")
ACC1.Jan12 <- read.csv(file = "../results/corACC1_Jan12.csv")
ACC1.Jan0112 <- read.csv(file = "../results/corACC1_Jan0112.csv")


JB1.all <- read.csv(file = "../results/corJB1_all.csv")
JB1.Jan01 <- read.csv(file = "../results/corJB1_Jan01.csv")
JB1.Jan12 <- read.csv(file = "../results/corJB1_Jan12.csv")
JB1.Jan0112 <- read.csv(file = "../results/corJB1_Jan0112.csv")

ACC1.Jan12[ACC1.Jan12$virus == "SARS-CoV-2" & ACC1.Jan12$animal == "Homo.sapiens", ]
ACC1.all[ACC1.all$virus == "Influenza A virus H3N2" & ACC1.all$animal == "Homo.sapiens", ]

# Load initial data
source("load_data.R")

# Plotting ####
## |- Function to plot the results by virus ####

wimage <- 4
himage <- 7
resimage <- 300
mai <- c(0.5, 0.1, 0.3, 1.2)

plotByVir <- function(vir, resdf, onlyMammals = TRUE, onlyBirds = FALSE, colSignif = colS, tit = "", latin = FALSE, png = FALSE, pdf = TRUE, filenamerq = ""){
  sub <- resdf[resdf$virus == vir, ]
  if(onlyMammals){
    sub <- sub[is.element(sub$animal, mammals), ]
  }
  if(onlyBirds){
    sub <- sub[is.element(sub$animal, birds), ]
  }
  sub <- sub[order(sub$commonName), ]
  
  if(png | pdf){
    if(onlyMammals) suffix <- "_mammals"
    if(onlyBirds) suffix <- "_birds"
    if(!onlyMammals & !onlyBirds) suffix <- "_all"
    
    fname <- paste0("../figs/correlations_", gsub(" ", "_", vir), suffix, "_", deparse(substitute(resdf)), filenamerq, "")
    if(png){
      png(width = wimage + mai[2] + mai[4], height = nrow(sub)*0.175 + mai[1] + mai[3], file = paste0(fname, ".png"), units = "in", res = resimage)
    }else{
      pdf(width = wimage + mai[2] + mai[4], height = nrow(sub)*0.175 + mai[1] + mai[3], file = paste0(fname, ".pdf"))
    }
  }
  
  par(las = 1, mai = mai)
  yy <- rev(seq_along(sub$commonName))
  
  par(mgp = c(1.25, 0.25, 0), tck = -0.02)
  plot(sub$cors, yy, xlim = c(-1, 1), 
       ylab = "", xlab = "correlation", axes = FALSE, type = "n")
  lines(c(0, 0), range(yy) + c(-1, 1), lty = 2)
  axis(1, lwd = 0, lwd.ticks = 1)
  par(xpd = TRUE)
  
  # Define colors
  cols <- rep(colSignif, length(yy))
  cols[sub$pval > 0.05] <- colNS # Non significant
  
  ii <- which(!is.na(sub$cors))
  if(latin){
    txt <- gsub("\\.", " ", sub$animal)
    ft <- 3
  }else{
    txt <- dicoSp[sub$animal]
    ft <- 1
  }
  text(sub$corsCI2[ii], yy[ii], labels = paste(" ", (txt)[ii]), adj = 0, cex = 0.9, 
       col = cols[ii], font = ft)  
  par(xpd = FALSE)
  segments(x0 = sub$corsCI1, x1 = sub$corsCI2, y0 = yy, y1 = yy, lwd = 2, 
           col = cols)
#  segments(x0 = sub$bootCI1, x1 = sub$bootCI2, y0 = yy, y1 = yy, lwd = 2, lty = 3, 
#           col = "red") # (to compare CIs)
  mtext(paste0(vir, tit), side = 3)
  
  points(sub$cors, yy, pch = 16, 
         col = cols)
  
  if(png | pdf){
    dev.off()
    system(paste0("open ", fname, "*"))
  }
  }; plotByVir("SARS-CoV-2", ACC1.all, tit = ", all dates, Mammals, ACC data", pdf = TRUE)

plotByVir("SARS-CoV-2", ACC1.Jan01, tit = ", 01 Jan 2020, Mammals, ACC data")
plotByVir("SARS-CoV-2", ACC1.Jan12, tit = ", 12 Jan 2020, Mammals, ACC data")
plotByVir("SARS-CoV-2", ACC1.Jan0112, tit = ", 01-12 Jan 2020, Mammals, ACC data")

plotByVir("SARS-CoV-2", ACC1.all, tit = ", all dates, ACC data", pdf = TRUE, onlyMammals = FALSE)
plotByVir("SARS-CoV-2", ACC1.Jan01, tit = ", 01 Jan 2020, ACC data", onlyMammals = FALSE)
plotByVir("SARS-CoV-2", ACC1.Jan12, tit = ", 12 Jan 2020, ACC data", onlyMammals = FALSE)
plotByVir("SARS-CoV-2", ACC1.Jan0112, tit = ", 01-12 Jan 2020, ACC data", onlyMammals = FALSE)

plotByVir("SARS-CoV-2", JB1.Jan01, onlyMammals = FALSE, latin = TRUE, tit = ", 01 Jan 2020, Bloom's data")
plotByVir("SARS-CoV-2", JB1.Jan12, onlyMammals = FALSE, latin = TRUE, tit = ", 12 Jan 2020, Bloom's data")
plotByVir("SARS-CoV-2", JB1.Jan0112, onlyMammals = FALSE, latin = TRUE, tit = ", 01 and 12 Jan 2020, Bloom's data")
plotByVir("SARS-CoV-2", JB1.all, onlyMammals = FALSE, latin = TRUE, tit = ", all dates, Bloom's data")

plotByVir("Influenza A virus H3N2", ACC1.all, tit = ", all dates")

plotByVir("Raccoon dog amdovirus", ACC1.all)
plotByVir("Civet kobuvirus", ACC1.all)
plotByVir("Bamboo rat coronavirus", ACC1.all)

plotByVir("Raccoon dog amdovirus", ACC1.Jan12, tit = ", 12 Jan 2020")
plotByVir("Civet kobuvirus", ACC1.Jan12, tit = ", 12 Jan 2020")
plotByVir("Bamboo rat coronavirus", ACC1.Jan12, tit = ", 12 Jan 2020")

plotByVir("Civet astrovirus", ACC1.all)
plotByVir("Civet astrovirus", ACC1.Jan12, tit = ", 12 Jan 2020")


