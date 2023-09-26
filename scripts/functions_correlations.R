# Define functions to compute correlations ####

# Compute correlation and CIs with bootstrap
# (Just for testing, but CIs are similar to the ones we obtain with cor.test)
corBoot <- function(x, y, method, nboot){
  # x, y: vectors on which correlation will be computed
  # method: "pearson" or "spearman" or "theil-sen"
  # nboot: number of resampling
  corResample <- function(x, y, method = method){
    out <- NA
    while(is.na(out)){ # Re-draw if leads to sd = 0
      # Draw positions with replacement
      ii <- sample(seq_along(x), size = length(x), replace = TRUE)
      if(method != "theil-sen"){
        # Compute correlation on the new vectors
        # Sometimes sd = 0; remove warnings
        out <- suppressWarnings(cor(x[ii], y[ii], method = method))
      }else{
        if((var(x[ii]) == 0) | var(y[ii]) == 0){
          out <- NA
        }else{
          out <- theilSen_r(x[ii], y[ii])
        }
      }
    }
    out
  }
  vout <- replicate(nboot, corResample(x, y, method = method))
  # Quantiles, ignoring NAs when sd = 0
  out <- quantile(vout, c(0.025, 0.975), na.rm = TRUE)
  out
}


# Theil-Sen correlation
# https://towardsdatascience.com/a-correlation-measure-based-on-theil-sen-regression-31b8b9ed64f1
theilSen_r <- function(x, y){
  m1 <- (deming::theilsen(y ~ x))$coefficients[2]
  m2 <- (deming::theilsen(x ~ y))$coefficients[2]
  out <- 0
  if (sign(m1 * m2) >=0){
    out <- sign(m1) * sqrt(m1 * m2)
  }
  out
}


# Compute correlation, with different options on the vectors
computeCor <- function(x, y, denom = NA, remove0.x = FALSE, remove0.y = FALSE, log = FALSE, minVal = NA, prop = TRUE, method, plot = FALSE, bootstrap = TRUE, nboot = 1000){
  # x, y vectors on which correlation will be computed
  # denom type of denominator
  # remove0: whether to include zeros or not
  # log: whether to log-transform
  # minVal: minimum value
  # prop: "none" if raw reads, "Mammals" if by mammals, "Chordates" if by chordates, "all" if by all reads
  # method: "pearson" or "spearman"
  # plot: whether to plot the vectors
  # bootstrap: whether to compute CIs by bootstrap
  # nboot: number of bootstrap replicates
  
  # Rename variables to modify them
  X <- x
  Y <- y
  
  # Remove zeros, in pair
  if(remove0.x){
    X <- X[X > 0]
    Y <- Y[X > 0]
  }
  
  if(remove0.y){
    X <- X[Y > 0]
    Y <- Y[Y > 0]
  }
  
  # Compute proportion
  if(prop){
    X <- X / denom
    Y <- Y / denom
  }
  
  # Log-transform
  if(log){
    #   Change zeros if needed
    if(!remove0.x){
      if(is.na(minVal)) stop("Need to enter minVal")
      X[X == 0] <- minVal
    }
    if(!remove0.y){
      if(is.na(minVal)) stop("Need to enter minVal")
      Y[Y == 0] <- minVal
    }
    # Log
    X <- log(X)
    Y <- log(Y)
  }
  
  if(plot){
    plot(X, Y, xlab = "", ylab = "")
  }
  
  # Compute correlation depending on the method
  if(method == "theil-sen"){
    out <- theilSen_r(X, Y)
  }else{
    # Pearson or Spearman
    out <- cor.test(X, Y, method = method, alternative = "t")
  }
  # If bootstrap
  if(bootstrap){
    bci <- corBoot(X, Y, method = method, nboot = nboot)
    out <- c(out, bootCI = bci)
  }
  out
}

# Compute denominators
data$denom.all <- data$Read_pairs_after_trimming
data$denom.chordates <- apply(data[, chordates], 1, sum)
data$denom.mammals <- apply(data[, mammals], 1, sum)


# Wrapper
wrapperCor <- function(dataset, virus, animal, remove0.x = FALSE, remove0.y = FALSE, log = TRUE, denomType = "all", prop = TRUE, method = "pearson", plot = FALSE, bootstrap = TRUE, nboot = 1000, printOutput = FALSE){
  # How we extract the vectors depends on the data source
  if(is.element("SARS2.reads", names(dataset))){
    # Source is ACC
    if(virus == "SARS-CoV-2"){
      dat <- dataset
      x <- dat$SARS2.reads
    }else{
      # Get virus results from the virus dataframe, which is treated as global variable
      # Make sure to only include samples that are in dataset
      tmp <- viral[which(viral$Virus_shorter == virus & is.element(viral$Sample.ID, dataset$Sample.ID)), ]
      dat <- merge(dataset, tmp[, c("Sample.ID", "read_count")], all = TRUE)
      x <- dat$read_count
      # Replace missing values by zeros
      x[is.na(x)] <- 0
    }
    y <- dat[, animal]
    
    denom <- 1
    if(denomType == "all"){
      denom <- dat$Read_pairs_after_trimming
    }
    if(denomType == "Chordates"){
      denom <- apply(dat[, chordates], 1, sum)
    }
    if(denomType == "Mammals"){
      denom <- apply(dat[, mammals], 1, sum)
    }
    
  }else{
    # Source is JB
    if(virus != "SARS-CoV-2") stop("Virus error! Only SARS-CoV-2 in JB's data")
    
    # Construct dataset
    dat <- dataset[which(dataset$species == "SARS-CoV-2"), ]
    names(dat)[names(dat) == "aligned_reads"] <- "SARS2.reads"
    
    # Add information about the species
    tmp <- dataset[which(dataset$species == gsub("\\.", " ", animal)), c("Sample.ID", "aligned_reads")]
    names(tmp)[names(tmp) == "aligned_reads"] <- animal
    
    dat <- merge(dat, tmp, all.x = TRUE)
    rm(tmp)
    # Replace NAs
    dat[is.na(dat[, animal]), animal] <- 0
    
    # Compute denominator
    denom <- 1
    if(denomType == "all"){
      denom <- dat$preprocessed_reads
    }
    # And I am not bothering to do the other types of denominators again, 
    # because it honestly makes no sense to restrict the denominator like this, 
    # especially when having SC2 at the numerator.
    
    x <- dat$SARS2.reads
    y <- dat[, animal]
    
  }  
  
  # Define a minimum value to replace zeros when plotting on log scale
  if(length(denom) == 1){
    minVal <- 0.5
  }else{
    minVal <- 0.5*min(1 / denom)
  }
  cr <- computeCor(x = x, y = y, denom = denom, remove0.x = remove0.x, remove0.y = remove0.y, log = log, minVal = minVal, prop = prop, method = method, plot = plot, bootstrap = bootstrap, nboot = nboot)
  txt <- ""
  nsignif <- 2
  
  if(method != "theil-sen"){
    txt <- paste0("cor = $", signif(cr$estimate, nsignif), 
                  "$ (95\\% CI: $", signif(cr$conf.int[1], nsignif), 
                  "$, $", signif(cr$conf.int[2], nsignif), 
                  "$), $p$ = $", signif(cr$p.value, nsignif), "$")
  }else{
    if(!bootstrap){
      txt <- paste0("cor = $", signif(cr, nsignif), "$")
    }else{
      txt <- paste0("cor = $", signif(cr[1], nsignif), "$", 
                    "$ (95\\% CI, bootstrap: $", signif(cr[2], nsignif), 
                    "$, $", signif(cr[3], nsignif), "$)")
    }
  }
  
  if(printOutput){
    print(txt)
    
  }
  list(cr = cr, txt = txt) 
}


# Another wrapper to be used with lapply
ww <- function(i, dd, parms){
  wrapperCor(data = dd, virus = parms$virus[i], animal = parms$animal[i], 
             remove0.x = FALSE, remove0.y = FALSE, log = TRUE, 
             prop = TRUE, denomType = "all", 
             method = "pearson", plot = FALSE, bootstrap = FALSE)$cr
}