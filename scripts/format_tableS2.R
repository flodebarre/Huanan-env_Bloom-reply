# Format Table S2 for the journal's requirements

# Load data
tS2 <- read.csv("../results/prop-chordates_human.csv", row.names = NULL)
# Round data for printing
tS2$prop.chordates <- round(tS2$prop.chordates, 4)
# Rename columns
names(tS2) <- c("Sample ID", "Lab code", "Address", "Sampling date", "Proportion human reads", "Number SARS-CoV-2 reads")

rownames(tS2) <- NULL


# Export the formatted table
if(!library(xtable, logical.return = TRUE, quietly = TRUE)){
  install.packages("xtable")
}

add.to.row <- list(pos = list(0), command = NULL)
command <- paste0("\\hline\n\\endhead\n",
                  "\\hline\n",
                  "\\multicolumn{", dim(tS2)[2] + 1, "}{l}", "{\\footnotesize Continued on next page}\n", "\\endfoot\n",
                  "\\endlastfoot\n")
add.to.row$command <- command

print(xtable::xtable(tS2[, c(1, 2, 5, 6)], type = "latex", digits = 4, caption = "Proportion of human reads among chordate reads, in samples with and without \\sct{} reads; data from \\citet{Bloom2023VE}. A machine-readable csv source is available on the projet's repository.", label = "tab:propHuman", hline.after = c(-1), add.to.row = add.to.row), file = "tmp.tex", include.rownames = FALSE, tabular.environment = "longtable", floating = FALSE)
# Remove table commands from the output
system("sed '/{table/d' tmp.tex > ../paper/tableS2.tex; rm tmp.tex")
