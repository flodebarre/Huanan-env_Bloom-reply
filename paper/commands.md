# Convert the ms to Word
pandoc -s Bloom-reply_final.tex --bibliography=biblio_clean.bib  --citeproc  -o Bloom-reply_final.docx
# (The tables are then extracted from the docx file)

# Convert figures to eps
pdftops -eps fig1.pdf

# Convert table to Word
pandoc -s tableS1.tex -o tableS1.docx

