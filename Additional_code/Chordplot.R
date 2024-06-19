library(readxl)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(tidyr)
library(circlize)

#table with miRDB score, miRNA, target mRNA and tartet's function
data_corr_file_path <- "Funciones_relacion_miRNAs_metabolitos.xlsx"

data <- read_xlsx(data_corr_file_path, col_names = TRUE) %>%
  filter(Score >= 80) 

colnames(data) <- c("Score", "miRNA", "mRNA", "Description")

to_chord = select(data, miRNA, mRNA)

pdf("chord_diagram-miRDB.pdf", 10, 10)

par(mar = c(2, 2, 2, 2))
chordDiagram(to_chord, 
             grid.col = c(brewer.pal(11, "Set3")[-10], rep("#B3B2AE", nrow(to_chord))), 
             annotationTrack = c("grid"), 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(to_chord))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.7)
}, bg.border = NA) # here set bg.border to NA is important

dev.off()

