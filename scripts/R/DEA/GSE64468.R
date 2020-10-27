#License

# © Copyright 2020 iGEM Concordia, Maher Hassanain, Benjamin Clark, Hajar Mouddene, Grecia Orlano
# This file is part of AstroBio.
# AstroBio is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any later version.
# AstroBio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with AstroBio.  If not, see <https://www.gnu.org/licenses/>.

###GSE64468 YEAST MICROARRAY
#BENJAMIN CLARK

library(GEOquery)
library(Biobase)
library(limma)
library(arrayQualityMetrics)
source("microarray_functions.R")

list.gse <- getGEO("GSE64468", GSEMatrix =TRUE, AnnotGPL=TRUE)

gse <- list.gse[[1]]

#take a look into the dataset
#head(gse)
View(gse$title)
#dim(gse)

#here we format the feature names. fvarLabels belongs the biobase package and is used to extract features from ExpressionSet Objects
fvarLabels(gse) <- make.names(fvarLabels(gse))

#WT groups
control <- c(1,2,3)
treatment <- c(4,5,6)
wt <- de.analysis(microgravity_group = treatment, ground_group = control, gse = gse)
filtered.tT <- remove.controls(wt$TopTable)


wt.filename <- "datasets/GSE64468_Scer/GSE64468.csv"
write.table(filtered.tT$TopTable, wt.filename, row.names = FALSE, sep = ",")


metaname <- "datasets/GSE64468_Scer/GSE64468_meta"
extractMetaData(gse_groups = list(wt), filename = metaname, microgravity_type = M.TYPE$SPACEFLOWN, metaLabels = c(""), strain = "BY4742")
