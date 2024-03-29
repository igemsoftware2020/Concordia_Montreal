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

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")



BiocManager::install()
BiocManager::install("GEOquery")
BiocManager::install("Biobase")
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("arrayQualityMetrics")


install.packages("devtools")
library(devtools)
devtools::install_github("aryoda/R_enumerations@v0.3.0-beta")
