install.packages('ggplot2')
install.packages('pheatmap')
install.packages('qqplotr')
install.packages('cowplot')
install.packages('ggthemes')

install.packages('limma')
install.packages('statmod')

install.packages('parallel')
#install.packages('pbmcapply')

install.packages('import')

source("https://bioconductor.org/biocLite.R")
biocLite("readat")
install.packages('gprofiler2')


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
install.packages('corrplot')
BiocManager::install("IHW")

install.packages('devtools')
Sys.unsetenv("GITHUB_PAT")

# Bioc version is too old, 2.1.0 required
# BiocManager::install("ComplexHeatmap")
devtools::install_github("jokergoo/ComplexHeatmap")

install.packages("pvclust")
