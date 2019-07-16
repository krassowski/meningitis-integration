install.packages('ggplot2', quiet=T)
install.packages('pheatmap', quiet=T)
install.packages('qqplotr', quiet=T)
install.packages('cowplot', quiet=T)
install.packages('ggthemes', quiet=T)

install.packages('limma', quiet=T)
install.packages('statmod', quiet=T)

install.packages('parallel', quiet=T)
#install.packages('pbmcapply')

install.packages('import', quiet=T)

source("https://bioconductor.org/biocLite.R")
biocLite("readat", quiet=T)
install.packages('gprofiler2', quiet=T)


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
install.packages('corrplot', quiet=T)
BiocManager::install("IHW")

install.packages('devtools')
Sys.unsetenv("GITHUB_PAT")

# Bioc version is too old, 2.1.0 required
# BiocManager::install("ComplexHeatmap")
devtools::install_github("jokergoo/ComplexHeatmap")

install.packages("pvclust", quiet=T)
install.packages("ggstatsplot", quiet=T)

BiocManager::install("mixOmics")
BiocManager::install("ropls")
