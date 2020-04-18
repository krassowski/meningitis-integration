if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", quiet=T)

# TODO:
# should be installed, but fails on Travis
# install.packages('Hmisc') # DESeq2 dependency
BiocManager::install("DESeq2", quiet=T)
BiocManager::install("edgeR")

install.packages('devtools', quiet=T)
install.packages('ggplot2', quiet=T)
install.packages('pheatmap', quiet=T)
install.packages('qqplotr', quiet=T)
install.packages('cowplot', quiet=T)
install.packages('ggthemes', quiet=T)

BiocManager::install('limma', quiet=T)
install.packages('statmod', quiet=T)

install.packages('parallel', quiet=T)
#install.packages('pbmcapply')

install.packages('import', quiet=T)

# pathological is required for readat and was removed from CRAN in the meantime
devtools::install_url('https://cran.r-project.org/src/contrib/Archive/pathological/pathological_0.1-2.tar.gz', quiet=T)
BiocManager::install("readat", quiet=T)

install.packages('gprofiler2', quiet=T)

install.packages('corrplot', quiet=T)
BiocManager::install("IHW", quiet=T)

# Bioc version is too old, 2.1.0 required
# BiocManager::install("ComplexHeatmap")
devtools::install_github("jokergoo/ComplexHeatmap@86907a77bfbb3e9325fcc677604f90be33ab574c")

install.packages("pvclust", quiet=T)
# TODO: remove the occurrences
#install.packages("ggstatsplot", quiet=T)

BiocManager::install("mixOmics", quiet=T)
BiocManager::install("ropls", quiet=T)

BiocManager::install("vsn", quiet=T)
BiocManager::install("TCC", quiet=T)

BiocManager::install("pROC", quiet=T)
BiocManager::install("ABSSeq")
# install.packages('ggforce', quiet=T) - maybe use it for graphics later, right now not needed
install.packages("cvAUC", quiet=T)
install.packages('hdi', quiet=T)

install.packages("ggnetwork", quiet=T)
install.packages('SNFtool', quiet=T)
install.packages('survminer')

install.packages('reshape')

devtools::install_github("krassowski/complex-upset")
install.packages('ggbeeswarm')
install.packages('imager')
install.packages('OmicsPLS')
