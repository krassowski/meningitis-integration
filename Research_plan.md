# Research plan

Here is the original proposed research plan from the project description as a reference:

> 1. Explore and QC the multi‐omics and clinical data using univariate and multivariate statistical tools. Initially each data set will be explored alone.
Data will be checked for outliers and expected characteristics. Both univariate (e.g. linear models) and multivariate (e.g. PCA, PLS) approaches will be used.
> 2. Integrate the multi‐omics data using multivariate techniques such as Orthogonal Partial Least Squares regression. Both conventional and O2PLS will
be explored to integrate the omics data sets. Further QC analysis on the joint data will be conducted. Common sources of variance will be identified.
> 3. The integrated multi‐omics data will be combined with the clinical parameters in a second stage data fusion, again using PLS and O2PLS techniques.
We will investigate the ability of the individual and integrated data sets to predict clinical outcomes. Predictive models will be interrogated to elucidate
proteins & transcripts linked to the disease diagnosis.
> 4. If time allows we will investigate the use of OnPLS, a multiblock extension of O2PLS which allows to
integrate data from more than two blocks. Results from this work will form the basis of future large‐scale clinical studies to further advance the diagnosis
of infectious meningitides and understanding its pathogenesis. The ultimate aim is to derive signatures that can be validated in other trials and applicable
for future clinical use.

# TODO (DRAF!)

## Data exploration:

### Molecular dataset
- [x] Basic exploration: extraction, reformatting, sanity checks ([notebook](data_exploration/Molecular_data_extraction.ipynb))
- [ ] Quality control (and maybe normalization) quantile-quantile plots (0.5 day) 
- [ ] Unsupervised exploration (PCA, hierarchical clustering) and correlations (1.5 days)

### Clinical dataset
- [x] Basic exploration: extraction, reformatting, sanity checks: ([notebook](data_exploration/Clinical_data_first_look.ipynb))
- [ ] Derived variables (e.g. age, survival) and correlations (e.g. CD4) (0.25 day)
  - which clinical variables can be used as covariates, which should be used (or, can be considered as) outcomes
  - descriptive statistics to summarize characteristics of the studied cohort

## Preliminary analyses

These analyses are proposed to be carried out in parallel to the PLS study.
I was aiming to chose techniques which are either established as standard in analysis of expression data, and/or easy to apply (so that the main focus is still on the PLS).

### Single omic
- [ ] Differential expression (DE) (1 day)
  - [ ] for RNA-seq data (performed by Dr Rachel, though I will re-analyse and plot the data)
  - [ ] for protein data
- [ ] Gene set enrichment analysis (GSEA) (1 day)
- [ ] Multiple regressions/ANOVA using top DEGs (1 day)

I do not have a strong intuition on use of multiple regressions in this setting:
 - while the simpler approach would be to apply the model to all the variables and correct for multiple testing,
 - I could also start with the specialised DE tools which would first highlight the genes/proteins which
 are differentially expressed between groups.

### Multiple omics
- [ ] Joint NMF clustering: a simple approach for unsupervised clustering (1 day)
- [ ] Joint pathways analysis (late integration): using combined GSEA results (ActivePathways + Cytoscape) (2 days)

## The actual study

### Preparation
- [ ] Review of available PLS implementations and recent developments (14 days)
  - [ ] Layman summary of PLS method with a graphical representation (3 days)
- [ ] Reading up on scaling, de-noising, transformations and normalization (4 days)
  - [ ] Application to the analyzed omics data (2-3 days)

### Comparison of patient groups (by condition)
- [ ] single-omic (O)PLS-DA discrimination analysis (X = [protein | RNA], Y=patient group)
- [ ] multi-omics (O)PLS-DA discrimination analysis (X = [protein, RNA], Y=patient group)

### Regression on clinical outcomes
Due to missing data, we may be a need to perform a few analyses:
  - comprehensive "high-level" analysis with the clinical variables available for all patients
    (e.g. survival, age, CSF color, CSF cloudiness, HIV status),
  - subgroup-specific analyses, including the variables which were measured only for these patients:
     - HIV positive patients, additional variables: CD4, ARV, UnARV
     - patients with history of TB, additional variables: PrevTBForm, PrevTBTreat, OnTBTreat,
     (potentially also: DateTBTreat and DateTBTreatStop, though these are not useful on their own as are mutually exclusive; however when combined with other dates it might of some use)

- [ ] (O2)-PLS: single-omic regression (X=[protein | RNA], Y=clinical outcomes)
- [ ] (O2)-PLS: multi-omics regression (X=[protein, RNA], Y=clinical outcomes)

### Multi-omics exploration?
- [ ] O2-PLS (X=protein, Y=RNA), O2 because the relation in either way is equally interesting.
[(Bylesjo et al, 2007)](https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1365-313X.2007.03293.x) has an example of such analysis (supervised by J. Trygg).

## Comparison to other multi-view approaches
- [ ] Compare to the joint NMF results
- [ ] Try applying alternative methods, two or more of:
  - PARADIGM, iCluster - promising approaches
  - CCA - similar to PLS
  - COCA - very basic approach
Beside comparing the 

## References:
- [Multi-omic and multi-view clustering algorithms: review and cancer benchmark](https://academic.oup.com/nar/article/46/20/10546/5123392)


## Side notes:
- As so far this is essentially a gene expression study, we could look at [R svapls package](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-236)
