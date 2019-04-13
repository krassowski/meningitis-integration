## Protein levels data

**Important notes from meeting:**
- During the study there were 186 admissions of meningitis patients in the study hospital. It is likely that many of the severe patients were not enrolled in the study, as it was not the priority of the clinician to ask for consent (or it would be impossible).
  - this might impact survival estimates, i.e. we may be missing some of the very severe patients who died too early for the necessary consents and biofluids to be collected
- Missing protein data → to little CSF, they need to decide whether these will be used to CSF or RNAseq. This may be random-like decision.

Note for the write-up: the term "protein abundance" might be more appropriate.

### Total CSF Protein

See [Total_CSF_protein.ipynb](Total_CSF_protein.ipynb)

Introduction (after consultation with Dr Rachel)
- Total CSF Protein was measured independently of the SOMAScan
- the signal is assumed to be dominated by albumin

#### Results

- There is a high correlation (0.91, p=2.9e-32) of Total CSF protein and sum of protein levels as measured by SOMAScan
  - which reassures us about the coherence of protein measurements (even though SOMAScan measures less than 10% of proteins)
- Total CSF protein correlates with the disease status (as expected)
- Pairwise comparison highlights outliers (24, 239) which are likely a result of a technical issue
- Total CSF protein correlates with SOMAScan measurements regardless of the disease status, though exact $\rho$ varies (generally in range of 0.72-0.9, Spearman).
  - the highest correlation for TBM (0.9) might mean that SOMAScan is better suited for analysis of TBM as it can explain more of the change in total CSF protein level variation (or, there might a difference in the distribution of proteins driving the change; it would mean that for TBM there are more proteins with consistent, but small contributions to the total protein CSF level)
- Albumin is highly correlated with the total CSF protein (0.89, p=9e-30)
  - but according to SOMAScan measurements accounts only for 0.5-0.6% of the proteins (while literature tells us about 35-80%!)
  - there are 112 proteins with higher correlation to the Total CSF protein level
  - and several proteins accounting for 1.5% of the proteins measured by SOMAScan each (so x3 the value of albumin)
  - another correlation of CSF proteins: albumin with IgG is not observed as well
  - **so there must be something wrong with albumin measurements:**
    - outside of the dynamic range?
    - filtered out in preprocessing?
    - a flaw in SOMAScan?
    - or (much less likely) the methodology that lead to the 35-80% figure is not sound
- IGF-I and Kallikrein 7 (among others) negatively correlate with the Total CSF protein (thus indirectly with the disease status). The former was [previously reported](https://journals.lww.com/infectdis/Abstract/2012/03000/Cerebrospinal_Fluid_and_Serum_Levels_of_Growth.9.aspx?trendmd-shared=0) as specific to TBM (as compared to BM and controls).

## Differential levels of proteins

See [Differential_levels.ipynb](Differential_levels.ipynb).

## Regression


## Classification


## Survival analysis

- has to account for right-censoring.
- left censoring is not obvious: we do not know when the patient became ill, but we do know when they were admitted
  - for patients with previous history of TB (having an end date for previous TB treatment) I could infer the range of dates...
  - though this would be overcomplicated...
- there are protein data 15 out of 22 patients who deceased before the end of 6 month follow-up period

See [Survival.ipynb](Survival.ipynb)