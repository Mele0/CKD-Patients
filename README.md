# Proteomics-Based Clustering for Endotype Discovery in Chronic Kidney Disease

**Alex Melendez, Chiao Yuen Lim, Joy Cheung, Parivrudh Sharma, Nicole Lim**  
Imperial College London, MSc Health Data Analytics and Machine Learning (2024–2025)

## Motivation

Chronic kidney disease (CKD) is a growing global health crisis, affecting more than 850 million individuals worldwide and ranking among the few non-communicable diseases whose mortality burden has risen steadily over the past three decades. Reliance on traditional markers such as estimated glomerular filtration rate and albuminuria provides limited insight into the complex molecular processes— inflammation, fibrosis, oxidative stress—that drive disease progression. A biologically grounded classification could enable earlier identification of high-risk patients, more precise risk stratification, and the development of targeted therapies.

## Overview

This project leverages high-dimensional plasma proteomics from the UK Biobank to uncover molecularly distinct CKD endotypes. We reduce the ~3 000 protein features to a two-dimensional embedding using Uniform Manifold Approximation and Projection (UMAP) and identify three stable clusters via Gaussian Mixture Models. Each endotype is characterised by its unique proteomic signature using tuned Random Forest classifiers and stability-enhanced LASSO, then linked to clinical outcomes—progression to end-stage renal disease (ESRD) and all-cause mortality—through Cox proportional hazards models. Pathway enrichment reveals divergent biological processes in each cluster, suggesting mechanistic subtypes that could inform precision nephrology.

## Data

Proteomic measurements were obtained from the UK Biobank, a population cohort of over half a million volunteers aged 40–69 at recruitment (2006–2010), which combines deep phenotyping, health records, and multimodal “omics.” We analysed 2 022 incident CKD cases drawn from a subset of 54 219 participants profiled with the Olink Explore 3072 proximity extension assay. Clinical covariates, comorbidities and laboratory data were extracted from baseline questionnaires, Hospital Episode Statistics and national death registries, with follow-up through late 2022. Raw proteomic NPX values underwent mixed-model adjustment for plate and processing date, removal of assays with > 50 percent missingness, and imputation of left-censored values via QRILC before standardisation to zero mean and unit variance.

## Methods

After quality control and standardisation, we applied UMAP with a cosine distance metric (30 nearest neighbours, minimum distance 0) to preserve local structure in the proteomic space. A grid search over Gaussian Mixture Models—exploring two to ten components and all covariance types—identified a three-component spherical model as optimal by silhouette score and Bayesian Information Criterion. Cluster-predictive proteins were then uncovered by fitting Random Forest classifiers (binary “cluster vs. non-cluster”) tuned for tree count, node size and feature sampling, ranking candidates by mean decrease in accuracy. Socio-demographic and clinical feature associations were evaluated through stability-enhanced LASSO logistic regression, ensuring robust variable selection across 100 subsamples. Differential expression and pathway enrichment employed limma contrasts against Reactome, while Cox models incorporating age, sex, creatinine and cluster membership quantified hazards for ESRD and mortality.

## Key Findings

A high-risk “orange” endotype exhibited a hazard ratio of 2.78 for ESRD (p 1.2 × 10⁻⁵) and 2.01 for mortality (p 1.4 × 10⁻¹³), with upregulated TNF-mediated inflammatory pathways. A moderate-risk “green” endotype was defined by Rho-GTPase signalling pathways linked to cytoskeletal maintenance. The lower-risk “blue” endotype showed downregulation of TP53 and SCF–KIT pathways. These molecular subgroups map to distinct clinical phenotypes and prognoses, demonstrating that proteomic endotyping can stratify CKD patients more precisely than conventional staging.

## Usage

To reproduce the analysis, clone this repository and install dependencies via the provided Conda environment. Launch the main Jupyter notebook to execute data loading, preprocessing, clustering, feature selection, pathway enrichment and survival modelling in sequence. Figures and tables will be generated automatically in the `Results` folder.

## Requirements

All Python and R dependencies are specified in `config/env.yml`. This includes `numpy`, `pandas`, `scikit-learn`, `umap-learn`, `lifelines`, `matplotlib`, `seaborn`, `limma`, `clusterProfiler` and `ReactomePA`. Create and activate the environment with  
```bash
conda env create -f config/env.yml  
conda activate ckd-proteomics  
