# Proteomics-Based Clustering for Endotype Discovery in Chronic Kidney Disease

## Motivation

Chronic kidney disease (CKD) is a growing global health crisis, affecting more than 850 million individuals worldwide and ranking among the few non-communicable diseases whose mortality burden has risen steadily over the past three decades. Reliance on traditional markers such as estimated glomerular filtration rate and albuminuria provides limited insight into the complex molecular processes— inflammation, fibrosis, oxidative stress—that drive disease progression. A biologically grounded classification could enable earlier identification of high-risk patients, more precise risk stratification, and the development of targeted therapies.

## Overview

This project leverages high-dimensional plasma proteomics from the UK Biobank to uncover molecularly distinct CKD endotypes. We reduce the ~3,000 protein features to a two-dimensional embedding using Uniform Manifold Approximation and Projection (UMAP) and identify three stable clusters via Gaussian Mixture Models. Each endotype is characterised by its unique proteomic signature using tuned Random Forest classifiers and stability-enhanced LASSO, then linked to clinical outcomes—progression to end-stage renal disease (ESRD) and all-cause mortality—through Cox proportional hazards models. Pathway enrichment reveals divergent biological processes in each cluster, suggesting mechanistic subtypes that could inform precision nephrology.

## Data

Proteomic measurements were obtained from the UK Biobank, a population cohort of over half a million volunteers aged 40–69 at recruitment (2006–2010), which combines deep phenotyping, health records, and multimodal “omics.” We analysed 2,022 incident CKD cases drawn from a subset of 54,219 participants profiled with the Olink Explore 3072 proximity extension assay. Clinical covariates, comorbidities, and laboratory data were extracted from baseline questionnaires, Hospital Episode Statistics, and national death registries, with follow-up through late 2022. 

Proteomic processing comprised:  
1. Adjustment for technical covariates (plate, processing date) via linear mixed models.  
2. Exclusion of assays with greater than 50 % missing values.  
3. Imputation of left-censored values by QRILC.  
4. Standardisation to zero mean and unit variance.

Clinical and demographic covariates, comorbidities (Hospital Episode Statistics ICD-10 codes), laboratory measurements and outcomes (end-stage renal disease, mortality) were obtained from baseline assessments, Hospital Episode Statistics and national death registries with follow-up through late 2022.

## Methods

After quality control and standardisation, we applied UMAP with a cosine distance metric (30 nearest neighbours, minimum distance 0) to preserve local structure in the proteomic space. A grid search over Gaussian Mixture Models—exploring two to ten components and all covariance types—identified a three-component spherical model as optimal by silhouette score and Bayesian Information Criterion. Cluster-predictive proteins were then uncovered by fitting Random Forest classifiers (binary “cluster vs. non-cluster”) tuned for tree count, node size, and feature sampling, ranking candidates by mean decrease in accuracy. Socio-demographic and clinical feature associations were evaluated through stability-enhanced LASSO logistic regression, ensuring robust variable selection across 100 subsamples. Differential expression and pathway enrichment employed limma contrasts against Reactome, while Cox models incorporating age, sex, creatinine, and cluster membership quantified hazards for ESRD and mortality.

## Results
Three distinct proteomic endotypes—designated Blue, Green and Orange—were identified, each exhibiting unique molecular signatures, pathway enrichments and clinical prognoses.

### Table 1. Differentially expressed genes and enriched pathways by cluster  
Upregulated gene sets are indicated by ▲; downregulated sets by ▼.

<table>
  <thead>
    <tr>
      <th></th>
      <th>Green (N = 646)</th>
      <th>Blue (N = 773)</th>
      <th>Orange (N = 603)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Differentially expressed genes</th>
      <td>Up ▲ 68<br>Down ▼ 0</td>
      <td>Up ▲ 0<br>Down ▼ 98</td>
      <td>Up ▲ 31<br>Down ▼ 0</td>
    </tr>
    <tr>
      <th>Significantly enriched gene sets</th>
      <td>Up ▲ 82<br>Down ▼ 0</td>
      <td>Up ▲ 0<br>Down ▼ 68</td>
      <td>Up ▲ 5<br>Down ▼ 0</td>
    </tr>
    <tr>
      <th>Top uniquely enriched gene sets</th>
      <td>
        ▲ Signaling by Rho GTPases<br>
        ▲ Signaling by Rho GTPases, Miro GTPases and RHOBTB3<br>
        ▲ RHO GTPase Cycle
      </td>
      <td>
        ▼ Signaling by SCF–KIT<br>
        ▼ Regulation of TP53 Activity
      </td>
      <td>
        ▲ TNFs Bind Their Physiological Receptors<br>
        ▲ TNFR2 Non-Canonical NF-κB Pathway<br>
        ▲ EPH–Ephrin Signaling<br>
        ▲ EPH-ephrin Mediated Repulsion of Cells
      </td>
    </tr>
  </tbody>
</table>

- The Orange endotype was associated with the highest risk: hazard ratio for progression to end-stage renal disease (HR 2.78, p 1.2 × 10⁻⁵) and for all-cause mortality (HR 2.01, p 1.4 × 10⁻¹³). Pathway analysis revealed marked upregulation of TNF-superfamily inflammatory signalling.
- The Green endotype was characterised by Rho-GTPase signalling pathways implicated in cytoskeletal maintenance.
- The Blue endotype displayed downregulation of TP53 and SCF–KIT pathways associated with fibrosis and cellular repair.

## Key Findings

A high-risk “orange” endotype exhibited a hazard ratio of 2.78 for ESRD (p 1.2 × 10⁻⁵) and 2.01 for mortality (p 1.4 × 10⁻¹³), with upregulated TNF-mediated inflammatory pathways. A moderate-risk “green” endotype was defined by Rho-GTPase signalling pathways linked to cytoskeletal maintenance. The lower-risk “blue” endotype showed downregulation of TP53 and SCF–KIT pathways. These molecular subgroups map to distinct clinical phenotypes and prognoses, demonstrating that proteomic endotyping can stratify CKD patients more precisely than conventional staging.


## References
1. Kovesdy CP. Epidemiology of chronic kidney disease: an update 2022. *Kidney International Supplements*. 2022;12:7–11.  
2. Naghavi M, Ong KL, Aali A, Ababneh HS, Abate YH, Abbafati C. Global burden of 288 causes of death and life expectancy decomposition in 204 countries and territories and 811 subnational locations, 1990–2021: a systematic analysis for the Global Burden of Disease Study 2021. *The Lancet*. 2024;403:2100–32.  
3. Ekart R, Ferjuc A, Furman B, Gerjevič Š, Hojs SB. Chronic Kidney Disease Progression to End Stage Renal Disease: A Single Center Experience of the Role of the Underlying Kidney Disease. *Therapeutic Apheresis and Dialysis*. 2013;17:363–7.  
4. Webster AC, Nagler EV, Morton RL, Masson P. Chronic Kidney Disease. *The Lancet*. 2017;389:1238–52.  
5. Delrue C, Speeckaert MM. Decoding Kidney Pathophysiology: Omics-Driven Approaches in Precision Medicine. *Journal of Personalized Medicine*. 2024;14:1157.  
6. Chen TK, Surapaneni AL, Schmidt IM, Waikar SS, Coresh J, Liu H, et al. Proteomics and Incident Kidney Failure in Individuals With CKD: The African American Study of Kidney Disease and Hypertension and the Boston Kidney Biopsy Cohort. *Kidney Medicine*. 2024;6:100921.  
7. Dubin RF, Deo R, Ren Y, Wang J, Zheng Z, Shou H, et al. Proteomics of CKD progression in the chronic renal insufficiency cohort. *Nature Communications*. 2023;14:6340.  
8. Govender MA, Brandenburg JT, Fabian J, Ramsay M. The Use of ‘Omics for Diagnosing and Predicting Progression of Chronic Kidney Disease: A Scoping Review. *Frontiers in Genetics*. 2021;12:12. 
9. Rroji M, Spasovski G. Omics Studies in CKD: Diagnostic Opportunities and Therapeutic Potential. *PROTEOMICS*. 2024;11.  
10. Ling XC, Kuo KL. Oxidative stress in chronic kidney disease. *Renal Replacement Therapy*. 2018;4:53.  

## Usage

To reproduce the analysis, clone this repository and install dependencies via the provided Conda environment. Launch the main Jupyter notebook to execute data loading, preprocessing, clustering, feature selection, pathway enrichment and survival modelling in sequence. Figures and tables will be generated automatically in the `Results` folder.

## Requirements

All Python and R dependencies are specified in `config/env.yml`. This includes `numpy`, `pandas`, `scikit-learn`, `umap-learn`, `lifelines`, `matplotlib`, `seaborn`, `limma`, `clusterProfiler` and `ReactomePA`. Create and activate the environment with  
```bash
conda env create -f config/env.yml  
conda activate ckd-proteomics


