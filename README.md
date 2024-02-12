# Breast Sensitivity Signature Collection (SSc breast)
The **Breast Sensitivity Signature Collection (SSc breast)** contains 1,372 gene signatures that reflect the **transcriptional differences between sensitive and resistant breast cancer cell lines before drug treatment.** The SSc breast collection has been used in the publication "Spatial Transcriptomics in Breast Cancer Reveals Tumour Microenvironment-Driven Drug Responses and Clonal Therapeutic Heterogeneity" and it is available to download on Zenodo (DOI: 10.5281/zenodo.10638906).

This is a fork of [cnio-bu/drug_susceptibility_collection](https://github.com/cnio-bu/drug_susceptibility_collection) repository, which is intended to compute a pan-cancer drug Sensitivity Signature Collection (SSc) for [Beyondcell](https://github.com/cnio-bu/beyondcell) [1].

## Description
In order to obtain these signatures, we performed a differential expression analysis against the area under the curve (AUC) with limma v3.54.0 [2] for all compounds tested in at least 10 different breast cancer cell lines. We selected the top 250 up- and down-regulated genes in sensitive versus resistant cancer cell lines, ranked by the t-statistic, to create bidirectional gene signatures of 500 genes each. The AUC was used to measure drug response because, contrary to IC50, it can always be estimated without extrapolation from the dose-response curve and has shown more accuracy in predicting drug response [3].

Expression and drug response data were retrieved from three independent pharmacogenomics assays: the Cancer Therapeutics Response Portal (CTRP) v2 [4–6], the Genomics of Drug Sensitivity in Cancer (GDSC) v2 [7–9] and the PRISM [10,11] repurposing compendium through the DepMap portal v22Q4 [12]. As these sources are independent, several signatures refer to the same compound. Consequently, the 1,372 transcriptomic signatures that form the SSc breast reflect the predicted response to >1,200 drugs.

To verify that the cancer cell lines included in the signature generation analyses were, in fact, representative of breast cancer patients, we used the corrected lineage reported in the Celligner project [13], which provides a framework to align cancer cell lines to human tumours from large cohorts of human patients such as the Cancer Genome Atlas (TCGA).

**References**

1. Fustero-Torre C, Jiménez-Santos MJ, García-Martín S, Carretero-Puche C, García-Jimeno L, Ivanchuk V, et al. [Beyondcell: targeting cancer therapeutic heterogeneity in single-cell RNA-seq data.](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-01001-x) *Genome Med.* 2021;**13**:187.
2. Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, et al. [limma powers differential expression analyses for RNA-sequencing and microarray studies.](https://academic.oup.com/nar/article/43/7/e47/2414268?login=false) *Nucleic Acids Res.* 2015;**43**:e47.
3. Jang IS, Neto EC, Guinney J, Friend SH, Margolin AA. [Systematic assessment of analytical methods for drug sensitivity prediction from cancer cell line data.](https://www.worldscientific.com/doi/abs/10.1142/9789814583220_0007) *Pac Symp Biocomput.* 2014;**63**–74.
4. Basu A, Bodycombe NE, Cheah JH, Price EV, Liu K, Schaefer GI, et al. [An interactive resource to identify cancer genetic and lineage dependencies targeted by small molecules.](https://www.cell.com/cell/fulltext/S0092-8674(13)00960-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867413009604%3Fshowall%3Dtrue) *Cell.* 2013;**154**:1151–61.
5. Seashore-Ludlow B, Rees MG, Cheah JH, Cokol M, Price EV, Coletti ME, et al. [Harnessing Connectivity in a Large-Scale Small-Molecule Sensitivity Dataset.](https://aacrjournals.org/cancerdiscovery/article/5/11/1210/4735/Harnessing-Connectivity-in-a-Large-Scale-Small) *Cancer Discov.* 2015;**5**:1210–23.
6. Rees MG, Seashore-Ludlow B, Cheah JH, Adams DJ, Price EV, Gill S, et al. [Correlating chemical sensitivity and basal gene expression reveals mechanism of action.](https://www.nature.com/articles/nchembio.1986) *Nat Chem Biol.* 2016;**12**:109–16.
7. Garnett MJ, Edelman EJ, Heidorn SJ, Greenman CD, Dastur A, Lau KW, et al. [Systematic identification of genomic markers of drug sensitivity in cancer cells.](https://www.nature.com/articles/nature11005) *Nature.* 2012;**483**:570–5.
8. Yang W, Soares J, Greninger P, Edelman EJ, Lightfoot H, Forbes S, et al. [Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells.](https://academic.oup.com/nar/article/41/D1/D955/1059448?login=false) *Nucleic Acids Res.* 2013;**41**:D955–61.
9. Iorio F, Knijnenburg TA, Vis DJ, Bignell GR, Menden MP, Schubert M, et al. [A Landscape of Pharmacogenomic Interactions in Cancer.](https://www.cell.com/cell/fulltext/S0092-8674(16)30746-2?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867416307462%3Fshowall%3Dtrue) *Cell.* 2016;**166**:740–54.
10. Yu C, Mannan AM, Yvone GM, Ross KN, Zhang Y-L, Marton MA, et al. [High-throughput identification of genotype-specific cancer vulnerabilities in mixtures of barcoded tumor cell lines.](https://www.nature.com/articles/nbt.3460) *Nat Biotechnol.* 2016;*34*:419–23.
11. Corsello SM, Nagari RT, Spangler RD, Rossen J, Kocak M, Bryan JG, et al. [Discovering the anti-cancer potential of non-oncology drugs by systematic viability profiling.](https://www.nature.com/articles/s43018-019-0018-6) *Nat Cancer.* 2020;**1**:235–48.
12. Ghandi M, Huang FW, Jané-Valbuena J, Kryukov GV, Lo CC, McDonald ER 3rd, et al. [Next-generation characterization of the Cancer Cell Line Encyclopedia.](https://www.nature.com/articles/s41586-019-1186-3) *Nature.* 2019;**569**:503–8.
13. Warren A, Chen Y, Jones A, Shibue T, Hahn WC, Boehm JS, et al. [Global computational alignment of tumor and cell line transcriptional profiles.](https://www.nature.com/articles/s41467-020-20294-x) *Nat Commun.* 2021;**12**:22.

## Installation
The code in this repository is distributed as a [snakemake](https://snakemake.readthedocs.io/en/stable/#) pipeline. It extensively uses Snakemake's integration with the [conda](https://docs.conda.io/en/latest/) package manager to take care of software requirements and dependencies automatically.

First, install conda by following the [installation instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Next, use the git clone command to create a local copy:

```
git clone -b breast https://github.com/cnio-bu/drug_susceptibility_collection
```

And create a snakemake environment:

```
conda env create -f envs/snakemake.yaml
```

## How to run

### Set up
You need to modify the configuration files:

* **`config.yaml`:** Contains all pipeline parameters. Please specify the location of the output and log files.
* **`datasets.csv`:** Please indicate the location of the datasets to be analysed.

All required datasets are public and can be downloaded from the indicated sites:

| Dataset                        | File                                                           | Download from                                                                                                 |
|--------------------------------|----------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------|
| raw_ccle_reads                 | OmicsExpressionGenesExpectedCountProfile.csv                   | [DepMap](https://depmap.org/portal/download/all/) (v22Q4)                                                     |
| celligner_data                 | Celligner_info.csv                                             | [figshare](https://figshare.com/articles/dataset/Celligner_data/11965269)                                     |
| ccle_sample_info               | sample_info.csv                                                | [DepMap](https://depmap.org/portal/download/all/) (v22Q2)                                                     |
| prism_response_curves          | secondary-screen-dose-response-curve-parameters.csv            | [DepMap](https://depmap.org/portal/download/all/) (v19Q4)                                                     |
| prism_treatment_info           | secondary-screen-replicate-collapsed-treatment-info.csv        | [DepMap](https://depmap.org/portal/download/all/) (v19Q4)                                                     |
| gdsc_response_curves           | GDSC2_fitted_dose_response_*.xlsx                              | [GDSC](https://www.cancerrxgene.org/downloads/bulk_download) (Accessed 25 February 2020)                      |
| gdsc_compound_meta             | gdsc_compound_meta.csv                                         | [GDSC](https://www.cancerrxgene.org/compounds) (Accessed 25 February 2020)                                                                                                           |
| crispr_gene_dependency_chronos | CRISPRGeneDependency.csv                                       | [DepMap](https://depmap.org/portal/download/all/) (v22Q4)                                                     |
| ctrp_response_curves           | CTRPv2.0_2015_ctd2_ExpandedDataset/v20.data.curves_post_qc.txt | [DepMap](https://depmap.org/portal/download/all/) (vCTRP CTD^2)                                               |
| ctrp_compound_meta             | CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt   | [DepMap](https://depmap.org/portal/download/all/) (vCTRP CTD^2)                                               |
| ctrp_cell_meta                 | CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_cell_line.txt  | [DepMap](https://depmap.org/portal/download/all/) (vCTRP CTD^2)                                               |
| ctrp_experiment_meta           | CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_experiment.txt | [DepMap](https://depmap.org/portal/download/all/) (vCTRP CTD^2)                                               |
| ccle_default_line              | OmicsDefaultModelProfiles.csv                                  | [DepMap](https://depmap.org/portal/download/all/) (v22Q4)                                                     |
| hgnc_protein_coding            | hgnc_gene_with_protein_product.tsv                             | [HGNC](https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt) (Accessed 22 March 2023) |
| hallmarks                      | h.all.v2023.1.Hs.symbols.gmt                                   | [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp)                                                   |

### Execution
Once the pipeline is configured, the user just needs to enter these commands:

```
conda activate snakemake
snakemake --use-conda -j 200
conda deactivate
```

The mandatory arguments are:

* **--use-conda**: To install and use the conda environments.
* **-j**: Number of threads/jobs provided to snakemake.

### Output
After all the jobs have been completed, the SSc breast collection can be found in `resultspath/drug_signatures_classic.gmt`.

## Authors

* Santiago García-Martín
* María José Jiménez-Santos

<!-- ## Citation -->

## Support
If you have any questions, feel free to submit an [issue](https://github.com/cnio-bu/SSc-breast/issues).
