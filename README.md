# CoBRAA: Conductance-Based Reciprocal Association Analysis

**Analysis code for:** *Reciprocal effective conductance-based inference of protein-protein relationships delineates protein network hubs and macromolecular properties*  
**Authors:** Rajendra K C, Sihang Zhou, Ruiying Cheng, Marco Andres Penaherrera, Simon Lizarazo, and Kevin Van Bortle.




## Overview

This repository contains the R scripts underlying the computational analyses presented in our study. We developed **CoBRAA (Conductance-Based Reciprocal Association Analysis)** to examine protein–protein relationships in the human interactome using effective conductance and reciprocal association analysis. The workflow starts with experimentally reported human physical protein–protein interactions (PPIs) from BioGRID, calculates effective conductance between protein pairs, and derives reciprocal association scores. These scores are used to examine retained and newly inferred protein relationships, protein-network properties, and interactions between multi-subunit protein complexes.

The repository also includes analyses of RNA polymerases I, II, and III, comparisons with pull-down mass spectrometry, protein structural features, simulations, and figure visualizations. **Newly inferred associations are network-based predictions, not direct experimental evidence of physical binding.**
<br>
<br>
<img width="1419" height="597" alt="Screenshot 2026-10-10 at 11 53 53 AM" src="https://github.com/user-attachments/assets/b5e25b77-559b-4cec-affa-9d1af62c0e76" />
<br>
<br>


## Analysis workflow and scripts

The scripts below are the original research scripts provided with the study. Their numbering reflects the working analysis organization; they are **not** a fully automated, sequential pipeline. Some files contain exploratory, alternative, or obsolete code blocks.

| Analysis | Description | R script |
| --- | --- | --- |
| **1. PPI network and effective conductance** | Construct the weighted human physical PPI matrix from BioGRID 4.4.243; retain the main connected network; compute effective conductance and reciprocal association matrices. | [`1_Effective_conductance.R`](scripts/1_Effective_conductance.R) |
| **2. Retained, gained, and lost interactors** | Compare observed PPIs with reciprocal-score-based relationships; examine interaction changes and Gene Ontology (GO) overlap using enrichment and permutation analyses. | [`2_About_interactors.R`](scripts/2_About_interactors.R) |
| **3. Protein complex interactions** | Process CORUM 5.3 annotations; calculate complex–complex interaction (CCI) scores and size-matched permutation backgrounds; examine custom-defined complexes, including RNA polymerases. | [`3_About_New_Complexes.R`](scripts/3_About_New_Complexes.R) |
| **4. RNA polymerases and mass spectrometry** | Examine within- and cross-polymerase interactions involving Pol I, II, and III subunits; compare network associations with pull-down mass-spectrometry data. | [`4_About_RNA_polymerases_and_MASSSPEC.R`](scripts/4_About_RNA_polymerases_and_MASSSPEC.R) |
| **5. NuRD, DNMT3B, and PRC2/3 example** | Exploratory visualization of interactions among selected chromatin-regulatory complexes after removing shared subunits. | [`4.1_About_Nurd.R`](scripts/4.1_About_Nurd.R) |
| **6. Network visualizations** | Explore heatmaps, network layouts, communities, embeddings, and protein-centered interaction plots. | [`5_Trying_to_plot_network.R`](scripts/5_Trying_to_plot_network.R) |
| **7. Protein features and connectivity** | Retrieve UniProtKB protein features through BioMart/UniProt and test their association with network connectivity using permutation-based scores, including TF and RBP subsets. | [`6_About_structural_features_of_hub_proteins.R`](scripts/6_About_structural_features_of_hub_proteins.R) |
| **8. Figure 1 analyses** | Generate PPI/conductance matrix views and reciprocal-plot and simulation illustrations used during Figure 1 preparation. | [`7_Figure_1_for_paper.R`](scripts/7_Figure_1_for_paper.R) |
| **9. Network simulations** | Construct simulated interaction networks and examine effective conductance and reciprocal analysis under varying protein popularity/experimental sampling. | [`8_SIMULATION.R`](scripts/8_SIMULATION.R) |


## Core computational methods

### Human PPI network and effective conductance

We downloaded BioGRID release **4.4.243**, retained physical interactions in which both interactors were human (NCBI taxonomy ID 9606), removed self-interactions, and represented each protein pair by its number of reported interaction observations. We retained the largest connected component of the network. The resulting weighted PPI matrix was treated as the conductance matrix of an electrical network. We computed its graph Laplacian, obtained the Moore–Penrose pseudoinverse, and calculated pairwise effective conductance values.

### Reciprocal association analysis

Effective conductance scores were standardized separately for each protein. For each focal protein, we fitted a smoothing spline to its reciprocal  plot and calculated deviations from the fitted relationship. We then combined the protein-centered scores from both directions to obtain a symmetric pairwise association matrix (`ECdz_min.Rdata`). The downstream analyses compare these scores with the original observed PPI network to identify retained, lost, and newly inferred protein relationships. A threshold of **2** is used in the principal retained/gained/lost interaction comparisons; individual exploratory plots may use other thresholds.

### Complex–complex interaction scores

Protein complexes were obtained from **CORUM version 5.3**. For each pair of complexes, we considered all possible inter-complex protein pairs and calculated the median of their underlying reciprocal association scores. For every combination of complex sizes, we generated a background distribution by randomly sampling protein sets of the corresponding sizes **10,000 times** and calculating the median score in each iteration. We converted each observed complex-pair median to a **Z-score** using the mean and standard deviation of its corresponding background distribution. These values were assembled into a complex–complex interaction matrix .

For RNA polymerase comparisons, the study focuses on enzyme-specific subunits, excluding shared subunits.

## Data sources and availability

| Source | Version / purpose | Access |
| --- | --- | --- |
| **BioGRID** | 4.4.243; human physical PPI observations | [Download ZIP](https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.243/BIOGRID-ALL-4.4.243.tab3.zip) |
| **CORUM** | 5.3; protein complex annotations | [CORUM website](https://mips.helmholtz-muenchen.de/corum/) |
| **UniProtKB** | Protein sequence and feature annotations | [UniProt](https://www.uniprot.org/) |
| **Ensembl BioMart** | Gene-symbol to UniProt mapping | [Ensembl BioMart](https://www.ensembl.org/info/data/biomart/index.html) |
| **Additional analysis inputs** | Pull-down MS tables, expression/proteomics reference data, and selected protein lists | See [`data/README.md`](data/README.md); not included in this package |

### Download BioGRID

From the repository root, run:

```bash
Rscript scripts/00_download_biogrid.R
```

This saves `BIOGRID-ALL-4.4.243.tab3.zip` in the repository root. **No manual extraction is needed**; the main analysis reads the ZIP directly:

```r
library(data.table)
ppi <- fread("BIOGRID-ALL-4.4.243.tab3.zip", sep = "\t")
```

Download the CORUM 5.3 table separately and save it as `CORUM_5.3.txt` in the repository root, matching the path expected by the complex-analysis script. For additional required inputs, see [`data/README.md`](data/README.md).

## Main intermediate files

The scripts write the following R objects (among others) to the current working directory:

| File | Description | Produced by |
| --- | --- | --- |
| `interaction_matrix.Rdata` | Weighted observed PPI matrix | Script 1 |
| `bias.Rdata` | Protein-level study/interaction frequency summary | Script 1 |
| `EC.Rdata` | Pairwise effective conductance matrix | Script 1 |
| `tem.Rdata`, `ECd.Rdata` | Intermediate reciprocal analysis matrices | Script 1 |
| `ECdz_min.Rdata` | Symmetric reciprocal association scores | Script 1 |
| `bias2.Rdata` | Counts of retained, gained, and lost partners | Script 2 |
| `complexinfo_CORUM_human.Rdata` | Processed CORUM complex annotations | Script 3 |
| `backg_mean29n.Rdata`, `backg_sd29n.Rdata` | Size-matched CCI background statistics | Script 3 |
| `CIM_corum.Rdata` | Complex–complex interaction score matrix | Script 3 |
| `h1_freq.Rdata`, `h1_length.Rdata` | Protein-feature summaries | Script 7 (protein features) |

Despite their `.Rdata` extension, these files are written with `saveRDS()` and should be loaded with **`readRDS()`**.

## Software and reproducibility notes

- Analyses were developed in **R**. Packages used across the scripts include `data.table`, `tidyverse`/`dplyr`, `igraph`, `MASS`, `ggplot2`, `clusterProfiler`, `org.Hs.eg.db`, `GO.db`, `biomaRt`, `IRanges`, `visNetwork`, `ggraph`, and `plotly`, with additional packages used in individual analyses.
- Several scripts call **custom helper functions** from `CODES/my_function.R` and `CODES/plotTriangle_ver2.R`. These helper files were not part of the supplied scripts and must be added before all analyses can be rerun; see [`CODES/README.md`](CODES/README.md).
- Some scripts reference **local file paths** (such as `~/Downloads/...`), existing R objects, or unpublished input tables. Update these paths and provide the corresponding inputs before running those sections.
- The core effective-conductance calculation and the complex-size permutation analysis can be computationally intensive.  Scripts are intended to be run **section by section**, rather than sourced end-to-end without review.

## Citation

If you use this code or derived resources, please cite the accompanying manuscript:

> K C, R., Zhou, S., Cheng, R., Penaherrera, M. A., Lizarazo, S., and Van Bortle, K. *Reciprocal effective conductance-based inference of protein-protein relationships delineates protein network hubs and macromolecular properties.* [Publication details / DOI to be added.]

**CoBRAA interactive tool (github) :** https://github.com/VanBortleLab/CoBRAA

**CoBRAA interactive tool (website):** http://35.222.228.120:3838/ppi-cci/ (Takes less than 1 minute to load)
