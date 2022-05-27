# Replication Instruction

Please follow below instructions in order to replicate both our simulation and real-data results. We manage to make people replicate the results with minimum effort.

If you have any questions, please contact Zeyun Lu at zeyunlu@usc.edu or Nicholas Mancuso at Nicholas.Mancuso@med.usc.edu.

```diff
- We detest usage of our software or scientific outcome to promote racial discrimination.
```

## Main Figures

You can find Figure 1 and Figure S1 in `/mainfigure/*`, and see below for other other figures.

## Simulation

We make bash script read in parameter files and then plug into python files with different scenarios.

You need to change the paths necessarily to make the flow run.

You also need to have 1000G plink data for EUR, AFR, and EAS.

### Running Scripts

#### Get regions

We need to sample 100 population-wise independent genomic regions resulted from LDetect. Two bash scripts are used and dependent on:

* `/simulation/codes/sample_regions/sample_genes.py`
* `/simulation/codes/sample_regions/grch37.eur.eas.afr.loci.bed` (For 3-pop MA-FOCUS)
* `/simulation/codes/sample_regions/grch37.eur.afr.loci.bed ` (For 2-pop MA-FOCUS)

Two bash scripts:

1. `/simulation/codes/sample_regions/mefocus_get_genes_eur_afr_eas.sh` (For 3-pop MA-FOCUS)
2. `/simulation/codes/sample_regions/mefocus_get_genes_eur_afr.sh` (For 2-pop MA-FOCUS)

#### Get parameter files

Specify and run `/simulation/generate_param.R` to make parameter files.

#### Run simulation for MA-FOCUS

There are 5 bash scripts `/simulation/codes/bash/me_focus_sim*.sh` for 5 scenarios:

1. eur_afr: to run general 2-pop MA-FOCUS working with parameter and python files:
	* `/simulation/codes/param/param_eur_afr_3pop.tsv`
	* `/simulation/codes/param/param_eur_afr.tsv`
	* `/simulation/codes/param/param_eur_afr_real.tsv`
	* `/simulation/codes/python/mefocus_run_eur_afr.py`
	* `/simulation/codes/python/mefocus_run_eur_real.py`

2. eur\_afr\_eas: to general run 3-pop MA-FOCUS working with parameter and python files:
	* `/simulation/codes/param/param_eur_afr_eas.tsv`
	* `/simulation/codes/python/mefocus_run_eur_afr_eas.py`

3. pop2h2ge: to run 2-pop MA-FOCUS when 2nd population has different h2ge working with parameter and python files:
	* `/simulation/codes/param/param_pop2h2ge.tsv`
	* `/simulation/codes/python/mefocus_run_pop2h2ge.py`

4. pop2weights: to run 2-pop MA-FOCUS when 2nd population eQTL weights are not available and substituted by 1st population eQTL weights working with parameter and python files:
	* `/simulation/codes/param/param_pop2weights.tsv`
	* `/simulation/codes/python/mefocus_run_pop2weights.py`

5. proxy: to run 2-pop MA-FOCUS when a proxy tissue is used for eQTL weights for 2nd population working with parameter and python files:
	* `/simulation/codes/param/param_proxy.tsv`
	* `/simulation/codes/python/mefocus_run_proxy.py`

There are three dependent files for the MA-FOCUS simulation python files:

1. `/simulation/codes/python/finemap_1.py`
2. `focus.functions_2pop.py`(For 2-pop MA-FOCUS)
3. `focus.functions_3pop.py`(For 3-pop MA-FOCUS)

### Data

Simulation data is in `/simulation/data/*.RDat` with above different scenario names.

### Analysis

#### Figures

The R codes that generate figures is in `/simulation/plot.R`, which depends on `/simulation/figure_style.R`

#### Manuscript

The result section 1 and 2 in manuscript contain the statistics computed in `/simulation/sec1-2.R`. The order of the statistics may not match that in manuscript, but all the numbers should be covered by this script. Sometimes, the P-value is divided by two because we want to test one tail while the R regression results give us two tailed P-value.


## Real-data

The scripts for the basic cleaning and manipulation of GWAS, TWAS, MA-FOCUS, and heritability results can be found in `/real-data/running/*`.

### Analysis

#### Figures

The R codes that generate figures is in `/real-data/plot.R`, which depends on `/real-data/figure_style.R`. All necessary data locate in `/real-data/data/`.

#### Manuscript

The R codes that generate result section 3 and 4 locates in `/real-data/sec3.R` and `/real-data/sec4.R`. All the required data locates in `/real-data/data/`.

The order of the statistics may not match that in manuscript, but all the numbers should be covered by this script. Sometimes, the P-value is divided by two because we want to test one tail while the R regression results give us two tailed P-value.

The entire codes for admixture analysis can be found in `real-data/admixture/*`
#### Table

The R codes that generate tables is in `/real-data/talbe.R`.
#### Data
The default path is in `/real-data/data/`
| Analysis Results | Paths |
| ------------- | ----- |
|  GWAS signal summary statistics for 115 regions that do not exhibit genome-wide signals | `115_gwas_signal.tsv` |
| Admixture 1000G sample size | `Adadmixture_sample_size.tsv` |
| TWAS | `twas.RData` |
| TWAS using GEUVADIS weights | `twas_geuvadis.RData` |
| MA-FOCUS all data (no matter TWAS significant or not) | `focus_all.RData` |
| MA-FOCUS data (no matter TWAS significant or not) using GEUVADIS weights | `focus_geuvadis.RData` |
| MA-FOCUS data (no matter TWAS significant or not) using maxcausal = 1 | `focus_maxgene1.RData` |
| MA-FOCUS data (no matter TWAS significant or not) using maxcausal = 5 | `focus_maxgene5.RData` |
| MA-FOCUS all data (TWAS significant for both ancestries; our analysis data) | `focus_analysis.RData` |
| Enrichment | `enrich.RData` |
| Silver | `silver.tsv` |
| DisGeNET Category | `DisGeNET_meta_categories.csv` |
| Traits used in enrichment analysis | `restricted_blood_traits.csv` |
| GENOA Heritability (all genes for ea) | `genoa_heritability_ea_all_genes.tsv` |
| GENOA Heritability (all genes for aa) | `genoa_heritability_aa_all_genes.tsv` |
| GENOA Heritability (analyzed genes) | `genoa_heritability_analyzed.tsv` |
| GEUVADIS Heritability | `geuvadis_her_total.tsv` |
| Genetic Variance CV R2 | `total_r2.tsv` |
| TWAS&GWAS Correlation | `twas_gwas_corr.tsv` |
| Gencode GRCh38 | `gencode.v26.GRCh38.genes.only.tsv` |
| GWAS signals on each LD block | `gwas_all_z2_signals.tsv` |
| 1000G subjects information | `igsr_samples.tsv` |
| Genomic regions with TWAS signals but without GWAS signals | `region_TWAS_nonGWAS.tsv` |
| LD blocks that doesn't have GWAS signals | `all_gwas_not_sig_ld.bed` |
| LD blocks of TWAS genes | `twas_all_ld_region.bed` |
| LD blocks that has GWAS signals | `gwas_ld_sig/*` |
| LD blocks that has TWAS signals | `twas_ld_sig` |