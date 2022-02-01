# Replication Instruction

Please follow below instructions in order to replicate both our simulation and real-data results. We manage to make people replicate the results with minimum effort.

If you have any questions, please contact Zeyun Lu at zeyunlu@usc.edu or Nicholas Mancuso at Nicholas.Mancuso@med.usc.edu. 

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

#### Run MA-FOCUS

There are 5 bash scripts `/simulation/codes/bash/me_focus_sim*.sh` for 5 scenarios:

1. eur_afr: to run general 2-pop MA-FOCUS working with parameter files:
  * `/simulation/codes/param/param_eur_afr_3pop.tsv`
  * `/simulation/codes/param/param_eur_afr.tsv`
  * `/simulation/codes/param/param_eur_afr_real.tsv`
  
  on python files:
  
  * `/simulation/codes/python/mefocus_run_eur_afr.py`
  * `/simulation/codes/python/mefocus_run_eur_real.py`

2. eur\_afr\_eas: to general run 3-pop MA-FOCUS working with parameter file:
  * `/simulation/codes/param/param_eur_afr_eas.tsv`

  on python file:
  
  * `/simulation/codes/python/mefocus_run_eur_afr_eas.py`

3. pop2h2ge: to run 2-pop MA-FOCUS when 2nd population has different h2ge working with parameter file:
  * `/simulation/codes/param/param_pop2h2ge.tsv`

  on python file:
  
  * `/simulation/codes/python/mefocus_run_pop2h2ge.py`

4. pop2weights: to run 2-pop MA-FOCUS when 2nd population eQTL weights are not available and substituted by 1st population eQTL weights working with parameter file:

  * `/simulation/codes/param/param_pop2weights.tsv`

  on python file:

  * `/simulation/codes/python/mefocus_run_pop2weights.py`

5. proxy: to run 2-pop MA-FOCUS when a proxy tissue is used for eQTL weights for 2nd population working with parameter file:

  * `/simulation/codes/param/param_proxy.tsv`

  on python files:

  * `/simulation/codes/python/mefocus_run_proxy.py`

There are three dependent files for the MA-FOCUS simulation python files:

1. `/simulation/codes/python/finemap_1.py`

2. `focus.functions_2pop.py`(For 2-pop MA-FOCUS)

3. `focus.functions_3pop.py`(For 3-pop MA-FOCUS)

### Data

Simulation data is in `/simulation/data/*.RDat` with above different scenario names.

### Analysis

#### Figures

The R codes that generate figures is in `/simulation/draw.R`, which depends on `/simulation/figure_style.R`

#### Manuscript

The result section 1 and 2 in manuscript contain the statistics computed in `/simulation/sec12analysis.R`. The order of the statistics may not match that in manuscript, but all the numbers should be covered by this script. Sometimes, the P-value is divided by two because we want to test one tail while the R regression results give us two tailed P-value.


## Real-data

We omit the scripts for the basic cleaning and manipulation of GWAS, TWAS, MA-FOCUS, and heritability results.

### Analysis

#### Figures

The R codes that generate figures is in `/real-data/draw.R`, which depends on `/real-data/figure_style.R`. All necessary data locate in `/real-data/data/`.

#### Manuscript

The R codes that generate result section 3 and 4locates in `/real-data/sec3.R` and `/real-data/sec4.R`. All the required data locates in `/real-data/data/`.

The order of the statistics may not match that in manuscript, but all the numbers should be covered by this script. Sometimes, the P-value is divided by two because we want to test one tail while the R regression results give us two tailed P-value.













