# Reassortment Dynamics

Supplementary code and data for inferring viral reassortment dynamics from phylogenetic data. This repository contains simulation frameworks for validating reassortment inference methods under different epidemic models, and applies these methods to analyze H5N1 avian influenza reassortment patterns in North American wild birds.

## Overview

The project addresses three main objectives:

1. **Simulation & Validation** — Forward simulations of reassortment under SIR, SIS, and structured population models, used to validate Bayesian phylogenetic inference of reassortment rates.
2. **Inference Models** — BEAST2 + [CoalRe](https://github.com/nicfel/CoalRe) configurations for estimating time-varying reassortment rates under constant, skygrowth, and Ne-dependent models.
3. **Application to H5N1** — Analysis of reassortment patterns in North American H5N1 influenza, testing for associations with migratory flyways and avian host orders.

## Repository Structure

```
├── Simulations/
│   ├── SIR/                  # SIR model simulations (with/without superspreading)
│   ├── SIR_SIS/              # Paired SIR vs. SIS model comparisons
│   ├── StructuredSIR/        # Structured population (50 demes) simulations
│   ├── WaitimeDistribution/  # Co-infection wait time analysis
│   └── EventRate/            # Reassortment rate vs. prevalence
├── Validation/               # BEAST2 MCMC validation (prior vs. posterior)
├── Applications/
│   ├── H5N1NorthAmerica/     # H5N1 sequences, metadata, and inference XMLs
│   └── NetworkViz/           # Reassortment network visualizations (BALTIC)
├── Scripts_and_XML/
│   └── code/                 # Statistical tests (observed/expected ratios,
│                             #   TMRCA extraction, detection counts)
├── Test/                     # MATLAB spline interpolation methods
└── Figures/                  # Publication-ready PDF figures
```

## Dependencies

- **[BEAST2](https://www.beast2.org/)** (v2.7+) with the [CoalRe](https://github.com/nicfel/CoalRe) package
- **[Seq-Gen](http://tree.bio.ed.ac.uk/software/seqgen/)** — sequence simulation along phylogenies
- **R** — `dplyr`, `ggplot2`, `tidyverse`, `coda`, `treeio`, `tidytree`, `patchwork`, `seqinr`
- **Python** — `pandas`, `matplotlib`, `seaborn`, `baltic`
- **MATLAB** (optional) — spline interpolation comparisons (dead code)

## Reproducing the Analyses

### Simulations

Each simulation directory contains an R script (`runSims.R` or similar) that generates BEAST2 XML files from templates with randomly sampled epidemiological parameters (transmission rate, recovery rate, sampling proportion, population size).

```bash
cd Simulations/SIR
Rscript runSims.R
```

This produces simulation XMLs in `xmls/`. Run them with BEAST2:

```bash
beast -seed 1234 xmls/SIR_simulations_1.xml
```

After inference completes, plot results comparing true vs. inferred parameters:

```bash
Rscript plotSIRResults.R
```

### H5N1 Application

1. **Build inference XMLs** from FASTA alignments and metadata:
   ```bash
   cd Applications/H5N1NorthAmerica
   Rscript buildXmls.R
   ```

2. **Run BEAST2 inference** on the generated XMLs.

3. **Statistical analysis** of genotype associations with flyways and host orders:
   ```bash
   cd Scripts_and_XML/code
   Rscript observed_expect_stattest_v3.R
   ```

4. **Extract TMRCAs** from MCC trees:
   ```bash
   Rscript MCC_targetedbeast_tmrca_geno_flyway.R
   Rscript MCC_targetedbeast_tmrca_geno_orders.R
   ```
### XML Generation

To generate coalRe XML files for a dataset/virus with a given number of segments see **buildCoalReXml.R** in Scripts_and_XML/code direcory. **README_buildCoalReXml.md** for more information on how to run the script. 
