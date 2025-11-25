# Cubature-scheme-for-likelihood-approximation-guidelines-for-its-setting-parameters
This repository contains codes and data to perform the simulation study and the application on real earthquake dataset related to the paper "Cubature scheme for likelihood approximation: guidelines for its setting parameters in three-dimensional Poisson point processes".

The material is provided for transparency and reproducibility.  
**Important:** running the full simulation study is computationally very expensive and is **not recommended** for routine use.

---

## Repository structure

### Root

- `README.md`  
  This file. Overview of the repository and its contents.

- `LICENSE`  
  Apache 2.0 license for the material in this repository.

---

### Folder `R/`

R scripts used for the simulation study and the real-data application.

- `R/README.md`  
  Short description of the R scripts contained in this folder.

- `R/codes_application_earthquake_dataset.R`  
  Script to reproduce the three-dimensional earthquake application presented in the paper.

- `R/codes_simulation_study_GET.R`  
  Script used in the simulation study involving the GET-based (global envelope test) diagnostics for the cubature schemes.

- `R/codes_simulation_study_mse.R`  
  Script used in the simulation study to compute mean square error (MSE) measures for different cubature configurations for all the simulated scenarios.

- `R/codes_simulation_study_validation_guidelines.R`  
  Script used to validate and summarise the proposed guidelines for choosing cubature parameters.

- `R/good_comb_grid_size.R`  
  Script to identify suitable combinations of cubature settings for the simulation scenarios.

> **Note:** The scripts in `R/` are designed to reproduce large-scale simulation experiments.  
> Running them end-to-end may require substantial computing time and memory (hours or days on a multi-core machine).

---

### Folder `RData/`

RData files containing data objects (datasets and functions) used by the scripts in `R/`.

- `RData/README.md`  
  Short description of the RData files contained in this folder.

- `RData/Libraries`  
  R script that contains all the R packages required by the other scripts.  

