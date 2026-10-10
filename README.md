# Kerner-Hybrid-Glucose-Models

## Overview

This repository contains MATLAB implementations of two glucose–insulin models:

1. Kerner model
2. Hybrid model

Both models were fitted to FSIGT (Frequently Sampled Intravenous Glucose Tolerance Test) data and evaluated using:

- RMSE
- R²
- AIC
- BIC
- Residual analysis
- Shapiro-Wilk normality test

---

## Repository Contents

Kerner_Model.m
    Complete implementation of the Kerner model.

Hybrid_Model.m
    Complete implementation of the Hybrid model.

swtest.m
    Third-party implementation of the Shapiro-Wilk normality test.

---

## Software Requirements

- MATLAB
- Optimization Toolbox

---

## Third-Party Function

This project uses:

swtest.m

for the Shapiro-Wilk normality test. The file is included in this repository to ensure reproducibility.

---

## Running the Code

1. Download or clone the repository.
2. Open MATLAB.
3. Set the Current Folder to the repository folder.
4. Run one of the following scripts:

Kerner_Model

or

Hybrid_Model

No additional data files are required because the FSIGT data are included directly within each script.

---

## Outputs

The scripts generate:

- Model simulations
- Parameter estimates
- Clinical indices (SI and SG)
- Goodness-of-fit metrics
- Information criteria (AIC and BIC)
- Residual analysis plots
- Shapiro-Wilk normality test results

---

## Authors

C. Maass

Z. de Villiers

Stellenbosch University

BSc Biomedical Mathematical Sciences
