# Circadian-Timing-of-24-Hour-Movement-Behaviors-and-Respiratory-Disease-Risk
This repository contains the analysis code for the study:

**Circadian Timing of 24-Hour Movement Behaviors and Respiratory Disease Risk: Insights from a Prospective Cohort Study of 92,326 Participants**

## Overview

This project investigates how the timing of 24-hour movement behaviors, including sedentary behavior (SB), light physical activity (LPA), moderate-to-vigorous physical activity (MVPA), and sleep, relates to respiratory disease risk.

The repository currently includes code for:

- mediation and multi-omics analysis
- behavioral recommendation modeling
- batch Mendelian randomization analysis
- model construction and validation in Python

## Requirements

### R
- R >= 4.3.0
- tidyverse
- survival
- data.table
- readxl
- TwoSampleMR
- purrr
- dplyr

### Python
- Python >= 3.9

Script Description
mediation_analysis.R

Performs the mediation and multi-omics analysis workflow, including biomarker screening, candidate mediator selection, and mediation analysis.

recommendation_analysis.R

Generates optimized 24-hour behavioral recommendations based on disease-specific risk patterns.

mr_analysis.R

Runs batch two-sample Mendelian randomization analyses between movement-related exposures and respiratory disease outcomes.

test_p.py

Python code for model construction.

test_val.py

Python code for model validation and evaluation.
