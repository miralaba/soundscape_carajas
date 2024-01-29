# Amazon Rainforest Soundscape Analysis


## Introduction

This GitHub repository hosts the data and R scripts for our study on the impact of environmental variables on soundscape indices in the Amazon rainforest's undegraded primary forests. The primary research questions investigate how temporal and frequency masking effects influence soundscape indices and whether these indices can be effectively used as a monitoring tool for environmental integrity. The project utilizes passive acoustic monitoring data and statistical analyses conducted in R.
Repository Structure

The repository is organized into three main folders:

    data/: Contains raw acoustic data files collected from 14 sites within two protected areas in the Amazon biome.
    scripts/: Includes two R scripts:
        _soundscape.R: For data preprocessing and extraction of soundscape indices.
        _analyses.R: Contains the statistical analysis workflow, as detailed in the methods section of the study.
    results/: Stores the main figures generated from our analyses, visualizing key findings.

## Prerequisites

To run the scripts, the following software and packages are required:

    R (Version 4.0 or later): Download R
    RStudio (Recommended for ease of use): Download RStudio
    Required R Packages: soundecology, tuneR, seewave, ggplot2, glmmTMB, MuMIn, DHARMa, rstatix, corrplot. Install these packages using install.packages("packageName").

Contact Information

For inquiries, please contact:

    Project Lead: Leonardo de Sousa Miranda
    Email: miralaba@gmail.com