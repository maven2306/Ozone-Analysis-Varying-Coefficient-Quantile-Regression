# Analysis on the Effect of Meteorological Factors on Ozone Concentration Using Varying-Coefficient Quantile Regression

> **Author:** Matteo Venturini
>
> **Course:** Non-parametric methods
> 
> **University:** Hasselt University
> 
> **Date:** June 8, 2025

### Full report

For a complete overview of the statistical models and in-depth discussion, please see the **[full project report](ozone_concentration_analysis.pdf)**.

## Project Overview

This project analyzes the dynamic relationship between meteorological factors and ozone concentration in New York from May to September 1973. The analysis employs varying-coefficient quantile regression to model how the effects of solar radiation, wind speed, and temperature on ozone levels change over time. The study highlights the seasonal patterns in ozone concentration and the varying influence of weather conditions, providing insights relevant to public health.

## Dataset

The analysis utilizes the `airquality` dataset from R, which contains 153 daily observations of air quality measurements. After removing rows with missing values, 111 complete observations remain.

The variables used in this analysis are:
*   **Ozone:** Mean ozone concentration in parts per billion (ppb) from 1300 to 1500 hours at Roosevelt Island.
*   **Solar.R:** Solar radiation in Langleys in the frequency band 4000-7700 Angstroms from 0800 to 1200 hours at Central Park.
*   **Wind:** Average wind speed in miles per hour (mph) at 0700 and 1000 hours at LaGuardia Airport.
*   **Temp:** Maximum daily temperature in degrees Fahrenheit (°F) at La Guardia Airport.

## Methods and Models

The analysis was conducted using R version 4.5.0 and involves several non-parametric techniques:

*   **Kernel Density Estimation:** To examine the distributions of the dependent and independent variables.
*   **Nadaraya-Watson and Local Polynomial Regression:** To model the bivariate relationships between ozone and each predictor.
*   **Varying-Coefficient Quantile Regression:** The primary model used to analyze the data, implemented with the `AHeVT` function from the `QRegVCM` package. This approach allows for the modeling of different quantiles of the ozone distribution and accounts for non-crossing quantile curves and heteroscedasticity.

The model equation is:
`Y_Ozone(t) = β₀(t) + β₁(t)X_Solar.R(t) + β₂(t)X_Wind(t) + β₃(t)X_Temp(t) + V(t)ε(t)`

## Key Findings

*   **Seasonal Trends:** Median ozone concentrations show a strong seasonal trend, peaking during the summer months (June-July).
*   **Solar Radiation:** Has a consistently positive but small effect on ozone levels, with the influence being most pronounced during the summer.
*   **Wind Speed:** Exhibits a notable negative effect on ozone, which becomes more accentuated from the onset of summer.
*   **Temperature:** Positively influences ozone levels, with its impact strengthening during the summer before slightly decreasing towards August.
*   **Public Health Implications:** Under average conditions, median ozone levels can approach the WHO's recommended 8-hour maximum of 50 ppb during summer peaks. In worst-case meteorological scenarios, even lower quantiles can significantly exceed this guideline.

## Repository Structure

*   **`Analysis.R`**: The R script containing all the code for data cleaning, analysis, and visualization.
*   **`Report.pdf`**: The full project report detailing the methodology, analysis, and conclusions.
*   **`README.md`**: This file.

## How to Run the Analysis

### Prerequisites

You will need to have R and RStudio installed on your system.

### R Packages

The following R packages are required to run the analysis. You can install them using the `install.packages()` function in R.

```R
install.packages(c("QRegVCM", "KernSmooth", "ggplot2", "tidyr", "dplyr"))
