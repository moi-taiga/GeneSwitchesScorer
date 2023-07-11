
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GeneSwitchesScorer

<!-- badges: start -->
<!-- badges: end -->

The goal of GeneSwitchesScorer is to identify the position of a sample
upon a trajectory.

## Installation

You can install the development version of GeneSwitchesScorer from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("moi-taiga/GeneSwitchesScorer")
```

## Example

[GSS
Workflow](https://moi-taiga.github.io/GeneSwitchesScorer/articles/GSS_Workflow.html)

## Limitations:

## Assumptions:

- Sample is found upon the chosen trajectory.
- Sample is from a distinct part of the trajectory. A sample with cells
  that are evenly distributed across the trajectory will have a
  predicted location of the centre of the trajectory.
