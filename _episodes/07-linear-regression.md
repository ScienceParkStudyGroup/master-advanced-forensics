---
title: "Linear, Multiple and Regularised Regressions"
teaching: 30
exercises: 30
questions:
- ""
objectives:
- ""
keypoints:
- ""

---

# Table of Contents
<!-- MarkdownTOC autolink="True" levels="1,2" -->

- [1. Introduction](#1-introduction)
  - [1.1 Definitions](#11-definitions)
- [Linear Regression](#linear-regression)
  - [Introduction](#introduction)
  - [Model performance](#model-performance)
  - [Issues with omics](#issues-with-omics)
- [The high-dimensionality curse of omic datasets](#the-high-dimensionality-curse-of-omic-datasets)
  - [Diagnostic plots](#diagnostic-plots)
  - [Making predictions](#making-predictions)
  - [Limitations and possible solutions](#limitations-and-possible-solutions)
- [Regularised Linear Regression](#regularised-linear-regression)
- [X. Going further](#x-going-further)
  - [8.1 Useful links](#81-useful-links)
  - [8.2. References](#82-references)

<!-- /MarkdownTOC -->

# 1. Introduction

## 1.1 Definitions

independent variable / explanatory variables
dependent variable / response variable


# Linear Regression

## Introduction 

Please read the [simple linear regression section](https://bradleyboehmke.github.io/HOML/linear-regression.html#simple-linear-regression) and the [multiple linear regression](https://bradleyboehmke.github.io/HOML/linear-regression.html#multi-lm) sections in [Bradley Boehmke "Hands-On Machine Learning with R" online book](https://bradleyboehmke.github.io/HOML/). 

## Model performance

Under the usual assumptions stated above, an unbiased estimate of the error variance is given as the sum of the squared residuals divided by $$n − p$$ (where $$p$$
is the number of variables in the model):

$$σ^2 = \frac{1}{n−p} \sum_{i=1}^{n} r_i^2$$  

where 

$$r_{i} = (Y_{i} - Y_{i})$$

$$r_{i}$$ is referred to as the _ith_ residual (i.e., the difference between the _ith_ observed and predicted response value). The quantity $$\sigma^{2}$$
is also referred to as the mean square error (MSE) and its square root is denoted __RMSE for Squared Root Mean Squared Error__. In R, the RMSE of a linear model can be extracted using the `sigma()` function.


## Issues with omics
- Linear relationship;
- There are more observations (n) than features (p) (n>p): very often not true since "omics" will measure thousands of variables (genes, CpG sites) at the same time.
- No or little multicollinearity: often not true since genes are co-regulated and methylation sites can be part of the same DNA region. 

# The high-dimensionality curse of omic datasets

High $$p$$ small $$n$$

High number of variables being measured (more than thousands usually, here 27,578) and the low number of samples (66 here). 

It is also called the residual sum of squares. In R, the `summary(fitted_regression_model)` will indicate this as the `Residual standard error`.

## Diagnostic plots
Residuals
Systematic errors versus no pattern/systematic pattern
Gain curve?

## Making predictions

Choose a few beta_values for sites X Y Z ... and report the predicted age. 

## Limitations and possible solutions
Colinearity of independent variables

# Regularised Linear Regression

# X. Going further 

## 8.1 Useful links
- [BiomartR](https://docs.ropensci.org/biomartr/)
- [Arabidopsis.org (TAIR) list of data mining tools](https://www.arabidopsis.org/portals/expression/microarray/microarrayExpressionV2.jsp)
- [ResearchGate related question](https://www.researchgate.net/post/How_can_I_analyze_a_set_of_DEGs_differentially_expressed_genes_to_obtain_information_from_them)	

## 8.2. References
* [The Cluster Profiler companion boo, a great place to start](https://yulab-smu.github.io/clusterProfiler-book/chapter2.html)
* Zhou et al. (2019). Metascape provides a biologist-oriented resource for the analysis of systems-level datasets. Nat Commun 10, 1523 (2019). [link](https://doi.org/10.1038/s41467-019-09234-6)
* Yates et al. (2020) Ensembl 2020, Nucleic Acids Research, Volume 48, Issue D1, 08 January 2020, Pages D682–D688, [Link](https://doi.org/10.1093/nar/gkz966)
* Tian et al. (2017) agriGO v2.0: a GO analysis toolkit for the agricultural community. _Nucleic Acids Research_, Volume 45, Issue W1, Pages W122–W129.[Link](https://doi.org/10.1093/nar/gkx382) 
* MapMan: [MapMan4: A Refined Protein Classification and Annotation Framework Applicable to Multi-Omics Data Analysis. Schwacke et al. _Molecular Plant_, 12(6):879-892](https://doi.org/10.1016/j.molp.2019.01.003)
* Drost et al. (2017) Biomartr: genomic data retrieval with R. _Bioinformatics_ 33(8): 1216-1217. [doi:10.1093/bioinformatics/btw821](https://academic.oup.com/bioinformatics/article/33/8/1216/2931816).
* Darzi et al. (2018) iPath3.0: interactive pathways explorer v3. _Nucleic Acids Research_, Volume 46, Issue W1, 2 July 2018, Pages W510–W513, [link](https://doi.org/10.1093/nar/gky299)


