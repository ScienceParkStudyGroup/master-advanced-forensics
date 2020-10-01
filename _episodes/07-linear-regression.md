---
title: "Regressions"
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
	- [Model performance](#model-performance)
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

## Model performance
R squared

## Diagnostic plots
Residuals
Systematic errors versus no pattern/systematic pattern
Gain curve?

## Making predictions

Choose a few beta_values for sites X Y Z ... and report the predicted age. 

## Limitations and possible solutions
Colinearity of independent variables

# Regularised Linear Regression

Taken from Bradley "Hands-On Machine Learning with R" book section on [regularised regressions](https://bradleyboehmke.github.io/HOML/regularized-regression.html): 	 
> Many real-life data sets, like those common to text mining and genomic studies are wide, meaning they contain a larger number of features ($$p > n$$). As $$p$$ increases, we’re more likely to violate some of the OLS [Ordinary Least Squares] assumptions and alternative approaches should be considered. [...] Having a large number of features invites additional issues in using classic regression models. For one, having a large number of features makes the model much less interpretable. Additionally, when $$p > n$$, there are many (in fact infinite) solutions to the OLS problem! 

This is our case as our dataset originally comprised 22,758 CpG sites (features) on 66 samples (observations). 


In such cases, it is useful (and practical) to assume that a smaller subset of the features exhibit the strongest effects (something called the bet on sparsity principle (see Hastie, Tibshirani, and Wainwright 2015, 2).). For this reason, we sometimes prefer estimation techniques that incorporate feature selection. One approach to this is called hard thresholding feature selection, which includes many of the traditional linear model selection approaches like forward selection and backward elimination. These procedures, however, can be computationally inefficient, do not scale well, and treat a feature as either in or out of the model (hence the name hard thresholding). In contrast, a more modern approach, called soft thresholding, slowly pushes the effects of irrelevant features toward zero, and in some cases, will zero out entire coefficients. As will be demonstrated, this can result in more accurate models that are also easier to interpret.

.




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


