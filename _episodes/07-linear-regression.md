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
- [2. Linear Regression](#2-linear-regression)
  - [2.1 Introduction](#21-introduction)
  - [2.2 Model performance](#22-model-performance)
  - [2.3 Diagnostic plots](#23-diagnostic-plots)
  - [2.4 Making predictions](#24-making-predictions)
  - [2.5 Limitations and possible solutions](#25-limitations-and-possible-solutions)
- [3. Regularised Linear Regression](#3-regularised-linear-regression)
  - [3.1 The high-dimensionality curse of omic datasets](#31-the-high-dimensionality-curse-of-omic-datasets)
  - [3.1 Ridge regression](#31-ridge-regression)
  - [3.2 LASSO regression](#32-lasso-regression)
- [8. References](#8-references)
  - [8.1 Online books and tutorials](#81-online-books-and-tutorials)
  - [8.2 Publications and books](#82-publications-and-books)

<!-- /MarkdownTOC -->

# 1. Introduction

independent variable / explanatory variables
dependent variable / response variable


# 2. Linear Regression

## 2.1 Introduction 

Please read the [simple linear regression section](https://bradleyboehmke.github.io/HOML/linear-regression.html#simple-linear-regression) and the [multiple linear regression](https://bradleyboehmke.github.io/HOML/linear-regression.html#multi-lm) sections in [Bradley Boehmke "Hands-On Machine Learning with R" online book](https://bradleyboehmke.github.io/HOML/). 

## 2.2 Model performance

Under the usual assumptions stated above, an unbiased estimate of the error variance is given as the sum of the squared residuals divided by $$n − p$$ (where $$p$$
is the number of variables in the model):

$$\widehat{σ}^2 = \frac{1}{n−p} \sum_{i=1}^{n} r_i^2$$  

where 

$$r_{i} = (Y_{i} - Y_{i})$$

$$r_{i}$$ is referred to as the _ith_ residual (i.e., the difference between the _ith_ observed and predicted response value). The quantity $$\widehat{\sigma}^{2}$$
is also referred to as the mean square error (MSE) and its square root is denoted __RMSE for Squared Root Mean Squared Error__. In R, the RMSE of a linear model can be extracted using the `sigma()` function.

It is also called the residual sum of squares. In R, the `summary(fitted_regression_model)` will indicate this as the `Residual standard error`.

$$\min(SSE = \sum_{i=1}^{n} (y_{i} - \widehat{y_{i}}))$$


## 2.3 Diagnostic plots
Residuals
Systematic errors versus no pattern/systematic pattern
Gain curve?

## 2.4 Making predictions

Choose a few beta_values for sites X Y Z ... and report the predicted age. 

## 2.5 Limitations and possible solutions
Colinearity of independent variables

# 3. Regularised Linear Regression

## 3.1 The high-dimensionality curse of omic datasets

There are at least three reasons why "omics" limit linear regression approaches: 
- __Linear relationship:__ "off/on" switch when a gene expression or CpG methylation reaches a certain level.  
- __There are more observations (n) than features (p) ($$n > p$$):__ very often not true since "omics" will measure thousands of variables (genes, CpG sites) at the same time.
- __No or little multicollinearity:__ often not true since genes are co-regulated and methylation sites can be part of the same DNA region. 


Taken from Bradley Boehmke "Hands-On Machine Learning with R" book section on [regularised regressions](https://bradleyboehmke.github.io/HOML/regularized-regression.html): 	 
> Many real-life data sets, like those common to text mining and genomic studies are wide, meaning they contain a larger number of features ($$p > n$$). As $$p$$ increases, we’re more likely to violate some of the OLS [Ordinary Least Squares] assumptions and alternative approaches should be considered. [...] Having a large number of features invites additional issues in using classic regression models. For one, having a large number of features makes the model much less interpretable. Additionally, when $$p > n$$, there are many (in fact infinite) solutions to the OLS problem! 

This is our case as our dataset originally comprised 22,758 CpG sites (features) on _only_ 66 samples (observations). Even after we filtered for non-differential CpG sites in relation to age, we are still left with 228 CpG sites which is more than 66 (but much closer).

What we should do in this case is called __feature selection__ that will select a set of features (here CpG sites) from a bigger ensemble assuming that only this small set plays an actual role in the studied phenomenon (age here).

To quote Bradley Boehmke again: 
> In such cases, __it is useful (and practical) to assume that a smaller subset of the features exhibit the strongest effects__ (something called the bet on sparsity principle[...]. For this reason, we sometimes prefer estimation techniques that incorporate __feature selection__. One approach to this is called hard thresholding feature selection, which includes many of the traditional linear model selection approaches like forward selection and backward elimination. These procedures, however, can be computationally inefficient, do not scale well, and treat a feature as either in or out of the model (hence the name hard thresholding). In contrast, a more modern approach, called soft thresholding, slowly pushes the effects of irrelevant features toward zero, and in some cases, will zero out entire coefficients. As will be demonstrated, this can result in more accurate models that are also easier to interpret.
{: .quotation}


$$\min(SSE = \sum_{i=1}^{n} (y_{i} - \widehat{y_{i}}) + P)$$

where $$P$$ is a penalty parameter. 

There are three common penalty parameters we can implement:
- Ridge.
- LASSO.
- Elastic net (or ENET), which is a combination of ridge and lasso.

## 3.1 Ridge regression

:construction_worker:
<center>Section under construction</center>
:construction_worker: 

## 3.2 LASSO regression

The LASSO (least absolute shrinkage and selection operator) penalty (Tibshirani 1996) is a true feature selection method since some variable coefficients will be "pushed" to zero thereby eliminating them from the model.

Mathematically: 

$$\min(SSE + \lambda \sum_{j=1}^j \vert\beta_{j}\vert)$$

The $$\lambda$$ is called the _tuning parameter_ and determines how strong the penalty will affect the coefficient estimates.  
The $$\vert\beta_{j}\vert$$ term increases the penalty term every time a coefficient is added to the model. 

Here is an illustration (Figure 6.3) based on a LASSO regression from housing data:


<figure>
  <img src="../img/07-lasso-regression.png" style="width:50%">
  <figcaption>Coefficients gradually decrease to 0 when lambda increases</figcaption>
</figure> 




# 8. References

## 8.1 Online books and tutorials
- Bradley Boehmke Hands-on Machine Learning with R book sections on regressions:
  - [Linear regression](https://bradleyboehmke.github.io/HOML/linear-regression.html)
  - [Regularised regression](https://bradleyboehmke.github.io/HOML/regularized-regression.html)

## 8.2 Publications and books
- LASSO regression: Hastie, T., R. Tibshirani, and M. Wainwright. 2015. Statistical Learning with Sparsity: The Lasso and Generalizations. Chapman & Hall/Crc Monographs on Statistics & Applied Probability. Taylor & Francis.
- Gareth James, Daniela Witten, Trevor Hastie and Robert Tibshirani. 2013. An Introduction to Statistical Learning [Link](https://faculty.marshall.usc.edu/gareth-james/ISL/)
