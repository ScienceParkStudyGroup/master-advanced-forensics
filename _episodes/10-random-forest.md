---
title: "Random Forest analysis"
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
- [2. Decision trees](#2-decision-trees)
- [2.1 Tree vocabulary](#21-tree-vocabulary)
- [3. From a single tree to a forest](#3-from-a-single-tree-to-a-forest)
	- [3.1 Bagged trees](#31-bagged-trees)
	- [3.2 Random Forests](#32-random-forests)
- [References](#references)

<!-- /MarkdownTOC -->

# 1. Introduction

## 1.1 Definitions

Random Forest

Forest comes from Decision Trees

# 2. Decision trees

# 2.1 Tree vocabulary
Node
Leaf
Terminal Node
Branch
Pruning = how far to you grow the tree?

Regression tree

Metric for split quality = Residual Sum of Squares = how homogeneous are the observations split at decision node X. 

# 3. From a single tree to a forest

## 3.1 Bagged trees
Bagging = bootstrap aggregation = reduce the variance due to random dataset split by combining the results of multiple decision trees.
You work on all predictors/variables at once (all CpG sites).
You generate N bootstrapped training datasets (using random sampling).
No tree pruning 


## 3.2 Random Forests
Similar to bagged trees but Random Forest = bagging on random subsets of variables and samples. 
You don't work on all predictors at once.

Since not all of toe observations/samples are used to build the tree, the Y predictions based on the train set are compared to their actual real values using the left-over data ("out of bag"). Formally, this is called the _Out-Of-Bag_ error estimate (OOB in short). 

# References
- [Data Camp decision tree guide](https://www.datacamp.com/community/tutorials/decision-trees-R)
- University of Cincinnati Business Analytics tutorials on [regression trees](https://uc-r.github.io/regression_trees) and [Random Forests](https://uc-r.github.io/random_forests).

