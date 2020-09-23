---
title: "Introduction to dataset #1 (gene expressions)"
teaching: 30
exercises: 0
questions:
- "Where does this dataset come from?"
- "How many genes will I consider in my analysis?"
- "How many tissues do I have in my dataset?"
- "In which unit are my gene expression measured?"
objectives:
- "Understand what this dataset comprises."
- "Be able to explain what the values mean and how they were obtained."
keypoints:
- "This dataset measures gene expression from various human tissues."
- "Gene expression is measured from the hybridization of mRNA molecules to microarray probes."
- "Some tissues have genes that are uniquely or strongly expressed in them which makes them gene-markers of that tissue."
- "Finding tissue-specific markers can be done through several methods: PCA, clustering or even custom ones."
---


## Table of Contents
1. [Overview](#overview)
2. [What you will learn](#what-you-will-learn)
3. [Dataset used](#dataset-used)


## What you will learn

1. **What are the important things to know before doing an RNA-Seq experiment** 
    - When should you perform a RNA-Seq experiment?  
    - RNA-Seq experiments have to comply with good experimental design practices just like any experiment.
    - What are biological replicates and why are they important?
2. **How can I assess the quality of my RNA-Seq sequencing results?**
    - FastQC.
    - PCA plot.
    - Sample clustering.
3. **How do I perform a differential expression analysis on RNA-Seq results using R?**
    - Raw and scaled counts: why do you need to scale counts between samples?
    - What are the gene expression units I need to know: RPKM, FPKM, TPM.
    - What are robust scaling/normalisation methods?
    - How does the DESeq method works?
4. **What are the plots that I can make from the differential analysis results?**
    - Heatmap coupled with gene and sample clustering.
    - Volcano plot.

## Dataset used 

The dataset used 

The original sequencing files can be found on the [Array Express database of the European Bioinformatic Institute](https://www.ebi.ac.uk/arrayexpress) by searching for the dataset accession number __E‐MTAB‐4683__.

> ## Question
> 1. How were gene expression measured? 
> HintL 
> > ## Solution
> > 1. A $$log2$$ equal to 1 means that gene X has a higher expression (x2, two-fold) in the DC3000 infected condition compared to the mock condition. 
> > 2. A $$log2$$ equal to -1 means that gene X has a smaller expression ($$\frac{1}{2}$$) in the DC3000 infected condition.   
> >  
> > ~~~
> > untreated = 230
> > treated = 750
> > log2(treated/untreated) # equals 1.705257
> > ~~~
> > {: .language-r}
> {: .solution}
{: .challenge}

## Credits

### Dataset
The original RNA-Seq dataset used comes from Vogel et al. 2016:  https://doi.org/10.1111/nph.14036.  

### Teaching materials
This lesson has been formatted according to the [Carpentries Foundation](https://carpentries.org/) lesson template and following their recommendations on how to teach researchers good practices in programming and data analysis.   

This material builds from a lot of fantastic materials developed by others in the open data science community. Most of the content derives from the [Harvard Chan Bioinformatics Core](https://github.com/hbctraining) which are greatly acknowledge for the quality of their teaching materials.

{% include links.md %}
