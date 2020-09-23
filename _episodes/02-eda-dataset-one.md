---
title: "Exploratory Data Analysis of dataset #1"
teaching: 5
exercises: 30
questions:
- "What are the main explorary data analysis step to perform?"
- "What would you do to compare the different tissues all at once?"
- "How would you already get an idea of the similarity between tissues?"
objectives:
- "Reckon the value of EDA to get a better understanding of a dataset."
- "Identify potential issues in this dataset."
keypoints:
- "See how EDA can help further downstream data analysis."
---


# Table of Contents
<!-- MarkdownTOC autolink="true" levels="1,2" -->

- [What is Exploratory Data Analysis?](#what-is-exploratory-data-analysis)
- [Data import sanity checks](#data-import-sanity-checks)
- [Getting Descriptive statistics](#getting-descriptive-statistics)
  - [Yet, eventually, EDA is about getting a better understanding of the studied dataset and preparing its downstream analysis \(e.g. fitting a regression model\).](#yet-eventually-eda-is-about-getting-a-better-understanding-of-the-studied-dataset-and-preparing-its-downstream-analysis-eg-fitting-a-regression-model)
- [Describing relationships between the different variables](#describing-relationships-between-the-different-variables)
- [References](#references)
- [Photo credits](#photo-credits)

<!-- /MarkdownTOC -->


# What is Exploratory Data Analysis?
EDA (short for Exploratory Data Analysis) can have different purposes:
* Handle Missing values.
* Removing duplicates.
* Outlier Treatment.
* Normalizing and Scaling( Numerical Variables).
* Encoding Categorical variables( Dummy Variables).
* Bivariate Analysis.


Here, we are not going to do all of this. Instead, we will perform some data import checks, plotting some statistics and get an idea of how each tissue relates to each other. 

# Data import sanity checks

~~~
df_expr <- read.csv("data/U133AGNF1B_gcrma_avg.csv", 
                    header = TRUE, 
                    stringsAsFactors = FALSE,
                    check.names = FALSE)
~~~
{: .language-r}

You can already check the number of rows (genes) and columns (tissues/samples) using simple R functions:

~~~
# number of rows / columns
nrow(df_expr)
ncol(df_expr)
~~~
{: .language-r}

This should give you __44,775 genes__ (probes) and __85 tissues__. 

How did R "understand" your different data types? 
~~~
glimpse(df)
~~~
{: .language-r}

This command will output this in your console:
~~~
Rows: 44,775
Columns: 85
$ probe                              <chr> "1007_s_at", "1053_at", "117_at", "121_at", "1255_g_at", "1294_a…
$ `721_B_lymphoblasts`               <dbl> 137.00, 81.75, 12.55, 13.00, 5.45, 30.80, 9.75, 7.10, 120.85, 8.…
$ Adipocyte                          <dbl> 31.10, 7.35, 10.60, 14.55, 4.50, 11.40, 9.15, 6.45, 17.00, 7.80,…
$ AdrenalCortex                      <dbl> 53.30, 8.90, 14.00, 14.10, 5.00, 13.30, 13.00, 8.10, 20.80, 10.9…
$ Adrenalgland                       <dbl> 81.20, 7.20, 9.85, 13.70, 3.95, 9.85, 7.95, 5.75, 21.70, 7.35, 3…
$ Amygdala                           <dbl> 489.35, 8.00, 10.15, 14.30, 4.85, 9.90, 9.65, 6.50, 19.95, 8.20,…
$ Appendix                           <dbl> 39.50, 13.15, 11.55, 16.45, 5.10, 11.60, 11.30, 7.70, 31.25, 10.…
$ AtrioventricularNode               <dbl> 20.45, 8.85, 33.60, 15.00, 3.75, 11.45, 11.35, 6.05, 18.95, 9.00…
$ `BDCA4+_DentriticCells`            <dbl> 17.10, 10.25, 36.15, 17.10, 4.90, 48.10, 9.10, 6.25, 45.80, 7.65…
...(more lines)
~~~
{: .output}

> ## Question
> 1. Did R convert the "probe" column to a suitable data type?
> 2. Did R convert all the other columns to a suitable data type?
>
> > ## Solution
> > 1. Yes, R has converted the "probe" column into the `chr` data type which stands for character. This is because we specify `stringsAsFactors = FALSE` in the 
> > `read.csv()` function at the beginning. 
> > 2. Yes, all tissue columns with gene expression measurements have been converted to the `dbl` data type which stands for "double class" (a double-precision floating point number). In short, it is a number with decimals. 
> {: .solution}
{: .challenge}

~~~
head(n = 5)
~~~
{: .language-r}

For big matrix, do `df_expr[1:5,1:5]` to show the first five lines and five columns for instance. 

# Getting Descriptive statistics


~~~
head(n = 5)
~~~
{: .language-r}

##
Yet, eventually, EDA is about getting a better understanding of the studied dataset and preparing its downstream analysis (e.g. fitting a regression model). 

> ## Question
> Calculate a series of descriptive metrics on the expression value for all tissues: 
  - minimumn 
   - maximum 
   - average 
   - median  
>  
> Hint: make your dataset tidy first and use the `%>%` operator to chain operations on your dataframe. 
>
> > ## Solution
> > ~~~
> > df_expr %>% 
> >   pivot_longer(cols = - probe, 
> >                names_to = "tissue", 
> >                values_to = "expression") %>% 
> >   summarise(max = max(expression), 
> >             min = min(expression), 
> >             average = mean(expression), 
> >             median = median(expression)
> >             )
> > ~~~~
> > {: .language-r}
> > 
> > This gives you the following results:
> > ~~~~
> > # A tibble: 1 x 4
> >      max   min average median
> >    <dbl> <dbl>   <dbl>  <dbl>
> > 1 46508.  0.55    80.0    6.9
> >~~ 
> > {: .output}
> {: .solution}
{: .challenge}



# Describing relationships between the different variables
    Analyzing relationships between variables


# References 
- [The power analysis section of the RNA-seq blog](https://www.rna-seqblog.com/tag/power-analysis/)
- [`pwr` R package vignette](https://cran.r-project.org/web/packages/pwr/vignettes/pwr-vignette.html)
- [The Scotty power analysis webtool](http://scotty.genetics.utah.edu/)
- [UCLA Stat consulting on power analysis](https://stats.idre.ucla.edu/r/dae/power-analysis-for-two-group-independent-sample-t-test/)
-  [Lei et al. (2015) Diminishing returns in next-generation sequencing (NGS) transcriptome data. _Gene_ 557(1):82-87](https://doi.org/10.1016/j.gene.2014.12.013)

# Photo credits

<a style="background-color:black;color:white;text-decoration:none;padding:4px 6px;font-family:-apple-system, BlinkMacSystemFont, &quot;San Francisco&quot;, &quot;Helvetica Neue&quot;, Helvetica, Ubuntu, Roboto, Noto, &quot;Segoe UI&quot;, Arial, sans-serif;font-size:12px;font-weight:bold;line-height:1.2;display:inline-block;border-radius:3px" href="https://unsplash.com/@ayahya09?utm_medium=referral&amp;utm_campaign=photographer-credit&amp;utm_content=creditBadge" target="_blank" rel="noopener noreferrer" title="Download free do whatever you want high-resolution photos from Ali Yahya"><span style="display:inline-block;padding:2px 3px"><svg xmlns="http://www.w3.org/2000/svg" style="height:12px;width:auto;position:relative;vertical-align:middle;top:-2px;fill:white" viewBox="0 0 32 32"><title>unsplash-logo</title><path d="M10 9V0h12v9H10zm12 5h10v18H0V14h10v9h12v-9z"></path></svg></span><span style="display:inline-block;padding:2px 3px">Ali Yahya</span></a>

