---
title: 
---


# Setup

<!-- MarkdownTOC autolink="True" levels="1,2,3" -->

- [Softwares \(R, RStudio and required R packages\)](#softwares-r-rstudio-and-required-r-packages)
	- [Option 1 \(preferred\): using a Docker image](#option-1-preferred-using-a-docker-image)
	- [Option 2: manual installation](#option-2-manual-installation)
- [Dataset #1 \(tissue gene expression\)](#dataset-1-tissue-gene-expression)
- [Dataset #2 \(DNA methylation and age\)](#dataset-2-dna-methylation-and-age)
	- [CpG site methylation values](#cpg-site-methylation-values)
	- [Sample to age correspondence](#sample-to-age-correspondence)
	- [CpG site to gene and related annotation](#cpg-site-to-gene-and-related-annotation)
- [References](#references)

<!-- /MarkdownTOC -->

# Softwares (R, RStudio and required R packages)

## Option 1 (preferred): using a Docker image

The preferred option to install all softwares and packages is to use a tailor-made Docker image. See [this nice introduction to Docker here](https://aws.amazon.com/docker/).   

There are two Docker images necessary to complete this RNA-seq lesson:
1. The command-line Docker `fastq-latest` image necessary to perform all bioinformatic analyses on the sequencing files: trimming, alignment and count table generation.
2. The RStudio Docker `rnaseq-latest` image necessary to perform all count-related analyses: EDA, differential expression and downstream functional analyses.   


So first thing first, we need to install Docker itself. 

> ## Install Docker
> Unfortunately, in many common situations installing Docker on your laptop will not straightforward if you do not have a large amount of technical experience. We have helpers on hand that have worked their way through the install process but be prepared for some troubleshooting.
> Please try to install the appropriate software from the list below depending on the operating system that your laptop is running:
> ### Microsoft Windows
> **You must have admin rights to run docker!** Some parts of the lesson will work without running as admin but if you are unable to `Run as admin` on your machine some of this workshop might not work easily.
> 
> If you have Windows 10 Pro Edition:
>  - First try to install the [Docker Desktop (Windows)](https://hub.docker.com/editions/community/docker-ce-desktop-windows), or **failing that**;
> - Install the [Docker Toolbox (Windows)](https://docs.docker.com/toolbox/toolbox_install_windows/).
>
> If you have Windows 10 Home Edition:
> - Install the [Docker Toolbox (Windows)](https://docs.docker.com/toolbox/toolbox_install_windows/).
>
> ### Apple macOS
> Either:
> - First, try to install the [Docker Desktop (Mac)](https://hub.docker.com/editions/community/docker-ce-desktop-mac), or **failing that**:
> - Install the [Docker Toolbox (Mac)](https://docs.docker.com/toolbox/toolbox_install_mac/).
> 
> ### Linux
> There are too many varieties of Linux to give precise instructions here, but hopefully you can locate documentation for getting Docker installed on your Linux distribution. It may already be installed. Note that Docker do list a number of versions of the Docker Engine for different Linux distributions [here](https://hub.docker.com/search/?type=edition&offering=community). 
>
> ### Troubleshooting
> Sometimes with git-bash and Windows, you can get issues listed here:   
> `the input device is not a TTY.  If you are using mintty, try prefixing the command with 'winpty'`. This can be troubleshooted following [this blog post](https://pitman.io/posts/tips-for-using-docker-from-git-bash-on-windows/).
{: .prereq}

> ## Before you start
>
> Before the training, please make sure you have done the following: 
>
> 1. First, install [Docker desktop](https://www.docker.com/products/docker-desktop) for your operating system. If not possible, install the Docker Toolbox (see above). 
> 2. If needed, install Shell Bash: [follow these instructions](http://swcarpentry.github.io/shell-novice/setup.html).
> 3. Open a new Shell Bash window and navigate to a folder that will be your workspace. For instance, you could create a folder named `rnaseq-tutorial/` on your Desktop and move inside with the Shell using `cd ~/Desktop/rnaseq-tutorial/`. 
> 4. In a Shell Bash window, type the following command: `docker run --rm --name rstudio_instance -v $PWD:/home/rstudio/ -e PASSWORD=mypwd -p 8787:8787 scienceparkstudygroup/master-advanced-forensics:latest`. This will download a Docker image for the course, create and run a container where RStudio will be running.   
> 4. Navigate to [http://localhost:8787](http://localhost:8787) in your web browser. You should have an RStudio session running. Type `rstudio` as the user name and `mypwd` as your password. 
> 5. To quit, close the web browser window where RStudio is running and exit the Shell too. 
{: .prereq}



> ## Important note
>
> You can save files to your disk when working inside the Docker-powered R session. You need to save them as you would normally. The files (e.g. `my_plot.png`) will be where you were working (the directory from which you launched the Docker container). 
>
{: .callout}


__Docker command-line explanations:__  
- The `--rm` removes the container when it has been run. No need to store it into your computer after use.      
- The `--name` gives a name to the running container for easy retrieval.  
- The `-p 8787:8787` follow the format `-p host_port:container_port`. Therefore the port 8787 inside the container will be exposed to the outside port on the host machine. That way, the running instance of RStudio can be access through the <IP address>:port format.

<br>

## Option 2: manual installation
This is the second way to install softwares and packages. It _should_ work but there is no guarantee that it _will_ work since R and packages versions on your machine might be different from the software and package versions used in this lesson. Thus, the preferred way is still to use the Docker image (option 1).  

> ## Before you start.
>
> Before the training, please make sure you have done the following: 
>
> 1. Download and install **up-to-date versions** of:
>    - R: [https://cloud.r-project.org](https://cloud.r-project.org).
>    - RStudio: [http://www.rstudio.com/download](http://www.rstudio.com/download). 
> 2. Read the workshop [Code of Conduct](https://docs.carpentries.org/topic_folders/policies/code-of-conduct.html) to make sure this workshop stays welcoming for everybody.
> 3. Get comfortable: if you're not in a physical workshop, be set up with two screens if possible. You will be following along in RStudio on your own computer while also following this tutorial on your own.
{: .prereq}

# Dataset #1 (tissue gene expression)

This dataset comes from the [Genotype-Tissue Expression (GTEx) database](https://www.gtexportal.org/home/) that gathers gene expression data from various human tissues. Specifically, the __GTEx Analysis v7__ version was used. 

One file contains the whole expression dataset while the other list 195 genes with an adipose-favored tissue expression. 

> ## Download the required datasets
> [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4068260.svg)](https://doi.org/10.5281/zenodo.4068260)
> 1. Go to the corresponding data record on [Zenodo](https://doi.org/10.5281/zenodo.4068260).
> 2. You will see two datasets that you need to download to your computer. 
> 3. Download the datasets to your machine and place it into a folder called `data/`. For instance, if you work on your Desktop, then place it in `~/Desktop/data/`
> 4. That's it!
{: .prereq}

It contains data from 53 different tissues for 56,202 human genes. The first 5 rows and 10 columns look like this:  

| gene_id           | Description  | Adipose - Subcutaneous | Adipose - Visceral (Omentum) | Adrenal Gland | Artery - Aorta | Artery - Coronary | Artery - Tibial | Bladder |
|-------------------|--------------|------------------------|------------------------------|---------------|----------------|-------------------|-----------------|---------|
| ENSG00000223972.4 | DDX11L1      | 0.056945               | 0.05054                      | 0.0746        | 0.03976        | 0.04386           | 0.04977         | 0.05878 |
| ENSG00000227232.4 | WASH7P       | 11.85                  | 9.753                        | 8.023         | 12.51          | 12.3              | 11.59           | 14.24   |
| ENSG00000243485.2 | MIR1302-11   | 0.06146                | 0.05959                      | 0.08179       | 0.04297        | 0.05848           | 0.05184         | 0.06097 |
| ENSG00000237613.2 | FAM138A      | 0.0386                 | 0.03245                      | 0.0405        | 0.02815        | 0.03678           | 0.03894         | 0.04113 |
| ENSG00000268020.2 | OR4G4P       | 0.035695               | 0                            | 0.03479       | 0              | 0                 | 0               | 0       |
| ENSG00000240361.1 | OR4G11P      | 0.04268                | 0.03988                      | 0.049065      | 0.03399        | 0                 | 0.04286         | 0       |
| ENSG00000186092.4 | OR4F5        | 0.05145                | 0.04558                      | 0.06136       | 0              | 0.04069           | 0.04669         | 0.05461 |
| ENSG00000238009.2 | RP11-34P13.7 | 0.1625                 | 0.1202                       | 0.087785      | 0.1351         | 0.1369            | 0.1472          | 0.143   |
| ENSG00000233750.3 | CICP27       | 0.1244                 | 0.1347                       | 0.1488        | 0.1026         | 0.1195            | 0.1145          | 0.0761  | 

<br>

# Dataset #2 (DNA methylation and age)

This dataset comes from the paper of [Bocklandt et al. 2011](https://doi.org/10.1371/journal.pone.0014821) on epigenetic prediction of human age.

> ## Download the required datasets
> [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4054842.svg)](https://doi.org/10.5281/zenodo.4054842)
> 1. Go to the corresponding data record on [Zenodo](https://doi.org/10.5281/zenodo.4054842).
> 2. You will see three datasets:
>    * __The methylation values of CpG sites:__        `cpg_methylation_beta_values.tsv`
>    * __The CpG site to gene functional annotation:__ `cpg_methylation_cpg_to_annotation.tsv`
>    * __The sample to age correspondence:__           `cpg_methylation_sample_age.tsv`
> 2. Download the three datasets to your machine and place them into a folder called `data/`. For instance, if you work on your Desktop, then place it in `~/Desktop/data/`. 
> 3. That's it!
{: .prereq}


## CpG site methylation values

This table contains the methylation so-called beta-values (between 0 and 1) of 27,578 CpG sites present on the [Illumina Methylation Assay](https://en.wikipedia.org/wiki/Illumina_Methylation_Assay) used to measure simultaneously 27,578 CpG dinucleotides from the Human genome. In total, 14,495 genes can have their CG site methylation state assayed.   

Here is a peek at the first 4 CpG sites and the first 10 samples. 

| CpG_id     | GSM712302  | GSM712303  | GSM712306  | GSM712307  | GSM712308  | GSM712309  | GSM712310  | GSM712311  | GSM712312 | GSM712313  |...|
|------------|------------|------------|------------|------------|------------|------------|------------|------------|-----------|------------|---|
| cg00000292 | 0.7901138  | 0.6837273  | 0.7662699  | 0.8209271  | 0.6333262  | 0.7478055  | 0.7912492  | 0.8109176  | 0.7846574 | 0.7107132  |...|
| cg00002426 | 0.6821175  | 0.4457058  | 0.5504405  | 0.7424242  | 0.1810345  | 0.524531   | 0.5488435  | 0.4585302  | 0.6980932 | 0.6316049  |...|
| cg00003994 | 0.07177814 | 0.07437458 | 0.07474366 | 0.07717327 | 0.07058448 | 0.06924484 | 0.07275518 | 0.07832465 | 0.0669353 | 0.07683486 |...|
| cg00005847 | 0.1930693  | 0.1557304  | 0.188825   | 0.1614361  | 0.1638656  | 0.2012987  | 0.1728916  | 0.1723153  | 0.1508261 | 0.1939234  |...|

... more lines ...

## Sample to age correspondence

This table contains the sample identifier (GSM712303) correspondence with the subject age (40). 
There are 66 samples with age ranging from 21 to 55. 

| sample    | title | source | pair | race   | age |
|-----------|-------|--------|------|--------|-----|
| GSM712302 | 111   | saliva | 1    | White  | 40  |
| GSM712303 | 112   | saliva | 1    | White  | 40  |
| GSM712306 | 811   | saliva | 8    | White  | 39  |
| GSM712307 | 812   | saliva | 8    | White  | 39  |
| GSM712308 | 911   | saliva | 9    | White  | 39  |
| GSM712309 | 912   | saliva | 9    | White  | 39  |
| GSM712310 | 1011  | saliva | 10   | White  | 52  |
| GSM712311 | 1012  | saliva | 10   | White  | 52  |
| GSM712312 | 1411  | saliva | 14   | White  | 27  |
| GSM712313 | 1412  | saliva | 14   | White  | 27  |
| GSM712314 | 2011  | saliva | 20   | White  | 35  |
| GSM712315 | 2012  | saliva | 20   | White  | 35  |
| GSM712316 | 2111  | saliva | 21   | Latino | 24  |
| GSM712317 | 2112  | saliva | 21   | Latino | 24  |

## CpG site to gene and related annotation 

This table links the CpG site identifier to a gene symbol and its functional annotation. It can be used to verify that CpG sites found to be 
important for age are indeed related to a known age-related gene. 

| CpG_id     | Symbol   | Description                                        |
|------------|----------|----------------------------------------------------|
| cg00000292 | ATP2A1   | ATPase; Ca++ transporting; fast twitch 1 isoform a |
| cg00002426 | SLMAP    | sarcolemma associated protein                      |
| cg00003994 | MEOX2    | mesenchyme homeo box 2                             |
| cg00005847 | HOXD3    | homeobox D3                                        |
| cg00006414 | ZNF398   | zinc finger 398 isoform b                          |
| cg00007981 | PANX1    | pannexin 1                                         |
| cg00008493 | COX8C    | cytochrome c oxidase subunit 8C                    |
| cg00008713 | IMPA2    | inositol(myo)-1(or 4)-monophosphatase 2            |
| cg00009407 | TTC8     | tetratricopeptide repeat domain 8 isoform B        |
| cg00010193 | FLJ35816 | hypothetical protein LOC401114                     |
| cg00011459 | PMM2     | phosphomannomutase 2                               |
| cg00012199 | RNASE4   | ribonuclease; RNase A family; 4 precursor          |
| cg00012386 | C1orf142 | hypothetical protein LOC116841                     |
| cg00012792 | TXNDC5   | thioredoxin domain containing 5 isoform 1          |

... more lines ...

# References

- Gene Expression Omnibus [accession GSE28746](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28746)
- Bocklandt S, Lin W, Sehl ME, Sánchez FJ, Sinsheimer JS, Horvath S, Vilain E. (2011) Epigenetic predictor of age. _PLoS One_. 6(6):e14821. [Link](https://doi.org/10.1371/journal.pone.0014821) 
