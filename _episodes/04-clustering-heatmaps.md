---
title: "Clustering analysis and heatmaps"
teaching: 30
exercises: 45
questions:
- "How can I group genes with similar profiles together?"
- "How do I represent the result of a clustering analysis?"
objectives:
- "Understand how distances are calculated between two tissues based on their gene expression profile."
- "Be able to name a few different clustering methods."
keypoints:
- ""
---

# Table of contents
<!-- TOC -->

- [1. Table of contents](#1-table-of-contents)
- [2. Introduction](#2-introduction)
    - [2.1. The fastq format](#21-the-fastq-format)
- [3. Quality control of FASTQ files](#3-quality-control-of-fastq-files)
    - [3.1. Running FastQC](#31-running-fastqc)
    - [3.2. Viewing the FastQC results](#32-viewing-the-fastqc-results)
        - [3.2.1. Decoding the FastQC outputs](#321-decoding-the-fastqc-outputs)
        - [3.2.2. Sequencing error profiles](#322-sequencing-error-profiles)
        - [3.2.3. Expected Errors](#323-expected-errors)
        - [3.2.4. Worrisome](#324-worrisome)
        - [3.2.5. Quality assessment](#325-quality-assessment)
        - [3.2.6. Per sequence quality scores](#326-per-sequence-quality-scores)
        - [3.2.7. Per base sequence content](#327-per-base-sequence-content)
        - [3.2.8. Per sequence GC content](#328-per-sequence-gc-content)
        - [3.2.9. Sequence duplication level](#329-sequence-duplication-level)
        - [3.2.10. Overrepresented sequences](#3210-overrepresented-sequences)
    - [3.3. Working with the FastQC text output](#33-working-with-the-fastqc-text-output)
    - [3.4. Documenting our work](#34-documenting-our-work)
- [4. Trimming and filtering](#4-trimming-and-filtering)
- [5. Alignment to a reference genome](#5-alignment-to-a-reference-genome)
    - [5.1. Index the reference genome](#51-index-the-reference-genome)
    - [5.2. Align reads to reference genome](#52-align-reads-to-reference-genome)
    - [5.3. The SAM/BAM format](#53-the-sambam-format)
- [6. Creating the counts file](#6-creating-the-counts-file)
- [7. Removal of Container and Image](#7-removal-of-container-and-image)
- [8. Keypoints](#8-keypoints)

<!-- /TOC -->

# References
- [Exploring gene expression patterns using clustering](https://tavareshugo.github.io/data-carpentry-rnaseq/04b_rnaseq_clustering.html)
