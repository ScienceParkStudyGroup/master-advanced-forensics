---
title: "Introduction"
teaching: 30
exercises: 0
questions:
- "What will I learn during this workshop?"
- "What are the tools that I will be using?"
- "What is Open Data Science?"
objectives:
- "Learn about Open Science."
- "Learn how to increase your data analysis efficacy"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---
# Overview {#overview}

Welcome. 

In this training you will learn R, RStudio, Git, and GitHub. It's going to be fun and empowering! You will learn a reproducible workflow that can be used in research and analyses of all kinds, including Ocean Health Index assessments. This is really powerful, cool stuff, and not just for data: this lesson website was actually made using some of these tools.

We will practice learning three main things all at the same time: coding with best practices (R/RStudio), collaborative version control (Git/GitHub), and communication/publishing (RMarkdown/GitHub). This training will teach these all together to reinforce skills and best practices, and get you comfortable with a workflow that you can use in your own projects. 

## What to expect

This is going to be a fun workshop. 

The plan is to expose you to a lot of great tools that you can have confidence using in your research. You'll be working hands-on and doing the same things on your own computer as we do live on up on the screen. We're going to go through a lot in these two days and it's less important that you remember it all. More importantly, you'll have experience with it and confidence that you can do it. The main thing to take away is that there *are* good ways to approach your analyses; we will teach you to expect that so you can find what you need and use it! A theme throughout is that tools exist and are being developed by real, and extraordinarily nice, people to meet you where you are and help you do what you need to do. If you expect and appreciate that, you will be more efficient in doing your awesome science.

You are all welcome here, please be respectful of one another. You are encouraged to help each other. 

Everyone in this workshop is coming from a different place with different experiences and expectations. But everyone will learn something new here, because there is so much innovation in the data science world. Instructors and helpers learn something new every time, from each other and from your questions. If you are already familiar with some of this material, focus on how we teach, and how you might teach it to others. Use these workshop materials not only as a reference in the future but also for talking points so you can communicate the importance of these tools to your communities. A big part of this training is not only for you to learn these skills, but for you to also teach others and increase the value and practice of open data science in science as a whole. 

## What you'll learn

- how to THINK about data 
    - how to think about data separately from your research questions
    - how and why to tidy data and analyze tidy data, rather than making your analyses accommodate messy data
    - how there is a lot of decision-making involved with data analysis, and a lot of creativity
- how to increase efficiency in your science
    - and increase reproducibility
    - and facilitate collaboration with others — especially Future You!
- how open science is a great benefit
    - find solutions faster
    - broaden the impact of your work
- how to learn with intention and community
    - think ahead instead of only to get a single job done now
    - the #rstats online community is fantastic. The tools we're using are developed by real people. Real, nice people. They are building powerful and empowering tools and are welcoming to all skill-levels


### Tidy data workflow

We will be learning about tidy data. And how to use a tidyverse suite of tools to work with tidy data.

[**Hadley Wickham**](http://hadley.nz/) and his team have developed a ton of the tools we'll use today. 
Here's an overview of techniques to be covered in Hadley Wickham and Garrett Grolemund of RStudio's book [R for Data Science](http://r4ds.had.co.nz/):

![](../img/r4ds_data-science.png)

We will be focusing on: 

- **Tidy**: `tidyr` to organize rows of data into unique values
- **Transform**: `dplyr` to manipulate/wrangle data based on subsetting by rows or columns, sorting and joining
- **Visualize**: 
    - `ggplot2` static plots, using grammar of graphics principles
- **Communicate**
    - dynamic documents with *R Markdown*
    
    
This is really critical. Instead of building your analyses around whatever (likely weird) format your data are in, take deliberate steps to make your data tidy. When your data are tidy, you can use a growing assortment of powerful analytical and visualization tools instead of inventing home-grown ways to accommodate your data. This will save you time since you aren't reinventing the wheel, and will make your work more clear and understandable to your collaborators (most importantly, Future You). 

    
## Learning with data that are not your own

One of the most important things you will learn is how to think about data separately from your own research context. Said in another way, you'll learn to distinguish your data questions from your research questions. Here, we are focusing on data questions, and we will use data that is not specific to your research.

We will be using several different data sets throughout this training, and will help you see the patterns and parallels to your own data, which will ultimately help you in your research.

## Emphasizing collaboration

Collaborating efficiently has historically been really hard to do. It's only been the last 20 years or so that we've moved beyond mailing things with the postal service. Being able to email and get feedback on files through track changes was a huge step forward, but it comes with a lot of bookkeeping and reproduciblity issues (did I do my analyses with `thesis_final_final.xls` or `thesis_final_usethisone.xls`?). But now, open tools make it much easier to collaborate. 

Working with collaborators in mind is critical for reproducibility. And, your most important collaborator is Future You. This training will introduce best practices using open tools, so that collaboration will become second nature to you!

## By the end of the course...

By the end of the course, you'll wrangle a few different data sets, and make your own graphics that you'll publish on webpages you've built collaboratively with GitHub and RMarkdown. Woop!

Here are some important things to keep in mind as you learn (these are joke book covers): 

![](../img/practical_dev_both.png)

## Prerequisites

Before the training, please make sure you have done the following: 

1. Download and install **up-to-date versions** of:
    - R: https://cloud.r-project.org
    - RStudio: http://www.rstudio.com/download 
    - Git: https://git-scm.com/downloads *Note: open the download and follow normal install procedures on your computer but you won't see any software installed when you're done*
1. Create a GitHub account: https://github.com *Note! Shorter names that kind of identify you are better, and use your work email!*
1. Get comfortable: if you're not in a physical workshop, be set up with two screens if possible. You will be following along in RStudio on your own computer while also watching a virtual training or following this tutorial on your own.

<!---
## Motivation 


More often than not, there are more than one way to do things. I'm going to focus mostly on what I have ended up using day-to-day; I try to incorporate better practices as I come upon them but that's not always the case. RStudio has some built-in redundancy too that I'll try to show you so that you can approach things in different ways and ease in.

- based on literature: best and good enough practices
- also based on our team's experience of how to do better science in less time




## Collaboration

Everything we learn today is to going to help you collaborate with your most important collaborator — YOU. Science is collaborative, starting with Future You, your current collaborators, and anyone wanting to build off your science later on. 

## Reproducibility

- record of your analyses. 
- rerun them!
- modify them, maybe change a threshold, try a different coefficient, etc, maybe today
- modify them, make a new figure, in 6 months! 

## Mindset

New but will become increasingly familiar. We’ll start you off with some momentum, like if you were going to learn to ride a bike or ...

Expect that there is a way to do what you want to do

- stop confounding data science with your science. Expect that someone has had your problem before or done what you want to do. 


If you plan to program mostly in one particular language on a single platform (such as Mac or Windows), you might try an integrated development environment (IDE). IDEs integrate text editing, syntax highlighting, version control, help, build tools, and debugging in one interface, simplifying development. 

http://r-bio.github.io/intro-git-rstudio/

## Data science is a discipline

It has theories, methods, and tools. 

Tidyverse and Hadley’s graphic. Tidy data.

Going to teach you how to think differently, get into some of the theory but in the context of hands-on work.


--->


## Credit

This material builds from a lot of fantastic materials developed by others in the open data science community. In particular, it pulls from the following resources, which are highly recommended for further learning and as resources later on. Specific lessons will also cite more resources.

- [R for Data Science](http://r4ds.had.co.nz/) by Hadley Wickham and Garrett Grolemund
- [STAT 545](http://stat545.com/) by Jenny Bryan
- [Happy Git with R](http://happygitwithr.com) by Jenny Bryan
- [Software Carpentry](https://software-carpentry.org/lessons/) by the Carpentries


{% include links.md %}
