# Master Advanced Forensics

This repository generates the lesson materials for the Master in Advanced Forensics based on the website template from [The Carpentries](https://carpentries.org/) Foundation. 

<!-- MarkdownTOC autolink="true" -->

- [Preview changes to the lesson locally](#preview-changes-to-the-lesson-locally)
- [Writing math symbols and formulas](#writing-math-symbols-and-formulas)
- [Credits](#credits)
	- [Maintainer\(s\)](#maintainers)
	- [Authors](#authors)
	- [Citation](#citation)
	- [Inspiration](#inspiration)

<!-- /MarkdownTOC -->

# Preview changes to the lesson locally
The lesson website is built through Github and Jekyll. 

__Option 1:__ follow the Carpentries setup: http://carpentries.github.io/lesson-example/setup.html 
The detailed instructions are listed in the "Jekyll Setup for Lesson Development" section.   

__Option 2:__ use a Docker container
1. Open a Shell window. 
2. Navigate to the `master-advanced-forensics/` folder using the `cd` command.
3. Since the lesson relies Jekyll 3.8.5, type within the Shell `export JEKYLL_VERSION=3.8.5`.
4. Make sure you have Docker for Windows or Mac installed: https://docs.docker.com/install/
5. With the Docker Desktop application running (you should see a little whale with containers at the top of your screen), type `docker run --rm --volume="$PWD:/srv/jekyll" -p 4000:4000 -it jekyll/jekyll:3.8.5 jekyll serve`  
6. Open a web browser and type `http://0.0.0.0:4000/` in the navigation bar. You should see the lesson website. Your changes should be automatically reflected online.  



# Writing math symbols and formulas
Sometimes, you need to render mathematical formula etc. GitHub pages which is using Jekyll to render html pages makes use of [MathJax])https://www.mathjax.org/). 

Visit the official website and consult [this thread](https://math.meta.stackexchange.com/questions/5020/mathjax-basic-tutorial-and-quick-reference) for quick ways to render math symbols.


# Credits

## Maintainer(s)

Current maintainers of this lesson are 

* Marc Galland, Data analyst and manager (University of Amsterdam, SILS, Plant Physiology Department).

## Authors

A list of contributors to the lesson can be found in [AUTHORS](AUTHORS)

## Citation

To cite this lesson, please consult with [CITATION](CITATION)

## Inspiration
...


