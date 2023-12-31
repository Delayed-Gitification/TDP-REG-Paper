# TDP-REG-Paper

This repository contains code and small data files relevant to our manuscript (Wilkins et al) on TDP-REG.

# R Markdown

An R Markdown file is provided which contains the code used to make almost every figure in the paper. Wherever feasible, this markdown contains the full code required to go from the output of processing pipelines to publication-ready figures, to help reproducibility and transparency.

Note that it is not feasible to store all the required output files on Github itself. Instead, a folder containing all the data is hosted on Google Drive, and is available on request.

# Scripts

This folder contains scripts and pipelilnes necessary for the processing of raw sequencing data. The outputs of these pipelines are subsequently analysed and plotted using R.

Pipelines are written in Snakemake. Individual software packages are either common published software packages (e.g. Minimap2), or are custom software packages which have been released on Github (e.g. QIAxcelR). In a couple of cases, custom scripts are included in the scripts folder itself.

## Cas9 library

This folder contains two custom scripts and a Snakemake pipeline for analysing the library of synonymous variants.

## Nanopore

This folder contains Snakemake pipelines (and accompanying config files, plus one custom pysam script) for analysing all the Nanopore targeted sequencing data used in this study.

# small_data_files

This, as the name suggests, contains small data files necessary for generating figures.

