
# colon_epimosaicism

## Embryonic origins of epigenetic mosaicism in the human colorectal epithelium

Hyein Won, Yunah Lee, Chang Hyun Nam, Beomki Lee, Joonoh Lim, Ji-Hyung Park, Jisong An, Soo A Oh, Myungsuk Choi, Jeong Yeon Kim, Jung Woo Park, Junehawk Lee, Jeonghwan Youk, Ji Won Park, Seung-Yong Jeong, Hyun Woo Kwon, Min Jung Kim, and Young Seok Ju.

On this page, you'll find R scripts used to generate the figures and perform the analyses presented in our paper.

## Requirements

- Operating system: Linux
- R version: 4.1.3
- Required R packages: `data.table`, `tidyverse`, `ggplot2`, `ape`, `phangorn`, `mclust`, `sva`, and other packages specified in individual scripts.

## Installation

The repository contains R scripts used for downstream analyses. No separate installation of custom software is required. Users should install R and the packages required by the relevant scripts.

## Input data

Input files for the analyses are provided in the `input/` directory. Additional sequencing data are available through the European Genome-phenome Archive (EGA; accession EGAS00001008522).

Raw sequencing data were processed using established bioinformatics tools, as described in the Methods section of the manuscript.

## Analysis scripts

The `script/` directory contains R scripts used for DNA methylation analysis, epigenotype assignment, phylogenetic inference, identification of embryonic epiallelic commitment (EEC) events, statistical analyses and figure generation.

## Demo

Input files are provided in the `input/` directory, with corresponding R scripts in the `script/` directory.

Expected outputs include processed data tables, statistical results and figures, depending on the analysis.

Execution times vary depending on the input dataset and available computing resources. Individual runtimes have not been separately benchmarked.

## Contact information

If you have any questions, please feel free to contact us at [hyein807@kaist.ac.kr](mailto:hyein807@kaist.ac.kr) or [yalee@kaist.ac.kr](mailto:yalee@kaist.ac.kr).
