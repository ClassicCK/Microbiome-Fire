# Microbiome-Fire

Repository for Chapter Two of Christopher Kilner's Dissertation: "***Forest fire drives reduction in soil microbial diversity in Appalachian Forests***"

Authors: Christopher L. Kilner, Jean P. Gibert

Affiliation: Duke University, Department of Biology

## Land Acknowledgement

*Great Smokly Mountains National Park is located on the traditional and sacred land of the S\'atsoyaha (Yuchi) and ᏣᎳᎫᏪᏘᏱ Tsalaguwetiyi (Cherokee, East) peoples, whose history, language, traditional ways of life, and culture continue to influence this community. This land acknowledgement is an intentional act to counter the erasure of indigenous people and to demonstrate respect for their sovereign rights.*

Our affiliation, Duke University, provides the following land acknowledgment:

*What is now Durham was originally the territory of several Native nations, including Tutelo (TOO-tee-lo) and Saponi (suh-POE-nee) -- speaking peoples. Many of their communities were displaced or killed through war, disease, and colonial expansion.  Today, the Triangle is surrounded by contemporary Native nations, the descendants of Tutelo, Saponi, and other Indigenous peoples who survived early colonization. These nations include the Haliwa-Saponi (HALL-i-wa suh-POE-nee), Sappony (suh-POE-nee), and Occaneechi (oh-kuh-NEE-chee) Band of Saponi. North Carolina\'s Research Triangle is also home to a thriving urban Native American community who represent Native nations from across the United States. Together, these Indigenous nations and communities contribute to North Carolina\'s ranking as the state with the largest Native American population east of Oklahoma.*

## Abstract

Fire directly and indirectly impacts the global carbon cycle that drives climate change by altering biogeochemical processes and restructuring biological communities. Soil microbiomes can be impacted by fire and this impact may determine ecosystems recovery to this increasingly intense disturbance. Here, we examine National Ecological Observatory Network (NEON) 16S amplicon sequencing data from a forest fire in the Great Smoky Mountains National Park to examine the impact of fire on the soil microbial community. We report drastic fire-induced shifts in bacterial composition post-fire, with a reduction in alpha- and beta-diversity and no significant recovery of the soil microbiome 1-year post-fire. We also show that certain bacterial clades are clear indicators of fire, with heat-tolerant taxa increasingly dominant post-fire within our study plots. These findings show how forest fires in landscapes adapted to infrequent fire-regimes---such as moist, montane communities---may contain soil microbiomes that are less likely to to recover, a concern as climate change alters the climate of many regions globally and large-scale forest fires become increasingly common in areas that have not historically experienced them.

## Structure

This GitHub repository contains all the code to replicate our analysis, as well as the .tex files of our manuscript (to be added upon publication). We have separated this repository into folders containing: code, data, and outputs---figures and tables.

All data is open source and pulled from the [National Ecological Observatory Network's](https://www.neonscience.org) online data repository through [R](https://www.r-project.org).

### Code

Code is broken down into 4 scripts:

-   00_Download_Raw_Data.R Sets up a file structure and pulls down metadata from NEON for amplicon sequence analysis

-   01_Process_16S.R Pulls down amplicon sequence and co-located environmental data and formats this into a phyloseq-object for downstream analysis

-   02_Analysis.R contains all analysis performed in the manuscript

-   Mapping.R contains the code for Figure 1 (map) and generates the information on burn status of the plots

### Data

Data is formatted with raw, mid-process (amplicon sequence) and final analysis data (phyloseq-object with environment and phylogeny information). Data will be made available upon publication of the manuscript.

### Outputs

To ensure reproducibility and transparency in our work, we have provided the code to generate all Figures and Tables above in our 02_Analysis.R script and Mapping.R script, as well as the outputs of this code. Any changes between the outputs in this repository and the figures/tables included in the manuscript are aesthetic in nature (no changes or alteration to data plotted or included) and were done in Adobe Illustrator or BioRender.

#### Figures

Figures are named to match their place in the manuscript.

#### Tables

Tables are provided as .tex files, as these were uploaded to the .tex manuscript file. Any alterations to the .tex file in the manuscript are aesthetic unless otherwise noted in the code 02_Analysis.R

## Manuscript

Manuscript .tex files will be made available upon publication. We will also make available all peer review feedback if permitted by the journal of publication.
