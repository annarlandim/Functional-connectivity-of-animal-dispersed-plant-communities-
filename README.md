# Functional-connectivity-of-animal-dispersed-plant-communities-
Data and code used in the paper: 

"Functional connectivity of animal-dispersed plant communities depends on the interacting effects of network specialization and resource diversity"

by: Anna R. Landim, Eike Lena Neuschulz, Isabel Donoso, Marjorie C. Sorensen, Thomas Mueller & Matthias Schleuning

# Code version 1

#########
SUMMARY
#########

Plant functional connectivity relies on seed dispersal between habitat patches, often performed by animals. We introduce a trait-based framework to estimate connectivity, integrating plant-frugivore networks, movement distance, and the selection of target patches by avian seed dispersers. We simulated how network specialization, dispersal distance, and resource diversity in a target patch affect the number and diversity of seeds dispersed. Results showed that connectivity generally decreased with network specialization, while resource diversity in the target patch increased both seed number and diversity, especially in specialized networks and at large distances. Our framework enhances the mechanistic understanding of plant functional connectivity and identifies processes limiting the connectivity of animal-dispersed plant communities.


#########
CONTENTS
#########

## Rscript
1) Functions_and_simulations.R - Basic R script for running the analysis.

## Rmarkdown
2) Functions_and_simulations.Rmd - R Markdown document detailing the analysis with explanations.

## HTML reports
3) Functions_and_simulations.html - HTML output of the R Markdown file for easy viewing in a web browser.

## CSV files
4) Manu_bird.community.csv - bird traits
5) Manu_plant.community.csv - plant traits

###########
DESCRIPTION
###########

Contained in this repository are all data files and computer code (in R language) used to reproduce the seed dispersal
simulations and calculations of plant functional connectivity. With the code it is also possible to reproduce the plots
in Fig.2 and Fig.3. 

An HTML report of the R script can be used to understand the steps of the simulations from used functions to plot making.

The data used is already publicly availabe at https://doi. org/10.5061/dryad.4f4qrfjdm.
The dataset included here is simply a streamlined version of the full dataset available in the link above.

######################
HOW TO RUN THE SCRIPTS
######################

To execute the analysis contained in this repository, begin by launching RStudio and opening the provided R script 
(Functions_and_simulations.R) and R Markdown file (Functions_and_simulations.Rmd). You can run the R script directly 
within RStudio by setting your working directory to the folder containing the script and data files, using the setwd()
function, and then executing the script using the source() command. This will perform the analysis as defined in the 
script. For a more detailed walkthrough, open the .Rmd file in RStudio and use the Knit button to compile it 
into an HTML document (analysis.html), which integrates both code and narrative text, providing a comprehensive view 
of the analysis steps and results. The HTML file is already compiled and can be viewed directly in any web browser, 
offering a convenient way to access the analyzed results without running the code.
