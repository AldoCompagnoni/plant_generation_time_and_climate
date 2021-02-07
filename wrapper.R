# required libraries
library(tidyverse)
library(gridExtra)
library(lme4)
library(nlme)
library(ggthemes)
library(testthat)
library(MASS)
library(ggpubr)
library(cowplot)
library(MuMIn)
library(bbmle)
library(lmerTest)

# create directories to store results in 
dir.create( 'results' )
dir.create( 'results/spp_spec' )
dir.create( 'results/spp_spec/prec' )
dir.create( 'results/spp_spec/temp' )

# run analysis script (this fills "results" with figures)
source( 'R/analysis_synthesis.R' )

# Write code for appendices
source( 'R/write_rmd_linear_figures.R' )

# "Knit" the Rmarkdown files for temperature figures
rmarkdown::render('R/tmp_bivariate_plots.Rmd',
                  output_dir  = getwd(),
                  output_file = 'tmp_bivariate_plots.pdf')

# "Knit" the Rmarkdown files for precipitation figures
rmarkdown::render('R/ppt_bivariate_plots.Rmd',
                  output_dir  = getwd(),
                  output_file = 'ppt_bivariate_plots.pdf')

