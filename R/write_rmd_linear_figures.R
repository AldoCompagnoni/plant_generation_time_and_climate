# Write down Rmarkdown files to produce linear plots

# read the files
ppt_plot_f <- list.files('results/spp_spec/prec')
tmp_plot_f <- list.files('results/spp_spec/temp')

# function to write the code for the Appendices with figures
write_it <- function( ii, file_v, clim_var ){
  
  paste0(
# '#',gsub('.png','',file_v[ii]) %>% gsub('_',' ',.),"\n",
"```{r, out.width='100%', fig.align='center', fig.cap='",
         gsub('.png','',file_v[ii]) %>% gsub('_',' ',.),
         "',echo=F}
knitr::include_graphics('C:/CODE/plant_generation_time_and_climate/results/spp_spec/",clim_var,"/",
file_v[ii],
          "')","\n",
          "```")
  
}

# headers ------------------------------------------------------------

ppt_header <- paste0(
"---
fontsize: 12pt
linespread: 2
mainfont: Arial
title: 'Graphs on linear responses to precipitation anomalies.'
knit: (function(inputFile, encoding) { 
      rmarkdown::render(inputFile,
                        encoding=encoding, 
                        output_file=file.path( paste0(getwd(),'ppt_linear_plots'))) })
output:
  pdf_document:
  fig_caption: true
number_sections: true
---
    
Each page of this appendix contains a scatter plot showing, for each one of the 162 populations analyzed in this study, log($\\lambda$) on the y-axis, and precipitation anomalies on the x-axis. The blue line is the mean prediction of the multiple linear model (Eq. 2 in the main text) calculated across the observed precipitation anomalies while keeping  temperature anomalies at the mean value observed during the study period. The shaded grey area is the 95% confidence interval around the mean prediction.
  
\\extrafloats{200}
\\maxdeadcycles=500 
\\tableofcontents
\\listoffigures
\\newpage",
"\n\n"
)

tmp_header <- paste0(
"---
fontsize: 12pt
linespread: 2
mainfont: Arial
title: 'Graphs on linear responses to temperature anomalies'
knit: (function(inputFile, encoding) { 
      rmarkdown::render(inputFile,
                        encoding=encoding, 
                        output_file=file.path( paste0(getwd(),'tmp_linear_plots'))) })
output:
  pdf_document:
  fig_caption: true
number_sections: true
---
  
  Each page of this appendix contains a scatter plot showing, for each one of the 162 populations analyzed in this study, log($\\lambda$) on the y-axis, and temperature anomalies on the x-axis. The blue line is the mean prediction of the multiple linear model (Eq. 2 in the main text) calculated across the observed temperature anomalies while keeping  precipitation anomalies at the mean value observed during the study period. The shaded grey area is the 95% confidence interval around the mean prediction.

\\extrafloats{200}
\\maxdeadcycles=500 
\\tableofcontents
\\listoffigures
\\newpage",
"\n\n"
)

# Write files -----------------------------------------------------------------------------

# Write markdown code
ppt_figs <- paste( lapply(1:162,write_it, ppt_plot_f, 'prec') %>% unlist, collapse = '\n\n')
tmp_figs <- paste( lapply(1:162,write_it, tmp_plot_f, 'temp') %>% unlist, collapse = '\n\n' )

# Rmarkdown files
writeLines( paste0(ppt_header, ppt_figs), con = 'R/ppt_appendix_figures.Rmd')
writeLines( paste0(tmp_header, tmp_figs), con = 'R/tmp_appendix_figures.Rmd')

