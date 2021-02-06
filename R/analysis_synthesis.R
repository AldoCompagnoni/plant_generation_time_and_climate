# Carry out plant review tests using Plant review + synthesis
rm(list=ls())
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
options(stringsAsFactors = F)
  

# read data --------------------------------------------------------------

# mean climate variables
clim_m <- read.csv( 'data/map_mat_wai.csv', 
                    stringsAsFactors = F) 

# demographic data and climatic anomalies
all_df <- read.csv( 'data/lambdas_vs_anomalies.csv',
                    stringsAsFactors = F ) 

# generation time for synthesis species
genT_S <- read.csv( 'data/generation_times.csv',
                    stringsAsFactors = F )
  

# H1: precipitation over temperature ------------------------------------

# fit the "big" model using all data
all_mod <- lmer( log_lambda ~ ppt_t0*tmp_t0 + (1 | spp_pop) + (0 + ppt_t0*tmp_t0 | spp_pop),
                 data = all_df )

# Fixed effects
fix_boxmod <- fixef( all_mod )
fix_ci     <- confint( all_mod )


# BIOME ANALYSES: Auto-regressive terms -------------------------------------------------------

# add squared terms
all_df        <- all_df %>% 
                   mutate( ppt_t02 = ppt_t0^2,
                           tmp_t02 = tmp_t0^2 )

# define populations
spp_pop_df    <- all_df %>% 
                  dplyr::select( SpeciesAuthor, 
                                 MatrixPopulation ) %>% 
                  unique %>% 
                  subset( !(SpeciesAuthor    == "Primula_farinosa_4" &
                            MatrixPopulation == "Dröstorp 2") )


mod_l <- list( NA, NA, NA, NA, NA ) %>% 
          setNames( c('pt0',    'pt', 
                      'pt0_p2', 'pt0_t2',
                      'pt0_p2t2' ) )


# mean function to estimate response to climate
coef_by_spp <- function( ii ){
  
  print(ii)
  
  # select the data
  an_df    <- all_df %>% 
        subset( SpeciesAuthor    == spp_pop_df$SpeciesAuthor[ii] ) %>% 
        subset( MatrixPopulation == spp_pop_df$MatrixPopulation[ii] )
  
  spp_name <- spp_pop_df$SpeciesAuthor[ii]
  spp_pop  <- paste0( unique(an_df$SpeciesAuthor), '_',
                      unique(an_df$MatrixPopulation) )
  pop      <- unique(an_df$MatrixPopulation)
  
  # print combination of spp and pop!
  print( unique(an_df$spp_pop) )
  
  yr_rep <- an_df$MatrixStartYear %>% unique %>% length
  
  if( nrow(an_df) == yr_rep | 
      # species with only one double year; not possible to fit RE model
      spp_pop_df$SpeciesAuthor[ii] %in% c('Agropyron_cristatum',
                                          'Cryptantha_flava') |
      unique(an_df$spp_pop) == 'Eryngium_alpinum_2_Fournel'
  ){
    
    mod_p      <- gls( log_lambda ~ ppt_t0,
                       correlation = corARMA( p=1, q=0),
                       data = an_df )
    
    mod_t      <- gls( log_lambda ~ tmp_t0,
                       correlation = corARMA( p=1, q=0),
                       data = an_df )
    
    mod_pt0    <- gls( log_lambda ~ ppt_t0 + tmp_t0,
                       correlation = corARMA( p=1, q=0),
                       data = an_df )
    
    # only fit nonlinear effects if yr_rep > 14
    if( yr_rep > 14 ){
    
      mod_pt     <- gls( log_lambda ~ ppt_t0 * tmp_t0,
                         correlation = corARMA( p=1, q=0),
                         data = an_df )
      
      mod_pt0_p2 <- gls( log_lambda ~ ppt_t0 + tmp_t0 + ppt_t02,
                         correlation = corARMA( p=1, q=0 ),
                         data = an_df )
      
      mod_pt0_t2 <- gls( log_lambda ~ ppt_t0 + tmp_t0 + tmp_t02,
                         correlation = corARMA( p=1, q=0 ),
                         data = an_df )
      
      mod_pt0_p2t2 <- gls( log_lambda ~ ppt_t0 + tmp_t0 + 
                             tmp_t02 + ppt_t02,
                           correlation = corARMA( p=1, q=0),
                           data = an_df )
    
      # model selection
      mod_l   <- list( mod_pt0,    mod_pt, 
                       mod_pt0_p2, mod_pt0_t2, 
                       mod_pt0_p2t2) %>% 
                  setNames( c('pt0',    'pt', 
                              'pt0_p2', 'pt0_t2',
                              'pt0_p2t2' ) )
      
      ms_df   <- bbmle::AICctab(mod_l, nobs=nrow(an_df), 
                                weights=T,
                                base = T,
                                delta = F,
                                sort = F) %>% 
                   as.data.frame %>% 
                   tibble::add_column( model = c('pt0',    'pt', 
                                                 'pt0_p2', 'pt0_t2',
                                                 'pt0_p2t2' ),
                                       .before=1 )
      
    }
    
    # Output file
    out <- data.frame( 
      
      # Coefficients
      a_pt0     = coef(mod_pt0)[1],
      ppt_pt0   = coef(mod_pt0)[2],
      tmp_pt0   = coef(mod_pt0)[3],
      
      # Coefficient standard deviations
      a_pt0_sd    = summary(mod_pt0)$tTable[2,'Std.Error'],
      ppt_pt0_sd  = summary(mod_pt0)$tTable[2,'Std.Error'],
      tmp_pt0_sd  = summary(mod_pt0)$tTable[3,'Std.Error'],
      
      # metadata
      SpeciesAuthor = spp_name,
      spp_pop       = unique(an_df$spp_pop),
      graminoid     = unique(an_df$graminoid),
      orgtype       = unique(an_df$OrganismType),
      yr_rep        = yr_rep,
      
      stringsAsFactors = F 
      
    )
    
  }else{
    
    mod_p   <- gls( log_lambda ~ ppt_t0,
                    correlation = corARMA( p=1, q=0,
                                           form = ~ 1 | MatrixReplicate),
                    data = an_df )
    
    mod_t   <- gls( log_lambda ~ tmp_t0,
                    correlation = corARMA( p=1, q=0,
                                           form = ~ 1 | MatrixReplicate),
                    data = an_df )
    
    mod_pt0 <- gls( log_lambda ~ ppt_t0 + tmp_t0,
                  correlation = corARMA( p=1, q=0,
                                         form = ~ 1 | MatrixReplicate),
                  data = an_df )
    
    # only fit nonlinear effects if yr_rep > 14
    if( yr_rep > 14){
      
      mod_pt  <- gls( log_lambda ~ ppt_t0 * tmp_t0,
                      correlation = corARMA( p=1, q=0,
                                             form = ~ 1 | MatrixReplicate),
                      data = an_df )
    
      mod_pt0_p2 <- gls( log_lambda ~ ppt_t0 + tmp_t0 + ppt_t02,
                         correlation = corARMA( p=1, q=0,
                                                form = ~ 1 | MatrixReplicate),
                         data = an_df )
      
      mod_pt0_t2 <- gls( log_lambda ~ ppt_t0 + tmp_t0 + tmp_t02,
                         correlation = corARMA( p=1, q=0,
                                                form = ~ 1 | MatrixReplicate),
                         data = an_df )
      
      mod_pt0_p2t2 <- gls( log_lambda ~ ppt_t0 + tmp_t0 + 
                                        tmp_t02 + ppt_t02,
                           correlation = corARMA( p=1, q=0,
                                                  form = ~ 1 | MatrixReplicate),
                           data = an_df )
      
      # model selection
      mod_l   <- list( mod_pt0,    mod_pt, 
                       mod_pt0_p2, mod_pt0_t2, 
                       mod_pt0_p2t2) %>% 
                  setNames( c('pt0',    'pt', 
                              'pt0_p2', 'pt0_t2',
                              'pt0_p2t2' ) )
      
      ms_df   <- bbmle::AICctab(mod_l, nobs=nrow(an_df), 
                                weights=T,
                                base = T,
                                delta = F,
                                sort = F) %>% 
                   as.data.frame %>% 
                   tibble::add_column( model = c('pt0',    'pt', 
                                                 'pt0_p2', 'pt0_t2',
                                                 'pt0_p2t2' ),
                                       .before=1 )
    }

    # output file
    out      <- data.frame(
      
      # Coefficients
      a_pt0     = coef(mod_pt0)[1],
      ppt_pt0   = coef(mod_pt0)[2],
      tmp_pt0   = coef(mod_pt0)[3],
      
      # Coefficient standard deviations
      a_pt0_sd    = summary(mod_pt0)$tTable[2,'Std.Error'],
      ppt_pt0_sd  = summary(mod_pt0)$tTable[2,'Std.Error'],
      tmp_pt0_sd  = summary(mod_pt0)$tTable[3,'Std.Error'],
      
      # metadata
      SpeciesAuthor = spp_name,
      spp_pop       = unique(an_df$spp_pop),
      graminoid     = unique(an_df$graminoid),
      orgtype       = unique(an_df$OrganismType),
      yr_rep        = yr_rep,
      
      stringsAsFactors = F 
    )
    
    if( yr_rep > 14){
      out <- out %>% 
        mutate( 
                # coefficients for NL models
                ppt_p2    = coef(mod_pt0_p2)[2],
                tmp_p2    = coef(mod_pt0_p2)[3],
                ppt2_p2   = coef(mod_pt0_p2)[4],
                
                ppt_t2    = coef(mod_pt0_t2)[2],
                tmp_t2    = coef(mod_pt0_t2)[3],
                tmp2_t2   = coef(mod_pt0_t2)[4], 
          
                ppt_p2t2  = coef(mod_pt0_p2t2)[2],
                tmp_p2t2  = coef(mod_pt0_p2t2)[3],
                tmp2_p2t2 = coef(mod_pt0_p2t2)[4],
                ppt2_p2t2 = coef(mod_pt0_p2t2)[5],
                
                # SE of coefficients 
                ppt_p2_sd     = summary(mod_pt0_p2)$tTable[2,'Std.Error'],
                tmp_p2_sd     = summary(mod_pt0_p2)$tTable[3,'Std.Error'],
                ppt2_p2_sd    = summary(mod_pt0_p2)$tTable[4,'Std.Error'],
                
                ppt_t2_sd     = summary(mod_pt0_t2)$tTable[2,'Std.Error'],
                tmp_t2_sd     = summary(mod_pt0_t2)$tTable[3,'Std.Error'],
                tmp2_t2_sd    = summary(mod_pt0_t2)$tTable[4,'Std.Error'],
                
                ppt_p2t2_sd   = summary(mod_pt0_p2t2)$tTable[2,'Std.Error'],
                tmp_p2t2_sd   = summary(mod_pt0_p2t2)$tTable[3,'Std.Error'],
                tmp2_p2t2_sd  = summary(mod_pt0_p2t2)$tTable[4,'Std.Error'],
                ppt2_p2t2_sd  = summary(mod_pt0_p2t2)$tTable[5,'Std.Error'],
                
                # model selection              
                aic_pt0       = subset(ms_df, model == 'pt0')$AICc,
                aic_pt        = subset(ms_df, model == 'pt')$AICc,
                aic_pt0_p2    = subset(ms_df, model == 'pt0_p2')$AICc,
                aic_pt0_t2    = subset(ms_df, model == 'pt0_t2')$AICc,
                aic_pt0_p2t2  = subset(ms_df, model == 'pt0_p2t2')$AICc,
                
                aic_pt0_w     = subset(ms_df, model == 'pt0')$weight,
                aic_pt_w      = subset(ms_df, model == 'pt')$weight,
                aic_pt0_p2_w  = subset(ms_df, model == 'pt0_p2')$weight,
                aic_pt0_t2_w  = subset(ms_df, model == 'pt0_t2')$weight,
                aic_pt0_p2t2_w= subset(ms_df, model == 'pt0_p2t2')$weight )
    }
    
  }
  
  # calculate confidence interval from scratch
  confint_calc <- function( pred_name, param_name ){
    
    # select the parameters we need!
    par_v  <- c('ppt_pt0', 'tmp_pt0')
    var_v  <- c('ppt_t0',  'tmp_t0')
    b0     <- dplyr::select(out, a_pt0 ) %>% as.numeric
    b1     <- dplyr::select(out, param_name ) %>% as.numeric
    b2     <- dplyr::select(out, setdiff(par_v, param_name) ) %>% as.numeric
    
    # other variables
    n      <- nrow(an_df)
    x      <- an_df[,pred_name]
    cov_x  <- an_df[,setdiff(var_v, pred_name)] %>% mean
    pred_x <- an_df[,pred_name] %>% sort
    
    # Predict y at the given value of x (argument pred.x); 
    # at the mean value of the covariate (cov_c * cov_x)
    pred_y <- b0 + (b1 * pred_x) + (b2 * cov_x)
    
    # Find SSE and MSE
    sse <- sum((an_df$log_lambda - mod_pt0$fitted)^2)
    mse <- sse / (n - 3)
    
    # Critical value of t
    t.val <- qt(0.975, n - 3) 
    
    # Standard error of the mean estimate  
    mean.se.fit <- (1 / n + (pred_x - mean(x))^2 / (sum((x - mean(x))^2))) 
    
    # Mean Estimate Upper and Lower Confidence limits at 95% Confidence
    mean.conf.upper <- pred_y + t.val * sqrt(mse * mean.se.fit)
    mean.conf.lower <- pred_y - t.val * sqrt(mse * mean.se.fit)
    
    # the confidence intervals
    data.frame( yhat  = pred_y,
                x_seq = pred_x,
                upr   = mean.conf.upper,
                lwr   = mean.conf.lower )
    
  }
  
  ci_p <- confint_calc( 'ppt_t0', "ppt_pt0")
  ci_t <- confint_calc( 'tmp_t0', "tmp_pt0")
  
  # Temperature plot
  ggplot(an_df) +
    geom_point(  aes( tmp_t0, log_lambda ) ) +
    geom_line( data = ci_t,
               aes( x = x_seq,
                    y = yhat ),
               color = '#0072B2',
               lwd   = 2 ) +
    geom_ribbon( data = ci_t,
                 aes( x    = x_seq,
                      ymin = lwr,
                      ymax = upr ),
                 alpha = 0.5,
                 fill = 'grey' ) +
    theme_minimal( ) +
    labs( x = "Temperature anomaly",
          y = expression('log('*lambda*')'),
          title = gsub('_',' ',spp_name),
          subtitle = paste0('Population: ',pop) ) +
    theme( axis.title = element_text( size = 18 ),
           axis.text  = element_text( size = 15 ),
           plot.title    = element_text( size = 20, hjust=0.5),
           plot.subtitle = element_text( size = 18, hjust=0.5) ) + 
    ggsave( paste0('results/spp_spec/temp/', spp_pop,'.png'),
            width=6.3, height=6.3 )
  
  # Precipitation plot
  ggplot(an_df) +
    geom_point(  aes( ppt_t0, log_lambda ) ) +
    geom_line( data = ci_p,
               aes( x = x_seq,
                    y = yhat ),
               color = '#0072B2',
               lwd   = 2 ) +
    geom_ribbon( data = ci_p,
                 aes( x    = x_seq,
                      ymin = lwr,
                      ymax = upr ),
                 alpha = 0.5,
                 fill = 'grey' ) +
    theme_minimal( ) +
    labs( x = "Precipitation anomaly",
          y = expression('log('*lambda*')'),
          title = gsub('_',' ',spp_name),
          subtitle = paste0('Population: ',pop) ) +
    theme( axis.title = element_text( size = 18 ),
           axis.text  = element_text( size = 15 ),
           plot.title    = element_text( size = 20, hjust=0.5),
           plot.subtitle = element_text( size = 18, hjust=0.5) ) + 
    ggsave( paste0('results/spp_spec/prec/', spp_pop,'.png'),
            width=6.3, height=6.3 )
  
  out 
  
}

# Store single pop coefficients 
ar_coef   <- lapply( 1:nrow(spp_pop_df), coef_by_spp) %>% 
               bind_rows %>% 
               left_join( clim_m ) 

# species for which the linear model has been selected
spp_linear <- ar_coef %>% 
                subset( yr_rep > 14 ) %>% 
                subset( aic_pt0_w < aic_pt0_p2_w |
                        aic_pt0_w < aic_pt0_t2_w |
                        aic_pt0_w < aic_pt0_p2t2_w ) %>% 
                .$SpeciesAuthor %>% 
                unique

# coefficients for the linear models
coef_pt0   <- ar_coef %>% 
                subset( !(SpeciesAuthor %in% spp_linear) ) %>% 
                # Temperature and precipitation coefficients of linear models
                mutate( ppt_b  = ppt_pt0,
                        tmp_b  = tmp_pt0,
                        ppt_sd = ppt_pt0_sd,
                        tmp_sd = tmp_pt0_sd,
                        ppt_w  = 1 / ppt_pt0_sd,
                        tmp_w  = 1 / tmp_pt0_sd )

# quadratic temperature models
coef_t2   <- ar_coef %>% 
               subset( SpeciesAuthor %in% c('Bouteloua_curtipendula',
                                            'Bouteloua_rothrockii',
                                            'Hesperostipa_comata',
                                            'Pseudoroegneria_spicata') ) %>% 
               mutate( ppt_b  = ppt_t2,
                       tmp_b  = tmp_t2 + tmp2_t2,
                       ppt_sd = ppt_t2_sd,
                       tmp_sd = tmp_t2_sd + tmp2_t2_sd,
                       ppt_w  = 1 / ppt_p2_sd ) %>% 
               mutate( tmp_w  = 1 / tmp_sd )

# quadratic precipitation models
coef_p2   <- ar_coef %>% 
              subset( SpeciesAuthor %in% c('Schizachyrium_scoparium',
                                           'Artemisia_tripartita') ) %>% 
              mutate( ppt_b  = ppt_p2 + ppt2_p2,
                      tmp_b  = tmp_p2,
                      ppt_sd = ppt_p2_sd + ppt2_p2_sd,
                      tmp_sd = tmp_p2_sd,
                      tmp_w  = 1 / tmp_p2_sd ) %>% 
              mutate( ppt_w  = 1 / ppt_sd )

# Double quadratic: temperature + precipitation models
coef_p2t2 <- ar_coef %>% 
              subset( SpeciesAuthor %in% c('Poa_secunda',
                                           'Bouteloua_eriopoda') ) %>% 
              mutate( ppt_b  = ppt_p2t2 + ppt2_p2t2,
                      tmp_b  = tmp_p2t2 + tmp2_p2t2,
                      ppt_sd = ppt_p2t2_sd + ppt2_p2t2_sd,
                      tmp_sd = tmp_p2t2_sd + tmp2_p2t2_sd ) %>% 
              mutate( tmp_w  = 1 / tmp_sd,  
                      ppt_w  = 1 / ppt_sd )

# coefficients for the hypothesis tests
coef_h <- list( coef_pt0, coef_p2,
                coef_t2,  coef_p2t2 ) %>% bind_rows
           

# Adjust for COVARIATES -------------------------------------

# Covariates BY POPULATION
cov_by_pop <- function(ii){
  
  print(ii)
  
  # select the data
  an_df    <- all_df %>% 
    subset( SpeciesAuthor    == spp_pop_df$SpeciesAuthor[ii] ) %>% 
    subset( MatrixPopulation == spp_pop_df$MatrixPopulation[ii] )
  
  spp_name <- spp_pop_df$SpeciesAuthor[ii]
  
  # print combination of spp and pop!
  print( unique(an_df$spp_pop) )
  
  yr_rep <- an_df$MatrixStartYear %>% unique %>% length
  
  if( nrow(an_df) <= (yr_rep*2) |
      spp_name == 'Eryngium_alpinum_2' 
  ){
    
    if( unique(all(an_df$cov_n > 0) & !is.na(an_df$cov_n)) ){
      
      # if predictor is continuous, transform cov in numeric
      if( all( grepl( 'continuous',an_df$type)) ){
        an_df$cov <- as.numeric(an_df$cov) 
        an_df$cov <- (an_df$cov - mean(an_df$cov,na.rm=T)) / (sd(an_df$cov,na.rm=T))
      }
      
      mod0   <- gls( log_lambda ~ ppt_t0 + tmp_t0, correlation = corARMA( p=1, q=0),
                     data = an_df )
      modi_t <- gls( log_lambda ~ ppt_t0 + tmp_t0*cov, correlation = corARMA( p=1, q=0),
                     data = an_df )
      modi_p <- gls( log_lambda ~ ppt_t0*cov + tmp_t0, correlation = corARMA( p=1, q=0),
                     data = an_df )
      
      # determine best model with covariate
      aic_i  <- ( AICc(modi_p, nobs=nrow(an_df)) < AICc(mod0,nobs=nrow(an_df)) |
                    AICc(modi_t, nobs=nrow(an_df)) < AICc(mod0,nobs=nrow(an_df)) )
      
      # is the best model interacting with precipitation or temperature?
      if( aic_i ){
        if( AICc(modi_p, nobs=nrow(an_df)) < AICc(mod0,nobs=nrow(an_df)) ) best_m <- 'p'
        if( AICc(modi_t, nobs=nrow(an_df)) < AICc(mod0,nobs=nrow(an_df)) ) best_m <- 't'
      }else{
        best_m <- NA
      }
      
      out      <- data.frame( # year 1 model
        best_aic   = aic_i %>% as.numeric,
        best_m     = best_m,
        levels     = an_df$levels %>% unique,
        SpeciesAuthor = spp_name,
        yr_rep        = yr_rep,
        type          = an_df$type %>% unique,
        id_pop        = ii,
        stringsAsFactors = F )
    }
    
  }else{
    
    if( all(an_df$cov_n > 0) & all(!is.na(an_df$cov_n)) ){
      
      # if predictor is continuous, transform cov in numeric
      if( all( an_df$type == 'continuous' ) ){
        an_df$cov <- as.numeric(an_df$cov) 
        an_df$cov <- (an_df$cov - mean(an_df$cov,na.rm=T)) / (sd(an_df$cov,na.rm=T)*2)
      }
      
      mod0     <- gls(log_lambda ~ ppt_t0 + tmp_t0,
                      correlation = corARMA( p=1, q=0, form = ~ 1 | MatrixReplicate), 
                      data = an_df )
      modi_p   <- gls(log_lambda ~ ppt_t0*cov + tmp_t0,
                      correlation = corARMA( p=1, q=0, form = ~ 1 | MatrixReplicate), 
                      data = an_df )
      modi_t   <- gls(log_lambda ~ ppt_t0 + tmp_t0*cov,
                      correlation = corARMA( p=1, q=0, form = ~ 1 | MatrixReplicate), 
                      data = an_df )
      
      # determine best model with covariate
      aic_i  <- ( AICc(modi_p, nobs=nrow(an_df)) < AICc(mod0,nobs=nrow(an_df)) |
                    AICc(modi_t, nobs=nrow(an_df)) < AICc(mod0,nobs=nrow(an_df)) )
      
      # is the best model interacting with precipitation or temperature?
      if( aic_i ){
        if( AICc(modi_p, nobs=nrow(an_df)) < AICc(mod0,nobs=nrow(an_df)) ) best_m <- 'p'
        if( AICc(modi_t, nobs=nrow(an_df)) < AICc(mod0,nobs=nrow(an_df)) ) best_m <- 't'
      }else{
        best_m <- NA
      }
      
      out      <- data.frame( # year 1 model
        best_aic   = aic_i %>% as.numeric,
        best_m     = best_m,
        levels     = an_df$levels %>% unique,
        SpeciesAuthor = spp_name,
        yr_rep        = yr_rep,
        type          = an_df$type %>% unique,
        id_pop        = ii,
        stringsAsFactors = F )
      
    }
    
  }    
  
  out %>% 
    mutate( pop_spp_i = an_df$spp_pop %>% unique )
  
}

# ids: covariates fit at the population level (covariate == MatrixTreatment)
cov_by_pop_i <- all_df %>% 
  subset( !is.na(cov_n) ) %>% 
  # Remove Dicerandra_frutescens' "no fire" plots
  subset( !(SpeciesAuthor == "Dicerandra_frutescens" & cov == 'no fire') ) %>%
  # Remove Dicerandra_frutescens' population which has missing data
  subset( !(SpeciesAuthor    == "Dicerandra_frutescens" & 
            MatrixPopulation == 'Oak-hickory scrub 1 (Pop 10)') ) %>%
  dplyr::select( SpeciesAuthor, MatrixPopulation, MatrixTreatment, cov ) %>% 
  unique %>% 
  count( SpeciesAuthor, MatrixPopulation ) %>% 
  as.data.frame %>% 
  # select only studies with more than one treatment
  subset( n > 1 ) %>% 
  right_join( spp_pop_df ) %>% 
  mutate( i = !is.na(n) ) %>% 
  .$i %>% 
  which

# model selection
cov_mod_sel_df <- lapply( cov_by_pop_i, cov_by_pop) %>% 
  bind_rows %>% 
  subset( best_aic == 1 )


# fit models separately (UGLY, but these are just 2 species)

# "Astragalus_cremnophylax_var._cremnophylax_Grand Canyon National Park"
ii <- 12

# select the data
an_df    <- all_df %>% subset( spp_pop == "Astragalus_cremnophylax_var._cremnophylax_Grand Canyon National Park")
 
# format covariate 
an_df$cov <- as.numeric(an_df$cov) 
an_df$cov <- (an_df$cov - mean(an_df$cov,na.rm=T)) / (sd(an_df$cov,na.rm=T))
modi_t    <- gls( log_lambda ~ ppt_t0 + tmp_t0*cov, 
                  correlation = corARMA( p=1, q=0),
                  data = an_df )  

# calculate differences (for percentage in the text)
abs(coef(modi_t)['tmp_t0:cov']) / abs(coef(modi_t)['tmp_t0']) 
abs(coef(modi_t)['tmp_t0:cov']) + abs(coef(modi_t)['tmp_t0']) 

# data frame to update coefficients
update_df <- data.frame( 
                ppt_b  = coef(modi_t)['ppt_t0'], 
                tmp_b  = coef(modi_t)['tmp_t0'] + coef(modi_t)['tmp_t0:cov'], 
                ppt_sd = summary(modi_t) %>% .$tTable %>% .[2,2], 
                tmp_sd = summary(modi_t) %>% .$tTable %>% .[3,2] +
                         summary(modi_t) %>% .$tTable %>% .[5,2] 
                )

# update coefficients for meta-analyses
row_id <- which( coef_h$SpeciesAuthor == "Astragalus_cremnophylax_var._cremnophylax" )
col_id <- which( names(coef_h) %in% c("ppt_b","tmp_b","ppt_sd","tmp_sd") )

# update 
coef_h[row_id,col_id] = update_df


# "Dicerandra_frutescens_Fireline Edge 1 (Pop 0)"
ii <- 62  # spp_pop_df[62,]

# select the data
an_df    <- all_df %>% subset( spp_pop == "Dicerandra_frutescens_Fireline Edge 1 (Pop 0)" )

# format covariate 
an_df$cov <- as.numeric(an_df$cov) 
an_df$cov <- (an_df$cov - mean(an_df$cov,na.rm=T)) / (sd(an_df$cov,na.rm=T))
modi_t    <- gls( log_lambda ~ ppt_t0 + tmp_t0*cov, 
                  correlation = corARMA( p=1, q=0),
                  data = an_df )  

# calculate differences (for percentage in the text)
abs(coef(modi_t)['tmp_t0:cov']) / abs(coef(modi_t)['tmp_t0']) 
abs(coef(modi_t)['tmp_t0:cov']) + abs(coef(modi_t)['tmp_t0']) 
abs(coef(modi_t)['tmp_t0']) 

# data frame to update coefficients
update_df <- data.frame( 
  ppt_b  = coef(modi_t)['ppt_t0'], 
  tmp_b  = coef(modi_t)['tmp_t0'] + coef(modi_t)['tmp_t0:cov'], 
  ppt_sd = summary(modi_t) %>% .$tTable %>% .[2,2], 
  tmp_sd = summary(modi_t) %>% .$tTable %>% .[3,2] +
           summary(modi_t) %>% .$tTable %>% .[5,2] 
)

# update coefficients for meta-analyses
row_id <- which( coef_h$SpeciesAuthor == "Dicerandra_frutescens_Fireline Edge 1 (Pop 0)" )
col_id <- which( names(coef_h) %in% c("ppt_b","tmp_b","ppt_sd","tmp_sd") )

# update 
coef_h[row_id,col_id] = update_df

        

# GRAMINOID vs. NONGRAMINOID --------------------------------

# update the organism
coef_h <- coef_h %>% 
            mutate( orgtype = replace(orgtype, 
                                      graminoid == 'yes',
                                      'Graminoid') ) %>% 
            mutate( orgtype = gsub(" perennial", '', orgtype) ) %>% 
            mutate( orgtype = gsub("Shrub", 'Woody', orgtype) ) %>% 
            mutate( orgtype = gsub("Tree",  'Woody', orgtype) ) 
  
# number of data points per organism type
org_n  <- count(coef_h, orgtype) %>% 
            mutate( y = -1) %>% 
            mutate( n = paste0('N=',n) )

# ORGANISM TYPE
ppt_orgtype <- coef_h %>% 
  ggplot( aes( x = orgtype, y = ppt_b) ) +
  geom_boxplot(  outlier.alpha = 0.00001 ) +
  geom_jitter( aes( x = orgtype, y = ppt_b), 
               width = 0.3,
               alpha = 0.3) +
  geom_text( data  = mutate(org_n,y = -0.5),
             aes(x = orgtype, y = y, label = n) ) +
  labs( y = expression("Precipitation effect"),
        x = NULL ) +
  ylim(-1,1) +
  theme_minimal() +
  theme( axis.title = element_text( size = 12),
         axis.text.x = element_text( size = 12, angle = 80) ) 

tmp_orgtype <- coef_h %>% 
  ggplot( aes( x = orgtype, y = tmp_b) ) +
  geom_boxplot(  outlier.alpha = 0.00001 ) +
  geom_jitter( aes( x = orgtype, y = tmp_b), 
               width = 0.3,
               alpha = 0.3) +
  geom_text( data = org_n,
             aes(x=orgtype, y=y, label = n) ) +
  labs( y = expression("Temperature effect"),
        x = NULL ) +
  theme_minimal() +
  theme( axis.title = element_text( size = 12),
         axis.text.x = element_text( size = 12, angle = 80) ) 
  
p_orgtype <- plot_grid( ppt_orgtype, tmp_orgtype,
                        labels = 'AUTO',
                        label_size = 18,
                        align='h', 
                        nrow = 2, ncol = 1 )

# clim effect versus LHS with A-D labels
ggsave('results/eff_size_by_orgtype.png', 
       p_orgtype,
       dpi=600, width=3.15, height=6.3)

# ANOVA tests
lm( ppt_b ~ orgtype, data = coef_h ) %>% anova
lm( tmp_b ~ orgtype, data = coef_h ) %>% anova

# post hoc comparisons
aov( ppt_b ~ orgtype, data = coef_h ) %>% 
  TukeyHSD %>% 
  .$orgtype %>% 
  as.data.frame() %>% 
  tibble::add_column( Comparison = row.names(.), .before = 1 ) %>% 
  write.csv('results/tukey_prec.csv',
            row.names=F)

aov( tmp_b ~ orgtype, data = coef_h ) %>% 
  TukeyHSD %>% 
  .$orgtype %>% 
  as.data.frame %>% 
  tibble::add_column( Comparison = row.names(.), .before = 1 ) %>% 
  write.csv('results/tukey_temp.csv',
            row.names=F)


# BOOTSTRAP H2/H3 ---------------------------------------------------------

# do all bootstrap analyses
bootstrap_linear <- function(par_df, clim_var ){
  
  if( clim_var == 'prec' ){
    par_df <- par_df %>% 
      dplyr::select( ppt_b, ppt_sd, wai) %>% 
      rename( b_t0    = ppt_b, 
              b_t0_sd = ppt_sd ) %>% 
      # dplyr::select( b_t0, b_t0_sd, wai) %>% 
      rename( clim_mean = wai ) %>% 
      drop_na
  }
  
  if( clim_var == 'tmp' ){
    par_df <- par_df %>% 
      dplyr::select( tmp_b, tmp_sd, mat) %>% 
      rename( b_t0    = tmp_b, 
              b_t0_sd = tmp_sd ) %>% 
      rename( clim_mean = mat ) %>% 
      drop_na
  }
  
  # number of populations and associated SD
  n_spp_pop     <- nrow(par_df)
  varcov        <- matrix(0, n_spp_pop, n_spp_pop)
  diag(varcov)  <- par_df$b_t0_sd^2
  
  set.seed( 1584 )
  
  # simulate beta values
  beta_sim      <- mvrnorm(1000, 
                           mu    = par_df$b_t0,
                           Sigma = varcov )
  
  # bootstrap linear regressions
  bootstrap_lm <- function(ii){
    
    mod <- lm(beta_sim[ii,] ~ par_df$clim_mean )
    
    data.frame( int  = coef(mod)[1],
                beta = coef(mod)[2],
                b_se = mod %>% summary %>% .$coefficients %>% .[2,2],
                p    = mod %>% summary %>% .$coefficients %>% .[2,4] )
    
  }
  
  # store bootstrapped betas
  beta_bstr    <- lapply(1:nrow(beta_sim), bootstrap_lm) %>% 
                    bind_rows
  
  # bootstrapped confidence intervals
  
  # sequence of x values (check if temp./precip. info)
  x_seq <- seq( min(par_df$clim_mean, na.rm=T),
                max(par_df$clim_mean, na.rm=T),
                length.out = 100 )
  
  # get bootstrappe confidence intervals
  bootstrap_ci <- function(ii){
    beta_bstr[ii,]$int + x_seq * beta_bstr[ii,]$beta
  }
  
  # bootstrapped regression lines
  bs_mat <- sapply(1:nrow(beta_bstr), bootstrap_ci)
  
  # bootstrapped confidence intervals
  ci_df  <- data.frame( 
    x      = x_seq,
    ci_upr = apply(bs_mat,1,quantile,prob=0.025),
    ci_lwr = apply(bs_mat,1,quantile,prob=0.975)
  )
  
  # meta-confidence intervals
  meta_ci <- function( ii ){
    set.seed( ii )
    rnorm(1000, beta_bstr$beta[ii],
                beta_bstr$b_se[ii] )
  }
  
  bs_met <- lapply(1:nrow(beta_bstr), meta_ci) %>% unlist
  
  list( beta_bstr  = beta_bstr,
        ci_df      = ci_df,
        ci_meta    = bs_met )
  
}

# bootstrapped values
bs_ppt     <- bootstrap_linear( coef_h, 'prec' )
bs_tm      <- bootstrap_linear( coef_h, 'tmp' )
bs_ppt_nog <- bootstrap_linear( subset(coef_h, orgtype != 'Graminoid' ), 
                               'prec' )
# mean coefficients
m_prec     <- lm( ppt_b ~ wai, data = coef_h)
m_tm       <- lm( tmp_b ~ mat, data = coef_h)
m_prec_nog <- lm( ppt_b ~ wai, data = subset(coef_h, orgtype != 'Graminoid' ) )


# 90% Marginally significant results for prec_eff ~ MAP
coef(m_prec)[2]
sum( bs_ppt$ci_meta < 0 ) / 1000000
bs_ppt$ci_meta %>% quantile(prob=c(0.025,0.975))

# Non-significant results for temp_eff ~ MAT
coef(m_tm)[2]
sum( bs_tm$ci_meta > 0 ) / 1000000
bs_tm$ci_meta %>% quantile(prob=c(0.025,0.975))
bs_tm$beta_bstr$beta %>% quantile(prob=c(0.025,0.975))

# No Graminoids: 74% betas are below 0
coef(m_prec_nog)[2]
sum( bs_ppt_nog$ci_meta < 0 ) / 1000000
bs_ppt_nog$ci_meta %>% quantile(prob=c(0.025,0.975))


# plot of all data
theme_set(theme_minimal()+
            theme(legend.title=element_blank(),
                  legend.text = element_text(size = 14),
                  legend.justification=c(1,1), legend.position=c(1,1),
                  text = element_text(size=25), #family = "Times"
                  axis.text.x  = element_text(size=25),
                  axis.text.y  = element_text(vjust=0.5, size=14),
                  axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))
            ))


# plots with all data 
p1 <- coef_h %>% 
  dplyr::select( ppt_b, wai, ppt_sd, orgtype) %>% 
  rename( b_t0         = ppt_b,
          b_t0_sd      = ppt_sd,
          `Plant type` = orgtype ) %>% 
  drop_na %>% 
  mutate(b_t0_w  = 1 / (b_t0_sd) ) %>% 
  ggplot() +
  geom_point( aes(wai, b_t0,
                  size  = b_t0_w,
                  color = `Plant type` ),
              alpha = 0.5 ) +
  geom_abline( intercept = coef(m_prec)[1],
               slope     = coef(m_prec)[2],
               lwd = 1.5,
               lty = 2 ) +
  geom_ribbon( data=bs_ppt$ci_df,
               aes(x, 
                   ymin = ci_lwr,
                   ymax = ci_upr),
               color = 'grey',
               fill  = 'grey',
               alpha = 0.5,
  ) +
  ylab( 'Effect of precipitation' ) +
  xlab( 'Water Availability Index (mm)') +
  geom_hline( yintercept = 0, lty=2) + 
  theme_minimal() + 
  guides( size = F ) +
  theme(axis.title = element_text(size=12) ) +
  scale_colour_colorblind() 

p2 <- coef_h %>% 
  dplyr::select( tmp_b, mat, tmp_sd, orgtype) %>% 
  rename( b_t0    = tmp_b,
          b_t0_sd = tmp_sd ) %>% 
  drop_na %>% 
  mutate(b_t0_w  = 1 / (b_t0_sd) ) %>% 
  ggplot() +
  geom_point( aes(x    = mat, 
                  y    = b_t0,
                  size = b_t0_w,
                  color = orgtype),
              alpha = 0.5 ) +
  geom_abline( intercept = coef(m_tm)[1],
               slope     = coef(m_tm)[2],
               lwd = 1.5,
               lty = 2 ) +
  geom_ribbon( data=bs_tm$ci_df,
               aes(x,
                   ymin = ci_lwr,
                   ymax = ci_upr),
               color = 'grey',
               fill  = 'grey',
               alpha = 0.5,
  ) +
  ylab( 'Effect of temperature' ) +
  xlab( 'Mean annual temperature (°C)') +
  geom_hline( yintercept = 0, lty=2) + 
  theme_minimal() + 
  guides( size = F ) +
  theme(axis.title = element_text(size=12) ) +
  scale_colour_colorblind() 

# TRICK to maintain color palette without graminoids
col_nog <- c("#000000", "#56B4E9", "#009E73", "#F0E442")

# make your own colour palette
scale_colour_discrete <- function(...) {
  scale_colour_manual(..., values = col_nog)
}

# plots without graminoids data 
p1_nog <- coef_h %>% 
  dplyr::select( ppt_b, wai, ppt_sd, orgtype) %>% 
  subset( orgtype != 'Graminoid' ) %>% 
  rename( b_t0         = ppt_b,
          b_t0_sd      = ppt_sd,
          `Plant type` = orgtype ) %>% 
  drop_na %>% 
  mutate(b_t0_w  = 1 / (b_t0_sd) ) %>% 
  ggplot() +
  geom_point( aes(wai, b_t0,
                  size  = b_t0_w,
                  color = `Plant type` ),
              alpha = 0.5 ) +
  geom_abline( intercept = coef(m_prec_nog)[1],
               slope     = coef(m_prec_nog)[2],
               lwd = 1.5,
               lty = 2 ) +
  geom_ribbon( data=bs_ppt_nog$ci_df,
               aes(x, 
                   ymin = ci_lwr,
                   ymax = ci_upr),
               color = 'grey',
               fill  = 'grey',
               alpha = 0.5,
  ) +
  ylab( 'Effect of precipitation' ) +
  xlab( 'Water Availability Index (mm)') +
  geom_hline( yintercept = 0, lty=2) + 
  theme_minimal() + 
  guides( size = F ) +
  theme(axis.title = element_text(size=12) ) +
  scale_colour_discrete() 

# put plots together
shared_legend_cowplot <- function( x, y ){
  
  h_gt     <- plot_grid( x + theme(legend.position="none"),
                         y + theme(legend.position="none"),
                         labels = 'AUTO',
                         label_size = 18,
                         align='h', 
                         nrow = 2, ncol = 1 )
  
  # extract the legend from one of the plots
  legend <- cowplot::get_legend(
    # create some space to the left of the legend
    x + theme(legend.box.margin = margin(0, 0, 0, 4))
  )
  
  # output the new plots
  plot_grid(h_gt, legend, rel_widths = c(3, 1) )
  
}

# 2 panel plot
h_12           <- shared_legend_cowplot(p1, p2)

# store graphs
ggsave('results/climEff_vs_meanClim_bs_WAI.png',
       h_12,
       dpi=600, width=4.5, height=6.3 )

ggsave('results/precEff_vs_meanClim_bs_WAI_no_graminoid.png',
       p1_nog,
       dpi=600, width=6.3, height=6.3 )


# H4: GENERATION-TIME ANALYSES -------------------------------------------------------


# precipitation-related data
par_prec <- coef_h %>% 
  mutate( b_t0    = ppt_b,
          b_t0_sd = ppt_sd ) %>% 
  mutate(b_t0_w  = 1 / (b_t0_sd) ) %>% 
  left_join( genT_S ) %>% 
  dplyr::select( b_t0, b_t0_w, b_t0_sd, GenT, graminoid, orgtype ) %>% 
  gather( lhs, lhs_val, GenT ) %>% 
  mutate( clim_var = 'Precipitation',
          b_t0_abs = abs(b_t0) )

# temperature-related data
par_tm <- coef_h %>% 
  # tmp_coef %>% 
  mutate( b_t0    = tmp_b,
          b_t0_sd = tmp_sd ) %>% 
  mutate(b_t0_w  = 1 / (b_t0_sd) ) %>% 
  left_join( genT_S ) %>% 
  dplyr::select( b_t0, b_t0_w, b_t0_sd, GenT, graminoid, orgtype ) %>% 
  gather( lhs, lhs_val, GenT ) %>% 
  mutate( clim_var = 'Temperature',
          b_t0_abs = abs(b_t0) )

# data frame for plotting
pars_all <- bind_rows(par_prec, par_tm ) %>% 
  mutate( lhs_val = log(lhs_val) ) %>% 
  mutate( lhs = replace(lhs,
                        lhs == 'GenT',
                        'Generation time') ) 

# do all bootstrap analyses
bootstrap_gamma <- function(pars_all, clim_type, lhs_type){
  
  # select relevant data
  par_df        <- pars_all %>% 
                    subset( clim_var == clim_type ) %>% 
                    subset( lhs == lhs_type ) %>% 
                    drop_na
  
  # number of populations and associated SD
  n_spp_pop     <- nrow(par_df)
  varcov        <- matrix(0, n_spp_pop, n_spp_pop)
  diag(varcov)  <- (par_df$b_t0_sd^2)
  
  set.seed( 1997 )
  
  # simulate beta values
  beta_sim      <- mvrnorm( 1000, 
                            mu    = par_df$b_t0,
                            Sigma = varcov ) %>% abs
  
  # mean prediction
  m_mod  <- glm(abs(par_df$b_t0) ~ par_df$lhs_val,
                family  = Gamma(link = "log") )
  
  # bootstrap linear regressions
  bootstrap_glm <- function(ii){
    
    # check if this is temp. or precip. info
    mod <- glm(beta_sim[ii,] ~ par_df$lhs_val,
               family  = Gamma(link = "log") )
    
    # here, use intercept of mean model!
    data.frame( int  = coef(mod)[1],
                beta = coef(mod)[2],
                b_se = mod %>% summary %>% .$coefficients %>% .[2,2],
                p    = mod %>% summary %>% .$coefficients %>% .[2,4] )
    
  }
  
  # store bootstrapped betas
  beta_bstr    <- lapply(1:nrow(beta_sim), 
                         bootstrap_glm) %>% bind_rows
  
  # bootstrapped confidence intervals
  
  # sequence of x values (check if temp./precip. info)
  x_seq <- seq( min(par_df$lhs_val, na.rm=T),
                max(par_df$lhs_val, na.rm=T),
                length.out = 100 )
  
  # get bootstrappe confidence intervals
  bootstrap_ci <- function(ii){
    exp( beta_bstr[ii,]$int + x_seq * beta_bstr[ii,]$beta )
  }
  
  # bootstrapped regression lines
  bs_mat  <- sapply(1:nrow(beta_bstr), bootstrap_ci)
  
  # get mean of bootstrapped values, 
  bs_mean <- exp( mean(beta_bstr$int) + x_seq * mean(beta_bstr$beta) )
  ci_upr  <- apply(bs_mat,1,quantile,prob=0.975)
  ci_lwr  <- apply(bs_mat,1,quantile,prob=0.025)
  
  # upper and lower deviation from the mean of the bootstrap
  upr_dev <- ci_upr - bs_mean
  lwr_dev <- bs_mean - ci_lwr
  pred_y  <- exp( coef(m_mod)[1] + coef(m_mod)[2] * x_seq )
  
  # bootstrapped confidence intervals
  ci_df  <- data.frame(
    x      = x_seq,
    # ci_upr = apply(bs_mat,1,quantile,prob=0.975),
    # ci_lwr = apply(bs_mat,1,quantile,prob=0.025),
    ci_upr = pred_y + upr_dev,
    ci_lwr = pred_y - lwr_dev,
    pred_y = pred_y,
    beta   = coef(m_mod)[2]
  )
  
  # meta-confidence intervals
  meta_ci <- function( ii ){
    set.seed( ii )
    rnorm(1000, beta_bstr$beta[ii],
                beta_bstr$b_se[ii] )
  }
  
  bs_met <- lapply(1:nrow(beta_bstr), meta_ci) %>% unlist
  
  list( beta_bstr  = beta_bstr,
        ci_df      = ci_df,
        ci_meta    = bs_met )
  
}


# bootstrapped values
bs_pr_T      <- bootstrap_gamma( pars_all, 'Precipitation', 'Generation time') 
bs_tm_T      <- bootstrap_gamma( pars_all, 'Temperature',   'Generation time') 
bs_pr_T_nog  <- bootstrap_gamma( subset(pars_all, orgtype != 'Graminoid'), 'Precipitation', 'Generation time') 

# CI 
bs_pr_T$beta_bstr$beta %>% quantile(prob=c(0.025,0.975) )
bs_tm_T$beta_bstr$beta %>% quantile(prob=c(0.025,0.975) )
bs_pr_T_nog$beta_bstr$beta %>% quantile(prob=c(0.025,0.975) )

# means 
bs_pr_T$beta_bstr$beta %>% mean
bs_tm_T$beta_bstr$beta %>% mean
bs_pr_T_nog$beta_bstr$beta %>% mean

# meta CI
bs_pr_T$ci_meta %>% quantile(prob=c(0.025,0.975) )
bs_tm_T$ci_meta %>% quantile(prob=c(0.025,0.975) )
bs_pr_T_nog$ci_meta %>% quantile(prob=c(0.025,0.975) )

# Proportion of slopes below 0
sum( bs_pr_T$ci_meta < 0 ) / 1000000
sum( bs_tm_T$ci_meta < 0 ) / 1000000
sum( bs_pr_T_nog$ci_meta < 0 ) / 1000000

# produce independent plots
graph_2x1 <- function( y_name, x_name, bs_df, lty_i,
                       outlier = 'yes',
                       graminoid = 'yes' ){
  
  xlab <- ifelse('Iteroparity' == x_name,
                 'log(S)','log(T)')
  
  if( outlier == 'no'){
    pars_graph <- pars_all %>% 
      subset( !(clim_var == 'Precipitation' & b_t0 > 0.8) )  %>% 
      subset( !(clim_var == 'Temperature' & b_t0 < -2) )
  }else{
    pars_graph <- pars_all
  }
  
  p1 <- subset( pars_graph ) %>% 
    subset( clim_var == y_name) %>% 
    subset( lhs == x_name ) %>% 
    mutate( b_t0 = abs(b_t0) ) %>%  
    rename( `Plant type` = orgtype ) %>% 
    drop_na %>% 
    ggplot() +
    geom_point( aes( x     = lhs_val,
                     y     = b_t0,
                     size  = b_t0_w,
                     color = `Plant type` ),
                alpha = 0.5 ) +
    geom_ribbon( data=bs_df$ci_df,
                 aes(x,
                     ymin = ci_lwr,
                     ymax = ci_upr),
                 color = 'grey',
                 fill  = 'grey',
                 alpha = 0.5,
    ) +
    geom_line( data = bs_df$ci_df,
               aes(x = x, 
                   y = pred_y),
               lwd = 1.5,
               lty = lty_i ) +
    theme_minimal( ) +
    guides( size =F ) +
    labs( y = paste0('|Effect of ',y_name,'|'),
          x = xlab ) +
    theme( axis.text  = element_text(size=10),
           axis.title = element_text(size=15),
           strip.text = element_text(size=10) ) +
    scale_color_colorblind()
  
  if( graminoid == 'no'){
    
    # remove graminoids
    pars_graph <- subset( pars_graph, orgtype != 'Graminoid' )
    
    # TRICK to maintain color palette without graminoids
    col_nog <- c("#000000", "#56B4E9", "#009E73", "#F0E442")
    
    # make your own colour palette
    scale_colour_discrete <- function(...) {
      scale_colour_manual(..., values = col_nog)
    }
    
    p1 <- subset( pars_graph ) %>% 
      subset( clim_var == y_name) %>% 
      subset( lhs == x_name ) %>% 
      mutate( b_t0 = abs(b_t0) ) %>%  
      rename( `Plant type` = orgtype ) %>% 
      drop_na %>% 
      ggplot() +
      geom_point( aes( x     = lhs_val,
                       y     = b_t0,
                       size  = b_t0_w,
                       color = `Plant type` ),
                  alpha = 0.5 ) +
      geom_ribbon( data=bs_df$ci_df,
                   aes(x,
                       ymin = ci_lwr,
                       ymax = ci_upr),
                   color = 'grey',
                   fill  = 'grey',
                   alpha = 0.5,
      ) +
      geom_line( data = bs_df$ci_df,
                 aes(x = x, 
                     y = pred_y),
                 lwd = 1.5,
                 lty = lty_i ) +
      theme_minimal( ) +
      guides( size =F ) +
      labs( y = paste0('|Effect of ',y_name,'|'),
            x = xlab ) +
      theme( axis.text  = element_text(size=10),
             axis.title = element_text(size=15),
             strip.text = element_text(size=10) ) +
      scale_colour_discrete()
  
  }
  
  p1
  
}

# make plots
p1           <- graph_2x1('Precipitation', 'Generation time', bs_pr_T,     1)
p2           <- graph_2x1('Temperature',   'Generation time', bs_tm_T,     1)
p1_nog       <- graph_2x1('Precipitation', 'Generation time', bs_pr_T_nog, 1, 
                          graminoid = 'no')

# put plots together
shared_legend_cowplot <- function( x, y ){

  h_gt     <- plot_grid( x + theme(legend.position="none"),
                         y + theme(legend.position="none"),
                         labels = 'AUTO',
                         label_size = 18,
                         align='h', 
                         nrow = 2, ncol = 1 )
  
  # extract the legend from one of the plots
  legend <- cowplot::get_legend(
    # create some space to the left of the legend
    x + theme(legend.box.margin = margin(0, 0, 0, 4))
  )
  
  # output the new plots
  plot_grid(h_gt, legend, rel_widths = c(3, 1) )

}

# Put plots together
h_gt           <- shared_legend_cowplot( p1, p2 )

# clim effect versus LHS with A-D labels
ggsave('results/clim_vs_GenT.png', 
       h_gt,
       dpi=600, width=4.5, height=6.3)

# clim effect versus LHS with A-D labels
ggsave('results/prec_vs_GenT_nograminoids.png', 
       p1_nog,
       dpi=600, width=6.3, height=6.3)
