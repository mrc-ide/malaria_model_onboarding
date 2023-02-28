################################################################################
##    title   single_country_model_launch
##    author  Lydia Haile
##    purpose launch all site files for a single country (no parameter changes)
################################################################################


# load packages ----------------------------------------------------------------
library(malariasimulation)
library(foresite)
library(data.table)
library(remotes)
library(didehpc)
library(conan)
library(drat)
library(vctrs)

# directories ------------------------------------------------------------------
setwd('Q:/')


# load packages you will need to run malariasimulation package  ----------------
packages<- c('dplyr', 'tidyr', 'data.table', 'malariasimulation')
src <- conan::conan_sources("github::mrc-ide/malariasimulation")

# save a context (working environment for your code) ---------------------------
# additional script contains helper functions for larger scale model runs
ctx <- context::context_save('pkgs', 
                             packages = packages, 
                             package_sources = src,
                             sources = '/run_malaria_model.R')


# load context into queue
obj <- didehpc::queue_didehpc(ctx)


# run bulk model runs ----------------------------------------------------------
# view sites associated with a country site file  ------------------------------
iso<- 'ETH'

# get site data
site_data<- foresite:::get_site(iso)



prep_inputs<- function(site_data){
  
  #' Prep inputs (without parameter updates)-- not explicitly passing site_data object in but this should be updated
  #'
  #' @param site_data dataset with site files for country
  #' output: list with site name, urban/rural grouping, iso code, and parameters to pass into cluster
  
  
  # how many sites in this country?
  jobs<- nrow(site_data$sites)
  
  message(paste0('prepping ', jobs, ' jobs for model launch'))
  
  prep_site_data<- function(num){
    site<- site::single_site(site_data, index= num) 
    
    ## get site info
    site_name<- site$sites$name_1
    ur<- site$sites$urban_rural
    iso<- site$sites$iso3c
    message(paste0('prepping inputs for site ', site_name, ' ', ur))
    
    # pull parameters for this site
    params<- site::site_parameters(
      interventions = site$interventions,
      demography = site$demography,
      vectors = site$vectors,
      seasonality = site$seasonality,
      eir= site$eir$eir[1],
      overrides = list(human_population= 1000)
    )
    
    inputs<- list('param_list'= params, 'site_name'= site_name, 'ur'= ur, 'iso'= iso)
    return(inputs)
  }
  output<- lapply(c(1:jobs), prep_site_data)
  

  
}

output<- prep_inputs(site_data)

# submit batch of jobs for single country run ----------------------------------
grp <- obj$lapply(output, run_malaria_model)

folder<- 'Q:/model_test_run/' # folder you would like to save outputs in

run_malaria_model<- function(output, folder){
  
  #' Prep inputs (without parameter updates)-- not explicitly passing site_data object in but this should be updated
  #' 
  #' @param output  list with site name, urban/rural grouping, and parameters to pass into cluster
  #'                generated from prep_inputs function
  #' output: model outputs for site with provided parameters (saved into pre-specified folder)
  
  # run the model
  message('running the model')
  model<- malariasimulation::run_simulation(timesteps = output$param_list$timesteps,
                                            parameters = output$param_list) 
  
  model<- data.table(model)
  model[, site_name:= output$site_name]
  model[, urban_rural:=output$ur]
  model[, iso:= output$iso]
  
  # save model runs somewhere
  message('saving the model')
  saveRDS(model, file= paste0(folder, 'raw_model_output_', output$site_name, '_', output$ur, '.RDS'))
}

