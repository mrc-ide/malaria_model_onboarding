################################################################################
##    title   specify_future_scenarios.R
##    purpose paramereterise two scenarios: one with intervention scale up and one without
##            started with MDA and RTS, S vaccinnation
##    author  Lydia Haile
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
setwd('Q:/')#

## view sites associated with a country site file  -----------------------------
iso<- 'ETH'

# get site data
site_data<- foresite:::get_site(iso)

prep_inputs<- function(site_data){
  
  #' Prep inputs (without parameter updates)
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

baseline<- output
intervention<- copy(output)
# change intervention coverage parameters --------------------------------------

specify_intervention_coverage<- function(output, mda= T, rtss= T){
  #' Parent function for modify_interventions function
  #' @param output parameter file created from prep_inputs
  #' @param change_mda if true, change MDA parameters
  #' @param change_rtss if true, change RTS,S parameters
  
  sites<- length(output) #number of sites you would like to modify coverage for

  modify_interventions<- function(index, mda= T, rtss= T){
   
    #' Specify intervention coverage for pre-existing site data prepped using prep_inputs
    #'
    #' @param index row of site file you will modify
    #' @param mda if true, change MDA parameters
    #' @param rtss if true, change RTS,S parameters
    
    
    params<- output[[index]]$param_list
    
    if (mda== T){
    # Add MDA strategy
    mda_events = data.frame(
      timestep = c(1, 2) * 365,
      name=c("MDA 1", "MDA 2")
    )  
    params <- set_mda(
      params,
      drug = 1,
      timesteps = mda_events$timestep,
      coverages = rep(.8, 2),
      min_ages = rep(0, length(mda_events$timestep)),
      max_ages = rep(200 * 365, length(mda_events$timestep))
    )
    
    }
    
    # different RTS,S scenario
    
    if (rtss== T){
    month<- 30
    params <- set_mass_rtss(
      params,
      timesteps = params$timesteps,
      coverages = rep(1, 2),
      min_wait = 0, # minimum acceptable time since the last vaccination
      min_ages = 5 * month, #target population minimum
      max_ages = 17 * month, #target population maximum
      boosters = 18 * month,
      booster_coverage = 0.95
    )
    }
    
    # reassign to input dataset
    output[[index]]$param_list<- params
    message(paste0('reassigned row ', index))
    
    return(output)
  }
  
  output<- lapply(c(1:sites), modify_interventions, mda= change_mda, rtss= change_rtss) 
    
  }
  
  
intervention<- specify_intervention_coverage(intervention)  


message(paste0('submitting ', length(output),  ' jobs'))

# load packages you will need to run malariasimulation package  ----------------
packages<- c('dplyr', 'tidyr', 'data.table', 'malariasimulation')
src <- conan::conan_sources("github::mrc-ide/malariasimulation")

# save a context (working environment for your code) ---------------------------
# additional script contains helper functions for larger scale model runs
ctx <- context::context_save('pkgs', 
                             packages = packages, 
                             package_sources = src,
                             sources = 'Q:/model_onboarding/run_malaria_model.R')


# load context into queue
obj <- didehpc::queue_didehpc(ctx)


# run baseline jobs
fold<- 'Q:/model_test_run/baseline/' # folder you would like to save outputs in
dir.create(fold)
grp1 <- obj$lapply(output, run_malaria_model, folder= fold)

# run intervention jobs
fold<- 'Q:/model_test_run/intervention/' # folder you would like to save outputs in
dir.create(fold)

grp2 <- obj$lapply(intervention, run_malaria_model, folder= fold)


