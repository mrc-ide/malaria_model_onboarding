################################################################################
##    title   specify_future_scenarios.R
##    purpose paramereterise two scenarios: one with intervention scale up and one without
##            started with MDA and RTS, S vaccination
##            helper functions for changing intervention parameters
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
options(repos = c(
  mrcide = "https://mrc-ide.r-universe.dev",
  CRAN = "https://cloud.r-project.org"))
drat::addRepo("malariaverse", "file:\\\\fi--didef3.dide.ic.ac.uk/malaria/malariaverse/drat")
remotes::install_github("mrc-ide/scene")
#remotes::install_github("mrc-ide/mvw")

library(scene)

# directories ------------------------------------------------------------------
setwd('Q:/')#


## view sites associated with a country site file  -----------------------------
iso<- 'ETH'

# get site data for a single country
site_data<- foresite:::get_site(iso)

# plot initial intervention coverage  ------------------------------------------
plot_interventions_combined(
  interventions = site_data$interventions,
  population = site_data$population,
  group_var = c("country", "name_1"),
  include = c("itn_use", "itn_input_dist", "tx_cov",
              "irs_cov", "rtss_cov", "smc_cov", "pmc_cov"),
  labels = c("ITN usage", "ITN model input", "Treatment",
             "IRS", "RTSS", "SMC", "PMC")
)


# use scene package to update site data
intervention <- copy(site_data)
group_var <- names(site_data$sites)
baseline<- copy(site_data)

# expand intervention years ----------------------------------------------------

# prep site data for model launch ----------------------------------------------
prep_inputs<- function(site_data){
  
  #' Prep inputs for batch launch
  #'
  #' @param site_data dataset with site files for country
  #' output: list with site name, urban/rural grouping, iso code, and parameters to pass into cluster
  
  
  # how many sites in this country?
  jobs<- nrow(site_data$sites)
  
  message(paste0('prepping ', jobs, ' jobs for model launch'))
  
  prep_site_data<- function(num){
    site<- site::single_site(site_file= site_data, index= num) 
    
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


# Expand the interventions for each site in the site file up to year 10
new_scenario$interventions <- new_scenario$interventions |>
  expand_interventions(max_year = 10, group_var = group_var)

# Add a target ITN usage of 60% in all sites by year 8
new_scenario$interventions <- new_scenario$interventions |>
  set_change_point(sites = new_scenario$sites, var = "itn_use", year = 8, target = 0.6)

# Add a target PMC coverage of 80% in site A
to_get_pmc <- new_scenario$sites[new_scenario$sites$site == "A", ]
new_scenario$interventions <- new_scenario$interventions |>
  set_change_point(sites = to_get_pmc, var = "pmc_cov", year = 10, target = 0.8)

# Add a target SMC coverage of 50% to any sites that have previously implemented SMC
to_get_smc <- ever_used(
  interventions = example_site$interventions,
  var = "smc_cov",
  group_var = group_var
)
new_scenario$interventions <- new_scenario$interventions |>
  set_change_point(sites = to_get_smc, var = "smc_cov", year = 10, target = 0.5)

# Linear scale up of coverage
new_scenario$interventions <- new_scenario$interventions |>
  linear_interpolate(vars = c("itn_use", "pmc_cov", "smc_cov"), group_var = group_var)

new_scenario$interventions <- new_scenario$interventions |>
  fill_extrapolate(group_var = group_var)

new_scenario$interventions <- new_scenario$interventions |>
  add_future_net_dist(group_var = group_var)

plot_interventions_combined(
  interventions = new_scenario$interventions,
  population = new_scenario$population,
  group_var = c("country", "site"),
  include = c("itn_use", "itn_input_dist", "fitted_usage", "tx_cov", "smc_cov", "pmc_cov"),
  labels = c("ITN usage", "ITN model input","ITN model usage", "Treatment","SMC", "PMC")
)



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
    #
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
  
  output<- lapply(c(1:sites), 
                  modify_interventions, 
                  mda= change_mda, 
                  rtss= change_rtss) 
    
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


