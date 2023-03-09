######################################################################################################################
##    title   specify_future_scenarios.R
##    purpose paramereterise two scenarios: one with intervention scale up and one without
##            started with MDA and RTS, S vaccination
##            helper functions for changing intervention parameters
##            code sourced from scene package vignette: https://mrc-ide.github.io/scene/articles/future-scenario.html
##    author  Lydia Haile (pulled from scene package vignette written by Pete Winskill)
######################################################################################################################

# load packages ----------------------------------------------------------------
options(repos = c(
  mrcide = "https://mrc-ide.r-universe.dev",
  CRAN = "https://cloud.r-project.org"))
drat::addRepo("malariaverse", "file:\\\\fi--didef3.dide.ic.ac.uk/malaria/malariaverse/drat")
#remotes::install_github("mrc-ide/scene")
#remotes::install_github("mrc-ide/mvw")

library(malariasimulation)
library(foresite)
library(data.table)
library(remotes)
library(didehpc)
library(conan)
library(drat)
library(vctrs)
library(scene)

# directories ------------------------------------------------------------------
setwd('Q:/')

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

# keep baseline data as is  ----------------------------------------------------
baseline<- copy(site_data)

# use scene package to update site data for intervention scenario
intvn <- copy(site_data)
group_var <- names(intvn$sites)

# settings you would like to modify for group of sites -------------------------
# I find it easier to keep track of changes in one place
# you can also hard-code them if you would like something more complex, this is a basic use case

expand_years<- 2050 # years you would like to expand intervention coverage out for. If you do not want to expand out to the future, set this to 0.

itn_change<- T    # do you want to modify ITN usage?
itn_target<- 0.6  # target for itn usage
itn_year<- 2045      # year you would like itn coverage to reach this target

pmc_change<- T    # do you want to modify PMC?
pmc_target<- 0.8  # target PMC
pmc_year<-   2040    # year for this target

rtss_change<- T   # do you want to modify RTSS?
rtss_target<- 0.5 # target for RTSS
rtss_year<-  2035    # year for target

smc_change<- T    # do you want to modify SMC?
smc_target<- 0.3  # target for SMC
smc_year<- 2025      # year of target

# expand intervention years ---------------------------------------------------
intvn$interventions <- intvn$interventions |>
  expand_interventions(max_year = expand_years,
                       group_var = group_var)

# ITN usage --------------------------------------------------------------------
if(itn_change== T){
  
intvn$interventions <- intvn$interventions |>
  set_change_point(sites = intvn$sites, 
                   var = "itn_use", 
                   year = itn_year, 
                   target = itn_target)
}
# PMC coverage  ----------------------------------------------------------------
if(pmc_change== T){
  
intvn$interventions <- intvn$interventions |>
  set_change_point(sites = intvn$sites, 
                   var = "pmc_cov", 
                   year = pmc_year, 
                   target = pmc_target)
}
# RTSS coverage ----------------------------------------------------------------
if (rtss_change== T){
  
intvn$interventions <- intvn$interventions |>
  set_change_point(sites = intvn$sites, 
                   var = "rtss_cov", 
                   year = rtss_year, 
                   target = rtss_target)
}

# SMC coverage  ----------------------------------------------------------------
if (smc_change== T){
  
intvn$interventions <- intvn$interventions |>
  set_change_point(sites = intvn$sites, 
                   var = "smc_cov", 
                   year = smc_year, 
                   target = smc_target)

}

# Linear scale up of coverage
intvn$interventions <- intvn$interventions |>
  linear_interpolate(vars = c("itn_use", "pmc_cov", "smc_cov", "rtss_cov"), 
                     group_var = group_var)

intvn$interventions <- intvn$interventions |>
  fill_extrapolate(group_var = group_var)

intvn$interventions <- intvn$interventions |>
  add_future_net_dist(group_var = group_var)


# plot the changes you made ----------------------------------------------------
plot_interventions_combined(
  interventions = intvn$interventions,
  population = intvn$population,
  group_var = c("country", "name_1"),
  include = c("itn_use", "itn_input_dist", "tx_cov", "smc_cov", "pmc_cov"),
  labels = c("ITN usage", "ITN model input", "Treatment","SMC", "PMC")
)


# plot baseline to make sure they look different  ------------------------------
plot_interventions_combined(
  interventions = baseline$interventions,
  population = baseline$population,
  group_var = c("country", "name_1"),
  include = c("itn_use", "itn_input_dist", "tx_cov", "smc_cov", "pmc_cov"),
  labels = c("ITN usage", "ITN model input", "Treatment","SMC", "PMC")
)

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


baseline<- prep_inputs(baseline)
intervention<- copy(intvn)


# submit jobs to cluster  ------------------------------------------------------
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


devtools::find_rtools()
