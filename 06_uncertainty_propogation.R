################################################################################
##  title   malariasimulation uncertainty propogation
##  author  Lydia Haile
##  purpose pulls development version of malariasimulation package to 
##          include model uncertainty in model runs
##          replicate 2 scenarios across 50 sets of model parameter posterior draws
##          info and example code sourced from PR: https://github.com/mrc-ide/malariasimulation/pull/219
################################################################################


# install development version of malariasimulation package
# before you do this you will need to detach the production version if you have it installed locally
# I also had to physically delete the package folder in my package directory
remove.packages("malariasimulation")

remotes::install_github('mrc-ide/malariasimulation@dev')

library(malariasimulation)
library(data.table)


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

# pick one site to obtain draws for so you can run this locally ----------------

baseline<- baseline[[1]]
intervention<- intervention[[1]]

# obtain 50 random draws from the model parameter posterior distribution
# made model run for a shorter time period so it would run more quickly

# Default pars
p <- baseline$param_list
sim <- run_simulation(timesteps= 1000, p)

plot(sim$timestep, sim$n_detect_730_3649 / sim$n_730_3649, t = "l", ylim = c(0, 1), ylab = "PfPr", xlab = "T")

cols<- rainbow(50)
dev.off()
for (i in 1:50){

message(paste0('obtaining parameter draw ', i))
param_draw<- baseline$param_list |>
             set_parameter_draw(sample(1:1000, 1)) |>
             set_equilibrium(init_EIR= 5) # appropriate init_EIR? Probably changes depending on the site

sim<- run_simulation(timesteps= 1000, param_draw)
lines(sim$timestep, sim$n_detect_730_3649 / sim$n_730_3649, col = cols[i])

}


# same process for intervention scenario  --------------------------------------
for (i in 1:50){
  param_draw<- intervention$param_list |>
    set_parameter_draw(sample(1:1000, 1)) |>
    set_equilibrium(init_EIR= 5) # appropriate init_EIR? Probably changes depending on the site
  
  sim<- run_simulation(timesteps= param_draw$timesteps, param_draw)
  lines(sim$timestep, sim$n_detect_730_3650 / sim$n_730_3650, col = cols[i])
  
}





#