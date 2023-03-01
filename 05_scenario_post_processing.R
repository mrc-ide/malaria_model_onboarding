################################################################################################
##  title   05_scenario_post_processing.R
##  author  Lydia Haile
##  purpose calculate cases, deaths, and DALYs averted (baseline versus intervention scenarios)
##          generally similar post processing to 03_model_post_processing.R
##          post-processes outputs from 04_specify_future_scenarios.R
################################################################################################


# load packages ----------------------------------------------------------------
library(malariasimulation)
library(data.table)
library(ggplot2)

# load in files  ------------------------------------------
dir<- '/model_test_run/' #directory where outputs are

intvn_files<- list.files(paste0(dir, 'intervention/'), full.names = T)
intvn<- rbindlist(lapply(files, readRDS), fill= T)
intvn[, scenario:= 'intervention']

baseline_files<- list.files(paste0(dir, 'baseline/'), full.names = T)
baseline<- rbindlist(lapply(files, readRDS))

# test on one output file
model<- data.table(readRDS('C:/Documents and Settings/lhaile/Documents/raw_model_output_Addis Abeba_rural.RDS'))


# functions for processing model outputs  --------------------------------------

aggregate_outputs<- function(dt, interval){
  
  #' aggregate incident cases based on a pre-determined time interval (expressed in days). 
  #' sum to the country level.
  #' 
  #' @param dt       raw model output from malariasimulation package
  #' @param interval time period you would like to calculate incidence over, expressed in days.
  #' 
  #' output: data table with summed cases and rates over specified time interval.
  
  # reformat case outputs to long
  # need clinical cases, severe cases, population, number treated
  
  message(paste0('aggregating outputs by time interval: ', interval, ' days'))
  dt <- dt |> 
    select(timestep, 
           contains("n_inc_clin"), contains("n_inc_sev"), contains("n_age"), 'n_treated') |>  
    pivot_longer(c(contains("n_inc_clin"), contains("n_inc_sev"), contains("n_age")),
                 names_to = "age", 
                 values_to = "value") |>
    mutate(type = ifelse(grepl('n_inc_clinical', age), 'clinical', 'severe')) |>
    mutate(type = ifelse(grepl('n_age', age), 'population', type))|>
    mutate(age = gsub('n_inc_clinical_', '', age),
           age = gsub('n_inc_severe_', '', age),
           age = gsub('n_age_', '', age)) |>
    spread(key= type, value= value) |>
    separate(age, into = c("age_years_start", "age_years_end"), sep = "_")
  
  
  # calculate incidence based on some time interval
  dt <- dt |>
    mutate(time = as.integer(timestep/ interval))
  
  
  # sum cases based on this interval, also by country
  dt<- dt |> 
    group_by(age_years_start, time, iso) |> 
    mutate(clinical = sum(clinical),
           severe = sum(severe),
           n_treated = sum(n_treated),
           population = round(mean(population))) |>
    select(-timestep) |>
    distinct()
  
  
  
  # calculate rates based on this interval
  dt<- dt |> 
    mutate(clin_rate = clinical/ population,
           severe_rate = severe/ population)
  
  return(dt)
  message('completed aggregation')
  
}
calculate_deaths_ylls<- function(dt, cfr= 0.215, treatment_scaler= 0.45, lifespan= 0.63){
  
  #' Calculate deaths + years of life lost (YLLs) per GTS method.
  #' Where severe cases= 0, deaths= 0. Additionally remove a proportion of the cases that have received treatment.
  #' 
  #' @param dt               malariasimulation model outputs with columns 'severe' for severe incidence and 'n_treated' for number of individuals treated
  #' @param cfr              per GTS method, deaths are calculated using a case fatality ratio (CFR) value that is applied to severe incidence.
  #'                         see World Malaria Report and Wilson et al. for more information.
  #' @param treatment_scaler we remove a proportion of cases that have received treatment, assuming that 45% of treated cases remit and
  #'                         are no longer susceptible to mortality.
  #' @param lifespan         expected lifespan used to calculate YLLs. 
  #'                         YLLs are calculated by multiplying deaths by the number of years an individual was expected to live past their year of death.
  #'                         when comparing YLLs across different locations, it is recommended to use the same lifespan across YLL calculations.
  #'                         Typically, you should use the highest observed life expectancy in the region/ location you are studying.
  #' Output: data table with columns titled 'deaths' and 'yll'
  
  
  message('calculating deaths and YLLs')
  
  dt<- dt |>
    mutate(deaths =  cfr * severe - (n_treated * treatment_scaler))
  
  
  dt<- data.table(dt)
  
  dt[deaths< 0 , deaths:= 0]
  dt[, deaths:= as.integer(deaths)]
  
  # calculate ylls
  dt[, `:=` (age_years_start= as.numeric(age_years_start),
             age_years_end= as.numeric(age_years_end))]
  
  dt <- dt |>
    mutate(yll= deaths * (lifespan - (age_years_end - age_years_start)/2))  
  
  message('completed calculation of deaths and YLLs')
  
  return(dt)
}
calculate_ylds_dalys<- function(dt, 
                                mild_dw= 0.006, 
                                moderate_dw= 0.051, 
                                severe_dw= 0.133,
                                clin_episode_length= 0.01375,
                                severe_episode_length= 0.04795){
  
  #' Calculate Years Lived with Disability (YLDs) and Disability-Adjusted Life-Years
  #' based on disability weights from the Global Burden of Disease study.
  #' 
  #' Disability weights sourced here:
  #' https://view.officeapps.live.com/op/view.aspx?src=https%3A%2F%2Fghdx.healthdata.org%2Fsites%2Fdefault%2Ffiles%2Frecord-attached-files%2FIHME_GBD_2017_DISABILITY_WEIGHTS_Y2018M11D08.XLSX&wdOrigin=BROWSELINK
  #' 
  #' Keep in mind this is an approximation of YLD estimation from the GBD study; disability due to comorbid conditions
  #' such as motor impairment and anemia are excluded.
  #' 
  #' For now, we assume that cases under 5 are moderate, due to the higher severity of malaria at younger ages.
  #' This assumption may be revisited in the future, 
  #' @param dt                    input dataset with columns 'clinical' for clinical incidence,
  #'                              'severe' for severe incidence, 'deaths', 'ylls', 'age_years_start', and 'age_years_end'
  #' @param mild_dw               disability weight for mild malaria
  #' @param moderate_dw           disability weight for moderate malaria
  #' @param severe_dw             disability weight for severe malaria
  #' @param clin_episode_length   length of an episode of clinical malaria
  #' @param severe_episode_length length of an episode of severe malaria 
  #' 
  #' Output: data table with columns 'ylds' and 'dalys'
  
  require(data.table)
  
  message('calculating YLDs and DALYs')
  
  # calculate YLDs  ---
  dt[ age_years_end < 5, yld:= severe * severe_dw * severe_episode_length + 
        clinical * moderate_dw * clin_episode_length]
  
  dt[ age_years_end >= 5, yld:= severe * severe_dw * severe_episode_length + 
        clinical * mild_dw * clin_episode_length]
  
  # calculate DALYs ---
  dt<- dt |>
    mutate(daly= yll+ yld)
  
  message('completed calculation of YLDs and DALYs')
  
  return(dt)
}


# process outputs --------------------------------------------------------------
bl<-aggregate_outputs(baseline, interval= 30)
bl<- calculate_deaths_ylls(bl)
bl<- calculate_ylds_dalys(bl)


intvn<- aggregate_outputs(intvn, interval= 30)
intvn<- calculate_deaths_ylls(intvn)
intvn<- calculate_ylds_dalys(intvn)

# calculate cases, deaths, and DALYs averted  ----------------------------------
setnames(intvn, c('clinical', 'severe', 'deaths', 'yll', 'yld','daly'),
                c('clinical_intvn', 'severe_intvn', 'deaths_intvn', 'yll_intvn', 'yld_intvn', 'daly_intvn'))


# merge dataframes by timesteps
dt<- data.table(merge(bl, intvn, by= c('iso', 'time', 'age_years_start', 'age_years_end')))

dt<- dt |>
mutate(cases_averted= clinical- clinical_intvn,
       severe_averted= severe- severe_intvn,
       deaths_averted= deaths- deaths_intvn,
       yll_averted= yll- yll_intvn,
       yld_averted= yld- yld_intvn,
       dalys_averted= daly- daly_intvn)

# save to output folder
output_folder<- 'Q:/model_test_run/final_output/'
dir.create(output_folder)

write.csv(dt, paste0(output_folder, 'processed_outputs.csv'))





