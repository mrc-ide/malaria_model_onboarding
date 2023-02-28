################################################################################
##    title     03_model_post_processing.R
##    author    Lydia Haile
##    purpose   post-process model outputs
##              aggregate cases, deaths, and DALYs
################################################################################

# load packages ----------------------------------------------------------------
library(malariasimulation)
library(data.table)
library(ggplot2)

# load in files  ------------------------------------------
dir<- '/model_test_run/' #directory where outputs are


files<- list.files(dir, full.names = T)
output<- rbindlist(lapply(files, readRDS), fill= T)

# test on one output file
model<- data.table(readRDS('C:/Documents and Settings/lhaile/Documents/raw_model_output_Addis Abeba_rural.RDS'))


# reformat and aggregate model outputs  ----------------------------------------
# reformat case outputs to long
raw_clin <- model |> 
  select(timestep, contains("n_inc_clin"), contains("n_inc_sev"), contains("n_age"), 'n_treated') |> # need clinical cases, severe cases, population, number treated 
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
# use D_count + U_count to calculate prevalence?
interval<- 30 # monthly

raw_clin <- raw_clin |>
  mutate(time = as.integer(timestep/ interval))


# sum cases based on this interval, also by country
raw_clin<- raw_clin |> 
  group_by(age_years_start, time, iso) |> 
  mutate(clinical = sum(clinical),
         severe = sum(severe),
         n_treated = sum(n_treated),
         population = round(mean(population))) |>
  select(-timestep) |>
  distinct()



# calculate rates based on this interval
raw_clin<- raw_clin |> 
  mutate(clin_rate = clinical/ population,
         severe_rate = severe/ population)


# calculate deaths -------------------------------------------------------------
# per GTS deaths are calculated using a scaling factor
# calculate deaths
treatment_scaler<- 0.45
cfr<- 0.215
mild_dw<- 0.006
moderate_dw<- 0.051
severe_dw<- 0.133
cfr<- 0.215
episode_length<-  0.01375
severe_episode_length<-  0.04795
lifespan<- 63

# in count space
# deaths= 0.215 * severe cases
# where severe cases= 0, deaths= 0
# remove a proportion of the cases that have received treatment
raw_clin<- raw_clin |>
  mutate(deaths =  cfr * severe - (n_treated * treatment_scaler))


raw_clin<- data.table(raw_clin)
raw_clin[deaths< 0 , deaths:= 0]
raw_clin[, deaths:= as.integer(deaths)]

# calculate ylls
raw_clin[, `:=` (age_years_start= as.numeric(age_years_start),
                 age_years_end= as.numeric(age_years_end))]

raw_clin <- raw_clin |>
  mutate(yll= deaths * (lifespan - (age_years_end - age_years_start)/2))

# calculate YLDs
raw_clin<- raw_clin |> 
  mutate(yld= severe * severe_dw * severe_episode_length + 
           clinical * mean(mild_dw+ moderate_dw) * episode_length)

# calculate DALYs
raw_clin<- raw_clin |>
  mutate(daly= yll+ yld)


# summarize
# raw_clin |> 
#  group_by(age_years_start) |> 
#  summarise(yll= sum(yll),
#            yld= sum(yld),
#            daly= sum(daly))
#

# plot outputs over time  ------------------------------------------------------

ggplot(data= raw_clin, mapping = aes(x= time, y= clin_rate))+
  geom_line(col = "grey80") +
  stat_smooth(col = "darkblue", se = FALSE) +
  facet_wrap(~age_years_start)

ggplot(data= raw_clin, mapping = aes(x= time, y= severe_rate))+
  geom_line(col = "grey80") +
  stat_smooth(col = "darkblue", se = FALSE) +
  facet_wrap(~age_years_start)

ggplot(data= raw_clin, mapping = aes(x= time, y= yll))+
  geom_line(col = "grey80") +
  stat_smooth(col = "darkblue", se = FALSE) +
  facet_wrap(~age_years_start)

ggplot(data= raw_clin, mapping = aes(x= time, y= deaths))+
  geom_line(col = "grey80") +
  stat_smooth(col = "darkblue", se = FALSE) +
  facet_wrap(~age_years_start)

ggplot(data= raw_clin, mapping = aes(x= time, y= daly))+
  geom_line(col = "grey80") +
  stat_smooth(col = "darkblue", se = FALSE) +
  facet_wrap(~age_years_start)
