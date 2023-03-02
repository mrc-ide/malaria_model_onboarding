################################################################################
##  title   01_cluster_setup.R
##  author  Lydia Haile
##  purpose run commands needed to run jobs on the cluster
##          can skip if you will only complete local runs
##  review this cluster documentation before running: 
##  https://mrc-ide.github.io/didehpc/articles/didehpc.html
################################################################################

# install packages  ------------------------------------------------------------
install.packages('drat')
install.packages('didehpc')
install.packages('pkgdepends')
install.packages('installr')
install.packages('conan')
install.packages('malariasimulation')
install_github('mrc-ide/malariasimulation')
install.packages('foresite')
install.packages('vctrs')
 
remotes::install_github("mrc-ide/individual")
options(repos = c(
mrcide = 'https://mrc-ide.r-universe.dev',
CRAN = 'https://cloud.r-project.org'))


drat::addRepo("malariaverse", "file:\\\\projects.dide.ic.ac.uk/malaria/malariaverse/drat")
install.packages("foresite", type = "source") # v0.1.0

# load packages ----------------------------------------------------------------
library(malariasimulation)
library(foresite)
library(data.table)
library(remotes)
library(didehpc)
library(conan)
library(drat)
library(vctrs)

# specify packages or functions you need for your job --------------------------
# you will need to have your code on a network drive to run it on the cluster
# I mapped my network drive to Q
setwd('Q:/')
options(didehpc.cluster = "fi--didemrchnb", didehpc.username = "lhaile")

# check the default configuration of your cluster setup ------------------------
didehpc::didehpc_config()

# web login to make sure you can get onto the cluster
# if this does not work, reach out to software dev team-- 
# you may need to be manually added to the permissions list
didehpc::web_login()

drat:::add("mrc-ide")
root<- "pkgs"

# load packages you will need to run malariasimulation package  ----------------
packages<- c('dplyr', 'tidyr', 'data.table', 'malariasimulation')
src <- conan::conan_sources("github::mrc-ide/malariasimulation")

# save a context (working environment for your code) ---------------------------
# additional script contains helper functions for larger scale model runs
ctx <- context::context_save('pkgs', 
                             packages = packages, 
                             package_sources = src,
                             sources = '/run_malaria_model.R')

# Manually install the dependencies, so that the pkgdepends does not
# get confused:  ---------------------------------------------------------------
obj <- didehpc::queue_didehpc(ctx, provision = "later")
obj$install_packages("mrc-ide/individual")
obj$install_packages("mrc-ide/malariaEquilibrium")
obj$install_packages("mrc-ide/malariasimulation")

# load context into queue
obj <- didehpc::queue_didehpc(ctx)

# test a run of the malariasimulation package on the cluster  ------------------

t <- obj$enqueue(run_simulation(timesteps= 100))

# check the status of the task
t$status()


# check the result of the task
t$result()


