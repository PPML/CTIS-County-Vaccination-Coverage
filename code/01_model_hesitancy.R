################################
## Model county-level hesitancy
################################

rm(list = ls())

## Load packages
library(data.table)
library(ggplot2)
library(tictoc)
library(lme4)
library(splitstackshape)
library(doParallel)
library(foreach)
library(lmeresampler)
library(survey)
library(fitdistrplus)
library(boot)

set.seed(123)

setwd(".")

group <- "child"

## County FIPS Lookup
fips <- as.data.table(tidycensus::fips_codes)
fips <- fips[, fips:=as.numeric(paste0(as.numeric(state_code), county_code))]
fips <- fips[, merge:=1]

## Prepare data
df <- fread("./data/child_race_ethnicity_wave11.csv")
df <- df[, wave:=11]
temp <- fread("./data/child_race_ethnicity_wave12.csv")
temp <- temp[, wave:=12]
df <- rbind(df, temp, fill = T)

df <- df[fips == 46113, fips:=46102] ## Oglala Lakota County County
df <- df[fips==2270, fips:=2158] ## Kusilvak Census Area
df <- df[fips==51515, fips:=51019] ## Bedford City Area
df <- df[, geo_code:=as.factor(fips)]

df <- df[raceethnicity=="Hispanic", race_grp:="Hispanic"]
df <- df[raceethnicity=="NonHispanicAmericanIndianAlaskaNative", race_grp:="American Indian or Alaska Native"]
df <- df[raceethnicity=="NonHispanicAsian", race_grp:="Asian"]
df <- df[raceethnicity=="NonHispanicBlackAfricanAmerican", race_grp:="Black"]
df <- df[raceethnicity=="NonHispanicMultipleOther", race_grp:="Other"]
df <- df[raceethnicity=="NonHispanicNativeHawaiianPacificIslander", race_grp:="Native Hawaiian or Other Pacific Islander"]
df <- df[raceethnicity=="NonHispanicWhite", race_grp:="White"]

df <- df[oldest_child_age=="5-11", child_age_indic:=1]
df <- df[is.na(child_age_indic), child_age_indic:=0]

## Create a smaller modeling dataset with just the variables of interest
mod_df <- df[, .(sex, age_cat_fb, education, race_grp, child_age_indic, geo_code, child_status, state_name, state_code, fips)]

## Generate the outcome indicator (Hesitancy = Def No or Prob No)
mod_df <- mod_df[child_status%in%c("Definitely No", "Probably No"), outcome:=1]
mod_df <- mod_df[child_status%in%c("Vaccinated", "Definitely Yes", "Probably Yes"), outcome:=0]

## Create grouped counties for counties with sample size smaller than 10
mod_df <- mod_df[, county_n:=.N, by = c("fips")]

## Sample characteristics for reporting
temp <- unique(mod_df[, .(state_name, fips, county_n)])
fips_temp <- as.data.table(tidycensus::fips_codes)
fips_temp <- fips_temp[, fips:=as.numeric(paste0(state_code, county_code))]
temp <- merge(temp, fips_temp[state_name%in%unique(temp$state_name)], by = c("fips", "state_name"), all=T)
temp <- temp[fips!=51515 & fips!=46113 & fips!= 2270]
nrow(temp[county_n>=100]) #1203
nrow(temp[county_n<10]) #261
nrow(temp[is.na(county_n)]) #23
rm(temp, fips_temp)

## Create grouped counties
mod_df <- mod_df[, group_county:=fips]
mod_df <- mod_df[county_n<=10, group_county:=state_code] ## If counties have sample size <= 10, then group them together within a state

## Fit mixed effect logistic regression
tic()
mod <- glmer(outcome ~ sex + age_cat_fb + education + race_grp + child_age_indic + (1 | state_name) + (1 | state_name:group_county),
             data = mod_df, family = "binomial", control = glmerControl(optimizer = "bobyqa",
                                                                        optCtrl=list(maxfun=2e9)))
toc()

## Save model object
save(mod, file = "./results/fitted_hesitancy_glmm_mod_childage_final.RData")
