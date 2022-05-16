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
mod_df <- df[, .(sex, age_cat_fb, education, race_grp, geo_code, child_age_indic, child_status, state_name, state_code, fips)]

## Generate the outcome indicator (Hesitancy = Def No or Prob No)
mod_df <- mod_df[child_status%in%c("Definitely No", "Probably No"), outcome:=1]
mod_df <- mod_df[child_status%in%c("Vaccinated", "Definitely Yes", "Probably Yes"), outcome:=0]

## Create grouped counties for counties with sample size smaller than 10
mod_df <- mod_df[, county_n:=.N, by = c("fips")]
mod_df <- mod_df[, group_county:=fips]
mod_df <- mod_df[county_n<=10, group_county:=state_code]

## Generate County-Level Hesitancy Draws
load("./results/fitted_hesitancy_glmm_mod_childage_final.RData")

parallel::detectCores()
n.cores <- 50

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

## START BOOTSTRAP LOOP ##
out <- foreach (big_loop = 1:1000, .combine = 'rbind',
                .packages = c("data.table", "tidycensus", "boot", "lme4")) %dopar% {
                  
                  ## Household analysis for children
                  pred_age_min <- 12
                  pred_age_max <- 17
                  
                  fips <- as.data.table(tidycensus::fips_codes)
                  fips <- fips[, fips:=as.numeric(paste0(as.numeric(state_code), county_code))]
                  fips <- fips[, merge:=1]
                  
                  pred_df <- as.data.table(expand.grid(merge = 1, sex = unique(df$sex),
                                                       age_cat_fb = unique(df$age_cat_fb), education = unique(df$education),
                                                       race_grp = unique(df$race_grp)))
                  
                  if (pred_age_max==11) {
                    pred_df <- pred_df[, child_age_indic:=1]
                  } else {
                    pred_df <- pred_df[, child_age_indic:=0]
                  }
                  pred_df <- merge(pred_df, fips, by = "merge", allow.cartesian = T)
                  pred_df <- pred_df[, merge:=NULL]
                  pred_df <- pred_df[, geo_code:=as.factor(fips)]
                  pred_df <- pred_df[, combine:=ifelse(!(fips%in%unique(mod_df$fips[mod_df$county_n>10])), 1, 0)]
                  pred_df <- pred_df[, county_merge:=paste0(state_name, ":", fips)]
                  pred_df <- pred_df[combine == 1, county_merge:=paste0(state_name, ":", as.numeric(state_code))]
                  
                  ## Take fitted model
                  n.sims = 1
                  set.seed(big_loop)
                  
                  ## Create sigma vector
                  sigmahat <- rep(1, n.sims)
                  
                  ## Prep RE output
                  re.xb <- vector(getME(mod, "n_rfacs"), mode = "list")
                  names(re.xb) <- names(ngrps(mod))
                  
                  ## Loop through REs and sample
                  for (j in rev(names(re.xb))) {
                    reMeans <- as.matrix(ranef(mod)[[j]])
                    reMatrix <- attr(ranef(mod, condVar = TRUE)[[j]], 
                                     which = "postVar")
                    
                    ## Keep observed levels that are needed for prediction
                    keep <- rownames(reMeans)
                    dimnames(reMatrix)[[3]] <- keep
                    reMatrix <- reMatrix[, , keep, drop = FALSE]
                    
                    ## Sample from mv distribution
                    tmpList <- vector(length = nrow(reMeans), mode = "list")
                    for (k in 1:nrow(reMeans)) {
                      meanTmp <- reMeans[k, ]
                      names(meanTmp) <- NULL
                      matrixTmp <- as.matrix(reMatrix[, , k])
                      tmpList[[k]] <- as.matrix(mvtnorm::rmvnorm(n = n.sims, 
                                                                 mean = meanTmp, sigma = matrixTmp, method = "chol"))
                    }
                    
                    out <- NULL
                    for (i in 1:nrow(reMeans)) {
                      temp <- data.table(tmp_geo = rownames(reMeans)[i], re_draws = as.numeric(tmpList[[i]]), sample = 1:n.sims)
                      setnames(temp, "tmp_geo", j)
                      out <- rbind(out, temp, fill = T)
                    }
                    if (j == "state_name:group_county") {
                      print("Merging Counties")
                      setnames(out, j, "county_merge")
                      setnames(out, "re_draws", "county_re")
                      pred_df <- merge(pred_df, out, by = c("county_merge", "sample"), all.x=T)
                    } else {
                      print("Merging States")
                      setnames(out, "re_draws", "state_re")
                      pred_df <- merge(pred_df, out, by = c("state_name"), all.x=T, allow.cartesian = T)
                    }
                  }
                  
                  var_re <- as.data.table(VarCorr(mod))
                  var_re <- as.numeric(var_re$vcov[var_re$grp=="state_name:group_county"])
                  
                  missing_re <- unique(pred_df$fips[is.na(pred_df$county_re) &
                                                      !pred_df$state_name%in%c("American Samoa", "Guam",
                                                                               "Northern Mariana Islands",
                                                                               "Puerto Rico", "U.S. Minor Outlying Islands",
                                                                               "U.S. Virgin Islands")])
                  
                  samp_re <- data.table(fips  = rep(missing_re, each = n.sims), temp_re = as.numeric(mvtnorm::rmvnorm(n = n.sims*length(missing_re), 
                                                                                                                      mean = as.vector(0), sigma = as.matrix(var_re), method = "chol")),
                                        sample = rep(c(1:n.sims), length(missing_re)))
                  pred_df <- merge(pred_df, samp_re, by = c("fips", "sample"), all.x=T)
                  pred_df <- pred_df[is.na(county_re), county_re:=temp_re]
                  pred_df <- pred_df[, temp_re:=NULL]  
                  
                  fe.tmp <- fixef(mod)
                  vcov.tmp <- as.matrix(vcov(mod))
                  
                  betaSim <- abind::abind(lapply(1:n.sims, function(x) mvtnorm::rmvnorm(n = 1, 
                                                                                        mean = fe.tmp, sigma = vcov.tmp, method = "chol")), 
                                          along = 1)
                  
                  betaSim <- as.data.table(betaSim)
                  betaSim <- betaSim[, sample:=.I]
                  setnames(betaSim, "(Intercept)", "Intercept")
                  ref_grps <- data.table(sexFemale = rep(0, n.sims),
                                         `age_cat_fb18-24` = rep(0, n.sims),
                                         `educationFour yr degree` = rep(0, n.sims),
                                         `race_grpAmerican Indian or Alaska Native` = rep(0, n.sims))
                  betaSim <- cbind(betaSim, ref_grps)
                  betaSim <- melt(betaSim, id.vars = c("sample", "Intercept", colnames(betaSim)[colnames(betaSim)%like%"age_cat" |
                                                                                                  colnames(betaSim)%like%"education" |
                                                                                                  colnames(betaSim)%like%"race_grp"],
                                                       "child_age_indic"))
                  setnames(betaSim, "variable", "sex")
                  setnames(betaSim, "value", "sex_fe")
                  betaSim <- melt(betaSim, id.vars = c("sample", "Intercept", "sex", "sex_fe", colnames(betaSim)[colnames(betaSim)%like%"age_cat" |
                                                                                                                   colnames(betaSim)%like%"education"],
                                                       "child_age_indic"))
                  setnames(betaSim, "variable", "race_grp")
                  setnames(betaSim, "value", "race_fe")
                  betaSim <- melt(betaSim, id.vars = c("sample", "Intercept","sex", "sex_fe", "race_grp", "race_fe", colnames(betaSim)[colnames(betaSim)%like%"age_cat"],
                                                       "child_age_indic"))
                  setnames(betaSim, "variable", "education")
                  setnames(betaSim, "value", "education_fe")
                  betaSim <- melt(betaSim, id.vars = c("sample", "Intercept", "sex", "sex_fe", "race_grp", "race_fe",
                                                       "education", "education_fe", "child_age_indic"))
                  setnames(betaSim, "variable", "age_cat_fb")
                  setnames(betaSim, "value", "age_fe")
                  setnames(betaSim, "child_age_indic", "child_age_indic_fe")
                  betaSim  <- betaSim[, child_age_indic:=1]
                  temp <- copy(betaSim)
                  temp <- temp[, child_age_indic:=0]
                  temp <- temp[, child_age_indic_fe:=0]
                  betaSim <- rbind(betaSim, temp)
                  betaSim <- betaSim[, sex:=gsub("sex", "", sex)]
                  betaSim <- betaSim[, race_grp:=gsub("race_grp", "", race_grp)]
                  betaSim <- betaSim[, education:=gsub("education", "", education)]
                  betaSim <- betaSim[, age_cat_fb:=gsub("age_cat_fb", "", age_cat_fb)]
                  
                  expand_df <- unique(pred_df[,.(fips)])
                  expand_df <- expand_df[, merge:=1]
                  
                  betaSim <- betaSim[, merge:=1]
                  betaSim <- merge(betaSim, expand_df, by = "merge", allow.cartesian = T)
                  betaSim <- betaSim[, merge:=NULL]
                  
                  pred_df <- merge(pred_df, betaSim,
                                   by = c("race_grp", "education", "age_cat_fb", "sex", "sample", "fips", "child_age_indic"))
                  
                  pred_df <- pred_df[, pred_prob:=boot::inv.logit(Intercept + state_re + county_re + sex_fe + race_fe + education_fe + age_fe + child_age_indic_fe)]
                  pred_df <- pred_df[, c("state_re", "county_re", "Intercept", "sex_fe", "race_fe", "education_fe",
                                         "age_fe", "county_merge", "child_age_indic_fe", "child_age_indic"):=NULL]
                  
                  pred_df <- pred_df[!(fips%in%c(46113, 2270, 51515))] ## Remove three counties that have been changed
                  
                  # ## Clean-up
                  # rm(df, expand_df, betaSim, matrixTmp, out, re.xb, ref_grps, reMeans, temp, tmpList, vcov.tmp,
                  #    fe.tmp, i, j, k, keep, meanTmp, reMatrix, sigmahat, mod_df)
                  
                  
                  ## Post-stratification
                  acs <- fread("./data/usa_00033.csv")
                  acs <- acs[GQ%in%c(1, 2, 5)] ## Remove institutionalized populations
                  
                  acs <- acs[HISPAN%in%c(1, 2, 3, 4), race_grp:="Hispanic"]
                  acs <- acs[is.na(race_grp) & RACE==1, race_grp:="White"]
                  acs <- acs[is.na(race_grp) & RACE==2, race_grp:="Black"]
                  acs <- acs[is.na(race_grp) & RACE==3, race_grp:="American Indian or Alaska Native"]
                  acs <- acs[is.na(race_grp) & RACE%in%c(4, 5), race_grp:="Asian"]
                  acs <- acs[is.na(race_grp) & RACE%in%c(7, 8, 9), race_grp:="Other"]
                  acs <- acs[is.na(race_grp) & RACED%in%c(630, 680, 681, 682, 685, 689, 690, 699), race_grp:="Native Hawaiian or Other Pacific Islander"]
                  acs <- acs[is.na(race_grp), race_grp:="Asian"]
                  
                  ## Education
                  acs <- acs[EDUC%in%c(0, 1, 2, 3, 4, 5, 6), education:="HS or less"]
                  acs <- acs[EDUCD%in%c(65, 71, 81), education:="Some college or 2 yr"]
                  acs <- acs[EDUCD%in%c(101), education:="Four yr degree"]
                  acs <- acs[EDUCD%in%c(114, 115, 116), education:="Graduate degree"]
                  
                  ## Age cat FB
                  acs <- acs[AGE<25 & AGE>=18, age_cat_fb:="18-24"]
                  acs <- acs[AGE<35 & AGE>=25, age_cat_fb:="25-34"]
                  acs <- acs[AGE<45 & AGE>=35, age_cat_fb:="35-44"]
                  acs <- acs[AGE<55 & AGE>=45, age_cat_fb:="45-54"]
                  acs <- acs[AGE<65 & AGE>=55, age_cat_fb:="55-64"]
                  acs <- acs[AGE>=65, age_cat_fb:="65+"]
                  
                  ## Sex
                  acs <- acs[SEX==1, sex:="Male"]
                  acs <- acs[SEX==2, sex:="Female"]
                  
                  ## Focus on parents of children < 18
                  acs <- acs[, hh_child:=ifelse(AGE<=pred_age_max & AGE>=pred_age_min, 1, 0)]
                  acs <- acs[, hh_child:=max(hh_child, na.rm=T), by = c("SERIAL", "STATEFIP", "PUMA", "FAMUNIT")]
                  acs_full <- copy(acs)
                  acs <- acs[hh_child==1]
                  
                  child_df <- acs[AGE <= pred_age_max & AGE >= pred_age_min]
                  child_df <- child_df[,.(SERIAL, FAMUNIT, PERWT, PUMA, STATEFIP)]
                  
                  ## First analyze linked parents in ACS
                  acs <- acs[AGE<=pred_age_max & AGE>=pred_age_min, mom_loc_hh:=min(MOMLOC, na.rm=T), by = c("SERIAL", "FAMUNIT", "STATEFIP", "PUMA")]
                  acs <- acs[, mom_loc_hh:=min(mom_loc_hh, na.rm=T), by = c("SERIAL", "FAMUNIT", "STATEFIP", "PUMA")]
                  
                  acs <- acs[AGE<=pred_age_max & AGE>=pred_age_min, mom_loc2_hh:=min(MOMLOC2, na.rm=T), by = c("SERIAL", "FAMUNIT", "STATEFIP", "PUMA")]
                  acs <- acs[, mom_loc2_hh:=min(mom_loc2_hh, na.rm=T), by = c("SERIAL", "FAMUNIT", "STATEFIP", "PUMA")]
                  
                  acs <- acs[AGE<=pred_age_max & AGE>=pred_age_min, pop_loc_hh:=min(POPLOC, na.rm=T), by = c("SERIAL", "FAMUNIT", "STATEFIP", "PUMA")]
                  acs <- acs[, pop_loc_hh:=min(pop_loc_hh, na.rm=T), by = c("SERIAL", "FAMUNIT", "STATEFIP", "PUMA")]
                  
                  acs <- acs[AGE<=pred_age_max & AGE>=pred_age_min, pop_loc2_hh:=min(POPLOC2, na.rm=T), by = c("SERIAL", "FAMUNIT", "STATEFIP", "PUMA")]
                  acs <- acs[, pop_loc2_hh:=min(pop_loc2_hh, na.rm=T), by = c("SERIAL", "FAMUNIT", "STATEFIP", "PUMA")]
                  
                  pred_df <- pred_df[, geo_code:=as.numeric(fips)]
                  pred_df <- pred_df[, state_code:=as.numeric(as.character(state_code))]
                  
                  puma_map <- fread("./data/geocorr2018.csv")
                  puma_map <- puma_map[!(county_fips%in%c(46113, 2270, 51515))] ## Remove three counties that have been changed
                  puma_map <- puma_map[, sum_afact:=sum(afact), by = c("PUMA", "state_code")]
                  puma_map <- puma_map[, afact:=afact/sum_afact]
                  
                  parent_df <- acs[PERNUM==mom_loc_hh | PERNUM==pop_loc_hh | PERNUM==mom_loc2_hh | PERNUM==pop_loc2_hh]
                  parent_df <- parent_df[AGE>=18]
                  
                  parent_df <- parent_df[,.(SERIAL, PERNUM, sex, AGE, age_cat_fb, education, race_grp, PUMA, PERWT, FAMUNIT, STATEFIP)]
                  
                  fips <- as.data.table(tidycensus::fips_codes)
                  fips <- unique(fips[,.(state_code, state_name)])
                  fips <- fips[,state_code:=as.numeric(state_code)]
                  setnames(parent_df, "STATEFIP", "state_code")
                  
                  parent_df <- merge(parent_df, fips, by = "state_code")
                  
                  # pred_probs <- pred_df[sample==i]
                  pred_probs <- copy(pred_df)
                  
                  pred_df_out <- merge(parent_df, puma_map[,.(county_fips, PUMA, state_code, afact)], by = c("PUMA", "state_code"), all.x=T, allow.cartesian = T)
                  pred_df_out <- pred_df_out[, geo_code:=county_fips]
                  
                  pred_df_out <- merge(pred_df_out, pred_probs[,.(sex, age_cat_fb, education, race_grp, county, geo_code, state_code, county_code, state_name, pred_prob)],
                                       by = c("state_name", "geo_code", "sex", "age_cat_fb", "education", "race_grp", "state_code"), all.x=T)
                  
                  pred_df_out <- pred_df_out[, avg_prob:=mean(pred_prob, na.rm=T), by = c("SERIAL", "FAMUNIT", "county", "state_code")] ## Average of parents' hesitancy
                  pred_df_out <- unique(pred_df_out[,.(state_name, geo_code, state_code, SERIAL, FAMUNIT, county_fips, county, avg_prob, afact)])
                  
                  pred_summaries <- merge(child_df, pred_df_out, by = c("SERIAL", "FAMUNIT"), all=T, allow.cartesian = T)
                  
                  unmerged <- pred_summaries[is.na(avg_prob)]
                  unmerged <- unique(unmerged[,.(SERIAL, FAMUNIT)])
                  
                  pred_summaries_parent <- pred_summaries[!is.na(avg_prob)]
                  
                  ## Then analyze linked grandparents in ACS
                  grandparent <- merge(acs, unmerged, by = c("SERIAL", "FAMUNIT"))
                  grandparent <- grandparent[GCRESPON>=1]
                  grandparent <- grandparent[,.(SERIAL, PERNUM, sex, AGE, age_cat_fb, education, race_grp, PUMA, PERWT, FAMUNIT, STATEFIP)]
                  setnames(grandparent, "STATEFIP", "state_code")
                  grandparent <- merge(grandparent, fips, by = "state_code")
                  
                  pred_df_out <- merge(grandparent, puma_map[,.(county_fips, PUMA, state_code, afact)], by = c("PUMA", "state_code"), all.x=T, allow.cartesian = T)
                  pred_df_out <- pred_df_out[, geo_code:=county_fips]
                  
                  pred_df_out <- merge(pred_df_out, pred_probs[,.(sex, age_cat_fb, education, race_grp, county, geo_code, state_code, county_code, state_name, pred_prob)],
                                       by = c("state_name", "geo_code", "sex", "age_cat_fb", "education", "race_grp", "state_code"), all.x=T)
                  
                  pred_df_out <- pred_df_out[, avg_prob:=mean(pred_prob, na.rm=T), by = c("state_code", "county", "SERIAL", "FAMUNIT")] ## Average of parents' hesitancy
                  pred_df_out <- unique(pred_df_out[,.(state_name, geo_code, state_code, SERIAL, FAMUNIT, county_fips, county, avg_prob, afact)])
                  
                  child_df_unmerged <- merge(child_df, unmerged, by = c("SERIAL", "FAMUNIT"))
                  pred_summaries_grandparent <- merge(child_df_unmerged, pred_df_out, by = c("SERIAL", "FAMUNIT"), all=T, allow.cartesian = T)
                  
                  unmerged <- pred_summaries_grandparent[is.na(avg_prob)]
                  unmerged <- unique(unmerged[,.(SERIAL, FAMUNIT)])
                  
                  pred_summaries_grandparent <- pred_summaries_grandparent[!is.na(avg_prob)]
                  
                  ## Then revert to average of all adults in family unit
                  other_unit <- merge(acs, unmerged, by = c("SERIAL", "FAMUNIT"))
                  other_unit <- other_unit[AGE>=18]
                  other_unit <- other_unit[,.(SERIAL, PERNUM, sex, AGE, age_cat_fb, education, race_grp, PUMA, PERWT, FAMUNIT, STATEFIP)]
                  setnames(other_unit, "STATEFIP", "state_code")
                  other_unit <- merge(other_unit, fips, by = "state_code")
                  
                  pred_df_out <- merge(other_unit, puma_map[,.(county_fips, PUMA, state_code, afact)], by = c("PUMA", "state_code"), all.x=T, allow.cartesian = T)
                  pred_df_out <- pred_df_out[, geo_code:=county_fips]
                  
                  pred_df_out <- merge(pred_df_out, pred_probs[,.(sex, age_cat_fb, education, race_grp, county, geo_code, state_code, county_code, state_name, pred_prob)],
                                       by = c("state_name", "geo_code", "sex", "age_cat_fb", "education", "race_grp", "state_code"), all.x=T)
                  
                  pred_df_out <- pred_df_out[, avg_prob:=mean(pred_prob, na.rm=T), by = c("state_code", "county", "SERIAL", "FAMUNIT")] ## Average of parents' hesitancy
                  pred_df_out <- unique(pred_df_out[,.(state_name, geo_code, state_code, SERIAL, FAMUNIT, county_fips, county, avg_prob, afact)])
                  
                  child_df_unmerged <- merge(child_df, unmerged, by = c("SERIAL", "FAMUNIT"))
                  pred_summaries_other_unit <- merge(child_df_unmerged, pred_df_out, by = c("SERIAL", "FAMUNIT"), all=T, allow.cartesian = T)
                  
                  unmerged <- pred_summaries_other_unit[is.na(avg_prob)]
                  unmerged <- unique(unmerged[,.(SERIAL)])
                  
                  unmerged_keep_unit <- pred_summaries_other_unit[is.na(avg_prob)]
                  unmerged_keep_unit <- unique(unmerged_keep_unit[,.(SERIAL, FAMUNIT)])
                  
                  pred_summaries_other_unit <- pred_summaries_other_unit[!is.na(avg_prob)]
                  
                  ## Then revert to average of all adults in a home
                  other_nonunit <- merge(acs_full, unmerged, by = c("SERIAL"))
                  other_nonunit <- other_nonunit[AGE>=18]
                  other_nonunit <- other_nonunit[,.(SERIAL, PERNUM, sex, AGE, age_cat_fb, education, race_grp, PUMA, PERWT, STATEFIP)]
                  setnames(other_nonunit, "STATEFIP", "state_code")
                  other_nonunit <- merge(other_nonunit, fips, by = "state_code")
                  
                  pred_df_out <- merge(other_nonunit, puma_map[,.(county_fips, PUMA, state_code, afact)], by = c("PUMA", "state_code"), all.x=T, allow.cartesian = T)
                  pred_df_out <- pred_df_out[, geo_code:=county_fips]
                  
                  pred_df_out <- merge(pred_df_out, pred_probs[,.(sex, age_cat_fb, education, race_grp, county, geo_code, state_code, county_code, state_name, pred_prob)],
                                       by = c("state_name", "geo_code", "sex", "age_cat_fb", "education", "race_grp", "state_code"), all.x=T)
                  
                  pred_df_out <- pred_df_out[, avg_prob:=mean(pred_prob, na.rm=T), by = c("state_code", "county", "SERIAL")] ## Average of parents' hesitancy
                  pred_df_out <- unique(pred_df_out[,.(state_name, geo_code, state_code, SERIAL, county_fips, county, avg_prob, afact)])
                  
                  child_df_unmerged <- merge(child_df, unmerged_keep_unit, by = c("SERIAL", "FAMUNIT"))
                  pred_summaries_other_nonunit <- merge(child_df_unmerged, pred_df_out, by = c("SERIAL"), all=T, allow.cartesian = T)
                  
                  unmerged <- pred_summaries_other_nonunit[is.na(avg_prob)]
                  unmerged <- unique(unmerged[,.(SERIAL, FAMUNIT)]) ## THESE ARE ENTIRELY HOUSEHOLDS WITH AGE <=17
                  
                  pred_summaries_other_nonunit <- pred_summaries_other_nonunit[!is.na(avg_prob)]
                  
                  pred_summaries <- rbind(pred_summaries_parent, pred_summaries_grandparent, pred_summaries_other_unit, pred_summaries_other_nonunit)
                  
                  pred_summaries <- pred_summaries[, fb_pred:=sum(PERWT*avg_prob*afact)/sum(PERWT*afact), by = c("state_name", "geo_code", "county")]
                  pred_summaries <- pred_summaries[, tot_pop:=sum(PERWT*afact), by = c("state_name", "geo_code", "county")]
                  pred_summaries <- unique(pred_summaries[,.(state_name, geo_code, county_fips, county, state_code, fb_pred, tot_pop)])
                  
                  pred_summaries <- pred_summaries[, bootstrap_i:=big_loop]  
                  
                  write.csv(pred_summaries, paste0("./results/bootstrap_hesitancy_cells_childage_", pred_age_min, "_", pred_age_max, "/pred_", big_loop, ".csv"), na = "", row.names = F)
                  rm(list=setdiff(ls(), c("df", "fips", "mod", "mod_df", "my.cluster", "n.cores")))
                }

parallel::stopCluster(cl = my.cluster)
