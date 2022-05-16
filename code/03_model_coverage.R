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

parallel::detectCores()
n.cores <- 40

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

## START BOOTSTRAP LOOP ##
out <- foreach (big_loop = 1:1000, .combine = 'rbind',
                .packages = c("data.table", "tidycensus", "boot", "lme4", "survey", "fitdistrplus", "fst")) %dopar% {
                  
                  set.seed(big_loop)
                  
                  preds <- fread(paste0("./results/bootstrap_hesitancy_cells_childage_12_17/pred_", big_loop, ".csv"))
                  preds_kids <- fread(paste0("./results/bootstrap_hesitancy_cells_childage_5_11/pred_", big_loop, ".csv"))
                  
                  ## FROM HESITANCY ESTIMATES, RUN SECOND REGRESSION TO COVERAGE ##
                  cdc <- fread("./data/COVID-19_Vaccinations_in_the_United_States_County.csv")
                  cdc <- cdc[,.(Date, FIPS, Series_Complete_12PlusPop_Pct, Series_Complete_18Plus, Series_Complete_18PlusPop_Pct, Series_Complete_12Plus, Completeness_pct)]
                  # cdc <- cdc[Date=="02/01/2022"] ## Average 9 months post 12-15 EUA / 16-17 eligibility
                  # cdc <- cdc[Date=="08/01/2021"] ## average 3 months post 12-15 EUA / 16-17 eligibility
                  cdc <- cdc[, vax_pct_12_17:=(Series_Complete_12Plus-Series_Complete_18Plus)/((Series_Complete_12Plus/Series_Complete_12PlusPop_Pct)-(Series_Complete_18Plus/Series_Complete_18PlusPop_Pct))]
                  cdc <- cdc[, cdc_pop:=(Series_Complete_12Plus-Series_Complete_18Plus)/(vax_pct_12_17/100)]
                  cdc <- cdc[, geo_code:=as.numeric(FIPS)]
                  
                  preds <- merge(preds, cdc, by = "geo_code")
                  
                  preds <- preds[, vax_pct_12_17:=vax_pct_12_17/100]
                  preds <- preds[is.na(cdc_pop) | is.infinite(cdc_pop), cdc_pop:=tot_pop]
                  
                  preds <- preds[county_fips%in%c(25019, 25001, 25007), Completeness_pct:=0]
                  preds <- preds[logit(vax_pct_12_17) < -4 | logit(vax_pct_12_17) > 3, Completeness_pct:=0]
                  preds <- preds[Completeness_pct>=80 & vax_pct_12_17<1 & vax_pct_12_17>0 & !is.na(vax_pct_12_17),
                                 county_n := .N, by = c("bootstrap_i", "state_name")]
                  preds <- preds[is.na(county_n), county_n := 0]
                  preds <- preds[, county_n:=max(county_n), by = c("bootstrap_i", "state_name")]
                  
                  regions <- fread("https://raw.githubusercontent.com/cphalpert/census-regions/master/us%20census%20bureau%20regions%20and%20divisions.csv")
                  regions <- regions[, geo_code:=tolower(`State Code`)]
                  
                  preds <- merge(preds, regions[,.(State, Division)], by.x = "state_name", by.y = "State", all.x=T)
                  
                  pred_coverage <- preds[,.(state_name, geo_code, county_fips, county, state_code,
                                            fb_pred, bootstrap_i, Completeness_pct, vax_pct_12_17,
                                            county_n, Division, cdc_pop)]
                  setnames(preds_kids, "fb_pred", "fb_pred_5_11")
                  pred_coverage <- merge(pred_coverage, preds_kids[,.(geo_code, fb_pred_5_11, bootstrap_i)], by = c("geo_code", "bootstrap_i"))
                  expand <- as.data.table(expand.grid(geo_code = unique(pred_coverage$geo_code), boot_i = 1:1000)) 
                  pred_coverage <- merge(pred_coverage, expand, by = "geo_code", allow.cartesian = T)
                  
                  state_preds <- unique(preds$state_name[preds$county_n>=10])
                  
                  design <- svydesign(id = ~1, weights = ~cdc_pop, data = preds[county_n >=10 & 
                                                                                  Completeness_pct>=80 &
                                                                                  vax_pct_12_17<1 & 
                                                                                  vax_pct_12_17>0 &
                                                                                  !is.na(vax_pct_12_17)])
                  mod <- svyglm(design = design,
                                formula = vax_pct_12_17 ~ logit(fb_pred)*state_name,
                                family = binomial(link = "logit"))
                  
                  betas = coef(mod)
                  vcov = vcov(mod)
                  betas_simulated = MASS::mvrnorm(1000, betas, vcov)
                  
                  betas_simulated <- as.data.table(betas_simulated)
                  ref_grp <- names(betas)[!(names(betas)%like%"logit")]
                  ref_grp <- ref_grp[ref_grp%like%"state_name"]
                  ref_grp <- gsub("state_name", "", ref_grp)
                  ref_grp <- state_preds[!(state_preds%in%ref_grp)]
                  betas_simulated <- betas_simulated[, paste0("logit(fb_pred):state_name", ref_grp):=0]
                  betas_simulated <- betas_simulated[, paste0("state_name", ref_grp):=0]
                  setnames(betas_simulated, "logit(fb_pred)", "fb_pred")
                  betas_simulated <- melt(betas_simulated, id.vars = c("(Intercept)", "fb_pred"))
                  betas_simulated <- betas_simulated[variable%like%"logit", term:="state_slope"]
                  betas_simulated <- betas_simulated[!(variable%like%"logit"), term:="state_intercept"]
                  betas_simulated <- betas_simulated[, variable:=gsub(pattern = "state_name", 
                                                                      replacement = "",
                                                                      x = variable)]
                  betas_simulated <- betas_simulated[, variable:=gsub(pattern = "logit(fb_pred):", 
                                                                      replacement = "",
                                                                      x = variable, fixed = TRUE)]
                  betas_simulated <- dcast(betas_simulated, formula = `(Intercept)` + fb_pred + variable ~ term, value.var="value")
                  setnames(betas_simulated, c("variable", "(Intercept)", "fb_pred"), c("state_name", "global_intercept", "global_slope"))
                  betas_simulated <- betas_simulated[, boot_i:=seq_len(.N), by = "state_name"]
                  
                  pred_coverage <- merge(pred_coverage, betas_simulated, by = c("state_name", "boot_i"), all.x=T)
                  
                  sd_dists <- cbind(preds[county_n >=10 & 
                                            Completeness_pct>=80 &
                                            vax_pct_12_17<1 & 
                                            vax_pct_12_17>0 &
                                            !is.na(vax_pct_12_17)], mod$residuals)
                  
                  for (s_small in state_preds) {
                    sd_dist_fit <- fitdist(sample(size = 10000, x = sd_dists$V2[!is.na(sd_dists$V2) & sd_dists$state_name==s_small],
                                                  prob = sd_dists$cdc_pop[!is.na(sd_dists$V2) & sd_dists$state_name==s_small],
                                                  replace = T), "norm")$estimate[["sd"]]
                    pred_coverage <- pred_coverage[state_name == s_small, error:=rnorm(.N, mean = 0, sd = sd_dist_fit)]
                    
                  }
                  
                  pred_coverage <- pred_coverage[state_name%in%state_preds, pred_12_17:=inv.logit(global_intercept + global_slope*logit(fb_pred) + state_intercept + state_slope*logit(fb_pred) + error)]
                  pred_coverage <- pred_coverage[state_name%in%state_preds, pred_5_11:=inv.logit(global_intercept + global_slope*logit(fb_pred_5_11) + state_intercept + state_slope*logit(fb_pred_5_11) + error)]
                  pred_coverage <- pred_coverage[, c("global_intercept", "global_slope", "state_intercept", "state_slope", "error"):=NULL]
                  
                  design <- svydesign(id = ~1, weights = ~cdc_pop, data = preds[Completeness_pct>=80 &
                                                                                  vax_pct_12_17<1 & 
                                                                                  vax_pct_12_17>0 & 
                                                                                  !is.na(vax_pct_12_17)])
                  mod <- svyglm(design = design,
                                formula = vax_pct_12_17 ~ logit(fb_pred)*Division,
                                family = binomial(link = "logit"))
                  
                  betas = coef(mod)
                  vcov = vcov(mod)
                  betas_simulated = MASS::mvrnorm(1000, betas, vcov)
                  
                  betas_simulated <- as.data.table(betas_simulated)
                  ref_grp <- names(betas)[!(names(betas)%like%"logit")]
                  ref_grp <- ref_grp[ref_grp%like%"Division"]
                  ref_grp <- gsub("Division", "", ref_grp)
                  ref_grp <- unique(regions$Division)[!(unique(regions$Division)%in%ref_grp)]
                  betas_simulated <- betas_simulated[, paste0("logit(fb_pred):Division", ref_grp):=0]
                  betas_simulated <- betas_simulated[, paste0("Division", ref_grp):=0]
                  setnames(betas_simulated, "logit(fb_pred)", "fb_pred")
                  betas_simulated <- melt(betas_simulated, id.vars = c("(Intercept)", "fb_pred"))
                  betas_simulated <- betas_simulated[variable%like%"logit", term:="state_slope"]
                  betas_simulated <- betas_simulated[!(variable%like%"logit"), term:="state_intercept"]
                  betas_simulated <- betas_simulated[, variable:=gsub(pattern = "Division", 
                                                                      replacement = "",
                                                                      x = variable)]
                  betas_simulated <- betas_simulated[, variable:=gsub(pattern = "logit(fb_pred):", 
                                                                      replacement = "",
                                                                      x = variable, fixed = TRUE)]
                  betas_simulated <- dcast(betas_simulated, formula = `(Intercept)` + fb_pred + variable ~ term, value.var="value")
                  setnames(betas_simulated, c("variable", "(Intercept)", "fb_pred"), c("Division", "global_intercept", "global_slope"))
                  betas_simulated <- betas_simulated[, boot_i:=seq_len(.N), by = "Division"]
                  
                  pred_coverage <- merge(pred_coverage, betas_simulated, by = c("Division", "boot_i"), all.x=T)
                  
                  sd_dists <- cbind(preds[Completeness_pct>=80 &
                                            vax_pct_12_17<1 & 
                                            vax_pct_12_17>0 &
                                            !is.na(vax_pct_12_17)], mod$residuals)
                  
                  for (d in unique(preds$Division)) {
                    sd_dist_fit <- fitdist(sample(size = 10000, x = sd_dists$V2[!is.na(sd_dists$V2) & sd_dists$Division==d],
                                                  prob = sd_dists$cdc_pop[!is.na(sd_dists$V2) & sd_dists$Division==d],
                                                  replace = T), "norm")$estimate[["sd"]]
                    pred_coverage <- pred_coverage[Division == d, error:=rnorm(.N, mean = 0, sd = sd_dist_fit)]
                    
                  }
                  
                  pred_coverage <- pred_coverage[!(state_name%in%state_preds), pred_12_17:=inv.logit(global_intercept + global_slope*logit(fb_pred) + state_intercept + state_slope*logit(fb_pred) + error)]
                  pred_coverage <- pred_coverage[!(state_name%in%state_preds), pred_5_11:=inv.logit(global_intercept + global_slope*logit(fb_pred_5_11) + state_intercept + state_slope*logit(fb_pred_5_11) + error)]
                  pred_coverage <- pred_coverage[, c("global_intercept", "global_slope", "state_intercept", "state_slope", "error"):=NULL]

                  write.fst(pred_coverage, paste0("./results/coverage_preds_rdata_4mo/pred_", big_loop, ".fst"))
                  
                  state_coverage <- copy(pred_coverage)
                  state_coverage <- state_coverage[, fb_pred:=sum(fb_pred*cdc_pop, na.rm=T)/sum(cdc_pop, na.rm=T), by = c("state_name", "boot_i")]
                  state_coverage <- state_coverage[, vax_pct_12_17 :=sum(vax_pct_12_17*cdc_pop, na.rm=T)/sum(cdc_pop, na.rm=T), by = c("state_name", "boot_i")]
                  state_coverage <- state_coverage[, fb_pred_5_11 :=sum(fb_pred_5_11*cdc_pop, na.rm=T)/sum(cdc_pop, na.rm=T), by = c("state_name", "boot_i")]
                  state_coverage <- state_coverage[, pred_12_17 :=sum(pred_12_17*cdc_pop, na.rm=T)/sum(cdc_pop, na.rm=T), by = c("state_name", "boot_i")]
                  state_coverage <- state_coverage[, pred_5_11  :=sum(pred_5_11*cdc_pop, na.rm=T)/sum(cdc_pop, na.rm=T), by = c("state_name", "boot_i")]
                  state_coverage <- state_coverage[, n_counties:=.N, by = c("state_name", "boot_i")]
                  
                  state_coverage <- unique(state_coverage[,.(state_name, state_code, fb_pred, vax_pct_12_17, fb_pred_5_11, pred_12_17, pred_5_11, boot_i, bootstrap_i, n_counties)])
                  
                  nat_coverage <- copy(pred_coverage)
                  nat_coverage <- nat_coverage[, fb_pred:=sum(fb_pred*cdc_pop, na.rm=T)/sum(cdc_pop, na.rm=T), by = c("boot_i")]
                  nat_coverage <- nat_coverage[, vax_pct_12_17 :=sum(vax_pct_12_17*cdc_pop, na.rm=T)/sum(cdc_pop, na.rm=T), by = c("boot_i")]
                  nat_coverage <- nat_coverage[, fb_pred_5_11 :=sum(fb_pred_5_11*cdc_pop, na.rm=T)/sum(cdc_pop, na.rm=T), by = c("boot_i")]
                  nat_coverage <- nat_coverage[, pred_12_17 :=sum(pred_12_17*cdc_pop, na.rm=T)/sum(cdc_pop, na.rm=T), by = c("boot_i")]
                  nat_coverage <- nat_coverage[, pred_5_11  :=sum(pred_5_11*cdc_pop, na.rm=T)/sum(cdc_pop, na.rm=T), by = c("boot_i")]
                  nat_coverage <- nat_coverage[, n_counties:=.N, by = c("boot_i")]
                  
                  nat_coverage <- unique(nat_coverage[,.(fb_pred, vax_pct_12_17, fb_pred_5_11, pred_12_17, pred_5_11, boot_i, bootstrap_i, n_counties)])
                  nat_coverage <- nat_coverage[, state_name:="United States"]
                  nat_coverage <- nat_coverage[, state_code:=0]
                  
                  state_coverage <- rbind(state_coverage, nat_coverage)
                  
                  write.csv(state_coverage, paste0("./results/coverage_preds_state_4mo/pred_", big_loop, ".csv"), na = "", row.names = F)
                  
                  
                  
                  rm(list=setdiff(ls(), c("df", "fips", "mod", "mod_df", "my.cluster", "n.cores")))
                }

parallel::stopCluster(cl = my.cluster)
