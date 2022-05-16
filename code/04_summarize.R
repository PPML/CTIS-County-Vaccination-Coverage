rm(list = ls())

library(data.table)
library(parallel)
library(foreach)
library(fst)

county_list <- read.fst("./results/coverage_preds_rdata_9mo/pred_1.fst", as.data.table = T)
county_list <- unique(county_list$state_name)

parallel::detectCores()
n.cores <- 24

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

for (m in c(3)) {
  
  ## START BOOTSTRAP LOOP ##
  out <- foreach (s = county_list, .combine = 'rbind',
                  .packages = c("data.table", "fst")) %dopar% {
                    
                      county_temp <- NULL
                      for (i in 1:1000) {
                        temp <- read.fst(paste0("./results/coverage_preds_rdata_", m, "mo/pred_", i, ".fst"), as.data.table = T)
                        temp <- temp[state_name == s]
                        county_temp <- rbind(county_temp, temp, fill = T)
                      }
                      
                      for (i in c(.025, .975)) {
                        county_temp <- county_temp[, paste0("percentile_", i*100, "_12_17"):=quantile(pred_12_17, p = i), by = "geo_code"]   
                        county_temp <- county_temp[, paste0("percentile_", i*100, "_5_11"):=quantile(pred_5_11, p = i), by = "geo_code"]   
                      }
                      county_temp <- county_temp[, mean_12_17:=mean(pred_12_17), by = "geo_code"]
                      county_temp <- county_temp[, mean_5_11:=mean(pred_5_11), by = "geo_code"]
                      keep_cols <- colnames(county_temp)
                      keep_cols <- keep_cols[keep_cols%in%c("state_name", "geo_code", "county") | keep_cols%like%"percentile" | keep_cols%like%"mean"]
                      county_temp <- unique(county_temp[,c(keep_cols), with = F])
                      
                      return(county_temp) 
                    
                  }
  
  write.csv(out, paste0("./results/summarized_coverage_estimates_", m, "mo.csv"), na = "", row.names = F)

state_temp <- NULL
for (i in 1:1000) {
  temp <- fread(paste0("./results/coverage_preds_state_", m, "mo/pred_", i, ".csv"))
  state_temp <- rbind(state_temp, temp, fill = T)
}

for (i in c(.025, .975)) {
  state_temp <- state_temp[, paste0("percentile_", i*100, "_12_17"):=quantile(pred_12_17, p = i), by = "state_code"]   
  state_temp <- state_temp[, paste0("percentile_", i*100, "_5_11"):=quantile(pred_5_11, p = i), by = "state_code"]   
}
state_temp <- state_temp[, mean_12_17:=mean(pred_12_17), by = "state_code"]
state_temp <- state_temp[, mean_5_11:=mean(pred_5_11), by = "state_code"]
keep_cols <- colnames(state_temp)
keep_cols <- keep_cols[keep_cols%in%c("state_name", "state_code") | keep_cols%like%"percentile" | keep_cols%like%"mean"]
state_temp <- unique(state_temp[,c(keep_cols), with = F])

write.csv(state_temp, paste0("./results/summarized_state_coverage_estimates_", m, "mo.csv"), na = "", row.names = F)

}

parallel::stopCluster(cl = my.cluster)



