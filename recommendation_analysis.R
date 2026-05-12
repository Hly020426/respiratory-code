suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

gc()

# =========================================================
# Multi-disease weighted recommendation
# =========================================================

# ====================== 1. Global configuration ======================
PROJECT_ROOT <- "."
DATA_DIR <- file.path(PROJECT_ROOT, "data")
LABEL_DIR <- file.path(DATA_DIR, "Moxingbiaoqian")

PATHS <- list(
  risk_scores = file.path(DATA_DIR, "J_risk_scores.csv"),
  movement_features = file.path(LABEL_DIR, "all_movement_features.csv"),
  group_file = file.path(LABEL_DIR, "group.xlsx"),
  output_file = file.path(DATA_DIR, "final_GAI_optimized.csv")
)

risk_percent <- 0.05
threshold <- 0.001
optimize_strength <- 0.20
hybrid_weight_opt <- 0.7
hybrid_weight_pop <- 0.3

night_hours <- c(22, 23, 0, 1, 2, 3, 4, 5, 6)
sleep_min_opt <- 0.7
mvpa_night_max_opt <- 0.1
sb_day_max_opt <- 0.8

# ====================== 2. Read data ======================
risk_data_raw <- fread(PATHS$risk_scores)
setnames(risk_data_raw, "Participant ID", "ID")

movement_data_raw <- fread(PATHS$movement_features)
setnames(movement_data_raw, names(movement_data_raw)[ncol(movement_data_raw)], "ID")

group_dt <- as.data.table(read_excel(PATHS$group_file))
setnames(group_dt, "Participant ID", "ID")

common_ids <- Reduce(intersect, list(movement_data_raw$ID, risk_data_raw$ID, group_dt$ID))
movement_data_raw <- movement_data_raw[ID %in% common_ids]
risk_data_raw <- risk_data_raw[ID %in% common_ids]
group_dt <- group_dt[ID %in% common_ids]

behavior_cols <- setdiff(names(movement_data_raw), "ID")
disease_cols <- grep("^J", names(risk_data_raw), value = TRUE)

# ====================== 3. Behavior column mapping ======================
map_dt <- data.table(original_col = behavior_cols)
map_dt[, idx := .I - 1]
map_dt[, activity := c("SB", "LPA", "MVPA", "Sleep")[idx %% 4 + 1]]
map_dt[, hour_raw := idx %/% 4]
map_dt[, day_type := ifelse(hour_raw < 24, "Weekday", "Weekend")]
map_dt[, hour := hour_raw %% 24]
map_dt[, idx := NULL]

# ====================== 4. Long-format movement data ======================
movement_long <- melt(movement_data_raw,
                      id.vars="ID",
                      measure.vars=behavior_cols,
                      variable.name="original_col",
                      value.name="value")

movement_long <- merge(movement_long, map_dt, by="original_col")
movement_long[, total := sum(value), by=.(ID,day_type,hour)]
movement_long[, prop := fifelse(total==0,0,value/total)]

movement_long <- merge(movement_long, map_dt, by = "original_col")
movement_long[, total := sum(value), by = .(ID, day_type, hour)]
movement_long[, prop := fifelse(total == 0, 0, value / total)]

# ====================== 5. Direction by disease ======================
get_direction <- function(disease_code){
  
  dt <- risk_data_raw[,.(ID, risk=get(disease_code))]
  dt <- dt[!is.na(risk)]
  
  low_t  <- quantile(dt$risk, risk_percent)
  high_t <- quantile(dt$risk, 1-risk_percent)
  dt[, grp := fifelse(risk<=low_t,"Low","High")]
  
  mv <- merge(movement_long, dt[,.(ID,grp)], by="ID")
  
  low_p  <- mv[grp=="Low",  .(P_low=mean(prop)),  by=.(day_type,hour,activity)]
  high_p <- mv[grp=="High", .(P_high=mean(prop)), by=.(day_type,hour,activity)]
  
  opt <- merge(low_p, high_p, by=c("day_type","hour","activity"))
  opt[, delta := P_low - P_high]
  
  opt[, dir :=
        fifelse(abs(delta)<threshold,"Retain",
                fifelse(delta>0,"Strengthen","Compress"))]
  
  return(opt[,.(day_type,hour,activity,dir)])
}

# ====================== 6. Optimize one disease ======================
optimize_one_disease <- function(dir_dt){
  
  movement_data <- copy(movement_data_raw)
  population_avg <- movement_data[, lapply(.SD, mean), .SDcols=behavior_cols]
  
  adjust_proportion <- function(p, dir){
    if(dir=="Strengthen") return(p*(1+optimize_strength))
    if(dir=="Compress")   return(p*(1-optimize_strength))
    return(p)
  }
  
  for(i in 1:nrow(map_dt)){
    colname <- map_dt$original_col[i]
    act  <- map_dt$activity[i]
    hr   <- map_dt$hour[i]
    dayt <- map_dt$day_type[i]
    
    dir_row <- dir_dt[day_type==dayt & hour==hr & activity==act]
    if(nrow(dir_row)==0) next
    
    movement_data[, (colname) := adjust_proportion(get(colname), dir_row$dir)]
    movement_data[, (colname) := hybrid_weight_opt*get(colname) +
                    hybrid_weight_pop*population_avg[[colname]]]
  }
  
  for(dayt in c("Weekday","Weekend")){
    for(hr in 0:23){
      
      cols_now <- map_dt[day_type==dayt & hour==hr, original_col]
      tmp <- movement_data[, ..cols_now]
      setnames(tmp, c("SB","LPA","MVPA","Sleep"))
      
      s <- rowSums(tmp)
      tmp <- tmp / ifelse(s==0,1,s)
      
      if(hr %in% night_hours){
        tmp$Sleep <- pmax(tmp$Sleep, sleep_min_opt)
        tmp$MVPA  <- pmin(tmp$MVPA, mvpa_night_max_opt)
        
        remain <- 1 - tmp$Sleep - tmp$MVPA
        remain[remain < 0] <- 0
        
        other_sum <- tmp$SB + tmp$LPA
        tmp$SB  <- remain * tmp$SB/ifelse(other_sum==0,1,other_sum)
        tmp$LPA <- remain * tmp$LPA/ifelse(other_sum==0,1,other_sum)
        
      } else {
        tmp$SB <- pmin(tmp$SB, sb_day_max_opt)
        
        remain <- 1 - tmp$SB
        other_sum <- tmp$LPA + tmp$MVPA + tmp$Sleep
        
        tmp$LPA   <- remain * tmp$LPA/ifelse(other_sum==0,1,other_sum)
        tmp$MVPA  <- remain * tmp$MVPA/ifelse(other_sum==0,1,other_sum)
        tmp$Sleep <- remain * tmp$Sleep/ifelse(other_sum==0,1,other_sum)
      }
      
      tmp <- tmp / rowSums(tmp)
      movement_data[, (cols_now) := tmp]
    }
  }
  
  return(movement_data)
}

# ====================== 7. Disease-specific recommendations ======================
cat("🔄 Generating disease-specific recommendations...\n")
disease_recommend_list <- list()

for(d in disease_cols){
  cat("   →", d, "\n")
  dir_dt <- get_direction(d)
  rec_dt <- optimize_one_disease(dir_dt)
  rec_dt[, disease := d]
  disease_recommend_list[[d]] <- rec_dt
}

all_recommend_dt <- rbindlist(disease_recommend_list)
# ====================== 8. Disease weights ======================
disease_weight_dt <- melt(group_dt,
                          id.vars="ID",
                          variable.name="disease",
                          value.name="has_disease")[has_disease==1]

disease_weight_dt <- disease_weight_dt[, .N, by=disease]
disease_weight_dt[, weight := N/sum(N)]

all_recommend_dt <- merge(all_recommend_dt,
                          disease_weight_dt[,.(disease,weight)],
                          by="disease")

all_recommend_dt[, weight_norm := weight / sum(weight), by=ID]

final_recommend <- all_recommend_dt[, lapply(.SD, weighted.mean, w=weight_norm),
                                    by=ID,
                                    .SDcols=behavior_cols]

# ====================== 9. Add healthy participants back ======================
healthy_ids <- setdiff(movement_data_raw$ID, final_recommend$ID)

if(length(healthy_ids) > 0){
  healthy_dt <- movement_data_raw[ID %in% healthy_ids]
  final_recommend <- rbind(final_recommend, healthy_dt, fill=TRUE)
}

# ====================== 10. Final normalization ======================
for(dayt in c("Weekday","Weekend")){
  for(hr in 0:23){
    cols_now <- map_dt[day_type==dayt & hour==hr, original_col]
    final_recommend[, (cols_now) := {
      tmp <- .SD
      tmp / rowSums(tmp)
    }, .SDcols=cols_now]
  }
}


setnames(final_recommend, "ID", "Participant ID")

all_cols <- names(final_recommend)
cols_except_id <- setdiff(all_cols, "Participant ID")
new_col_order <- c(cols_except_id, "Participant ID")

setcolorder(final_recommend, new_col_order)

# ====================== 11. Export ======================
fwrite(final_recommend, PATHS$output_file, na = "0")
cat("Multi-disease weighted recommendation completed. Output saved.\n")
