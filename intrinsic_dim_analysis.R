
insula_hunger_paths <- c("Y:\\livneh\\itayta\\data\\IC19\\day_150911\\",
                         "Y:\\livneh\\itayta\\data\\IC17\\day_150615\\",
                         "Y:\\livneh\\itayta\\data\\IC13\\day_150406\\",
                         "Y:\\livneh\\itayta\\data\\IC13\\day_150407\\",
                         "Y:\\livneh\\itayta\\data\\IC32\\day_161214\\",
                         "Y:\\livneh\\itayta\\data\\IC42\\day_161117\\")

all <- get_all_paths("Y:\\livneh\\itayta\\data")

#IC44_days <- c("170518", "170523", "170519", "170524")
IC44_paths <- all[grepl("IC44", all) | grepl("IC52", all)]

#insula_thirst_paths <- IC44_paths[which(rowSums(sapply(IC44_days, function(d) {as.numeric(grepl(d, IC44_paths))})) > 0)]
insula_thirst_paths <- IC44_paths

v1_paths <- c("Y:\\livneh\\itayta\\v1_controls\\fov1\\day_140524\\",
              "Y:\\livneh\\itayta\\u1_controls\\fov3\\day_140920\\",
              "Y:\\livneh\\itayta\\v1_controls\\fov3\\day_140921\\",
              "Y:\\livneh\\itayta\\v1_controls\\fov5\\day_150723\\")

por_paths <- c("Y:\\livneh\\itayta\\por_controls\\fov1\\day_141023\\",
               "Y:\\livneh\\itayta\\por_controls\\fov2\\day_140805\\",
               "Y:\\livneh\\itayta\\por_controls\\fov3\\day_150411\\")
#"Y:\\livneh\\itayta\\por_controls\\fov4\\day_150112\\")


paths_all <- c(insula_hunger_paths,
               insula_thirst_paths,
               v1_paths,
               por_paths)


path_categories <- c(rep("Hunger", times=len(insula_hunger_paths)),
                     rep("Thirst", times=len(insula_thirst_paths)),
                     rep("V1", times=len(v1_paths)),
                     rep("POR", times=len(por_paths)))

categories <- path_categories
categories[categories %in% c("V1", "POR")] <- "Control"

id_analysis <- function(path) {
  
  
  id_dim_df <- c()
  ext=""
  for (idx in 1:len(paths_all)) {
    print(idx)
    
    
    work_path <- paths_all[idx]
    control <- (categories[idx] %in% c("Control"))#
    activity_thres <- ifelse(categories[idx] %in% c("Hunger", "Control"),
                              0.25,
                              0.2)
    
    # neur_mat_func <- ifelse(control,
    #                         get_reduced_mat_full_day_control,
    #                         get_reduced_mat_full_day)
    
    idims <- c()
    for (preset in c("hummus", "skippy", "jiffie", "nutella", "hashahar")) {
      
      fmat <- get_mat_with_preset(work_path, preset, oldscope = control, activity_threshold=activity_thres)
      
      trim_005 <-  twonn(fmat, method = "linfit", c_trimmed = 0.05)$est[2]
      trim_001 <-  twonn(fmat, method = "linfit", c_trimmed = 0.01)$est[2]
      trim_none <-  twonn(fmat, method = "linfit", c_trimmed = 0)$est[2]
      
      idims <- c(idims, c(trim_005, trim_001, trim_none))
    }
    
    
    id_dim_df <- rbind(id_dim_df, idims)
  }
}

