library(RColorBrewer)
library(intRinsic)
insula_hunger_paths <- c("Y:\\livneh\\itayta\\data\\IC19\\day_150911\\",
                         "Y:\\livneh\\itayta\\data\\IC17\\day_150615\\",
                         "Y:\\livneh\\itayta\\data\\IC13\\day_150406\\",
                         "Y:\\livneh\\itayta\\data\\IC13\\day_150407\\",
                         "Y:\\livneh\\itayta\\data\\IC32\\day_161214\\",
                         "Y:\\livneh\\itayta\\data\\IC42\\day_161117\\")

all <- get_all_paths("Y:\\livneh\\itayta\\data")

#IC44_days <- c("170518", "170523", "170519", "170524")
IC44_paths <- all[grepl("IC44", all) | grepl("IC52", all) | grepl("IC47", all)]

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

paths_all <- paths_all[-15]
path_categories <- path_categories[-15]

categories <- path_categories
categories[categories %in% c("V1", "POR")] <- "Control"

id_analysis <- function(path) {
  
  
  id_dim_df <- c()
  ext=""
  for (preset in c("skippy", "dalnatran", "jiffie", "bnd", "shufersal",   "hummus", "nitzat", "taaman")){ #,"milka", "hashahar")) {
    idims <- c()
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
    
    mat <- get_mat_with_preset(work_path, "hummus", oldscope=control, activity_threshold=activity_thres)
    
    if (control) {
      stim_master_mat <- 
        get_color_palettes_control(work_path, just_mat = T, window_size = 15)
    } else {
      stim_master_mat <- 
        get_stim_mat(work_path,just_mat = T, window_size = 15)             
    }
    
    spec_cg <- colorRampPalette(rev(brewer.pal(n = 11,  name = "Spectral")))
    behav <- get_behavior_variables(stim_master_mat, mat_size = nrow(mat), window_size = 15, seconds_after=1)
    cp_base <- viridis(nrow(mat))
    reward_col <- adjustcolor(cp_base, alpha=0.1); reward_col[behav[,"Reward"] == 1] <- "red"
    tt3_col <- adjustcolor(cp_base, alpha=0.1); tt3_col[behav[,"TrialType 3"] == 1] <- "red"
    tt4_col <- adjustcolor(cp_base, alpha=0.1); tt4_col[behav[,"TrialType 4"] == 1] <- "red"
    tt5_col <- adjustcolor(cp_base, alpha=0.1); tt5_col[behav[,"TrialType 5"] == 1] <- "red"
    cumreward_col <- adjustcolor(rev(spec_cg(len(unique(behav[,"Cumulative reward"]))))[behav[,"Cumulative reward"]],
                                 0.1)
    
    cols = list(Time=cp_base,
             Reward=reward_col,
             TrialType3=tt3_col,
             TrialType4=tt4_col,
             TrialType5=tt5_col,
             CumulReward=cumreward_col)
    
    
    gen_pdf=F
    

    
    dir.create("~/IntrinsicDim")
      
      fmat <- get_mat_with_preset(work_path, preset, oldscope = control, activity_threshold=activity_thres)
      
      
      if (path_categories[idx] %in% c("V1", "POR")) {
        pattern = "fov"
      } else {
        pattern = "IC"
      }
      
      split_r <- str_split(str_split(work_path, pattern)[[1]][[2]], "\\\\")
      
      trim_005 <-  twonn(fmat, method = "linfit", c_trimmed = 0.05)$est[2]
      trim_001 <-  twonn(fmat, method = "linfit", c_trimmed = 0.01)$est[2]
      trim_none <-  twonn(fmat, method = "linfit", c_trimmed = 0)$est[2]
      
      idims <- rbind(idims, c(trim_005, trim_001, trim_none))
      
      dir.create(sprintf("~/IntrinsicDim/%s%s", pattern, split_r[[1]][1]) )
      day_path <- sprintf("~/IntrinsicDim/%s%s/%s", pattern, split_r[[1]][1], split_r[[1]][[2]])
      dir.create(day_path)
      
      if (preset != "hummus") {
        for (var_name in names(cols)) {
          final_path <- sprintf("%s/%s", day_path, var_name)
          print(sprintf("Creating %s!", final_path))
          dir.create(final_path)
          var_col <- cols[[var_name]]
          
          
          png(file=sprintf("%s/%s.png", final_path, preset),height=10,width=10,
              res=300, units="in")
          pairs(fmat, col=var_col, pch=19)
          dev.off()
          
          png(file=sprintf("%s/%s_small.png", final_path, preset),height=5,width=5,
              res=300, units="in")
          pairs(fmat, col=var_col, pch=19)
          dev.off()
          
          pdf(sprintf("%s/%s.pdf", final_path, preset),height=10,width=10)
          pairs(fmat, col=var_col, pch=19)
          dev.off()
          
          downs <- sample(size = 3000, x=1:nrow(fmat))
          
          pdf(sprintf("%s/%s_downsampled.pdf", final_path, preset),height=6,width=6, pointsize = 0.5)
          pairs(fmat[downs,], col=var_col[downs], pch=19)
          dev.off()
        }
      }
    }
    
    
    id_dim_df <- cbind(id_dim_df, idims)
  }
}

