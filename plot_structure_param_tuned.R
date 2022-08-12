library("tidyverse")

mat_frames = 57600

plot_structure <- function(path, window_size = 30) {
  
  paths <- list.dirs(path, recursive=F)[grep("IC", list.dirs(path, recursive=F))]
  path_pairs <- list()
  
  # Get all pairs of possible paths
  for (p in paths) {
    mice_num <- str_split(p, "IC")[[1]][[2]]
    runs <- list.dirs(p, recursive = F)
    days <- sapply(str_split(runs, "day_"), function(l) {return(l[[2]])})
    
    pairs_mat <- apply(combn(days,2), 2, 
                       function(p) {
                         return(sprintf("%s\\IC%s\\day_%s", path, mice_num, p))
                       })
    
    path_pairs <-
      append(path_pairs,
             lapply(seq_len(ncol(pairs_mat)), function(i) pairs_mat[,i]))
  }
  
  paths_all <- unique(unlist(path_pairs))
  
  
  knn1_range <- c(seq(0.2, 0.275, by=0.025), seq(0.325, 0.425, by=0.05))
  knn2_range <- seq(0.05, 0.075, by=0.025)
  print(knn1_range)
  print(knn2_range)
  
  
      for (path_to_tune in paths_all) {
        
        params_mat <- cross_df(list(knn1=knn1_range, knn2=knn2_range))
        params_mat <- params_mat[-8,]
        
        mat <- get_reduced_mat_full_day(path_to_tune, 
                                        knn1 = params_mat[1,]$knn1, 
                                        knn2 = params_mat[1,]$knn2,
                                        window_size = window_size)
        
        session_size <- mat_frames / window_size
        runs <- dim(mat)[1] / (session_size)
        
        behav_files <- list.files(sprintf("%s\\behavior\\",path_to_tune))
        behavior_mat_list <- lapply(behav_files,
                                    function(tv_mat) {
                                      return(readMat(sprintf("%s\\behavior\\%s",
                                                             path_to_tune,
                                                             tv_mat))$Stim.Master)
                                    })
        
        # Good luck reading this lol
        behavior_indices <- 
          unlist(lapply(unlist(lapply(str_split(behav_files, ".TrialVar"), 
                                      function(sp) {
                                        sp[[1]][1]
                                       })),
                        function(p) {
                            tmp <- strsplit(p, "")[[1]]; 
                            return(as.numeric(paste(tmp[(length(tmp) - 1):(length(tmp))], collapse="")))
                            }))
        
        
        cp_list <- 
          get_fixed_color_palette(behav_mat_list = behavior_mat_list,
                                  behavior_indices = behavior_indices,
                                  runs = runs,
                                  mat_frames = mat_frames,
                                  session_size=nrow(mat))

        
        use_path <- sprintf("%s\\%s", path_to_tune, "tune_pairs")
        dir.create(use_path)
        for (i in 1:nrow(params_mat)) {
          
          reduced_mat <- 
          get_reduced_mat_full_day(path_to_tune,
                                   knn1=params_mat[i,]$knn1,
                                   knn2=params_mat[i,]$knn2,
                                   window_size = window_size)
          
          png(sprintf("%s\\lem_knn1%.3f_knn2%.3f.png",
                      use_path,
                      round(params_mat[i,]$knn1 ,digits=3),
                      round(params_mat[i,]$knn2, digits=3)),
              width=2500, height=2500)
          
          pairs(reduced_mat, col=cp_list$time_cp, pch=19, cex=3.2)
          dev.off()
          
          png(sprintf("%s\\lem_knn1%.3f_knn2%.3f_response.png",
                      use_path,
                      round(params_mat[i,]$knn1 ,digits=3),
                      round(params_mat[i,]$knn2, digits=3)),
              width=2500, height=2500)
          
          pairs(reduced_mat, col=cp_list$response_cp, pch=19, cex=3.2)
          dev.off()          
          
          
        }
        

      }
}

get_fixed_color_palette <- function(behav_mat_list, mat_frames, behavior_indices, runs, session_size) {
  
  fm <- behav_mat_list
  ret_list <- list()
  
  
  for (i in 1:runs) {
    if (!i %in% behavior_indices) {
      print("Hello!")
      ret_list <- append(ret_list, list(NA))
      next
    }
    
    j <- which(i == behavior_indices)
    
    fm[[j]][,1] <- fm[[j]][,1] + (i-1) * mat_frames
    fm[[j]][,6] <- fm[[j]][,6] + (i-1) * mat_frames
    ret_list <- append(ret_list, list(fm[[j]]))
  }

  
  
  response_ind <- unlist(lapply(ret_list, function(sdf) {
              if (!all(is.na(sdf))) { return(sdf[,6])}
              }))
  
  response_ind <- response_ind[!is.nan(response_ind)]
  binned_trials <- get_binned_index(response_ind, window_size)
  
  time_cp=viridis(session_size)
  response_cp = adjustcolor(time_cp, alpha=0.2)
  response_cp[binned_trials] <- "red"
  return(list(time_cp=time_cp, response_cp=response_cp))
}


tune_similarities <- function(path, window_size = 30, nclusters=100) {
  
  paths <- list.dirs(path, recursive=F)[grep("IC", list.dirs(path, recursive=F))]
  path_pairs <- list()
  
  # Get all pairs of possible paths
  for (p in paths) {
    mice_num <- str_split(p, "IC")[[1]][[2]]
    runs <- list.dirs(p, recursive = F)
    days <- sapply(str_split(runs, "day_"), function(l) {return(l[[2]])})
    
    pairs_mat <- apply(combn(days,2), 2, 
                       function(p) {
                         return(sprintf("%s\\IC%s\\day_%s", path, mice_num, p))
                       })
    
    path_pairs <-
      append(path_pairs,
             lapply(seq_len(ncol(pairs_mat)), function(i) pairs_mat[,i]))
  }
  
  paths_all <- unique(unlist(path_pairs))
  
  
  knn1_range <- c(seq(0.225, 0.275, by=0.025), seq(0.325, 0.425, by=0.05))
  knn2_range <-c(0.075)
  print(knn1_range)
  print(knn2_range)
  
  
  structure_sim_mat <- c()
  
  for (path_to_tune in paths_all) {
    
    params_mat <- cross_df(list(knn1=knn1_range, knn2=knn2_range))

    for (i in 1:nrow(params_mat)) {
      
      reduced_mat <- 
        get_reduced_mat_full_day(path_to_tune,
                                 knn1=params_mat[i,]$knn1,
                                 knn2=params_mat[i,]$knn2,
                                 window_size = window_size)
      
      
      
      
      for (path_to_tune2 in paths_all) {
        
        params_mat <- cross_df(list(knn1=knn1_range, knn2=knn2_range))
        params_mat <- params_mat[-8,]
        
        if (path_to_tune2 == path) { 
          
        }
        
        for (j in 1:nrow(params_mat)) {
          
          if (path_to_tune == path_to_tune2) {
            cor_level = 0
          } else {
            reduced_mat2 <- 
              get_reduced_mat_full_day(path_to_tune2,
                                       knn1=params_mat[j,]$knn1,
                                       knn2=params_mat[j,]$knn2,
                                       window_size = window_size)
            
           
            km_day1 <- kmeans(apply(reduced_mat, 2, scale), nclusters, iter.max=300)
            km_day2 <- kmeans(apply(reduced_mat2, 2, scale), nclusters, iter.max=300)
  
            dist_matrix_day1 <- as.matrix(dist(km_day1$centers))
            dist_matrix_day2 <- as.matrix(dist(km_day2$centers))
            
            h1 <- hclust(dist(dist_matrix_day1))
            h2 <- hclust(dist(dist_matrix_day2))
            
            cor_level <- cor(c(dist_matrix_day1[h1$order, h1$order]), 
                             c(dist_matrix_day2[h2$order, h2$order]))
              
  
          }
          
          structure_sim_mat <- 
            rbind(structure_sim_mat,
                  data.frame(path1 = split_name(path_to_tune),
                             path2 = split_name(path_to_tune2),
                             knn1_path1 = params_mat[i,]$knn1,
                             knn2_path1 = params_mat[i,]$knn2,
                             knn1_path2 = params_mat[j,]$knn1,
                             knn2_path2 = params_mat[j,]$knn2,
                             cor = cor_level))
        }
      }
    }
  }
}

split_name <- function(path_to_tune) {
  return(sprintf("IC%s", str_split(path_to_tune, "IC")[[1]][2]))
}




build_sim_mat <- function(structure_sim_mat, params_mat) {

sim_mat <- matrix(rep(0, times=132**2), nrow = 132, ncol=132)

counter = 1
up <- unique(structure_sim_mat$path1)
uk <- unique(structure_sim_mat$knn1_path1)
for (p1 in up) {
  for (i in 1:nrow(params_mat)) {
    for (p2 in up) {
      for (j in 1:nrow(params_mat)) {
       
        c1 <- structure_sim_mat$path1 == p1
        c2 <- structure_sim_mat$path2 == p2
        c3 <- structure_sim_mat$knn1_path1 == params_mat[i,]$knn1
        c4 <- structure_sim_mat$knn2_path1 == params_mat[i,]$knn2
        c5 <- structure_sim_mat$knn1_path2 == params_mat[j,]$knn1
        c6 <- structure_sim_mat$knn2_path2 == params_mat[j,]$knn2
        
        r <- which(c1 & c2 & c3 & c4 & c5 & c6)
        if(len(r) != 1) {
          print("WHaT went wrong????")
        }
        
        p1_idx <- which(p1 == up)
        p2_idx <- which(p2 == up)
        
        kp1_idx <- which(structure_sim_mat[r,]$knn1_path1 == uk)
        kp2_idx <- which(structure_sim_mat[r,]$knn1_path2 == uk)
        
        row_idx <- (p1_idx - 1) * len(uk) + kp1_idx
        col_idx <- (p2_idx - 1) * len(uk) + kp2_idx
        print(col_idx)
        print(row_idx)
        
        sim_mat[row_idx, col_idx] <- structure_sim_mat[r,]$cor
        #print(counter)
        counter <- counter + 1
        
      }
    }
  }
}
}

closest <- apply(sim_mat, 1, function(r) {(which.max(r) - (which.max(r) %% 6)) / 6 + ifelse(which.max(r) %% 6 == 0, 0, 1)})
params <- apply(sim_mat, 1, function(r) {ifelse(which.max(r) %% 6 == 0, 6, which.max(r) %% 6)})
params_mat[]


show <- function(path_to_tune, i) {
  
  mat <- get_reduced_mat_full_day(path_to_tune, 
                                  knn1 = params_mat[1,]$knn1, 
                                  knn2 = params_mat[1,]$knn2,
                                  window_size = window_size)
  
  session_size <- mat_frames / window_size
  runs <- dim(mat)[1] / (session_size)
  
  behav_files <- list.files(sprintf("%s\\behavior\\",path_to_tune))
  behavior_mat_list <- lapply(behav_files,
                              function(tv_mat) {
                                return(readMat(sprintf("%s\\behavior\\%s",
                                                       path_to_tune,
                                                       tv_mat))$Stim.Master)
                              })
  
  # Good luck reading this lol
  behavior_indices <- 
    unlist(lapply(unlist(lapply(str_split(behav_files, ".TrialVar"), 
                                function(sp) {
                                  sp[[1]][1]
                                })),
                  function(p) {
                    tmp <- strsplit(p, "")[[1]]; 
                    return(as.numeric(paste(tmp[(length(tmp) - 1):(length(tmp))], collapse="")))
                  }))
  
  
  cp_list <- 
    get_fixed_color_palette(behav_mat_list = behavior_mat_list,
                            behavior_indices = behavior_indices,
                            runs = runs,
                            mat_frames = mat_frames,
                            session_size=nrow(mat))
  
  rmt <- get_reduced_mat_full_day(path_to_tune,
                           knn1=params_mat[i,]$knn1,
                           knn2=params_mat[i,]$knn2,
                           window_size = window_size)
  # plot3d(rmt, 
  #        col = cp_list$response_cp, 
  #        size=5)
  
  colnames(rmt) <- c("x","y","z")
  rmt$group <- 1:nrow(rmt)
  fig <- plotly::plot_ly(rmt, x = ~x, 
                 y = ~y, 
                 z = ~z, 
                 color=~group
                 colors = cp_list$response_cp)
  fig <- fig %>% plotly::add_trace(mode="line+markers", marker=list(size=2.5))
  fig <- fig %>% plotly::layout(scene = list(xaxis = list(title = 'Dim 1'),
                                     yaxis = list(title = 'Dim 2'),
                                     zaxis = list(title = 'Dim 3')))
  
  fig
  
}

show_friend <- function(i) {
  show(paths_all[closest[i]], params[i])
}

show_me <- function(i) {
    show(paths_all[((i - (i %% 6)) / 6 + ifelse(i %% 6 ==0,0,1))], 
         ifelse(i%%6 ==0, 6, i %%6))
}
