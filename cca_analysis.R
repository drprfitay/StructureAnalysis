library("zoo")
library("reshape2")

variables_to_smooth <- c("Cumulative reward", "Time", "Step")
smoothed_variables_names <- c("Smoothed cumulative", "Smoothed time", "Smoothed step")
smoothing_reps <- c(10, 33)
smoothed_variables_names <- c(sapply(smoothed_variables_names, function(vn) {paste(vn, as.character(smoothing_reps))}))

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
              "Y:\\livneh\\itayta\\v1_controls\\fov3\\day_140920\\",
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

dir.create("..//Desktop//cca_analysis")
dir.create("..//Desktop//cca_analysis//structure_vs_cca")
dir.create("..//Desktop//cca_analysis//figures_cca")
dir.create("..//Desktop//cca_analysis//modes_comp")


smooth_projection_to_max <- function(projection, size) {
  smoothed_proj <- rollmean(projection, size)
  smoothed_proj[(len(projection) - size + 1):len(projection)] <-
    smoothed_proj[len(smoothed_proj)]
  
  return(smoothed_proj)
}

get_neur_mat_and_stim_mat <- function(path, activity, control, window_size=30, mat_frames=57600) {
  neur_mat_func <- ifelse(control,
                          get_reduced_mat_full_day_control,
                          get_reduced_mat_full_day)
  
  mat <- t(neur_mat_func(path, 
                         normed=F,
                         activity_threshold = activity,
                         window_size=window_size,
                         just_original_mat = T))
  
  
  if (control) {
    stim_master_mat <- 
      get_color_palettes_control(path, just_mat = T, window_size = window_size)
  } else {
    stim_master_mat <- 
      get_color_palettes(path,
                         runs=nrow(mat) / (mat_frames / window_size),
                         just_mat = T, 
                         window_size = window_size)             
  }
  
  return(list(neur_mat=mat, stim_mat=stim_master_mat))
}

get_elbow <- function(vec_to_use_elbow_method) {

    elbow_vec <- c(scale(-vec_to_use_elbow_method))
  elbow_mat <- 
    rbind(elbow_vec, (1:len(elbow_vec) / len(elbow_vec)) * max(elbow_vec))
  
  elbow_dist <- sapply(1:ncol(elbow_mat),
                       function(pt_idx) {euc_dist(elbow_mat[,pt_idx], c(0,0))})
  
  elbow_point <- which.min(elbow_dist)
  
  return(elbow_point)
}

get_behavior_variables <- function(stim_master_mat, mat_size, window_size=30, seconds_after=2) {
  trials_ind <- sapply(stim_master_mat[,"Frames"], function(i) {get_binned_index(i, window_size)})

  non_na_ind <- !is.na(trials_ind)
  stim_master_mat <- stim_master_mat[non_na_ind,]
  
  cor_var_base <- rep(0, times=mat_size)
  rewards <- !is.nan(stim_master_mat[,"Reward"]) & stim_master_mat[,"TrialType"] %in% c(3)
  rewards_ind <- sapply(stim_master_mat[rewards,"Reward"],
                        function(i) {get_binned_index(i, window_size)})
  
  tmp <- rep(1:(len(rewards_ind) - 1), times=diff(rewards_ind))
  
  cumulative_reward <- cor_var_base
  cumulative_reward[1:len(tmp)] <- tmp
  cumulative_reward[(len(tmp) + 1):len(cumulative_reward)] <- max(tmp)
  
  
  rand_ind <-  trials_ind[which(stim_master_mat[,"TrialType"] == 4 | is.nan(stim_master_mat[,"Reward"]))]
  tmp_rand <- rep(1:(len(rand_ind) - 1), times=diff(rand_ind))
  
  random_cumulative <- cor_var_base
  random_cumulative[1:len(tmp_rand)] <- tmp_rand
  random_cumulative[(len(tmp_rand) + 1):len(random_cumulative)] <- max(tmp_rand)
  
  step_of_cumulative <- rep(0, times=len(cumulative_reward))
  step_of_cumulative[get_elbow(cumulative_reward):len(cumulative_reward)] <- max(cumulative_reward)

  
  cont_rewards_ind <- sapply(stim_master_mat[rewards,"Reward"],
                                         function(i) {get_binned_index(i:(i + seconds_after * window_size), window_size)})
  cont_rewards_ind <- unique(cont_rewards_ind)
  
  rewards <- cor_var_base
  rewards[cont_rewards_ind] <- 1
  
  time_var <- 1:len(cor_var_base)
 
  variables_mat <- cbind(cumulative_reward,
                         rewards,
                         time_var,
                         step_of_cumulative,
                         -cumulative_reward,
                         -random_cumulative) 
  
  for (ttype in c(3,4,5)) {
    
    presented_stimulus <- cor_var_base; 
    t_ind <- unique(c(sapply(trials_ind[which(stim_master_mat[,"TrialType"] == ttype)], function(i) {i:(i+seconds_after)})))
    presented_stimulus[t_ind] <- 1
    variables_mat <- cbind(variables_mat, presented_stimulus)
  }
  
  colnames(variables_mat) <- c(c("Cumulative reward",
                            "Reward",
                            "Time",
                            "Step",
                            "Inverse cumulative",
                            "Inverse random"),
                          sprintf("TrialType %d", c(3,4,5)))  
  
  return(variables_mat)
}

theta <- function(a,b, degrees=F) {rad = acos(sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ));
                                   return(ifelse(degrees, (rad/pi) * 180, rad))}

cancor_analysis <- function(neur_mat, stim_master_mat, window_size=30) {
  
  variables_mat <- get_behavior_variables(stim_master_mat,
                   mat_size=nrow(neur_mat),
                   window_size=window_size)
  
  cancor_list <- 
    lapply(1:ncol(variables_mat), 
           function(i) {return(cancor(neur_mat, 
                                      as.matrix(variables_mat[,i])))})
  
  modes_angel_mat <- matrix(rep(0, times=len(cancor_list) ** 2),
                            nrow=len(cancor_list))
  
  for (i in 1:len(cancor_list)) {
    for(j in 1:len(cancor_list)) {
      
      if (i == j) {
        modes_angel_mat[i,j] <- 0
        next
      }
    
      orth <- theta(cancor_list[[i]]$xcoef[,1],
                   cancor_list[[j]]$xcoef[,1], degrees=T)
      modes_angel_mat[i,j] <- orth
    }
  }
  
  correlations_vec <- sapply(cancor_list,
                             function(cc)
                             {cc$cor[1]})

  projected_mat <- 
  do.call(cbind,
          lapply(cancor_list,
                 function(cc) {neur_mat %*% cc$xcoef[,1]}))
  
  spearman_vec <- 
    sapply(1:ncol(variables_mat), 
           function(i) {abs(cor(projected_mat[,i], variables_mat[,i], method="spearman"))})
  
  
  scaled_dist_vec <- 
    sapply(1:ncol(projected_mat),
           function(proj_idx)
           { return(euc_dist(scale(projected_mat[,proj_idx]), scale(variables_mat[,proj_idx])))})
  
  stretched_dist_vec <- sapply(1:ncol(projected_mat),
                     function(proj_idx)
                       { return(euc_dist(stretch_to_edges(projected_mat[,proj_idx]), stretch_to_edges(variables_mat[,proj_idx])))})
  modes_mat <- 
    do.call(cbind,
            lapply(cancor_list,
                   function(cc) {cc$xcoef[,1]}))
  
  
  colnames(projected_mat) <- colnames(variables_mat)
  
  for (idx in 1:len(variables_to_smooth)) {
  
    var_to_smooth <- variables_to_smooth[idx]
    for (stime in smoothing_reps * 60 * (30 / window_size)) {
      smoothed_variable <- smooth_projection_to_max(projected_mat[,var_to_smooth], stime)
      projected_mat <- cbind(projected_mat, smoothed_variable)
      variables_mat <- cbind(variables_mat, variables_mat[,var_to_smooth])
      
      correlations_vec <- c(correlations_vec, 
                            cor(smoothed_variable, variables_mat[,var_to_smooth]))
      
      spearman_vec <- c(spearman_vec, 
                        cor(smoothed_variable, variables_mat[,var_to_smooth], method="spearman"))
      
      stretched_dist_vec <- c(stretched_dist_vec,
                              euc_dist(stretch_to_edges(smoothed_variable), 
                                       stretch_to_edges(variables_mat[,var_to_smooth])))
      
      scaled_dist_vec <- c(scaled_dist_vec,
                           euc_dist(scale(smoothed_variable), 
                                    scale(variables_mat[,var_to_smooth])))
    }
  }
  
  
  scaled_projection_mat <- apply(projected_mat, 2, scale)
  stretched_projection_mat <- apply(projected_mat, 2, stretch_to_edges)
  scaled_variables_mat <- apply(variables_mat, 2, scale)
  stretched_variables_mat <- apply(variables_mat, 2, stretch_to_edges)
  
  
  colnames(variables_mat)[(ncol(variables_mat) - len(smoothed_variables_names) + 1):ncol(variables_mat)] <- 
    smoothed_variables_names
  

  mode_names <- colnames(variables_mat)

  colnames(projected_mat) <- mode_names
  colnames(scaled_projection_mat) <- mode_names
  colnames(stretched_projection_mat) <- mode_names
  
  names(correlations_vec) <- mode_names
  names(spearman_vec) <- mode_names
  names(scaled_dist_vec) <- mode_names
  names(stretched_dist_vec) <- mode_names
  
  rownames(modes_angel_mat) <- mode_names[1:(len(mode_names) - len(smoothed_variables_names))]
  colnames(modes_angel_mat) <- mode_names[1:(len(mode_names) - len(smoothed_variables_names))]
  colnames(modes_mat) <- mode_names[1:(len(mode_names) - len(smoothed_variables_names))]
  
  return(list(angels_mat=modes_angel_mat,
              projection_mat=projected_mat,
              correlations_vec=correlations_vec,
              spearman_vec=spearman_vec,
              modes_mat=modes_mat,
              variables_mat=variables_mat,
              scaled_dist_vec=scaled_dist_vec,
              stretched_dist_vec=stretched_dist_vec,
              scaled_projection_mat=scaled_projection_mat,
              scaled_variables_mat=scaled_variables_mat,
              stretched_projection_mat=stretched_projection_mat,
              stretched_variables_mat=stretched_variables_mat))
}

stretch_to_edges <- function(vec) {
  return((vec - vec[1]) / (vec[len(vec)] - vec[1]))
}

get_heatmap_plot <- function(angles_mat) {
  plot_df <- expand.grid(M1=rownames(angles_mat),
                         M2=rownames(angles_mat))
  plot_df$Angle = c(angles_mat)
  
  hmp <- 
    ggplot(plot_df,aes(x=M2, y=M1, fill=Angle))+ 
    geom_tile() + 
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    xlab("Mode 1") +
    ylab("Mode 2") + 
    scale_fill_distiller(palette = "Spectral")
  return(hmp)
}


get_mode_plot <- function(mode_obj, i) {
  plot_df  <- data.frame(`Scaled variable`=mode_obj$cancor_res$scaled_variables_mat[,i],
                         `Stretched variable`=mode_obj$cancor_res$stretched_variables_mat[,i],
                         Time=1:nrow(mode_obj$cancor_res$projection_mat),
                         `Scaled mode`=mode_obj$cancor_res$scaled_projection_mat[,i],
                         `Stretched mode`=mode_obj$cancor_res$stretched_projection_mat[,i])
  
  if (colnames(mode_obj$cancor_res$projection_mat)[i] %in% 
      c("Cumulative reward", "Time", "Step", "Inverse cumulative", "Inverse random",
        smoothed_variables_names)){
    mode_plot <- ggplot(plot_df, aes(x=Time, y=`Stretched.mode`)) +
      geom_line(color="royalblue4", size=0.5) +
      geom_line(aes(y=`Stretched.variable`), col="red", linetype="longdash", alpha=0.8,
                size=2) + 
      theme_light() +     
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(color="black"),
            panel.border = element_blank(),
            panel.background = element_blank()) + 
      xlab("Time (frames)") + 
      ylab(sprintf("Mode - %s (AU)", colnames(mode_obj$variables_mat)[i])) +
      ggtitle(sprintf("Projection on %s, dist = %f, pearson = %f",
                      colnames(mode_obj$variables_mat)[i],
                      mode_obj$cancor_res$dist_vec[i],
                      mode_obj$cancor_res$correlations_vec[i]))
    
    
  } else {
    mode_plot <- ggplot(plot_df, aes(x=Time, y=`Scaled.variable`)) +
      geom_line(col="red",  size=1) + 
      geom_line(aes(y=`Scaled.mode`), alpha=0.8, color="royalblue4", size=1) +
      theme_light() +     
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(color="black"),
            panel.border = element_blank(),
            panel.background = element_blank()) + 
      xlab("Time (frames)") + 
      ylab(sprintf("Mode - %s (AU)", colnames(mode_obj$variables_mat)[i])) +
      ggtitle(sprintf("Projection on %s, dist = %f, pearson = %f",
                      colnames(mode_obj$variables_mat)[i],
                      mode_obj$cancor_res$dist_vec[i],
                      mode_obj$cancor_res$correlations_vec[i]))
  }
  
  return(mode_plot)
}


analyse_modes <- function(window_size=30, mat_frames=57600) {

  
  
  modes_analysis_list <- 
  lapply(1:len(paths_all),
         function(idx) {
           
           work_path <- paths_all[idx]
           control <- (categories[idx] %in% c("Control"))
           activity_thres <- ifelse(categories[idx] %in% c("Hunger", "Control"),
                                    0.25,
                                    0.2)
           
           print(sprintf("%d %f %s", 
                         idx,
                         activity_thres, 
                         categories[idx]))
           
          res <- get_neur_mat_and_stim_mat(work_path, activity_thres, control)
           

           cancor_res=cancor_analysis(res$neur_mat , res$stim_mat, window_size=window_size)

          return(list(cancor_res=cancor_res))
         })
  
      comparsion_methods_vec <- c("Pearson"="correlations_vec",
                                  "Spearman"="spearman_vec",
                                  "Stretched distance"="stretched_dist_vec",
                                  "Scaled distance"="scaled_dist_vec")
      
      comparsion_methods_axes <- c("Pearson" = "Correlation (pearson r)",
                                   "Spearman" = "Correlation (spearman r)",
                                   "Stretched distance" = "Euclidean distance - stretched (AU)",
                                   "Scaled distance" = "Euclidean distance - Scaled (AU)")
      
      for (compare_method_name in names(comparsion_methods_vec)) {
        compare_method <- comparsion_methods_vec[compare_method_name]
        
        compare_df <- 
          do.call(rbind, lapply(modes_analysis_list,
                                function(m) {return(m$cancor_res[[compare_method]])}))
      
        compare_df <- as.data.frame(compare_df)
        compare_df <- cbind(compare_df, categories)
        compare_df <- melt(compare_df)[,-4]
        
        compare_modes_graph <-  ggplot(compare_df, aes(x=variable, y=value)) + 
                            geom_bar(stat="summary", 
                                     aes(group=categories, fill=categories), 
                                     color="NA", 
                                     position = "dodge") + 
                            geom_point(aes(x=variable, y=value, group=factor(categories)), 
                                       position=position_dodge(1), 
                                       color="gray60") + 
                            scale_fill_manual(name="Task",
                                              breaks=c("Control", "Thirst", "Hunger"), 
                                              values=c("goldenrod2", "royalblue4", "darksalmon" ),
                                              labels=c("V1/POR", "Insula - thirst", "Insula - hunger")) + 
                            theme_light() +     
                            theme(panel.grid.major = element_blank(), 
                                  panel.grid.minor = element_blank(),
                                  axis.line = element_line(colour = "black"),
                                  panel.border = element_blank(),
                                  panel.background = element_blank()) + 
                            ylab(comparsion_methods_axes[compare_method_name]) +
                            xlab("Modes") + 
                            ggtitle(sprintf("Modes comparsion - %s", compare_method_name))
        
        pdf(sprintf("../Desktop/modes_comp/Modes_%s.pdf", tolower(str_replace(compare_method_name, " ", "_"))),
            height=8,
            width=30)
        
        plot(compare_modes_graph)
        dev.off()
        
      }

      
      all_mode_plots <- 
      lapply(modes_analysis_list,
             function(mode_obj)
               {
                hmp_plot <- get_heatmap_plot(mode_obj$cancor_res$angels_mat)
                mode_plot_list <- 
                lapply(1:ncol(mode_obj$cancor_res$projection_mat),
                       
                       function(i) {
                         return(get_mode_plot(mode_obj, i))
                       })
                
                mode_plot_list$hmp_plot <- hmp_plot
                mode_plot_list$nrow = 3
                final <- 
                do.call(grid.arrange,mode_plot_list)
               
                return(final)
             })
      
      for (i in 1:len(paths_all)) {
        path <- paths_all[i]
        cat = path_categories[i]
        
        if (cat %in% c("V1", "POR")) {
          pattern = "fov"
        } else {
          pattern = "IC"
        }
      
        split_r <- str_split(str_split(path, pattern)[[1]][[2]], "\\\\")
      
        pdf(sprintf("../Desktop//figures_cca//%s_%s%s_%s.pdf",
                    cat,
                    pattern,
                    split_r[[1]][1],
                    split_r[[1]][2]),
            height=30,
            width=40)
        
        plot(all_mode_plots[[i]])
        dev.off()
    }
}



mode_contribution_analysis <- function(neur_mat, mode_weight_vector, mode_variable) {
  
  
  individual_correlation <- sapply(1:ncol(neur_mat),
                                   function(neur_idx) {
                                      projection <- as.matrix(neur_mat[,neur_idx]) %*% 
                                                          mode_weight_vector[neur_idx]
                                      
                                      return(cor(projection, mode_variable))
                                   })
  
  
  cumulative_shape_dist <- 
  sapply(1:ncol(neur_mat), 
         function(len_of_subsample) {
            idx_to_use <- order(individual_correlation, decreasing = T)[1:len_of_subsample]
            projection <- as.matrix(neur_mat[,idx_to_use]) %*% mode_weight_vector[idx_to_use]
            return(euc_dist(scale(smooth_projection_to_max(projection, 600)), scale(mode_variable)))
         })
  
  
  cumulative_cor <- 
    sapply(1:ncol(neur_mat), 
           function(len_of_subsample) {
             idx_to_use <- order(individual_correlation, decreasing = T)[1:len_of_subsample]
             projection <- as.matrix(neur_mat[,idx_to_use]) %*% mode_weight_vector[idx_to_use]
             return(cor(projection, mode_variable))
           })

  
  # Find elbow point
  elbow_point <- get_elbow(cumulative_cor)

  
  interesting_cells <- order(individual_correlation, decreasing = T)[1:elbow_point]
  projection <- as.matrix(neur_mat[,interesting_cells]) %*% mode_weight_vector[interesting_cells]
  
}

cool_pairs <- function(reduced_mat, modes_n, modes_labs) {
  axes_names <- c("x", "y", "z")
  axes_labs <- c("Component 1", "Component 2", "Component 3")
  colnames(reduced_mat) <- axes_names
  reduced_mat$ind <- 1:nrow(reduced_mat)
  
  red_modes_df <- reduced_mat[(nrow(reduced_mat) - modes_n + 1):nrow(reduced_mat),]
  red_to_use <- reduced_mat[1:(nrow(reduced_mat) - modes_n),]
  
  plt_theme <- 
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "gray30", size=1),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = "NA")
  
  filler_theme <-       
    theme_light() +     
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          #axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = "NA",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  plot_list <- list()
  
  for (i in 1:3) {
    for(j in 1:3) {
      
      if (i == j) {
        df <- data.frame(x=0,
                         y=0,
                         text=axes_labs[i])
        p <- ggplot(df, aes(x,y)) + filler_theme + 
          geom_text(aes(label=text)) +
          xlab("") + 
          ylab("")
      } else {
        p <- ggplot(red_to_use) +
          geom_point(aes(x=!!sym(axes_names[i]),
                         y=!!sym(axes_names[j]),
                         col=ind)) +
          plt_theme +
          scale_colour_viridis() +
          xlab("") + 
          ylab("")
        
        
        for (mode_i in 1:nrow(red_modes_df)) {
          p_mode_df <- data.frame(x=c(0, red_modes_df[mode_i, "x"]),
                                  y=c(0, red_modes_df[mode_i, "y"]),
                                  z=c(0, red_modes_df[mode_i, "z"]),
                                  label=modes_to_use[mode_i])
          
          p <- 
            p + 
            geom_line(data=p_mode_df,
                      aes(x=!!sym(axes_names[i]),
                          y=!!sym(axes_names[j])),
                      color="black",
                      size=1.5)
          
          if (i > j) {
            
            p_label_df <- data.frame(x=c(red_modes_df[mode_i, "x"]),
                                     y=c(red_modes_df[mode_i, "y"]),
                                     z=c(red_modes_df[mode_i, "z"]),
                                     label=modes_labs[mode_i])
            
            p <- 
              p + 
              geom_label(data=p_label_df,
                         aes(x=!!sym(axes_names[i]),
                             y=!!sym(axes_names[j]),
                             label=label))
            
          }
        }
      }
      
      
      plot_list <- append(plot_list, list(p))
    }
  }
  
  plot_list$nrow=3
  
  gr <- do.call(grid.arrange, plot_list)
  aligned <- align_plots(gr, align = "v")
  
  return(aligned[[1]])
}


cca_and_structure_analysis <- function(path, activity_thres, control) {
  
  ext=""
  for (idx in 11:17) {
    work_path <- paths_all[idx]
    control <- (categories[idx] %in% c("Control"))
    activity_thres <- ifelse(categories[idx] %in% c("Hunger", "Control"),
                             0.25,
                             0.2)
    
    modes_to_use <- c("Cumulative reward", "Reward",  "TrialType 3")
    res <- get_neur_mat_and_stim_mat(work_path, activity_thres, control)
    cca_res <-  cancor_analysis(res$neur_mat, res$stim_mat)
    
    modes_mat_scaled <- cca_res$modes_mat[,modes_to_use] * 1000
    max_val <- colSums(modes_mat_scaled)[which.max(abs(colSums(modes_mat_scaled)))]
    modes_mat_scaled_2 <- apply(modes_mat_scaled, 2, function(c) {c * abs((max_val / sum(c)))})
    
    joint_matrices <- rbind(res$neur_mat, t(modes_mat_scaled_2))
    joint_reduced_mat <- reduce_dim(joint_matrices, "lem2", ndim=3, knn1=0.35, knn2=0.125)
    
    # Guys, this is a chunk of code for the purposes of research.
    # I'm not developing something that's going to NASDAQ here, stop judging me about the naming convention
    aligned_plot <- 
     cool_pairs(joint_reduced_mat, 
                len(modes_to_use), 
              modes_labs=c("Cumulative rewrad", "Reward",
                           "Reward cue"))
    
    cat = path_categories[idx]
    
    if (cat %in% c("V1", "POR")) {
      pattern = "fov"
    } else {
      pattern = "IC"
    }
    
    split_r <- str_split(str_split(work_path, pattern)[[1]][[2]], "\\\\")
    
    
    pdf(sprintf("..//Desktop//cca_analysis//structure_vs_cca//%s_%s%s_%s%s.pdf",
                cat,
                pattern,
                split_r[[1]][1],
                split_r[[1]][2],
                ext),
        height=10,
        width=10)
    
    plot(aligned_plot)
    dev.off()
  }
}






mode_variability_analysis <- function(path, activity_thres, control) {
  
  ext=""
  for (idx in 1:len(paths_all)) {

    work_path <- paths_all[idx]
    control <- (categories[idx] %in% c("Control"))
    activity_thres <- ifelse(categories[idx] %in% c("Hunger", "Control"),
                             0.25,
                             0.2)
    
    modes_to_use <- c("Cumulative reward", "Reward",  "TrialType 3")
    res <- get_neur_mat_and_stim_mat(work_path, activity_thres, control, window_size = 15)
    cca_res <-  cancor_analysis(res$neur_mat, res$stim_mat, window_size = 1)
    
    rss_all <- apply(res$neur_mat, 2, function(c) {var(c)})
    rss_projections <- apply(cca_res$projection_mat, 2, function(c) {var(c)})

    
    cat = path_categories[idx]
    
    if (cat %in% c("V1", "POR")) {
      pattern = "fov"
    } else {
      pattern = "IC"
    }
    
    split_r <- str_split(str_split(work_path, pattern)[[1]][[2]], "\\\\")
    
    
    pdf(sprintf("..//Desktop//cca_analysis//structure_vs_cca//%s_%s%s_%s%s.pdf",
                cat,
                pattern,
                split_r[[1]][1],
                split_r[[1]][2],
                ext),
        height=10,
        width=10)
    
    plot(aligned_plot)
    dev.off()
  }
}
