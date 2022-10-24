library("zoo")
library("RColorBrewer")
library("pheatmap")
library("tidyverse")

hmp_col <- colorRampPalette(rev(brewer.pal(n = 11,  name = "Spectral")))

slow_time_scale_mat <- function(reduced, nclusters=30, ntimes=50, roll_window=1, use_hcluster=F, scaled=F) {
  
  if (use_hcluster) {
    hc <- hclust(dist(reduced), method='ward.D')
    km <- list()
    km$cluster <- cutree(hc, nclusters)
  } else {
    km <- kmeans(reduced, nclusters, iter.max=500)
  }
  
  
  
  cluster_mat <- matrix(rep(0, times=(nrow(reduced) * nclusters)), nrow=nclusters)
  

  cluster_mat <-   
    lapply(sort(unique(km$cluster)),  
           function(clust_idx) {
             cluster_mat[clust_idx, which(km$cluster == clust_idx)] <- 1;
             return(cluster_mat[clust_idx,])
             })
  
  cluster_mat <- do.call(rbind, cluster_mat)
  chunk_size <- floor(nrow(reduced) / ntimes) 
  chunk_indices <- lapply(1:ntimes, 
                          function(i) {((i - 1) * chunk_size + 1):(i * chunk_size)})
  
  slow_time_scale_mat  <- t(apply(cluster_mat,
                                1,
                                function(cluster_row)
                                  {rate_vec <- unlist(lapply(chunk_indices, function(ind) {mean(cluster_row[ind])}));
                                   return(rollmean(rate_vec, roll_window))}))
  
  if (scaled) {
    return(as.matrix(scale(slow_time_scale_mat)[,]))
  }
  
  return(slow_time_scale_mat)
}


transition_prob_mat <- function(reduced, nclusters=30) {
  km <- kmeans(reduced, nclusters, iter.max=500)
  cluster_time_series <-  c(km$cluster, -1) # -1 for last time frame
  tabled <-  table(km$cluster)
 
  prob_mat <- matrix(rep(0, times=nclusters ** 2), nrow=nclusters) 
  
  clusters = sort(unique(km$cluster))
  

  for (cluster_idx in clusters) {
    tabulated_transitions = table(cluster_time_series[which(cluster_time_series == cluster_idx) + 1])
    
    # Last time frame has no transition
    if ("-1" %in% names(tabulated_transitions)) {
      print(" found last one!")
      tabulated_transitions = tabulated_transitions[-which(names(tabulated_transitions) == "-1")]
      
      # One no transition
      tabled[cluster_idx] <- tabled[cluster_idx] - 1
    }
    
    # In case there's a state we're not transiting to, set probability as 0
    no_transitions <- which(!as.character(1:nclusters) %in% names(tabulated_transitions))
    
    if (len(no_transitions) > 0) {
      print(sprintf("No transitions! length of (%d)", len(no_transitions)))
      new_tabulated_transitions <- c(tabulated_transitions, rep(0, times=len(no_transitions)))
      names(new_tabulated_transitions) <- c(names(tabulated_transitions), as.character(1:nclusters)[no_transitions])
      
      #Resort 
      tabulated_transitions <- new_tabulated_transitions[order(as.numeric(names(new_tabulated_transitions)))]
    }
    
    assert(sum(tabulated_transitions) == tabled[cluster_idx])
    prob_mat[cluster_idx, ] <- tabulated_transitions / tabled[cluster_idx]
  }
  
  return(prob_mat)
}



get_slow_time_cor_dist <- function(ins_dirs, por_dirs, v1_dirs, metric=euc_dist, nclusters=50, ntimes=100, roll_window=20) {
ins_trans <- lapply(ins_dirs, function(v)  {sorted_mat(v, knn1=0.325, knn2=0.1, nclusters=nclusters, ntimes=ntimes, trans=F, roll_window=roll_window)})
por_trans <- lapply(por_dirs, function(v)  {sorted_mat(v, knn1=0.325, knn2=0.1, nclusters=nclusters, ntimes=ntimes, trans=F, roll_window=roll_window, control = T)})
v1_trans <- lapply(v1_dirs, function(v)  {sorted_mat(v, knn1=0.325, knn2=0.1, nclusters=nclusters, ntimes=ntimes, trans=F, roll_window=roll_window, control=T)})


  por_ins <- c() 
  for (i in 1:len(por_trans)) {
    for (j in 1:len(ins_trans)) {
      por_ins <- c(por_ins, metric(c(por_trans[[i]]), 
                                   c(ins_trans[[j]])))
    }
  } 
  
  v1_ins <- c()
  for (i in 1:len(v1_trans)) {
    for (j in 1:len(ins_trans)) {
      v1_ins <- c(v1_ins, metric(c(v1_trans[[i]]), c(ins_trans[[j]])))
    }
  }
  
  
  por_v1 <- c()
  
  for (i in 1:len(v1_trans)) {
    for (j in 1:len(por_trans)) {
      por_v1 <- c(por_v1, metric(c(v1_trans[[i]]), c(por_trans[[j]])))
    }
  }
  
  
  ins_ins <- c()
  for (i in 1:len(ins_trans)) {
    for (j in 1:len(ins_trans)) {
      if (i != j) {
        ins_ins <- c(ins_ins, metric(c(ins_trans[[i]]), c(ins_trans[[j]])))
      }
    }
  }
  
  v1_v1 <- c()
  for (i in 1:len(v1_trans)) {
    for (j in 1:len((v1_trans))) {
      if (i != j) {
        v1_v1 <- c(v1_v1, metric(c(v1_trans[[i]]), c(v1_trans[[j]])))
      }
    }
  }
  
  por_por <- c()
  for (i in 1:len(por_trans)) {
    for (j in 1:len(por_trans)) {
      if (i != j) {
        por_por <- c(por_por, metric(c(por_trans[[i]]), c(por_trans[[j]])))
      }
    }
  }
  
  par(mfrow=c(3,2))
  
  plot(density(ins_ins), lwd=2, main=sprintf("Insula vs Insula: %f", mean(ins_ins)), xlim=c(0,1))
  abline(v=mean(ins_ins), col="red", lty=2)
  
  plot(density(v1_ins), lwd=2, main=sprintf("V1 vs Insula: %f", mean(v1_ins)), xlim=c(0,1))
  abline(v=mean(v1_ins), col="red", lty=2)

  plot(density(por_ins), lwd=2, main=sprintf("POR vs Insula: %f", mean(por_ins)), xlim=c(0,1))
  abline(v=mean(por_ins), col="red", lty=2)
    
  plot(density(por_v1), lwd=2, main=sprintf("POR vs V1: %f", mean(por_v1)), xlim=c(0,1))
  abline(v=mean(por_v1), col="red", lty=2)
  
  plot(density(por_por), lwd=2, main=sprintf("POR vs POR: %f", mean(por_por)), xlim=c(0,1))
  abline(v=mean(por_por), col="red", lty=2)
  
  plot(density(v1_v1), lwd=2, main=sprintf("V1 vs V1: %f", mean(v1_v1)), xlim=c(0,1))
  abline(v=mean(v1_v1), col="red", lty=2)
}



get_prob_trans_cor_dist_no_diag <- 
  function(ins_dirs, por_dirs, v1_dirs, metric=euc_dist, nclusters=50) {
  ins_trans <- lapply(ins_dirs, function(v)  {sorted_mat(v, knn1=0.325, knn2=0.1, nclusters=nclusters, trans=T, control=F)})
  por_trans <- lapply(por_dirs, function(v)  {sorted_mat(v, knn1=0.325, knn2=0.1, nclusters=nclusters, trans=T, control=T)})
  v1_trans <- lapply(v1_dirs, function(v)  {sorted_mat(v, knn1=0.325, knn2=0.1, nclusters=nclusters, trans=T, control=T)})
  
  
  por_ins <- c() 
  for (i in 1:len(por_trans)) {
    for (j in 1:len(ins_trans)) {
      por_ins <- c(por_ins, metric(c(por_trans[[i]])[-(((1:nclusters)- 1) * nclusters + 1:nclusters)], 
                                   c(ins_trans[[j]])[-(((1:nclusters)- 1) * nclusters + 1:nclusters)]))
    }
  } 
  
  v1_ins <- c()
  for (i in 1:len(v1_trans)) {
    for (j in 1:len(ins_trans)) {
      v1_ins <- c(v1_ins, metric(c(v1_trans[[i]])[-(((1:nclusters)- 1) * nclusters + 1:nclusters)], 
                                 c(ins_trans[[j]])[-(((1:nclusters)- 1) * nclusters + 1:nclusters)]))
    }
  }
  
  
  por_v1 <- c()
  
  for (i in 1:len(v1_trans)) {
    for (j in 1:len(por_trans)) {
      por_v1 <- c(por_v1, metric(c(v1_trans[[i]])[-(((1:nclusters)- 1) * nclusters + 1:nclusters)], 
                                 c(por_trans[[j]])[-(((1:nclusters)- 1) * nclusters + 1:nclusters)]))
    }
  }
  
  
  ins_ins <- c()
  for (i in 1:len(ins_trans)) {
    for (j in 1:len(ins_trans)) {
      if (i != j) {
        ins_ins <- c(ins_ins, metric(c(ins_trans[[i]])[-(((1:nclusters)- 1) * nclusters + 1:nclusters)], 
                                     c(ins_trans[[j]])[-(((1:nclusters)- 1) * nclusters + 1:nclusters)]))
      }
    }
  }
  
  v1_v1 <- c()
  for (i in 1:len(v1_trans)) {
    for (j in 1:len((v1_trans))) {
      if (i != j) {
        v1_v1 <- c(v1_v1, metric(c(v1_trans[[i]])[-(((1:nclusters)- 1) * nclusters + 1:nclusters)],
                                 c(v1_trans[[j]])[-(((1:nclusters)- 1) * nclusters + 1:nclusters)]))
      }
    }
  }
  
  por_por <- c()
  for (i in 1:len(por_trans)) {
    for (j in 1:len(por_trans)) {
      if (i != j) {
        por_por <- c(por_por, metric(c(por_trans[[i]])[-(((1:nclusters)- 1) * nclusters + 1:nclusters)], 
                                     c(por_trans[[j]])[-(((1:nclusters)- 1) * nclusters + 1:nclusters)]))
      }
    }
  }
  
  par(mfrow=c(3,2))
  
  plot(density(ins_ins), lwd=2, main=sprintf("Insula vs Insula: %f", mean(ins_ins)), xlim=c(0,1))
  abline(v=mean(ins_ins), col="red", lty=2)
  
  plot(density(v1_ins), lwd=2, main=sprintf("V1 vs Insula: %f", mean(v1_ins)), xlim=c(0,1))
  abline(v=mean(v1_ins), col="red", lty=2)
  
  plot(density(por_ins), lwd=2, main=sprintf("POR vs Insula: %f", mean(por_ins)), xlim=c(0,1))
  abline(v=mean(por_ins), col="red", lty=2)
  
  plot(density(por_v1), lwd=2, main=sprintf("POR vs V1: %f", mean(por_v1)), xlim=c(0,1))
  abline(v=mean(por_v1), col="red", lty=2)
  
  plot(density(por_por), lwd=2, main=sprintf("POR vs POR: %f", mean(por_por)), xlim=c(0,1))
  abline(v=mean(por_por), col="red", lty=2)
  
  plot(density(v1_v1), lwd=2, main=sprintf("V1 vs V1: %f", mean(v1_v1)), xlim=c(0,1))
  abline(v=mean(v1_v1), col="red", lty=2)
}

sorted_mat <- function(v, knn1, knn2, nclusters, ntimes=100, trans=T, control=F, roll_window=1) {
  if (control) {
    red <- get_reduced_mat_full_day_control(v, knn1=knn1, knn2=knn2)
  } else {
    red <- get_reduced_mat_full_day(v, knn1=knn1, knn2=knn2)
  }
  
  if (trans) {
    
    mat <- transition_prob_mat(red, nclusters=nclusters)
  } else {
    mat <- slow_time_scale_mat(red, nclusters=nclusters, ntimes=ntimes, roll_window=roll_window)
  }
  hc <- hclust(dist(mat))
  
  if (trans) {
    return(mat[hc$order, hc$order])
  } else {
    return(mat[hc$order,])
  }
}

get_prob_trans_cor_dist <- function(ins_dirs, por_dirs, v1_dirs, metric=euc_dist, nclusters=50) {
  ins_trans <- lapply(ins_dirs, function(v)  {sorted_mat(v, knn1=0.325, knn2=0.1, nclusters=nclusters, trans=T, control=F)})
  por_trans <- lapply(por_dirs, function(v)  {sorted_mat(v, knn1=0.325, knn2=0.1, nclusters=nclusters, trans=T, control=T)})
  v1_trans <- lapply(v1_dirs, function(v)  {sorted_mat(v, knn1=0.325, knn2=0.1, nclusters=nclusters, trans=T, control=T)})
  
  
  por_ins <- c() 
  for (i in 1:len(por_trans)) {
    for (j in 1:len(ins_trans)) {
      por_ins <- c(por_ins, metric(c(por_trans[[i]]), c(ins_trans[[j]])))
    }
  } 
  
  v1_ins <- c()
  for (i in 1:len(v1_trans)) {
    for (j in 1:len(ins_trans)) {
      v1_ins <- c(v1_ins, metric(c(v1_trans[[i]]), c(ins_trans[[j]])))
    }
  }
  
  
  por_v1 <- c()
  
  for (i in 1:len(v1_trans)) {
    for (j in 1:len(por_trans)) {
      por_v1 <- c(por_v1, metric(c(v1_trans[[i]]), c(por_trans[[j]])))
    }
  }
  
  
  ins_ins <- c()
  for (i in 1:len(ins_trans)) {
    for (j in 1:len(ins_trans)) {
      if (i != j) {
        ins_ins <- c(ins_ins, metric(c(ins_trans[[i]]), c(ins_trans[[j]])))
      }
    }
  }
  
  v1_v1 <- c()
  for (i in 1:len(v1_trans)) {
    for (j in 1:len((v1_trans))) {
      if (i != j) {
        v1_v1 <- c(v1_v1, metric(c(v1_trans[[i]]), c(v1_trans[[j]])))
      }
    }
  }
  
  por_por <- c()
  for (i in 1:len(por_trans)) {
    for (j in 1:len(por_trans)) {
      if (i != j) {
        por_por <- c(por_por, metric(c(por_trans[[i]]), c(por_trans[[j]])))
      }
    }
  }
  
  #tiff("C:\\Users\\itayta.WISMAIN\\Desktop\\struct_plots\\15c_compare.tiff", res=400, units="in", heigh=10, width=7)
  par(mfrow=c(3,2))
  
  plot(density(ins_ins), lwd=2, main=sprintf("Insula vs Insula: %f", mean(ins_ins)), xlim=c(0,1))
  abline(v=mean(ins_ins), col="red", lty=2)
  
  plot(density(v1_ins), lwd=2, main=sprintf("V1 vs Insula: %f", mean(v1_ins)), xlim=c(0,1))
  abline(v=mean(v1_ins), col="red", lty=2)
  
  plot(density(por_ins), lwd=2, main=sprintf("POR vs Insula: %f", mean(por_ins)), xlim=c(0,1))
  abline(v=mean(por_ins), col="red", lty=2)
  
  plot(density(por_v1), lwd=2, main=sprintf("POR vs V1: %f", mean(por_v1)), xlim=c(0,1))
  abline(v=mean(por_v1), col="red", lty=2)
  
  plot(density(por_por), lwd=2, main=sprintf("POR vs POR: %f", mean(por_por)), xlim=c(0,1))
  abline(v=mean(por_por), col="red", lty=2)
  
  plot(density(v1_v1), lwd=2, main=sprintf("V1 vs V1: %f", mean(v1_v1)), xlim=c(0,1))
  abline(v=mean(v1_v1), col="red", lty=2)
  
  #dev.off()
}

fix_missing <- function(tabled, nclusters) {
  missing <- which(!1:nclusters %in% names(tabled))
  tabled_fixed <- c(tabled, rep(0, times=len(missing)))
  names(tabled_fixed) <- c(names(tabled), missing)
  return(tabled_fixed[order(as.numeric(names(tabled_fixed)))])
}


plot_differentially_visited_regions <- function(path, knn1=0.325, knn2=0.1, nclusters=50, q, activity, control=F, plot=T, use_hclust=T, threshold=0) {
  
  mat_aq_func <- ifelse(control, get_reduced_mat_full_day_control, get_reduced_mat_full_day)
  mat <- mat_aq_func(path, knn1=knn1, knn2=knn2, activity_threshold = activity)  
  
  if (use_hclust) {
    hc <- hclust(dist(mat))
    labels <- cutree(hc, nclusters)
  } else {
    km <- kmeans(mat, centers = nclusters, iter.max=500)
    labels <- km$cluster
  }

  ind_first_portion <- 1:(len(labels) * q) 
  ind_second_portion <- (len(labels) * (1-q) + 1):(len(labels))
  tabled_first_portion <- fix_missing(table(labels[ind_first_portion]), nclusters)
  tabled_second_portion <- fix_missing(table(labels[ind_second_portion]), nclusters)
  
  delta_portion <- log(tabled_first_portion + 10^-30, 2) - log(tabled_second_portion + 10^-30,2)
  if (plot) {
    barplot(sort(delta_portion), col="royalblue4")
  }
  
  return(delta_portion[delta_portion > threshold])
}

# 
# 
# plot_differential_visitation_distributions() {
#   params_df <- list(nclusters=c(20,50,75, 100, 150, 200), q=seq(0.05,0.5, by=0.05))
#   params_df <- as.matrix(cross_df(params_df))
#   
#   
#   for (idx in 1:nrow(params_df)) {
#    
#     nclusters = params_df[idx, 1]
#     q = params_df[idx,2] 
#     
#     v1_diff <- c(unlist(sapply(1:4,
#                                function(i) {plot_differentially_visited_regions(v1_paths[i], q=q, activity=0.25, nclusters=nclusters)})))
#     ins_thirst_diff <- c(unlist(sapply(1:4,
#                                        function(i) {plot_differentially_visited_regions(insula_thirst_paths[i], q=q, activity=0.2, nclusters=nclusters)})))
#     ins_hunger_diff <- c(unlist(sapply(1:6,
#                                        function(i) {plot_differentially_visited_regions(insula_paths[i], q=q, activity=0.25, nclusters=nclusters)})))
# 
#     png(sprintf("%s\\fig_2_nc%d_q%.3f.png",
#               "Y:\\livneh\\itayta\\figures\\figure_2_differential_time_logfc",
#               nclusters, q),
#         units="in",
#         height=18,
#         width=8, 
#         res=300)
#     
#     par(mfrow=c(3,1))
# 
#     for (thr in c(0.01,0.02,0.03)) {
#       
#       if (thr > max(c(v1_diff[v1_diff > thr],
#               ins_thirst_diff[ins_thirst_diff > thr],
#               ins_hunger_diff[ins_hunger_diff > thr]))) {
#         next
#       }
#         
# 
#       boxplot(v1_diff[v1_diff > thr],
#               ins_thirst_diff[ins_thirst_diff > thr],
#               ins_hunger_diff[ins_hunger_diff > thr],
#               col=c("lightpink1", "royalblue1", "green"),
#               names=c("V1", "Thirst", "Hunger"),
#               main=sprintf("NC: %d, Q: %.3f, Threshold: %.3f",
#                             nclusters, q, thr))
#     }
#     
#     dev.off()
#   }
# }
# 




plot_differential_visitation_distributions_via_sd <- function() {
  params_df <- list(nclusters=c(20,30,50,75,100, 200), ntimes=seq(10,100, by=10))
  params_df <- as.matrix(cross_df(params_df))
  sd_threshold = 1

  insula_paths <- c("Y:\\livneh\\itayta\\data\\IC19\\day_150911\\",
                    "Y:\\livneh\\itayta\\data\\IC17\\day_150615\\",
                    "Y:\\livneh\\itayta\\data\\IC13\\day_150406\\",
                    "Y:\\livneh\\itayta\\data\\IC13\\day_150407\\",
                    "Y:\\livneh\\itayta\\data\\IC32\\day_161214\\",
                    "Y:\\livneh\\itayta\\data\\IC42\\day_161117\\")
  
  
  v1_paths <- c("Y:\\livneh\\itayta\\v1_controls\\fov1\\day_140524\\",
                "Y:\\livneh\\itayta\\v1_controls\\fov3\\day_140920\\",
                "Y:\\livneh\\itayta\\v1_controls\\fov3\\day_140921\\",
                "Y:\\livneh\\itayta\\v1_controls\\fov5\\day_150723\\")
  
  por_paths <- c("Y:\\livneh\\itayta\\por_controls\\fov1\\day_141023\\",
                 "Y:\\livneh\\itayta\\por_controls\\fov2\\day_140805\\",
                 "Y:\\livneh\\itayta\\por_controls\\fov3\\day_150411\\",
                 "Y:\\livneh\\itayta\\por_controls\\fov4\\day_150112\\")
  
  all <- get_all_paths("Y:\\livneh\\itayta\\data")
  
  IC44_days <- c("170518", "170523", "170519", "170524")
  IC44_paths <- all[grep("IC44", all)]
  insula_thirst_paths <- IC44_paths[which(rowSums(sapply(IC44_days, function(d) {as.numeric(grepl(d, IC44_paths))})) > 0)]
  
  for (i in 1:nrow(params_df)) { 
    
    nc <- params_df[i,1]
    ntimes <- params_df[i,2]
    
    for (wind in c(2,5,10,15)){
      
      if (wind >= ntimes) {
        next
      }
      
    print(sprintf("%d %d %d", nc, ntimes, wind))
    v1_sd <- lapply(1:4,
                    function(i) {slow_time_scale_mat(get_reduced_mat_full_day_control(v1_paths[i], knn1=0.325, knn2=0.1, activity=0.25), nc, ntimes, wind, T, T)})
    
    por_sd <- lapply(1:4,
                    function(i) {slow_time_scale_mat(get_reduced_mat_full_day_control(por_paths[i], knn1=0.325, knn2=0.1, activity=0.25), nc, ntimes, wind, T, T)})    
    
    ins_thirst_sd <-lapply(1:4,
                               function(i) {slow_time_scale_mat(get_reduced_mat_full_day(insula_thirst_paths[i], knn1=0.325, knn2=0.1, activity=0.2), nc, ntimes, wind, T, T)})
    ins_hunger_sd <-lapply(1:6,
                               function(i) {slow_time_scale_mat(get_reduced_mat_full_day(insula_paths[i], knn1=0.325, knn2=0.1, activity=0.25), nc, ntimes, wind, T, T)})
      
    
    v1_sd <- unlist(lapply(v1_sd, function(mt) {apply(mt, 1, sd)}))
    ins_thirst_sd <- unlist(lapply(ins_thirst_sd, function(mt) {apply(mt, 1, sd)}))
    ins_hunger_sd <- unlist(lapply(ins_hunger_sd, function(mt) {apply(mt, 1, sd)}))
    por_sd <- unlist(lapply(por_sd, function(mt) {apply(mt, 1, sd)}))
    
    v1_sd_filt <- v1_sd[v1_sd >= sd_threshold]
    ins_thirst_sd_filt <- ins_thirst_sd[ins_thirst_sd >= sd_threshold]
    ins_hunger_sd_filt <- ins_hunger_sd[ins_hunger_sd >= sd_threshold]
    por_sd_filt <- por_sd[por_sd >= sd_threshold]
    
    if (len(c(v1_sd_filt, ins_thirst_sd_filt, ins_hunger_sd_filt, por_sd_filt)) == 0) {
      next
    }
    
    png(sprintf("%s\\fig_2_occupancy_sd_nc%d_ntimes%d_wind%d.png",
                "Y:\\livneh\\itayta\\figures\\figures_2_slow_time_diff",
                nc, ntimes,  wind),
        units="in",
        height=5,
        width=5, 
        res=300)
    
        plot(y=c(v1_sd_filt,
              ins_thirst_sd_filt,
              ins_hunger_sd_filt,
              por_sd_filt),
             
             x=c(rep(1, times=len(v1_sd_filt)),
                 rep(2, times=len(ins_thirst_sd_filt)),
                 rep(3, times=len(ins_hunger_sd_filt)),
                 rep(4, times=len(por_sd_filt))),
             
             bg=c(rep("darksalmon", times=len(v1_sd_filt)),
                  rep("royalblue4", times=len(ins_thirst_sd_filt)),
                  rep("darkseagreen2", times=len(ins_hunger_sd_filt)),
                  rep("goldenrod2", times=len(por_sd_filt))),
             xlim=c(0,5),
             xlab="Group",
             ylab="Manifold occupancy SD",
             xaxt="n",
              main=sprintf("NC: %d, Ntimes: %d, Roll: %d", nc, ntimes, wind),
             cex=2,
             pch=21)
        
        axis(1, at=c(1,2,3,4), labels=c("V1", "Thirst", "Hunger", "POR"))
        
        dev.off()
        
    }
  }
}




plot_differential_visitation_distributions_via_cor <- function() {
  params_df <- list(nclusters=c(5, 10, 20,30,50,75,100, 200, 300), ntimes=seq(10,100, by=10))
  params_df <- as.matrix(cross_df(params_df))
  cor_threshold = 0.8
  
  insula_paths <- c("Y:\\livneh\\itayta\\data\\IC19\\day_150911\\",
                    "Y:\\livneh\\itayta\\data\\IC17\\day_150615\\",
                    "Y:\\livneh\\itayta\\data\\IC13\\day_150406\\",
                    "Y:\\livneh\\itayta\\data\\IC13\\day_150407\\",
                    "Y:\\livneh\\itayta\\data\\IC32\\day_161214\\",
                    "Y:\\livneh\\itayta\\data\\IC42\\day_161117\\")
  
  
  v1_paths <- c("Y:\\livneh\\itayta\\v1_controls\\fov1\\day_140524\\",
                "Y:\\livneh\\itayta\\v1_controls\\fov3\\day_140920\\",
                "Y:\\livneh\\itayta\\v1_controls\\fov3\\day_140921\\",
                "Y:\\livneh\\itayta\\v1_controls\\fov5\\day_150723\\")
  
  por_paths <- c("Y:\\livneh\\itayta\\por_controls\\fov1\\day_141023\\",
                 "Y:\\livneh\\itayta\\por_controls\\fov2\\day_140805\\",
                 "Y:\\livneh\\itayta\\por_controls\\fov3\\day_150411\\",
                 "Y:\\livneh\\itayta\\por_controls\\fov4\\day_150112\\")
  
  all <- get_all_paths("Y:\\livneh\\itayta\\data")
  
  IC44_days <- c("170518", "170523", "170519", "170524")
  IC44_paths <- all[grep("IC44", all)]
  insula_thirst_paths <- IC44_paths[which(rowSums(sapply(IC44_days, function(d) {as.numeric(grepl(d, IC44_paths))})) > 0)]
  
  for (i in 1:nrow(params_df)) { 
    
    nc <- params_df[i,1]
    ntimes <- params_df[i,2]
    
    for (wind in c(1, 2, 5, 10, 20)){
      
      if (wind >= ntimes) {
        next
      }
      
      print(sprintf("%d %d %d", nc, ntimes, wind))
      v1_mats <- lapply(1:4,
                      function(i) {slow_time_scale_mat(get_reduced_mat_full_day_control(v1_paths[i], knn1=0.325, knn2=0.1, activity=0.25), nc, ntimes, wind, T, F)})
      
      por_mats <- lapply(1:4,
                       function(i) {slow_time_scale_mat(get_reduced_mat_full_day_control(por_paths[i], knn1=0.325, knn2=0.1, activity=0.25), nc, ntimes, wind, T, F)})    
      
      ins_thirst_mats <-lapply(1:4,
                             function(i) {slow_time_scale_mat(get_reduced_mat_full_day(insula_thirst_paths[i], knn1=0.325, knn2=0.1, activity=0.2), nc, ntimes, wind, T, F)})
      ins_hunger_mats <-lapply(1:6,
                             function(i) {slow_time_scale_mat(get_reduced_mat_full_day(insula_paths[i], knn1=0.325, knn2=0.1, activity=0.25), nc, ntimes, wind, T, F)})


      
      for (i in 1:len(v1_mats)) {
        pdf(sprintf("../Desktop//slow_heatmaps//v1//heatmap_V1_%d_%d_ntimes%d_wind%d.pdf", i, nc, ntimes,  wind),
            height=8,
            width=8)
        pheatmap(v1_mats[[i]], cluster_cols=F, cluster_rows=F, 
                 main=sprintf("NC: %d, Ntimes: %d, Roll: %d", nc, ntimes, wind), color=hmp_col(200))
        dev.off()
      }
      
      for (i in 1:len(ins_thirst_mats)) {
        pdf(sprintf("../Desktop//slow_heatmaps//thirst//heatmap_thirst_%d_%d_ntimes%d_wind%d.pdf", i, nc, ntimes,  wind),
            height=8,
            width=8)
        pheatmap(ins_thirst_mats[[i]], cluster_cols=F, cluster_rows=F, 
                 main=sprintf("NC: %d, Ntimes: %d, Roll: %d", nc, ntimes, wind), color=hmp_col(200))
        dev.off()
      }

      for (i in 1:len(ins_hunger_mats)) {
        pdf(sprintf("../Desktop//slow_heatmaps//hunger//heatmap_thirst_%d_%d_ntimes%d_wind%d.pdf", i, nc, ntimes,  wind),
            height=8,
            width=8)
        pheatmap(ins_hunger_mats[[i]], cluster_cols=F, cluster_rows=F, 
                 main=sprintf("NC: %d, Ntimes: %d, Roll: %d", nc, ntimes, wind), color=hmp_col(200))
        dev.off()
      }      
      
      for (i in 1:len(por_mats)) {
        pdf(sprintf("../Desktop//slow_heatmaps//por//heatmap_POR_%d_%d_ntimes%d_wind%d.pdf", i, nc, ntimes,  wind),
            height=8,
            width=8)
        pheatmap(por_mats[[i]], cluster_cols=F, cluster_rows=F, 
                 main=sprintf("NC: %d, Ntimes: %d, Roll: %d", nc, ntimes, wind), color=hmp_col(200))
        dev.off()
      }
      
      for (cor_t in c("spearman", "pearson", "kendall",  "sd")) {
      print(cor_t)
      if (cor_t %in% c("spearman", "pearson", "kendall")) {
        v1_cor <- (lapply(v1_mats, function(mt) {apply(mt, 1, function(occ) {abs(cor(occ, 1:((ntimes - wind) + 1), method=cor_t))})}))
        ins_thirst_cor <- (lapply(ins_thirst_mats, function(mt) {apply(mt, 1, function(occ) {abs(cor(occ, 1:((ntimes - wind) + 1),  method=cor_t))})}))
        ins_hunger_cor <- (lapply(ins_hunger_mats, function(mt) {apply(mt, 1, function(occ) {abs(cor(occ, 1:((ntimes - wind) + 1),  method=cor_t))})}))
        por_cor <- (lapply(por_mats, function(mt) {apply(mt, 1, function(occ) {abs(cor(occ, 1:((ntimes - wind) + 1),  method=cor_t))})}))
      } else {
        v1_cor <- lapply(v1_mats, function(mt) {apply(mt, 1, function(occ) {sd(occ)})})
        ins_thirst_cor <- (lapply(ins_thirst_mats, function(mt) {apply(mt, 1, function(occ) {sd(occ)})}))
        ins_hunger_cor <- (lapply(ins_hunger_mats, function(mt) {apply(mt, 1, function(occ) {sd(occ)})}))
        por_cor <- (lapply(por_mats, function(mt) {apply(mt, 1, function(occ) {sd(occ)})}))
      }
      
      # v1_cor <- unlist(lapply(v1_cor, function(cor_vec) {sum(cor_vec > cor_threshold) / len(cor_vec)}))
      # ins_thirst_cor <- unlist(lapply(ins_thirst_cor, function(cor_vec) {sum(cor_vec > cor_threshold) / len(cor_vec)}))
      # ins_hunger_cor <- unlist(lapply(ins_hunger_cor, function(cor_vec) {sum(cor_vec > cor_threshold) / len(cor_vec)}))
      # por_cor <- unlist(lapply(por_cor, function(cor_vec) {sum(cor_vec > cor_threshold) / len(cor_vec)}))


      
      png(sprintf("Y:\\livneh\\itayta\\figures\\figures_2_slow_time_diff_%s\\fig_2_occupancy_cor_nc%d_ntimes%d_wind%d.png",
                  cor_t,
                  nc, ntimes,  wind),
          units="in",
          height=5,
          width=12, 
          res=300)
      par(mfrow=c(1,2))
      plot(y=c(unlist(v1_cor),
               unlist(ins_thirst_cor),
               unlist(ins_hunger_cor),
               unlist(por_cor)),
           
           x=c(rep(1:len(v1_cor), each=nc),
               rep((len(v1_cor) + 1):(len(v1_cor) + len(ins_thirst_cor)), each=nc),
               rep((len(v1_cor) + len(ins_thirst_cor) + 1):(len(v1_cor) + len(ins_thirst_cor) + len(ins_hunger_cor)), each=nc),
               rep((len(v1_cor) + len(ins_thirst_cor) + len(ins_hunger_cor) + 1):(len(v1_cor) + len(ins_thirst_cor) + len(ins_hunger_cor) + len(por_cor)), each=nc)),
           
           bg=c(rep("darksalmon", times=nc * len(v1_cor)),
                rep("royalblue4", times=nc * len(ins_thirst_cor)),
                rep("darkseagreen2", times=nc * len(ins_hunger_cor)),
                rep("goldenrod2", times=nc * len(por_cor))),
           xlim=c(0,20),
           xlab="Group",
           ylab=ifelse(cor_t == "sd", "SD",
                       sprintf("Cluster occupancy and time correlation (%s)", cor_t)),
           xaxt="n",
           main=sprintf("NC: %d, Ntimes: %d, Roll: %d", nc, ntimes, wind),
           cex=2,
           pch=21)
      
      axis(1, at=c(1:18), labels=c(rep("V1", times=len(v1_cor)),
                                   rep("Thirst", times=len(ins_thirst_cor)),
                                   rep("Hunger", times=len(ins_hunger_cor)),
                                   rep("POR", times=len(por_cor))))
      
      boxplot(append(append(append(v1_cor, ins_thirst_cor), ins_hunger_cor), por_cor),
              col=c(rep("darksalmon", times=len(v1_cor)),
                        rep("royalblue4", times=len(ins_thirst_cor)),
                        rep("darkseagreen2", times=len(ins_hunger_cor)),
                        rep("goldenrod2", times=len(por_cor))),
              names=c(rep("V1", times=len(v1_cor)),
                      rep("Thirst", times=len(ins_thirst_cor)),
                      rep("Hunger", times=len(ins_hunger_cor)),
                      rep("POR", times=len(por_cor))))
              
      
      dev.off()
      }
      
    }
  }
  
}
