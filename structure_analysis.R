library("reticulate")
library("testit")
library("dplyr")
library("plyr")
library("OneR")
library("stringr")
library("R.matlab")
library("ggplot2")
library('latex2exp')
library("gridExtra") 
library("dimRed")
library("scatterplot3d")
library("Rtsne")
library("viridis")
library("umap")

list.of.packages <- c("rgl","ggplot2","knitr","rglwidget")

len <- length
#
#You should be able to simply reuse the following lines of code as is
#
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#
if(length(new.packages)) install.packages(new.packages)
#
# By now we have installed the requisite packages. Time to load them .
#
lapply(list.of.packages,function(x){library(x,character.only=TRUE)})


matrices_path <- "C:\\Users\\97254\\Documents\\yoav_livneh_lab\\temp\\"
sample_path <- "C:\\Users\\97254\\Documents\\yoav_livneh_lab\\sample_data\\IC44_170518_run1\\dFF_mat2.mat"
sample_path2 <- "C:\\Users\\97254\\Documents\\yoav_livneh_lab\\sample_data\\IC44_170518_run2\\dFF_mat2.mat"
sample_path3 <- "C:\\Users\\97254\\Documents\\yoav_livneh_lab\\sample_data\\IC44_170518_run3\\dFF_mat2.mat"

stim_trace_path <- "C:\\Users\\97254\\Documents\\yoav_livneh_lab\\sample_data\\IC44_170518_run1\\Stim_Master.mat"
stim_trace_path2 <- "C:\\Users\\97254\\Documents\\yoav_livneh_lab\\sample_data\\IC44_170518_run2\\Stim_Master.mat"
stim_trace_path3 <- "C:\\Users\\97254\\Documents\\yoav_livneh_lab\\sample_data\\IC44_170518_run3\\Stim_Master.mat"

lick_path <- "C:\\Users\\97254\\Documents\\yoav_livneh_lab\\sample_data\\IC44_170518_run1\\licking.mat"
bl_lick_path <- "C:\\Users\\97254\\Documents\\yoav_livneh_lab\\sample_data\\IC44_170518_run1\\bl.mat"
lick_path2 <- "C:\\Users\\97254\\Documents\\yoav_livneh_lab\\sample_data\\IC44_170518_run2\\LickingFrames.mat"
lick_path3 <- "C:\\Users\\97254\\Documents\\yoav_livneh_lab\\sample_data\\IC44_170518_run3\\LickingFrames.mat"

result_path <- "C:\\Users\\97254\\Documents\\yoav_livneh_lab\\results\\Individual mice\\"

trial_type_name <- c("Pavlov", "Reward", "Aversive", "Neutral", "Blank")
names(trial_type_name) <- c("1","3","4","5","6")

color_scheme <- list(c("deeppink", "chocolate1", "cadetblue1", "chartreuse3"),
                     c("darkseagreen1", "chocolate1", "cadetblue1", "darkseagreen3"),
                     c("dodgerblue1", "chocolate1", "cadetblue1", "coral3"),
                     c("darkorchid2", "chocolate1", "cadetblue1", "aquamarine3"),
                     c("yellow1", "chocolate1", "cadetblue1", "yellow1"),
                     c("black", "chocolate1", "cadetblue1", "black"))

color_scheme <- do.call(rbind, color_scheme)
rownames(color_scheme) <- c("3","4","1","5","6", "2")
colnames(color_scheme) <- c("type", "lick", "nolick", "iti")

# mat <- readMat(sample_path)
# mat2 <- readMat(sample_path2)
# mat3 <- readMat(sample_path3)
# stim <- readMat(stim_trace_path)
# stim2 <- readMat(stim_trace_path2)
# stim3 <- readMat(stim_trace_path3)



euc_dist <- function(a, b) sqrt(sum((a - b)^2))

mat_stability <- function(mt)
{return(c(0,sapply(2:ncol(mt), function(i) {euc_dist(mt[,i], mt[,i-1])})))}

get_binned_index <- function(i, window_size) {
  rst <- floor((i)/window_size) + ifelse(i %% window_size == 0, 0,1)
}


# Return cells in which activity in above 3sigma 
get_outliers <- function(spike_train, plot=F, threshold=3) {
  avg_activity <- rowMeans(spike_train, na.rm = T) / dt
  outliers <- c(which(is.nan(avg_activity)))
  
  avg_activity <- avg_activity[!is.nan(avg_activity)]
  sigma <- (avg_activity - mean(avg_activity)) / sd(avg_activity)
  
  outliers <- c(outliers, which(abs(sigma) > threshold))
  

  
  if (plot) {
    ind <- 1:len(avg_activity)
    ind <- which(!ind %in% outliers)
    
    hist(avg_activity[ind], breaks=50, col="lightblue", xlab="Average df/f",
        main="Activity histogram")
  }
  
  return(outliers)
}

get_final_spike_train <- function(spike_train, outliers) {
  return(spike_train[!1:nrow(spike_train) %in% outliers,])
}

time_bin_average_vec <- function(vec, window_size) {

  indices <- 1:(len(vec) / window_size)
  
  final_vec <- c()
  for (ind in indices) {
    subset_ind <- ((ind - 1 ) * window_size + 1):(ind * window_size)
    final_vec <-  c(final_vec,
                       mean(vec[subset_ind]))
  }
  
  return(final_vec)
}

time_bin_average <- function(spike_train, window_size) {
  indices <- 1:(ncol(spike_train) / window_size)
  
  final_spike_train <- c()
  for (ind in indices) {
    subset_ind <- ((ind - 1 ) * window_size + 1):(ind * window_size)
    final_spike_train <- cbind(final_spike_train,
                               rowMeans(spike_train[,subset_ind]))
  }
  
  return(final_spike_train)
}




plot_3d_pc <- function(df, extra_colors=F, extra_ind=NA, extra_col=NA, plotly=F, lty="p", mcol_pal=NA, ylim=NA, xlim=NA, zlim=NA, xlab=NA, ylab=NA, zlab=NA, screen=NA) {
  library(plotly)
  
  #col_pal <- plasma(nrow(df))
  col_pal <- c(viridis(nrow(df)))
  #col_pal <- c(col_pal, inferno(nrow(df) / 2))

  
  if (extra_colors && len(extra_col) == len(extra_ind)) {
    print("Changing col")
    col_pal[extra_ind] <- extra_col
    print(col_pal[extra_ind])
  }
  
  df$group <- 1:nrow(df)
  
  if (sum(!is.na(mcol_pal)) > 0) {
    col_pal = mcol_pal
  }
  
  if (plotly) {
  fig <- plot_ly(df, x = ~x, 
                     y = ~y, 
                     z = ~z, 
                     color = ~group,
                     colors = col_pal)
  fig <- fig %>% add_trace(mode="line+markers", marker=list(size=1.5))
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'Dim 1'),
                                     yaxis = list(title = 'Dim 2'),
                                     zaxis = list(title = 'Dim 3')))
  
  fig
  
  } else {
    col_pal[1] <- "gray50"
    plot3d(df$x, df$y, df$z,col=col_pal, 
           pch=19, size=5, type=lty, screen=list(x=60, y=-30, z=-100),
           xlim=xlim, ylim=ylim, zlim=zlim,
           xlab=xlab,zlab=zlab,ylab=ylab)
  }
}

get_stim_indices <- function(stim_mat, allowed_stim=c(1:6), window_size, response=F) {
  stim_indices <- stim_mat[,1]
  c_to_use <- ifelse(response, 8, 7)
  stim_types <- stim_mat[,c_to_use]
  
  
  ind_to_use <- which(stim_types %in% allowed_stim)
  stim_indices <- stim_indices[ind_to_use]
  stim_types <- stim_types[ind_to_use]
  fixed_indices <- sapply(stim_indices, function(i) {get_binned_index(i, window_size)})
  
  return(list(ind=fixed_indices, type=stim_types))
}


reduce_dim <- function(data, method, ndim=3, knn1=0.275, knn2=0.075) {
  print(sprintf("performing dim reduction on method! %s, data dimensonality (%d x %d)", 
         method,
         dim(data)[1],
         dim(data)[2]
         ))
  if (method == "tsne") {
    red <- Rtsne(data, dims = ndim, perplexity=100, check_duplicates = F)
    red_df <- do.call(cbind, lapply(1:ndim, function(i) {red$Y[,i]}))
    red_df <- as.data.frame(red_df)
    
  } else if (method == "lem") {
    red <- dimRed::embed(data, "LaplacianEigenmaps", ndim=ndim)
    red_df <- do.call(cbind, lapply(1:ndim, function(i) {red@data@data[,i]}))
    red_df <- as.data.frame(red_df)
  } else if (method == "isomap") {
    red <- dimRed::embed(data, "Isomap", ndim=ndim)
    red_df <- do.call(cbind, lapply(1:ndim, function(i) {red@data@data[,i]}))
    red_df <- as.data.frame(red_df)
  } else if (method == "lem2") {

    print(sprintf("Using KNN1=%d, KNN2=%d",
          floor(nrow(data) * knn1),
          floor(nrow(data) * knn2)))
    
    red <- dimRed::embed(data, 
                 "LaplacianEigenmaps", 
                 ndim=15,
                 knn=floor(nrow(data) * knn1),
                 t=Inf)

    red <- do.call(cbind, lapply(1:15, function(i) {red@data@data[,i]}))
    red <- as.data.frame(red)    
    
    red <- dimRed::embed(red, 
                 "LaplacianEigenmaps", 
                 ndim=ndim,
                 knn=floor(nrow(data) * knn2),
                 t=Inf)
    
    red_df <- do.call(cbind, lapply(1:ndim, function(i) {red@data@data[,i]}))
    red_df <- as.data.frame(red_df)
  } else if (method == "umap") {
    red <- umap(data, 
                n_components=ndim,
                metric="cosine",
                #min_dist=0.05,
                n_neighbors=50)
    red_df <- do.call(cbind, lapply(1:ndim, function(i) {red$layout[,i]}))
    red_df <- as.data.frame(red_df)
  }  else if (method == "umap_denoised") {
    pc <- prcomp(data)
    
    if (ncol(pc$x) < 100) {
      red <- umap(pc$x, 
                  n_components=ndim,
                  metric="cosine",
                  #min_dist=0.05,
                  n_neighbors=50)  
    } else {
    red <- umap(pc$x[,1:100], 
                n_components=ndim,
                metric="cosine",
                #min_dist=0.05,
                n_neighbors=50)
    }
    red_df <- do.call(cbind, lapply(1:ndim, function(i) {red$layout[,i]}))
    red_df <- as.data.frame(red_df)
    
  } else if (method == "mds") {
    red <- cmdscale(dist(data),eig=TRUE, k=ndim)
    red_df <- do.call(cbind, lapply(1:ndim, function(i) {red$points[,i]}))
    red_df <- as.data.frame(red_df)
    
  } else {
   return(data) 
  }
  
  colnames(red_df) <- c("x","y","z")[1:ndim]
  return(red_df)
}


split <- function(binned_st, factor, alg) {
  
  return(lapply(lapply(1:factor,function(i) {((i - 1) * (ncol(binned_st) / factor) + 1):(i * (ncol(binned_st) / factor))}),
       function(ind) {
         print("Reducing!")
         return(reduce_dim(t(binned_st)[ind,], alg))
       }))
  
}


calc_lick_rate <- function(lick_vec, mat, window_size) {
  binary_lick <- rep(0, ncol(mat))
  binary_lick[lick_vec] <- 1
  avg_rate <- time_bin_average_vec(binary_lick, window_size)
  return(avg_rate)
}


create_event_matrix <- function(ca_mat, mad_threshold=4, remove_previous=F) {
  
  emat <- matrix(rep(0, times=ncol(ca_mat) * nrow(ca_mat)), 
                 nrow = nrow(ca_mat),
                 ncol = ncol(ca_mat))
  
  for(i in 1:nrow(ca_mat)) {
    
    diffr <- diff(ca_mat[i,])
    diffr[which(diffr < 0)] <- -1
    diffr[which(diffr > 0)] <- 1
    diffr[which(diffr == 0)] <- 0
    
    local_maxima <- which(diff(diffr) < 0) + 1
    #e_ind <- abs(ca_mat[i,] - median(ca_mat[i,])) / 
    #                median(abs(ca_mat[i,] - median(ca_mat[i,]))) > mad_threshold
    local_maxima_mad <- ca_mat[i,local_maxima] / 
                        median(abs(ca_mat[i,] - median(ca_mat[i,])))    
    
    # Only frames with calcium transients higher than previous frame
    #e_ind <- e_ind & (diff(ca_mat[1,]) < 0)
    print(sum(local_maxima_mad > mad_threshold))
    emat[i,local_maxima] = as.numeric(local_maxima_mad > mad_threshold)
    
    
    
    if (verbose) {
      print(i)
    }
  }
  
  return(emat)
}

preprocess_neur_activity <- function(ca_mat, window_size, scale=F) {
   
  if(scale) {
    print("SCALING$$$$$$$$$$$$$$$$$$$$$$$$$$$!")
    ca_mat <- smooth_ca_trace(ca_mat)
  }
  
  emat <- create_event_matrix(ca_mat, mad_threshold = 4)
  
  emat_clean <- get_final_spike_train(emat, get_outliers(emat))
  ca_mat_clean <- get_final_spike_train(ca_mat, get_outliers(ca_mat))
  
  emat_binned <- time_bin_average(emat_clean, window_size)
  ca_mat_binned <- time_bin_average(ca_mat_clean, window_size)
  
  return(list(event=emat_binned, ca=ca_mat_binned))
}



plot_rand_trace <- function(ca_mat, emat, trace_size) {
  n  <- floor(runif(1, min=1, max=nrow(emat)))
  
  
  s <- floor(runif(1, min=1, max=(ncol(emat) - trace_size)))
  
  plot(ca_mat[n,s:(s+trace_size)], col="black", type="l")
  lines(emat[n,s:(s+trace_size)], col="red", type="l")
}

median_filter_ca_mat <- function(ca_mat, hz, filter_size) {
  return(apply(ca_mat, 1, function(trace) { trace - runmed(trace, hz * filter_size)}))
}

reduce_mat <- function(mat_list, alg, ndim) {
  return(lapply(mat_list, function(m) {reduce_dim(t(m), alg,ndim)}))
}

calc_var <- function(pc) {
  egs <- pc$sdev ** 2
  return(egs / sum(egs))
}


smooth_ca_trace <- function(ca_mat) {
  return(t(apply(ca_mat, 1, FUN=scale)))
}


generate_response_colors <- function(nframes, stim_master, window_size) {
  
  col_pal <- rep("gray60", times=nframes)
  stim <- stim_master[,7]
  stim_indices <- stim_master[,1]
  stim_indices_fixed <- sapply(stim_indices, 
                               function(i) {get_binned_index(i, window_size)})
  stim_responses <- stim_master[,8]
  lick_nolick <- c("cadetblue1", "chocolate1")
  names(lick_nolick) <- c(2,1)
  # # Trial structure, 2s window response
  # col_pal[stim_indices_fixed] <- color_scheme[,"type"][stim]
  # col_pal[stim_indices_fixed + 1] <- color_scheme[,"type"][stim]
  # col_pal[stim_indices_fixed + 2] <- color_scheme[,"type"][stim]
  # col_pal[stim_indices_fixed + 3] <- color_scheme[,"type"][stim]
  
  # 2s response window
  col_pal[stim_indices_fixed + 4] <- lick_nolick[stim_responses %% 2 + 1]
  col_pal[stim_indices_fixed + 5] <- lick_nolick[stim_responses %% 2 + 1]
  
  # 6s ITI  
  col_pal[stim_indices_fixed + 8] <- color_scheme[,"iti"][stim]
  col_pal[stim_indices_fixed + 9] <- color_scheme[,"iti"][stim] # 1s
  col_pal[stim_indices_fixed + 10] <- color_scheme[,"iti"][stim]
  col_pal[stim_indices_fixed + 11] <- color_scheme[,"iti"][stim] # 2s
  col_pal[stim_indices_fixed + 12] <- color_scheme[,"iti"][stim]
  col_pal[stim_indices_fixed + 13] <- color_scheme[,"iti"][stim] # 3s
  col_pal[stim_indices_fixed + 14] <- color_scheme[,"iti"][stim]
  col_pal[stim_indices_fixed + 15] <- color_scheme[,"iti"][stim] # 4s
  
  return(col_pal)
}

get_ind_from_stim <- function(stim_master, b, a, r, window_size) {

  ind <- stim_master[which(stim_master[,7] == r),1]
  
  ind <- sapply(ind, function(i) {get_binned_index(i, window_size)})
  
  return(as.vector(sapply(ind, function(idx) {(idx - b):(idx + a)})))
}


get_time_binned_mat <- function(mat, window_size, clean=T) {
  if (clean) {
    return(time_bin_average(get_final_spike_train(mat, get_outliers(mat)),window_size))
  } 
  
  return(time_bin_average(mat,window_size))
}

get_binary_licking <- function(lick_mat, nframes) {
  bl <- rep(0, times <- nframes)
  bl[lick_mat[1,]] <- 1
  return(bl)
}


tuning_correlation <- function(mat, 
                               stim_master, 
                               window_size=30,
                               trials_to_use  =c(3),
                               seconds_before=0,
                               seconds_after=8,
                               outout_path="",
                               trial_str="NA",
                               fname_ext="") {
  
  f_stim_ind <- c() # Sorted indices by trial
  n_trials <- c() # Vector containing how many trials for each trial type
  lick_ind <- c() # Response vector
  
  for (st in stim_to_use) {
    # Ind of trial
    stim_ind <- stim_master[which(stim_master[,7] %in% c(st)),1]
    f_stim_ind <- c(f_stim_ind, stim_ind)
    
    # Num of each trial
    n_trials <- c(n_trials, len(stim_ind))
    
    # Trial response
    l_ind <- stim_master[which(stim_master[,7] %in% c(st)),8]
    l_ind <- as.character(as.numeric(l_ind %% 2 == 1))
    l_ind[which(l_ind == "1")] <- "No lick"
    l_ind[which(l_ind == "0")] <- "Lick"
    lick_ind <- c(lick_ind, l_ind)
  }
  

  smpl <- data.frame(Trials = rep(trial_type_name[as.character(stim_to_use)], 
                                  n_trials),
                     Response=lick_ind)
  if (sort_trials) {
    f_stim_ind <- sort(f_stim_ind)
    
    lick_ind <- stim_master[which(stim_master[,7] %in% stim_to_use),8]
    lick_ind <- as.character(as.numeric(lick_ind %% 2 == 1))
    lick_ind[which(lick_ind == "1")] <- "No lick"
    lick_ind[which(lick_ind == "0")] <- "Lick"
    
    smpl <- data.frame(Trials = trial_type_name[as.character(stim_master[which(stim_master[,7] %in% stim_to_use),7])],
                       Response = lick_ind)
  }
  
  f_out_path <- sprintf("%s\\Individual Neurons", out_path)
  stim_mat_lst <- lapply(f_stim_ind,
                       function(si) {
                         get_time_binned_mat(mat[,(si - seconds_before * 30):(si + seconds_after * 30)], 
                                             8 * 30,
                                             clean=F)
                       })
  
  
  mean_cor <- c()
  median_cor <- c()

  sym = F
  
  for (neur in 1:nrow(mat)) {
    
    
    
    print(sprintf("Analyzing neuron %d", neur))
    
    neur_tuning_per_trial_list <- 
      lapply(1:len(stim_mat_lst),
             function(i) {
               return(stim_mat_lst[[i]][neur,])
             })
    
    neur_tuning_mat <- do.call(cbind, neur_tuning_per_trial_list)
    neur_tuning_cor_mat <- cor(neur_tuning_mat, neur_tuning_mat)
    dir.create(sprintf("%s\\%d", f_out_path, neur))

    
    colnames(neur_tuning_cor_mat) <- 1:nrow(neur_tuning_cor_mat)
    rownames(neur_tuning_cor_mat) <- 1:nrow(neur_tuning_cor_mat)
    rownames(smpl) <-  1:nrow(neur_tuning_cor_mat)
    
    cl  <- color_scheme[,"type"][as.character(stim_to_use)]
    names(cl) <- trial_type_name[as.character(stim_to_use)]
    cl_rs <- c("lightblue", "lightpink1")
    names(cl_rs) <- c("Lick", "No lick")
    
    clf=list(Trials=cl, Response=cl_rs)
    
    # heatmap(neur_tunng_cor_mat,
    #         col=inferno(200),
    #         Colv=l,
    #         Rowv=NA,
    #         symm=T,
    #         xlab = "",
    #         ylab= "",
    #         main=sprintf("Trials correlation matrix for neuron %d (%s)",
    #                      neur,
    #                      trial_str))
    
    png(filename=sprintf("%s\\%d\\Clustered_corr_matrix_%s.png", 
                         f_out_path, 
                         neur,
                         paste(trial_type_name[as.character(stim_to_use)], collapse=".")),
        width=600,
        height=600)
    
    
    pheatmap(neur_tuning_cor_mat, 
             col=turbo(200), 
             annotation_col=smpl, 
             annotation_row=smpl,
             cluster_cols=T, 
             cluster_rows=T,
             show_rownames=F,
             show_colnames=F,
             annotation_names_col = F,
             annotation_names_row = F,
             annotation_colors=clf,
             main=sprintf("Trials correlation matrix neuron %d", neur))    
    
    dev.off()
    
    
    png(filename=sprintf("%s\\%d\\Unclustered_corr_matrix_%s.png", 
                         f_out_path, 
                         neur,
                         paste(trial_type_name[as.character(stim_to_use)], collapse=".")),
        width=600,
        height=600)
    
    pheatmap(neur_tuning_cor_mat, 
             col=turbo(200), 
             annotation_col=smpl,
             annotation_row=smpl,
             cluster_cols=F, 
             cluster_rows=F,
             show_rownames=F,
             show_colnames=F,
             annotation_names_col = F,
             annotation_names_row=F,
             annotation_colors =clf,
             main=sprintf("Trials correlation matrix neuron %d", neur))    
    
    dev.off()
    
    cmbnation <- combn(1:51,2)
    
    all_cors <- c()
    for (mat_c in 1:ncol(cmbnation)) {
      cmb <- cmbnation[,mat_c]
      
      
      all_cors <- c(all_cors,
                    cor(neur_tuning_mat[,cmb[1]], neur_tuning_mat[,cmb[2]]))
      
    }
    
    png(filename=sprintf("%s\\%d\\Pairwise_cor_hist_%s.png", 
                         f_out_path, 
                         neur,
                         paste(trial_type_name[as.character(stim_to_use)], collapse=".")),
        width= 800)
    
    hist(all_cors, 
         breaks = 100,
         col="lightblue",
         main=sprintf("Pairwise correlation distribution of all trials for neuron %d (%s)",
                      neur,
                      trial_str),
         xlab="Pearson correlation",
         plot=T)
    
    medc <- median(all_cors)
    avgc <- mean(all_cors)
    abline(v=medc, col="red", lty=2,lwd=2)
    abline(v=avgc, col="purple", lty=2,lwd=2)
    dev.off()
    
    median_cor <- c(median_cor, medc)
    mean_cor <- c(mean_cor, avgc)
  }
}




pv_correlation <- function(mat, 
                           stim_master, 
                           window_size=30,
                           trials_to_use=c(3),
                           seconds_before=0,
                           seconds_after=8,
                           outout_path="",
                           trial_str="NA",
                           fname_ext="") {
  
  
  f_stim_ind <- c() # Sorted indices by trial
  n_trials <- c() # Vector containing how many trials for each trial type
  
  for (st in trials_to_use) {
    # Ind of trial
    stim_ind <- stim_master[which(stim_master[,7] %in% c(st)),1]
    f_stim_ind <- c(f_stim_ind, stim_ind)
    
    # Num of each trial
    n_trials <- c(n_trials, len(stim_ind))
  }
  
  
  smpl <- data.frame(Trials = rep(trial_type_name[as.character(trials_to_use)], 
                                  n_trials))
  if (sort_trials) {
    f_stim_ind <- sort(f_stim_ind)
    smpl <- data.frame(Trials = trial_type_name[as.character(stim_master[which(stim_master[,7] %in% stim_to_use),7])])
  }
  
  f_out_path <- sprintf("%s\\Individual Neurons", out_path)
  stim_mat_lst <- lapply(f_stim_ind,
                         function(si) {
                           get_time_binned_mat(mat[,(si - seconds_before * window_size):(si + seconds_after * window_size)], 
                                               8 * window_size,
                                               clean=F)
                         })
  
  pv_mat <- do.call(cbind, stim_mat_lst)
  dir.create(sprintf("%s\\%s", f_out_path, "pv_matrices\\"))
  
  cor_mat <- cor(pv_mat, pv_mat)
  dist_mat <- as.matrix(dist(t(pv_mat)))
  
  colnames(cor_mat) <- 1:nrow(cor_mat)
  rownames(cor_mat) <- 1:nrow(cor_mat)
  rownames(dist_mat) <- 1:nrow(dist_mat)
  colnames(dist_mat) <- 1:nrow(dist_mat)
  rownames(smpl) <-  1:nrow(cor_mat)
    
  cl  <- color_scheme[,"type"][as.character(trials_to_use)]
  names(cl) <- trial_type_name[as.character(trials_to_use)]
  
  clf=list(Trials=cl)
    
    
    png(filename=sprintf("%s\\%d\\Clustered_corr_matrix_%s.png", 
                         f_out_path, 
                         neur,
                         paste(trial_type_name[as.character(stim_to_use)], collapse=".")),
        width=600,
        height=600)
    
    
    pheatmap(cor_mat, 
             col=turbo(200), 
             annotation_col=smpl, 
             annotation_row=smpl,
             cluster_cols=F, 
             cluster_rows=F,
             show_rownames=F,
             show_colnames=F,
             annotation_names_col = F,
             annotation_names_row = F,
             annotation_colors=clf)    
    
    dev.off()
    
    
    png(filename=sprintf("%s\\%d\\Unclustered_corr_matrix_%s.png", 
                         f_out_path, 
                         neur,
                         paste(trial_type_name[as.character(stim_to_use)], collapse=".")),
        width=600,
        height=600)
    
    pheatmap(neur_tuning_cor_mat, 
             col=turbo(200), 
             annotation_col=smpl,
             annotation_row=smpl,
             cluster_cols=F, 
             cluster_rows=F,
             show_rownames=F,
             show_colnames=F,
             annotation_names_col = F,
             annotation_names_row=F,
             annotation_colors =clf,
             main=sprintf("Trials correlation matrix neuron %d", neur))    
    
    dev.off()
}



mult_mat <- function(mat_list, 
                     window_size, 
                     activity_threshold=0.2, 
                     fnames=F,
                     norm=F) {
  
  if(fnames) {
   mat_list <- lapply(mat_list, function(f) {load(f);return(fmat)}) 
  }
  
  avgd_mat_list <- lapply(mat_list, function(mat) {time_bin_average(mat, window_size)})
  
  mean_avgd <- lapply(mat_list, function(mat) {rowMeans(mat)})
  
  outliers  <- lapply(mean_avgd, function(avg) {which(abs(avg) > activity_threshold)})
  outliers <- unique(unlist(outliers))
  
  to_keep <- 1:nrow(avgd_mat_list[[1]])
  to_keep <- to_keep[which(!to_keep %in% outliers)]
  
  final_mat <- do.call(cbind,
                       lapply(avgd_mat_list, 
                              function(mat) {return(mat[to_keep,])}))
  
  
  outliers <- which(apply(final_mat, 1, function(r) {sum(is.nan(r)) > 0}))
  to_keep <- 1:nrow(final_mat)
  to_keep <- to_keep[which(!to_keep %in% outliers)]
  final_mat <- final_mat[to_keep,]
  
  if (norm) {
    print("Z-scoring!!!")
    final_mat2  <- smooth_ca_trace(final_mat)
    print(final_mat2[1:10,1])
    print(final_mat[1:10,1])
    final_mat <- final_mat2
  }
  
  
  return(final_mat)
}

reduce_final_mat <- function(final_mat,
                             threeD=F,
                             use_isomap=F,
                             use_lem_extra=F,
                             umap_4=F,
                             lem_4=F,
                             denoised=F) {
  
  res <- list(umap_2=NA,
              lem_2=NA,
              lem_ex=NA,
              iso_2=NA,
              umap_3=NA,
              lem_3=NA,
              iso_3=NA,
              lem_4=NA,
              umap_4=NA)
  
  if (denoised) {
    
    res$umap_3 <- reduce_dim(t(final_mat), "umap_denoised", 3)
    res$lem_3 <- reduce_dim(t(final_mat), "lem2", 3)
    return(res)
  }
  
  res$umap_2 <- reduce_dim(t(final_mat), "umap", 2)
  res$lem_2 <- reduce_dim(t(final_mat), "lem", 2)
  
  
  if (use_isomap) {
    res$iso_2 <- reduce_dim(t(final_mat), "isomap", 2)
  }
  
  if (use_lem_extra) {
    res$lem_ex <- reduce_dim(t(final_mat), "lem2", 2)
  }
  
  if (threeD) {
    res$umap_3 <- reduce_dim(t(final_mat), "umap", 3)
    res$lem_3 <- reduce_dim(t(final_mat), "lem", 3)
    
    if (use_isomap) {
      res$iso_3 <- reduce_dim(t(final_mat), "isomap", 3)
    }
  }
  
  if (lem_4) {
    res$lem_4 <- reduce_dim(t(final_mat), "lem", 4)
  }
  
  if (umap_4) {
    res$umap_4 <- reduce_dim(t(final_mat), "umap", 4)
  }
  
  return(res)
}


get_reduced_outliers <- function(reduced_df, sigma_t=3, ndim=3) {
  
  exc <- lapply(1:ndim, function(d){which(abs(scale(reduced_df[,d])) > sigma_t)})
  
  include <- which(!1:nrow(reduced_df) %in% unique(unlist(exc)))
  
  return(include)
}



run_analysis_by_mice_path <- function(path) {
  
  runs <- list.dirs(path, recursive = F)
  days <- sapply(str_split(runs, "day_"), function(l) {return(l[[2]])})
  
  for (day in days) {
    print(sprintf("Running day! %s", day))
    run_analysis_by_day_path(path, day)
  }
}


run_analysis_by_day_path <- function(path, day) {
  base_path <- path
  full_path <- sprintf("%s\\day_%s", base_path, day) 
  
  runs <- list.files(full_path, recursive=F)
  runs <- runs[which(grepl("\\.R", runs))]
  runs <- sapply(str_split(runs, ".R"), function(l){return(l[[1]])})
  
  # Load all matrices for all runs
  matrices_list <- lapply(runs, 
                          function(r) {load(sprintf("%s\\%s.R", full_path, r)); 
                                       return(fmat)})
  
  # nr <- nrow(matrices_list[[1]])
  # 
  # if (nr > 190) {
  #   print("Uniforming matrices!!!")
  #   ind <- sort(sample(1:nr, 190))
  #   
  #   matrices_list <- lapply(matrices_list,
  #                           function(mat) {return(mat[ind,])})
  # }
  
  
  
  behavior_mat_list <- lapply(list.files(sprintf("%s\\behavior\\",full_path)),
                              function(tv_mat)
                                {
                                  return(readMat(sprintf("%s\\behavior\\%s",
                                                         full_path,
                                                         tv_mat))$Stim.Master)
                                })
  
  names(matrices_list) <- runs

  if (len(behavior_mat_list) == len(matrices_list)) {
    names(behavior_mat_list) <- runs
  }
  
  analyse_matrices_all(matrices_list, behavior_mat_list, full_path, shuffled=F)
  #analyse_matrices_combined(matrices_list, behavior_mat_list, full_path)
  #analyse_matrices_solo(matrices_list, NA, full_path)
}


analyse_matrices_all <- function(matrices_list, behavior_mat_list, full_path, just_pairs=F, shuffled=F) {

  mat_list <- matrices_list
  behav_list <- behavior_mat_list
  output_path <- full_path
  
  res_path <- sprintf("%s\\%s", output_path, "results")
  dir.create(res_path)
  run_output_path <- sprintf("%s\\%s", res_path, ifelse(shuffled, "all_runs_shuffled","all_runs_cleaned"))
  dir.create(run_output_path)
  
  # fm <- mult_mat(mat_list,
  #                 window_size=30,
  #                 norm=F)
  # 
  # 
  # res <- reduce_final_mat(fm,
  #                         use_lem_extra = F,
  #                         use_isomap = F,
  #                         threeD = T,
  #                         denoised=T)
  
  #res$umap_3[,4] <- 1:nrow(res$umap_3)
  # umap_2_time <- res$umap_2
  # lem_2_time <- res$lem_2
  # umap_2_time[,3] <- 1:nrow(res$umap_2)
  # lem_2_time[,3] <- 1:nrow(res$lem_2)
  
  res = list()
  
  res$lem_3 <- get_reduced_mat_full_day(full_path, override=F, shuffled = shuffled)
  res$lem_3[,4] <- 1:nrow(res$lem_3)

  
  if (!shuffled && len(behav_list) == len(mat_list)) {
  corrected_behav <- get_corrected_behavior_matrices(behav_list[c(1:len(behav_list))],
                                                     ncol(mat_list[[1]]))
  
  color_scheme <- get_color_palettes(corrected_behav,
                                     window_size = 30,
                                     session_size = nrow(res$lem_3))
  
  color_scheme_2 <- get_color_palettes(corrected_behav,
                                       window_size = 30,
                                       session_size = nrow(res$lem_3),
                                       before=3,
                                       after=5)
  
  color_scheme$time <- viridis(nrow(res$lem_3))
  color_scheme_2$time <- viridis(nrow(res$lem_3))
  
  if ("2" %in% names(color_scheme)){
    names(color_scheme) <- c("Pavlov", "Unknown", "Reward", "Aversive", "Neutral", "Blank", "Response", "Time")
    names(color_scheme_2) <- c("Pavlov", "Unknown", "Reward", "Aversive", "Neutral", "Blank", "Response", "Time")
  } else {
    names(color_scheme) <- c("Pavlov", "Reward", "Aversive", "Neutral", "Blank", "Response", "Time")
    names(color_scheme_2) <- c("Pavlov", "Reward", "Aversive", "Neutral", "Blank", "Response", "Time")
  }
  
  assert(all(unlist(lapply(color_scheme_2, function(it) {len(it)}))))
  all(len(color_scheme_2$Pavlov) == nrow(res$lem_3))
  
  } else {
    color_scheme <- list(viridis(nrow(res$lem_3)))
    color_scheme_2 <- list(viridis(nrow(res$lem_3)))
    names(color_scheme) <- c("Time")
    names(color_scheme_2) <- c("Time")
    
  }
  # Color structure in all possible ways
  for (cs_name in names(color_scheme)) {
    
    cs <- unlist(color_scheme[cs_name])
    cs2 <- unlist(color_scheme_2[cs_name])
    
    # umap_path <- sprintf("%s\\%s", run_output_path, "umap_2d")
    # dir.create(umap_path)
    # save_3d_reduced_matrix(res$umap_2,
    #                        umap_path,
    #                        color_scheme = cs,
    #                        plot_movie = F,
    #                        prefix=sprintf("umap_2d_by_%s_first_last", cs_name))
    # 
    # save_3d_reduced_matrix(res$umap_2,
    #                        umap_path,
    #                        color_scheme = cs2,
    #                        plot_movie = F,
    #                        prefix=sprintf("umap_2d_by_%s_first_last_trail", cs_name))  
    # 
    # 
    # lem_path <- sprintf("%s\\%s", run_output_path, "lem_2d")
    # dir.create(lem_path)
    # save_3d_reduced_matrix(res$lem_2,
    #                        lem_path,
    #                        color_scheme = cs,
    #                        plot_movie = F,
    #                        prefix=sprintf("lem_2d_by_%s_first_last", cs_name))
    # 
    # save_3d_reduced_matrix(res$lem_2,
    #                        lem_path,
    #                        color_scheme = cs2,
    #                        plot_movie = F,
    #                        prefix=sprintf("lem_2d_by_%s_first_last_trail", cs_name))                         

    
    pairs_path <- sprintf("%s\\%s", run_output_path, "pairs")
    dir.create(pairs_path)

    png(sprintf("%s\\pairs_lem_3_%s_first_last.png",
                pairs_path,
                cs_name),
        width=2500, height=2500)
    pairs(res$lem_3, col=cs, pch=19, cex=3.2)
    dev.off()
    
    png(sprintf("%s\\pairs_lem_3_%s_first_last_trail.png",
                pairs_path,
                cs_name),
        width=2500, height=2500)
    pairs(res$lem_3, col=cs2, pch=19, cex=3.2)
    dev.off() 
    
    
    
    # 
    # png(sprintf("%s\\pairs_lem_2_%s_first_last.png",
    #             pairs_path,
    #             cs_name),
    #     width=2500, height=2500)
    # pairs(lem_2_time, col=cs, pch=19, cex=3.2)
    # dev.off()
    # 
    # png(sprintf("%s\\pairs_lem_2_%s_first_last_trail.png",
    #             pairs_path,
    #             cs_name),
    #     width=2500, height=2500)
    # pairs(lem_2_time, col=cs2, pch=19, cex=3.2)
    # dev.off()  
    # 
    # png(sprintf("%s\\pairs_umap_2_%s_first_last.png",
    #             pairs_path,
    #             cs_name),
    #     width=2500, height=2500)
    # pairs(umap_2_time, col=cs, pch=19, cex=3.2)
    # dev.off()
    # 
    # png(sprintf("%s\\pairs_umap_2_%s_first_last_trail.png",
    #             pairs_path,
    #             cs_name),
    #     width=2500, height=2500)
    # pairs(umap_2_time, col=cs2, pch=19, cex=3.2)
    # dev.off() 
    
  }
}

analyse_matrices_combined <- function(mat_list, behav_list, output_path) {
  
  last_fourth_matrices <- mat_list[names(mat_list)[c(1,len(mat_list))]]
  runs  <- names(mat_list)
  res_path <- sprintf("%s\\%s", output_path, "results")
  dir.create(res_path)
  run_output_path <- sprintf("%s\\%s", res_path, "multiple_runs")
  dir.create(run_output_path)
  
  res <- mult_mat(last_fourth_matrices,
                  window_size=15,
                  use_lem_extra = T,
                  use_isomap = F,
                  threeD = T,
                  umap_4 = T,
                  lem_4 = T)
  
  res$umap_4[,5] <- 1:nrow(res$umap_4)
  res$lem_4[,5] <- 1:nrow(res$lem_4)
  res$umap_3[,4] <- 1:nrow(res$umap_3)
  res$lem_3[,4] <- 1:nrow(res$lem_3)
  umap_2_time <- res$umap_2
  lem_2_time <- res$lem_2
  umap_2_time[,3] <- 1:nrow(res$umap_2)
  lem_2_time[,3] <- 1:nrow(res$lem_2)
  
  colnames(res$umap_4) <- c("Component 1", "Component 2", "Component 3", "Component 4", "Time")
  colnames(res$lem_4) <- c("Component 1", "Component 2", "Component 3", "Component 4", "Time")
  
  corrected_behav <- get_corrected_behavior_matrices(behav_list[c(1,len(behav_list))],
                                                     ncol(mat_list[[1]]))
  
  color_scheme <- get_color_palettes(corrected_behav,
                                     window_size = 15,
                                     session_size = nrow(res$umap_2))
  
  color_scheme_2 <- get_color_palettes(corrected_behav,
                                     window_size = 15,
                                     session_size = nrow(res$umap_2),
                                     before=4,
                                     after=8)
  
  color_scheme$time <- viridis(nrow(res$umap_2))
  color_scheme_2$time <- viridis(nrow(res$umap_2))
  
  if ("2" %in% names(color_scheme)){
    names(color_scheme) <- c("Pavlov", "Unknown", "Reward", "Aversive", "Neutral", "Blank", "Response", "Time")
    names(color_scheme_2) <- c("Pavlov", "Unknown", "Reward", "Aversive", "Neutral", "Blank", "Response", "Time")
  } else {
    names(color_scheme) <- c("Pavlov", "Reward", "Aversive", "Neutral", "Blank", "Response", "Time")
    names(color_scheme_2) <- c("Pavlov", "Reward", "Aversive", "Neutral", "Blank", "Response", "Time")
  }
  
  assert(all(unlist(lapply(color_scheme_2, function(it) {len(it)}))))
  all(len(color_scheme_2$Pavlov) == nrow(res$umap_2))
  # Color structure in all possible ways
  for (cs_name in names(color_scheme)) {
    
  cs <- unlist(color_scheme[cs_name])
  cs2 <- unlist(color_scheme_2[cs_name])
    
  umap_path <- sprintf("%s\\%s", run_output_path, "umap_2d")
  dir.create(umap_path)
  save_3d_reduced_matrix(res$umap_2,
                         umap_path,
                         color_scheme = cs,
                         plot_movie = F,
                         prefix=sprintf("umap_2d_by_%s_first_last", cs_name))
  
  save_3d_reduced_matrix(res$umap_2,
                         umap_path,
                         color_scheme = cs2,
                         plot_movie = F,
                         prefix=sprintf("umap_2d_by_%s_first_last_trail", cs_name))  
  
  
  lem_path <- sprintf("%s\\%s", run_output_path, "lem_2d")
  dir.create(lem_path)
  save_3d_reduced_matrix(res$lem_2,
                         lem_path,
                         color_scheme = cs,
                         plot_movie = F,
                         prefix=sprintf("lem_2d_by_%s_first_last", cs_name))
  
  save_3d_reduced_matrix(res$lem_2,
                         lem_path,
                         color_scheme = cs2,
                         plot_movie = F,
                         prefix=sprintf("lem_2d_by_%s_first_last_trail", cs_name))                         
  
  
  lem_path <- sprintf("%s\\%s", run_output_path, "lem_ex_2d")
  dir.create(lem_path)
  save_3d_reduced_matrix(res$lem_ex,
                         lem_path,
                         color_scheme = cs,
                         plot_movie = F,
                         prefix=sprintf("lem_ex_by_%s_first_last", cs_name))

  save_3d_reduced_matrix(res$lem_ex,
                         lem_path,
                         color_scheme = cs2,
                         plot_movie = F,
                         prefix=sprintf("lem_ex_by_%s_first_last_trail", cs_name))  
  
  
  pairs_path <- sprintf("%s\\%s", run_output_path, "pairs")
  dir.create(pairs_path)

    png(sprintf("%s\\pairs_lem_4_%s_first_last.png",
              pairs_path,
              cs_name),
      width=2500, height=2500)
  pairs(res$lem_4, col=cs, pch=19, cex=3.2)
  dev.off()
  
  png(sprintf("%s\\pairs_lem_4_%s_first_last_trail.png",
              pairs_path,
              cs_name),
      width=2500, height=2500)
  pairs(res$lem_4, col=cs2, pch=19, cex=3.2)
  dev.off()  
  
  png(sprintf("%s\\pairs_umap_4_%s_first_last.png",
              pairs_path,
              cs_name),
      width=2500, height=2500)
  pairs(res$umap_4, col=cs, pch=19, cex=3.2)
  dev.off()

  png(sprintf("%s\\pairs_umap_4_%s_first_last_trail.png",
              pairs_path,
              cs_name),
      width=2500, height=2500)
  pairs(res$umap_4, col=cs2, pch=19, cex=3.2)
  dev.off() 
  
  
  png(sprintf("%s\\pairs_lem_3_%s_first_last.png",
              pairs_path,
              cs_name),
      width=2500, height=2500)
  pairs(res$lem_3, col=cs, pch=19, cex=3.2)
  dev.off()
  
  png(sprintf("%s\\pairs_lem_3_%s_first_last_trail.png",
              pairs_path,
              cs_name),
      width=2500, height=2500)
  pairs(res$lem_3, col=cs2, pch=19, cex=3.2)
  dev.off()  
  
  png(sprintf("%s\\pairs_umap_3_%s_first_last.png",
              pairs_path,
              cs_name),
      width=2500, height=2500)
  pairs(res$umap_3, col=cs, pch=19, cex=3.2)
  dev.off()
  
  png(sprintf("%s\\pairs_umap_3_%s_first_last_trail.png",
              pairs_path,
              cs_name),
      width=2500, height=2500)
  pairs(res$umap_3, col=cs2, pch=19, cex=3.2)
  dev.off() 
  
  

  
  png(sprintf("%s\\pairs_lem_2_%s_first_last.png",
              pairs_path,
              cs_name),
      width=2500, height=2500)
  pairs(lem_2_time, col=cs, pch=19, cex=3.2)
  dev.off()
  
  png(sprintf("%s\\pairs_lem_2_%s_first_last_trail.png",
              pairs_path,
              cs_name),
      width=2500, height=2500)
  pairs(lem_2_time, col=cs2, pch=19, cex=3.2)
  dev.off()  
  
  png(sprintf("%s\\pairs_umap_2_%s_first_last.png",
              pairs_path,
              cs_name),
      width=2500, height=2500)
  pairs(umap_2_time, col=cs, pch=19, cex=3.2)
  dev.off()
  
  png(sprintf("%s\\pairs_umap_2_%s_first_last_trail.png",
              pairs_path,
              cs_name),
      width=2500, height=2500)
  pairs(umap_2_time, col=cs2, pch=19, cex=3.2)
  dev.off() 
  
  }
}


analyse_matrices_solo <- function(mat_list, behav_list, output_path) {
  
  runs  <- names(mat_list)
  res_path <- sprintf("%s\\%s", output_path, "results")
  dir.create(res_path)
  res_path <- sprintf("%s\\%s", res_path, "single_runs")
  dir.create(res_path)
  for (run in runs) {
    print(sprintf("Analysis of run %s!", run))
    
    res <- mult_mat(mat_list[run],
                    window_size=30,
                    use_lem_extra = T,
                    use_isomap = T,
                    threeD = F)
    

    run_output_path <- sprintf("%s\\%s", res_path, run)
    dir.create(run_output_path)
    
    
    corrected_behav <- get_corrected_behavior_matrices(behavior_mat_list,
                                                         ncol(mat_list[[1]]))
    
    color_scheme <- get_color_palettes(list(corrected_behav[[run]]),
                                       window_size = 30,
                                       session_size = nrow(res$umap_2))
    
    
    umap_path <- sprintf("%s\\%s", run_output_path, "umap_2d")
    dir.create(umap_path)
    save_3d_reduced_matrix(res$umap_2,
                           umap_path,
                           color_scheme = viridis(nrow(res$umap_2)),
                           plot_movie = F,
                           prefix="umap_2d_by_time")
    
    lem_path <- sprintf("%s\\%s", run_output_path, "lem_2d")
    dir.create(lem_path)
    save_3d_reduced_matrix(res$lem_2,
                           lem_path,
                           color_scheme = viridis(nrow(res$lem_2)),
                           plot_movie = F,
                           prefix="lem_2d_by_time")
    
    
    lem_path <- sprintf("%s\\%s", run_output_path, "lem_ex_2d")
    dir.create(lem_path)
    save_3d_reduced_matrix(res$lem_ex,
                           lem_path,
                           color_scheme = viridis(nrow(res$lem_ex)),
                           plot_movie = F,
                           prefix="lem_ex_2d_by_time")
    
    iso_path <- sprintf("%s\\%s", run_output_path, "iso_2d")
    dir.create(iso_path)
    save_3d_reduced_matrix(res$iso_2,
                           iso_path,
                           color_scheme = viridis(nrow(res$iso_2)),
                           plot_movie = F,
                           prefix="iso_2d_by_time")
  }
  }


save_3d_reduced_matrix <- function(df, 
                                   out_path, 
                                   color_scheme,
                                   plot_movie=F,
                                   prefix="",
                                   main="") {
  
  time_view_matrix_1 <- list(c( 0.996785522, -0.080081321, -0.002386985,   0),
                             c( 0.002787631,  0.004891566,  0.999984145,   0),
                             c(-0.080068380, -0.996776283,  0.005099070,   0),
                             c( 0.000000000,  0.000000000,  0.000000000,   1))
  
  side_view_matrix_1 <- list(c( 0.5039489, -0.8542808, -0.1274353,    0),
                             c( 0.3166651,  0.0454708,  0.9474469,    0),
                             c(-0.8035911, -0.5178189,  0.2934358,    0),
                             c( 0.0000000,  0.0000000,  0.0000000,    1))
  
  
  time_view_matrix_2 <-  list(c(-0.99928027,  0.02641503, -0.02722498,   0),
                              c(-0.03079081, -0.14563999,  0.98885846,   0),
                              c( 0.02215547,  0.98898482,  0.14634864,   0),
                              c( 0.00000000,  0.00000000,  0.00000000,   1))
                            
  
  side_view_matrix_2 <- list(c(-0.5639776,  0.7865990, 0.2513785,    0),
                             c(-0.2832605, -0.4702138, 0.8358603,    0),
                             c( 0.7756886,  0.4002008, 0.4880027,    0),
                             c( 0.0000000,  0.0000000, 0.0000000,    1))
  
  bottom_view_matrix <- list(c(-0.7773179, 0.6162394,  0.1265937,    0),
                             c( 0.4466669, 0.3989020,  0.8008534,    0),
                             c( 0.4430187, 0.6790625, -0.5853276,    0),
                             c( 0.0000000, 0.0000000,  0.0000000,    1))
  
  time_view_matrix_1 <- do.call(rbind, time_view_matrix_1)
  time_view_matrix_2 <- do.call(rbind, time_view_matrix_2)
  side_view_matrix_1 <- do.call(rbind, side_view_matrix_1)
  side_view_matrix_2 <- do.call(rbind, side_view_matrix_2)
  bottom_view_matrix <- do.call(rbind, bottom_view_matrix)
  
  views <- list(Time_1=time_view_matrix_1,
         Time_2=time_view_matrix_2,
         Side_1=side_view_matrix_1,
         side_2=side_view_matrix_2,
         Bottom=bottom_view_matrix)  
  
  for (line_type in c("l", "p")) {
    
  line_str <- ifelse(line_type == "l", "lines", "points")
  if (ncol(df) == 2) {
    
    outliers <- get_reduced_outliers(df, 3,2)
    plot3d(y=df[,1],
           ylab="Component 1",
           ylim=c(min(df[outliers,1]) * 1.2, 
                  max(df[outliers,1]) * 1.2),
           z=df[,2],
           zlab="Component 2",
           zlim=c(min(df[outliers,2]) * 1.2, 
                  max(df[outliers,2]) * 1.2),
           x=c(1:nrow(df)),
           xlab="Time",
           col=color_scheme,
           type=line_type,
           size=10)
    
    par3d(zoom=0.8)
    par3d(windowRect = c(0, 0, 1000, 1000))
    rgl.snapshot(sprintf("%s\\%s_%s", out_path, prefix, 
                         sprintf("%s_%s", line_str, 'default.png')), 
                 fmt = 'png')
    
    if (plot_movie) {
      movie3d(spin3d(axis = c(0, 0, 1), rpm=8), duration = 8,  dir = out_path, 
              movie=sprintf("%s_%s_animated", prefix, line_str))
    
    }
    
    lapply(names(views),
           function(view_name) {
             print(sprintf("Saving for view %s", view_name))
             par3d(userMatrix=views[[view_name]])
             rgl.snapshot(sprintf("%s\\%s_%s", out_path, prefix, 
                          sprintf('%s_%s.png', line_str, view_name)), 
                          fmt = 'png')
           })
    
    while (rgl.cur() > 0) { rgl.close() }
    
  } else if (ncol(df) == 3) {
    outliers <- get_reduced_outliers(df, 3,2)
    plot3d(y=df[,1],
           ylab="Component 1",
           ylim=c(min(df[outliers,1]) * 1.2, 
                  max(df[outliers,1]) * 1.2),
           z=df[,2],
           zlab="Component 2",
           zlim=c(min(df[outliers,2]) * 1.2, 
                  max(df[outliers,2]) * 1.2),
           x=df[,3],
           xlim=c(min(df[outliers,3]) * 1.2, 
                  max(df[outliers,3]) * 1.2),
           xlab="Component 3",
           col=color_scheme,
           type=line_type)
    
    par3d(zoom=0.8)
    par3d(windowRect = c(0, 0, 1000, 1000))
    rgl.snapshot(sprintf("%s\\%s_%s", out_path, prefix, 'default.png'), 
                 fmt = 'png')
    
    movie3d(spin3d(axis = c(0,0,1), rpm = 8), 
            duration = 15, 
            fps = 50)

    if (plot_movie) {
      movie3d(spin3d(axis = c(0, 0, 1), rpm=8), duration = 8,  dir = out_path, 
              movie=sprintf("%s_%s_animated", prefix, line_str))
      
    }
    
    lapply(names(views),
           function(view_name) {
             print(sprintf("Saving for view %s", view_name))
             par3d(userMatrix=views[[view_name]])
             rgl.snapshot(sprintf("%s\\%s_%s", out_path, prefix, 
                                  sprintf('%s.png', view_name)), 
                          fmt = 'png')
           })
    
    while (rgl.cur() > 0) { rgl.close() }
  }
  }
} 


get_corrected_behavior_matrices <- function(behav_mat_list, session_size) {
  
  fm <- behav_mat_list
  
  for (i in 1:len(fm)) {
    fm[[i]][,1] <- fm[[i]][,1] + (i-1) * session_size
    fm[[i]][,6] <- fm[[i]][,6] + (i-1) * session_size
  }
  
  return(fm)
}


get_color_palettes <- function(behav_mat_list,
                               session_size,
                               window_size,
                               before=0,
                               after=0) {
  df <- behav_mat_list # Just for easy naming
  color_scheme <- list(c("deeppink", "chocolate1", "cadetblue1", "chartreuse3"),
                       c("darkseagreen1", "chocolate1", "cadetblue1", "darkseagreen3"),
                       c("dodgerblue1", "chocolate1", "cadetblue1", "coral3"),
                       c("darkorchid2", "chocolate1", "cadetblue1", "aquamarine3"),
                       c("yellow1", "chocolate1", "cadetblue1", "yellow1"),
                       c("black", "chocolate1", "cadetblue1", "black"))
  color_scheme <- list(c("deeppink", "chocolate1", "cadetblue1", "chartreuse3"),
                       c("darkseagreen1", "chocolate1", "cadetblue1", "darkseagreen3"),
                       c("dodgerblue1", "chocolate1", "cadetblue1", "coral3"),
                       c("darkorchid2", "chocolate1", "cadetblue1", "aquamarine3"),
                       c("yellow1", "chocolate1", "cadetblue1", "yellow1"),
                       c("black", "chocolate1", "cadetblue1", "black"))
  
  color_scheme <- do.call(rbind, color_scheme)
  rownames(color_scheme) <- c("3","4","1","5","6", "2")
  colnames(color_scheme) <- c("type", "lick", "nolick", "iti")
  
  trials_type <- unlist(lapply(as.list(df), function(sdf) {sdf[,7]}))
  trials_ind <- unlist(lapply(as.list(df), function(sdf) {sdf[,1]}))
  response_ind <- unlist(lapply(as.list(df), function(sdf) {sdf[,6]}))
   
  trials <- sort(unique(trials_type))
  col_scheme_list <- list()
  
  for (t in trials) {
    print(sprintf("Using %d, %s", t, color_scheme[as.character(t), "type"]))
    specific_t_ind <- which(trials_type == t)
    
    binned_trials <- get_binned_index(trials_ind[specific_t_ind], window_size)
    binned_trials <- as.vector(sapply(binned_trials, function(i) {c((i-before):(i+after))}))
    
    #col_f <- rep(adjustcolor("gray60", alpha=0.3), times=session_size)
    #col_f[binned_trials] <- color_scheme[as.character(t), "type"]
    
    col_f <- adjustcolor((viridis(session_size)), alpha=0.1)[1:session_size]
    col_f[binned_trials] <- "red" #color_scheme[as.character(t), "type"]
    col_scheme_list <- append(col_scheme_list, list(col_f))
  }
  
  response_ind <- response_ind[!is.nan(response_ind)]
  binned_trials <- get_binned_index(response_ind, window_size)
  binned_trials <- as.vector(sapply(binned_trials, function(i) {c((i-before):(i+after))}))
  
  # col_f <- rep(adjustcolor("gray60", alpha=0.3), times=session_size)
  # col_f[binned_trials] <- "cadetblue1"
  
  col_f <- adjustcolor((viridis(session_size)), alpha=0.1)[1:session_size]
  col_f[binned_trials] <- "red"
  col_scheme_list <- append(col_scheme_list, list(col_f))
  names(col_scheme_list) <- c(trials,"res")
  
  
  return(col_scheme_list)
}
  




pl <- function(i) {
  mcol <- adjustcolor(rbPal(100)[as.numeric(cut(final_mat[,i],breaks = 100))], alpha=0.15)
  
  
  plot(res$lem_2[,1], res$lem_2[,2], col=mcol, pch=19)
}



rate_correlation_mat <- function(matrices_list, behavior_mat_list) {
  trials_rate_mat <- c()
  for (i in 1:len(matrices_list)) {
  
      mt <- matrices_list[[i]]
      #mt <- t(apply(mt, 1, scale))
      bhv <- behavior_mat_list[[i]]
      
      for (t in unique(bhv[,7])) {
        print(sprintf("Computing rate for trial %d, session %d", t,i))
        ind <- unlist(lapply(bhv[(bhv[,7] == t),1], function(i) {return(c((i - 1 * 31):(i + 6*31)))}))
        rate <- rowMeans(mt[,ind])
        trials_rate_mat <- cbind(trials_rate_mat, rate)
      }
  }
}



run_all <- function(path, indices=NA) {
  paths <- list.dirs(path, recursive=F)[grep("IC", list.dirs(path, recursive=F))]
  
  ind_to_use <- 1:len(paths)
  if(sum(!is.na(indices) > 0)){
    ind_to_use <- indices
    print(sprintf("Using indices: %s", paste(ind_to_use, collapse= " ")))
  }
  
   
  for(p in paths[ind_to_use]) {
    mice_num <- str_split(p, "IC")[[1]][[2]]
    print(sprintf("Running mice %s", mice_num))
    run_analysis_by_mice_path(sprintf("%s\\IC%s", path, mice_num))
   
  }
}





pairs_analysis <- function(path, 
                           across_mice=F, 
                           nclusters=20,
                           scale_cluster_dist_matrix=T,
                           compare_shuffle=F,
                           whole_matrix=T) {
  
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
  
  
  # Pairs across mice, and remove pairs from within mice  
  if (across_mice) {
    all_days <- c()
    
    for (p in paths) {
      mice_num <- str_split(p, "IC")[[1]][[2]]
      runs <- list.dirs(p, recursive = F)
      days <- sapply(str_split(runs, "day_"), function(l) {return(l[[2]])})
      all_days <- c(all_days, sprintf("%s\\IC%s\\day_%s", path, mice_num, days))
    }
  
    # All possible combinations of all days
    pairs_mat <- combn(all_days,2)
    across <- lapply(seq_len(ncol(pairs_mat)), function(i) pairs_mat[,i])
    
    # Remove days which are within mice
    path_pairs <- across[which(!across %in% path_pairs)]
  }
  
  
  
  structure_correlations <- c()
  
  for (pair in path_pairs) {
    day1 <- pair[1]
    day2 <- pair[2]
    
    s1 <- F
    s2 <- F
    
    if (compare_shuffle) {
      # Randomly select a shuffle comparsion
      s1 <- sample(c(T,F),1)
      s2 <- !s1
      print(s1)
      print(s2)
    }

    day1_mat <- get_reduced_mat_full_day(day1, shuffled=s1)
    day2_mat <- get_reduced_mat_full_day(day2, shuffled=s2)
    
    if (scale_cluster_dist_matrix) {
      km_day1 <- kmeans(apply(day1_mat, 2, scale), nclusters, iter.max=300)
      km_day2 <- kmeans(apply(day2_mat, 2, scale), nclusters, iter.max=300)
    } else {
      km_day1 <- kmeans(day1_mat, nclusters, iter.max=300)
      km_day2 <- kmeans(day2_mat, nclusters, iter.max=300)
    }
    
    dist_matrix_day1 <- as.matrix(dist(km_day1$centers))
    dist_matrix_day2 <- as.matrix(dist(km_day2$centers))
    
    h1 <- heatmap(dist_matrix_day1, plot=F)
    h2 <- heatmap(dist_matrix_day2, plot=F)
    
    if (whole_matrix) {
      central_tendency_cor <- (cor(c(dist_matrix_day1[h1$rowInd, h1$colInd]), 
                                         c(dist_matrix_day2[h2$rowInd, h2$colInd])))
    
    } else {
      central_tendency_cor <- median(cor(dist_matrix_day1[h1$rowInd, h1$colInd], 
                                         dist_matrix_day2[h2$rowInd, h2$colInd]))      
    }
    
    print(sprintf("Got mean cor(%f) between %s to %s", 
                  central_tendency_cor, 
                  str_split(day1, "data\\\\")[[1]][2], 
                  str_split(day2, "data\\\\")[[1]][2]))
    
    structure_correlations <- c(structure_correlations, central_tendency_cor)
  }
  
  return(structure_correlations)
}


 get_reduced_mat_full_day <- function(day_path, 
                                     type="lem2", 
                                     ndim=3, 
                                     window_size=30,
                                     normed=T,
                                     shuffled=F,
                                     time_shuffled=T,
                                     matrix_subpath="reduced_matrices_full_day",
                                     override=F,
                                     knn1=0.275,
                                     knn2=0.075) {
  
  knn1 <- round(knn1, digits=3)
  knn2 <- round(knn2, digits=3)
  
  reduced_mat_name <- sprintf("t%s_d%d_w%d_knnf%.3f_knns%.3f_zs%s%s",
                              type,
                              ndim,
                              window_size,
                              knn1,
                              knn2,
                              ifelse(normed, "1", "0"),
                              ifelse(shuffled, ifelse(time_shuffled, "_time_shuffled", "_cell_shuffled"), ""))
  
  # Firstly check whether there exists a matrice for that day
  if (!override && matrix_subpath %in% list.dirs(day_path, recursive = F, full.names = F)) {
    
    if (sprintf("%s.R", reduced_mat_name) %in% 
        list.files(sprintf("%s\\%s\\", day_path, matrix_subpath), full.names = F)) {
      
      print(sprintf("Found reduced matrix %s, loading", reduced_mat_name))
      load(sprintf("%s\\%s\\%s.R", day_path, matrix_subpath, reduced_mat_name))
      return(reduced)
    }
  }
  
  
  print(sprintf("Reduced matrix %s does not exist! creating!", reduced_mat_name))
  
  runs <- list.files(day_path, recursive=F)
  runs <- runs[which(grepl("\\.R", runs))]
  runs <- sapply(str_split(runs, ".R"), function(l){return(l[[1]])})
  
  # Load all matrices for all runs
  matrices_list <- lapply(runs, 
                          function(r) {load(sprintf("%s\\%s.R", day_path, r)); 
                            return(fmat)})
  
  
  final_mat <- mult_mat(matrices_list, window_size=window_size, norm=normed)
  
  if (shuffled) {
    if (time_shuffled) {
      print("Time shuffling!")
      final_mat <- generate_shuffled_matrices(final_mat, time_shuffled=T)
    } else {
      
      print("Cell shuffling!")
      final_mat <- generate_shuffled_matrices(final_mat, time_shuffled=F)
    }
  }
  
  reduced <- reduce_dim(t(final_mat), type, ndim, knn1=knn1,knn2=knn2)
  
  
  dir.create(sprintf("%s\\%s", day_path, matrix_subpath))
  save(reduced, file=sprintf("%s\\%s\\%s.R", day_path, matrix_subpath, reduced_mat_name))
  print(sprintf("Saving mat! %s", reduced_mat_name))
  return(reduced)
}

generate_shuffled_matrices <- function(final_mat, time_shuffled=T) {
  
  if (time_shuffled) {
    shuffled_mat <- t(apply(final_mat, 1, function(row) {return(row[sample(1:len(row), len(row))])}))
    return(shuffled_mat)
  }
  
  
  shuffled_mat_cells <-  apply(final_mat, 2, function(col) {return(col[sample(1:len(col), len(col))])})
  return(shuffled_mat_cells)
}


analyse_structure <- function(day_path, nclusters=NA, centroid_threshold=NA, shuffle=F, hdbscan_pts=NA, dbscan_eps = NA, dbscan_minpts=NA) {
  mat <- get_reduced_mat_full_day(day_path, shuffled = shuffle)
  eucl_dist <- function(x1, x2){
    return(sqrt(sum((x1 - x2)^2)))
  }
  
  cp = viridis(nrow(mat))
  
  if(!is.na(nclusters)) {
    km <- kmeans(mat, nclusters, iter.max=300)
    cp = viridis(nclusters)[km$cluster]
    barplot(table(km$cluster))
    
    
    if (!is.na(centroid_threshold)) {
      
      keep <- c()
      for (cluster_id in unique(km$cluster)) {
        centroid <- km$centers[cluster_id,]
        temp_mat <- mat[which(km$cluster == cluster_id),]
        distances <- apply(temp_mat, 1, function(pt) {eucl_dist(pt, centroid)})
        scaled_distances <- distances / sd(distances)

        
        cluster_keep <- which(abs(scaled_distances) < centroid_threshold)
        print(sprintf("Removing %d points from cluster %d", 
                      nrow(temp_mat) - len(cluster_keep),
                      cluster_id))
        
        keep <- c(keep, rownames(temp_mat[cluster_keep,]))
      }
      
      print(sprintf("Totally removed %d points", nrow(mat) - len(keep)))
      
      mat <- mat[keep,]
      cp <- cp[keep]
    }
  } else if (!is.na(hdbscan_pts)) {
    hc <- hdbscan(mat, hdbscan_pts)
    
    mat <- mat[which(hc$cluster != 0),]
    
    print(sprintf("Using hdbscan! Removing %d", sum(hc$cluster ==0)))
    cp <- viridis(max(hc$cluster))[hc$cluster[which(hc$cluster != 0)]]
  } else if (!is.na(dbscan_eps) && !is.na(dbscan_minpts)) {
    hc <- dbscan(mat, dbscan_eps, dbscan_minpts)
    
    mat <- mat[which(hc$cluster != 0),]
    
    print(sprintf("Using hdbscan! Removing %d", sum(hc$cluster ==0)))
    cp <- viridis(max(hc$cluster))[hc$cluster[which(hc$cluster != 0)]]
  }
  
  
  plot_3d_pc(mat, mcol_pal=cp) 
  
}


tune_knn_parameters_lem <- function(path) {
  
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
  
  
  knn1_range <- seq(0.275,0.525, by=0.05)
  knn2_range <- seq(0.05, 0.075, by=0.025)
  print(knn1_range)
  print(knn2_range)
  
  
  for (knn1 in knn1_range) {
    for (knn2 in knn2_range)
      for (path_to_tune in paths_all) {
        
        print(sprintf("Generating lem2 reduction knn1(%f) knn2(%f) for path (%s)",
                      knn1,
                      knn2,
                      path_to_tune))
        
        get_reduced_mat_full_day(path_to_tune,
                                 knn1 = knn1,
                                 knn2 = knn2,
                                 override = F)
        
        get_reduced_mat_full_day(path_to_tune,
                                 knn1 = knn1,
                                 knn2 = knn2,
                                 override = F,
                                 shuffled = T)              
    }
  }
}


