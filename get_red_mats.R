
source(sprintf ("%s/matrix_presets.R", str_split(rstudioapi::getSourceEditorContext()$path, "get_red_mats.R")[[1]][1]))
                 
len <- length
euc_dist <- function(a, b) sqrt(sum((a - b)^2))

get_binned_index <- function(i, window_size) {
  rst <- floor((i)/window_size) + ifelse(i %% window_size == 0, 0,1)
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
  
  if (window_size == 1) {
    return(spike_train)  
  }
  
  indices <- 1:(ncol(spike_train) / window_size)
  
  
  
  final_spike_train <- c()
  for (ind in indices) {
    subset_ind <- ((ind - 1 ) * window_size + 1):(ind * window_size)
    final_spike_train <- cbind(final_spike_train,
                               rowMeans(spike_train[,subset_ind]))
  }
  
  return(final_spike_train)
}



reduce_dim <- function(data, method, ndim=3, knn1=0.275, knn2=0.075) {
  print(sprintf("performing dim reduction on method! %s, data dimensonality (%d x %d)", 
                method,
                dim(data)[1],
                dim(data)[2]
  ))
  
  nn1 = floor(nrow(data) * knn1)
  nn2 = floor(nrow(data) * knn2)
  
  if (method == "tsne") {
    red <- Rtsne(data, dims = ndim, perplexity=nrow(data) * 0.3, check_duplicates = F)
    red_df <- do.call(cbind, lapply(1:ndim, function(i) {red$Y[,i]}))
    red_df <- as.data.frame(red_df)
    
  } else if (method == "lem") {
    red <- dimRed::embed(data, "LaplacianEigenmaps", ndim=ndim)
    red_df <- do.call(cbind, lapply(1:ndim, function(i) {red@data@data[,i]}))
    red_df <- as.data.frame(red_df)
  } else if (method == "isomap") {
    print(sprintf("Using KNN1=%d",
                  floor(nrow(data) * knn1)))
    red <- dimRed::embed(data, "Isomap", ndim=ndim, k=floor(nrow(data) * knn1))
    red_df <- do.call(cbind, lapply(1:ndim, function(i) {red@data@data[,i]}))
    red_df <- as.data.frame(red_df)
  } else if (method == "lem2") {
    
    print(sprintf("Using KNN1=%d, KNN2=%d",nn1, nn2))
    
    red <- dimRed::embed(data, "LaplacianEigenmaps", ndim=15, nn1, t=Inf)
    
    red <- do.call(cbind, lapply(1:15, function(i) {red@data@data[,i]}))
    red <- as.data.frame(red)    
    
    red <- dimRed::embed(red, "LaplacianEigenmaps",ndim=ndim, nn2, t=Inf)
    
    red_df <- do.call(cbind, lapply(1:ndim, function(i) {red@data@data[,i]}))
    red_df <- as.data.frame(red_df)
  } else if (method == "umap") {
    
    print(sprintf("Using KNN1=%d, min_dist=%f",
                  nn1,
                  knn2))
    red <- umap(data, 
                n_components=ndim,
                metric="cosine",
                min_dist=knn2,
                n_neighbors=nn1)
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
                  n_neighbors=nn1)
    }
    red_df <- do.call(cbind, lapply(1:ndim, function(i) {red$layout[,i]}))
    red_df <- as.data.frame(red_df)
    
  } else if (method == "mds") {
    red <- cmdscale(dist(data),eig=TRUE, k=ndim)
    red_df <- do.call(cbind, lapply(1:ndim, function(i) {red$points[,i]}))
    red_df <- as.data.frame(red_df)
    
  } else if (method == "dmaps") {
    red <- dimRed::embed(data, "DiffusionMaps", ndim=ndim)
    red_df <- do.call(cbind, lapply(1:ndim, function(i) {red@data@data[,i]}))
    red_df <- as.data.frame(red_df)
  } else if (method == "hlle") {
    
    if (knn1 != 0) {
      print(sprintf("Using KNN1=%d", nn2))
      red <- dimRed::embed(data, "HLLE", ndim=ndim, knn=nn1)
    } else {
      red <- dimRed::embed(data, "HLLE", ndim=ndim)
    }
    
    red_df <- do.call(cbind, lapply(1:ndim, function(i) {red@data@data[,i]}))
    red_df <- as.data.frame(red_df)
  } else {
    return(data) 
  }
  
  colnames(red_df) <- c("x","y","z")[1:ndim]
  
  if (ndim > 3){
    colnames(red_df)[4:ndim] <- 4:(ndim)
  }
  return(red_df)
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
  
  rm <- rowMeans(final_mat)
  outliers <- which(is.nan(rm))
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

smooth_ca_trace <- function(ca_mat) {
  return(t(apply(ca_mat, 1, FUN=scale)))
}


get_reduced_mat_full_day_control <- 
  function(day_path, 
           type="lem2", 
           ndim=3, 
           window_size=30,
           normed=T,
           shuffled=F,
           time_shuffled=T,
           matrix_subpath="reduced_matrices_full_day",
           override=F,
           knn1=0.275,
           knn2=0.075, 
           just_original_mat=F,
           in_matlab=T,
           activity_threshold=0.2,
           first_q=0,
           last_q=0,
           chunk=0) {
    
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
    
    
    if (activity_threshold != 0.2) {
      reduced_mat_name <- sprintf("%s_act%.3f", reduced_mat_name, activity_threshold)
    }
    
    if (all(chunk != 0)) {
      reduced_mat_name <- sprintf("%s_%s", reduced_mat_name, paste(chunk, collapse="_"))
    }
    
    if ((first_q != 0) || (last_q != 0)) {
      if (first_q != 0) {
        reduced_mat_name <- sprintf("%s_fq%.3f", reduced_mat_name, first_q)
      } else {
        reduced_mat_name <- sprintf("%s_lq%.3f", reduced_mat_name, first_q)
      }
    }
    
    # Firstly check whether there exists a matrice for that day
    if (!just_original_mat && !override && matrix_subpath %in% list.dirs(day_path, recursive = F, full.names = F)) {
      
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
    
    
    final_mat <- mult_mat(list(merged_mat), window_size=window_size, norm=normed)
    
    if ((first_q != 0) || (last_q != 0)) {
      n_frames <- ncol(final_mat)
      
      if (first_q != 0) {
        ind <- 1:(n_frames * first_q)
      } else {
        ind <- (n_frames - last_q * n_frames):n_frames
      }
      
      
      print(sprintf("Subsetting mat (%d) frames", len(ind)))
      final_mat <- final_mat[,ind]
    }
    
    if (just_original_mat) {
      print("Returning original matrix")
      return(final_mat)
    }
    
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
                                     knn2=0.075, 
                                     just_original_mat=F,
                                     activity_threshold=0.2,
                                     first_q=0,
                                     last_q=0,
                                     chunk=0,
                                     control=F,
                                     preset=list()) {
  
  
  if (len(preset) != 0 ) {
    reduced_mat_name = preset$name
  } else if (just_original_mat) { 
    reduced_mat_name <- sprintf("original_w%d_act%.3f_zs%s%s",
                                window_size,
                                activity_threshold, 
                                ifelse(normed, "1", "0"),
                                ifelse(shuffled, ifelse(time_shuffled, "_time_shuffled", "_cell_shuffled"), "")                                )
  } else {
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
    
    if (activity_threshold != 0.2) {
      reduced_mat_name <- sprintf("%s_act%.3f", reduced_mat_name, activity_threshold)
    }
    
    if (all(chunk != 0)) {
      reduced_mat_name <- sprintf("%s_%s", reduced_mat_name, paste(chunk, collapse="_"))
    }
    
    
    if ((first_q != 0) || (last_q != 0)) {
      if (first_q != 0) {
        reduced_mat_name <- sprintf("%s_fq%.3f", reduced_mat_name, first_q)
      } else {
        reduced_mat_name <- sprintf("%s_lq%.3f", reduced_mat_name, first_q)
      }
    }
  }
  
  
  # Firstly check whether there exists a matrice for that day
  if (!override && matrix_subpath %in% list.dirs(day_path, recursive = F, full.names = F)) {
    
    if (sprintf("%s.R", reduced_mat_name) %in% 
        list.files(sprintf("%s\\%s\\", day_path, matrix_subpath), full.names = F)) {
      
      print(sprintf("Found reduced matrix %s, loading", reduced_mat_name))
      load(sprintf("%s\\%s\\%s.R", day_path, matrix_subpath, reduced_mat_name))
      return(reduced)
    }
  }
  
  ## Use a preset matrix (for two dimension reductions e.g.)
  if (len(preset) != 0 ) {
    arguments = preset[-which(names(preset) == "name")]
    arguments$day_path = day_path
    
    if (control) {
      arguments$control <- control
    }
    
    final_mat <- t(do.call(get_reduced_mat_full_day, arguments))
#    final_mat = t(preset$preset_matrix)
    print(sprintf("Using preset matrix size of %d, %d", nrow(final_mat), ncol(final_mat)))
  } else {
    print(sprintf("Reduced matrix %s does not exist! creating!", reduced_mat_name))
    runs <- list.files(day_path, recursive=F)
    runs <- runs[which(grepl("\\.R", runs))]
    runs <- sapply(str_split(runs, ".R"), function(l){return(l[[1]])})
    
    
    ###### In case using control (Old scope mats)
    if (control) {
      matlab_mat <- readMat(sprintf("%s\\%s.mat", day_path, "dff"))
      merged_mat <- matlab_mat[[1]]
      
      if (all(chunk != 0)) {
        merged_mat <- merged_mat[,unlist(lapply(chunk, 
                                                function(chnk_i) {((chnk_i - 1) * 57600 + 1):(chnk_i * 57600)}))]
      }
      
      matrices_list = list(merged_mat)
    } 
    ####### In case using new matrices (New scope mats)
    else { 
      
      if (all(chunk != 0)) {
        runs <- sort(runs)[chunk]
      }
      
      # Load all matrices for all runs
      matrices_list <- lapply(runs, 
                              function(r) {load(sprintf("%s\\%s.R", day_path, r)); 
                                return(fmat)})
      
      
      
      are_na_matrices <- unlist(lapply(matrices_list, function(mt) {all(is.na(mt))}))
      
      if (sum(are_na_matrices) > 0) {
        print("Removing matrix at indexes:")
        print(which(are_na_matrices))
        matrices_list = matrices_list[which(!are_na_matrices)]
      }
      
      are_same_dim_1 <- unlist(lapply(matrices_list, function(mt) {dim(mt)[1]}))
      are_same_dim_2 <- unlist(lapply(matrices_list, function(mt) {dim(mt)[2]}))
      
      if (!all(are_same_dim_1 == are_same_dim_1[1])) {
        print("Incompatible dimensions!")
        return(list())
      }
    }
    
    final_mat <- mult_mat(matrices_list, 
                          window_size=window_size, 
                          norm=normed, 
                          activity_threshold = activity_threshold)
  }
  
  if ((first_q != 0) || (last_q != 0)) {
    n_frames <- ncol(final_mat)
    
    if (first_q != 0) {
      ind <- 1:(n_frames * first_q)
    } else {
      ind <- (n_frames - last_q * n_frames):n_frames
    }
    
    print(sprintf("Subsetting mat (%d) frames", len(ind)))
    final_mat <- final_mat[,ind]
  }
  
  
  if (just_original_mat) {
    print("Returning original matrix")
    reduced <- final_mat
  } else {
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
  }
  
  dir.create(sprintf("%s\\%s", day_path, matrix_subpath))
  save(reduced, file=sprintf("%s\\%s\\%s.R", day_path, matrix_subpath, reduced_mat_name))
  print(sprintf("Saving mat! %s", reduced_mat_name))
  
  return(reduced)
}


get_mat_with_preset <- function(path, preset_name, activity_threshold=0.2, oldscope=F) {
  preset <- get_preset_of_choice(preset_name)
  preset$activity_threshold <- activity_threshold
  
  if (len(preset$preset) > 0) {
    preset$preset$activity_threshold <- activity_threshold
  }
  
  preset$day_path <- path
  preset$control <-  oldscope
  
  mat <- do.call(get_reduced_mat_full_day, preset)
  
  return(mat)
}
