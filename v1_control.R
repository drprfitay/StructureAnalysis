

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

  if (in_matlab) {
    matlab_mat <- readMat(sprintf("%s\\%s.mat", day_path, "dff"))
    merged_mat <- matlab_mat[[1]]
    
    if (all(chunk != 0)) {
      merged_mat <- merged_mat[,unlist(lapply(chunk, function(chnk_i) {((chnk_i - 1) * 57600 + 1):(chnk_i * 57600)}))]
    }
    
  }   else {
    
  }
  

  
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
