

paths = list("Y:\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_4\\IC47\\171121_IC47\\171121_IC47_run1\\IC47_171121_001.mat2",
         "Y:\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_4\\IC47\\171121_IC47\\171121_IC47_run2\\IC47_171121_002.mat2",
         "Y:\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_4\\IC47\\171121_IC47\\171121_IC47_run3\\IC47_171121_003.mat2",
         "Y:\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_4\\IC47\\171121_IC47\\171121_IC47_run4\\IC47_171121_004.mat2",
         "Y:\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_4\\IC47\\171121_IC47\\171121_IC47_run5\\IC47_171121_005.mat2",
         "Y:\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_4\\IC47\\171122_IC47\\171122_IC47_run1\\IC47_171122_001.mat2",
         "Y:\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_4\\IC47\\171122_IC47\\171122_IC47_run2\\IC47_171122_002.mat2",
         "Y:\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_4\\IC47\\171122_IC47\\171122_IC47_run3\\IC47_171122_003.mat2",
         "Y:\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_4\\IC47\\171122_IC47\\171122_IC47_run4\\IC47_171122_004.mat2",
         "Y:\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_4\\IC47\\171122_IC47\\171122_IC47_run5\\IC47_171122_005.mat2",
         "Y:\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_4\\IC47\\180102_IC47\\180102_IC47_run1\\IC47_180102_001.mat2",
         "Y:\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_4\\IC47\\180102_IC47\\180102_IC47_run2\\IC47_180102_002.mat2",
         "Y:\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_4\\IC47\\180102_IC47\\180102_IC47_run3\\IC47_180102_003.mat2",
         "Y:\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_4\\IC47\\180102_IC47\\180102_IC47_run4\\IC47_180102_004.mat2",
         "Y:\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_4\\IC47\\180102_IC47\\180102_IC47_run5\\IC47_180102_005.mat2",
         "Y:\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_4\\IC47\\180102_IC47\\180102_IC47_run6\\IC47_180102_006.mat2",
         "Y:\\yoavlivn\\YoavLivneh\\2p_data\\yoav\\2p_data_4\\IC47\\180102_IC47\\180102_IC47_run7\\IC47_180102_007.mat2")

preprocess_matrices <- function(path){
  
  output_path <- "C:\\Users\\97254\\Documents\\yoav_livneh_lab\\temp\\"
  mt <- readMat(path)
  cellsort <- mt$cellsort
  fmat <- c()
 
  if (dim(cellsort)[1] == 1) {
    s <- which(dimnames(mt$cellsort[1,1,1]$timecourse)[[1]] == "dff.axon")  
  } else {
    s <- which(dimnames(mt$cellsort[12,1,1]$timecourse)[[1]] == "dff.axon")  
  }
  
  for (i in 1:dim(cellsort)[3]) {
    if (dim(cellsort)[1] == 1) {
      fmat <- rbind(fmat, as.vector(cellsort[1,1,i]$timecourse[[s]]))
    } else {
      fmat <- rbind(fmat, as.vector(cellsort[12,1,i]$timecourse[[s]]))
    }
    print(sprintf("%d",i))
  }
  
  
  rm(mt)
  
  fpath <- paste(output_path, str_split(path, "\\\\")[[1]][10], sep="\\")
  final_path <- paste(strsplit(fpath, ".mat")[[1]][1], ".R", sep="")
  
  save(fmat, file=final_path[[1]])
  
  return(final_path)
}

get_path_by_prefix_vector <- function(prefix_vec) {
  return(sprintf("%s\\%s.R", output_path, prefix_vec))
}




