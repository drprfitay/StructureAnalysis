

hummus_preset <- list(just_original_mat=T,
                      window_size=15,
                      preset=list())

skippy_preset <- list(type="lem",
                      knn1=0.075,
                      knn2=0,
                      ndim=6,
                      preset=list(name="Skippy", 
                                  type = "lem", 
                                  ndim = 20, 
                                  knn1=0.25, 
                                  knn2=0,
                                  window_size = 15))

jiffie_preset <- list(type="lem",
                      knn1=0.15,
                      knn2=0,
                      ndim=8,
                      preset=list(name="Jiffie", 
                                  type = "lem", 
                                  ndim = 20, 
                                  knn1=0.25, 
                                  knn2=0,
                                  window_size = 15))

bnd_preset <- list(type="lem",
                   knn1=0.1,
                   knn2=0,
                   ndim=8,
                   preset=list(name="BnD", 
                               type = "lem", 
                               ndim = 20, 
                               knn1=0.25, 
                               knn2=0,
                               window_size = 15))


shufersal_preset <- list(type="lem",
                         knn1=0.1,
                         knn2=0,
                         ndim=6,
                         preset=list(name="Shufersal", 
                                     type = "lem", 
                                     ndim = 20, 
                                     knn1=0.25, 
                                     knn2=0,
                                     window_size = 15))
        

nutella_preset <- list(type="dmaps",
                       ndim=8,
                       knn1=0,
                       knn2=0,
                       preset=list(name="Nutella", 
                                   type = "dmaps", 
                                   knn1=0,
                                   knn2=0,
                                   ndim = 20, 
                                   window_size = 15))


milka_choc_preset <- list(type="dmaps",
                          ndim=6,
                          knn1=0,
                          knn2=0,
                          preset=list(name="MilkaChoc", 
                                      type = "dmaps", 
                                      ndim = 20,
                                      knn1=0,
                                      knn2=0,
                                      window_size = 15))


hashahar_preset <- list(type="hlle",
                        ndim=8,
                        knn1=0.05,
                        knn2=0,
                        preset=list(name="HaShahar", 
                                    type = "dmaps", 
                                    ndim = 20, 
                                    knn1=0,
                                    knn2=0,
                                    window_size = 15))



preset_list=list("hashahar"=hashahar_preset,
                 "milka"=milka_choc_preset,
                 "milka_choc"=milka_choc_preset,
                 "nutella"=nutella_preset,
                 "shufersal"=shufersal_preset,
                 "bnd"=bnd_preset,
                 "jiffie"=jiffie_preset,
                 "jiffy"=jiffie_preset,
                 "skippy"=skippy_preset,
                 "hummus"=hummus_preset,
                 "humus"=hummus_preset,
                 "original"==hummus_preset)

get_preset_of_choice <- function(preset_name) {
  lowered_name <- as.character(tolower(preset_name))
  return(preset_list[[lowered_name]])
}
