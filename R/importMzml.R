#' mzml file import
#'
#' Read centroided spectra from .mzml files
#' @param filePath The path to the folder containing the .mzml files
#' @param SN Signal to noise ratio (noise level will be determined automatically)
#' @param parallel Parallel processing (TRUE or FALSE, default: TRUE)
#' @return A list containing for each sample the spectra, noise level, and sample name
#' @export

importMzml <- 
function(filePath,
         SN,
         parallel = TRUE){
  `%>%` <- magrittr::`%>%`
  if (parallel) {
    require(future)
    future::plan(multisession) #use plan(multicore) under linux
  } else {
    future::plan(sequential)
  }
  filePaths <- list.files(filePath, full.names = T, pattern = '.mzML$')
  SN_import <- SN
  time1 <- Sys.time()
  progressr::handlers(progressr::handler_progress(format="[:bar] :percent :eta :message")) #formatting of the progress bar
  progressr::with_progress({
    p <- progressr::progressor(along = 1:length(filePaths))
    import <- future.apply::future_lapply(1:length(filePaths), FUN = function(x) {
      p("importing spectra...")
      
      # load an averaged and peak picked mzml file
      mzint <- as.data.frame(
        MALDIquant::as.matrix(
          MALDIquantForeign::importMzMl(filePaths[x],
                                             centroided = TRUE)[[1]]
          )
        ) %>%
        suppressWarnings()
      noise_sample_specific <- MFAssignR::KMDNoise(mzint)$Noise
      #detect (mono)isotope peaks
      MonoIso <- try(MFAssignR::IsoFiltR(mzint,
                                         SN = noise_sample_specific  * SN_import),
                     silent = TRUE)
      
      # if no isotope peaks are detected, we set them to
      # "none" and continue without them   
      peakList <- list()
      if (class(MonoIso) == "try-error"){
        MonoIso <- list()
        MonoIso$Mono <- mzint
        MonoIso$Iso <- "none" 
      }
      
      tmp <- list(data = MonoIso,
                  noiselevel = noise_sample_specific,
                  samplename = fs::path_ext_remove(basename(filePaths[x])))
      return(tmp)
    }#,
    #future.chunk.size = 5  #during parallel processing, this controls how many tasks(samples) are sent to each worker at a time
    )
  })
  time2 <- Sys.time()
  print(round(difftime(time2,time1, units = "mins"),1))
  return(import)
}
