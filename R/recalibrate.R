#' Mass spectra recalibration
#'
#' Recalibrates the input mass spectra using function MFAssignR::Recal. Automatically selects homologous series for each sample based on function MFAssignR::RecalList.
#' @param import The result of the importMzml function
#' @param ionizationMode "neg" or "pos"
#' @param LowMassLimit Lower end of the measured mass range in Da (default 100)
#' @param HighMassLimit Upper end of the measured mass range in Da (default 1000)
#' @param CHO_assign_error Error tolerance for preliminary sum formula assignment in ppm (default 3)
#' @param SN_recal Minimum signal-to-noise ratio for peaks to be used as potential recalibrants (default 10)
#' @param RecalRangeCheck TRUE/FALSE  Require that the recalibration covers a specified mass range (set below, default TRUE)
#' @param recalRange TRUE/FALSE  Sample must be recalibrated within this entire range, else the sample is reported as having failed recalibration (default c(50, 500))
#' @param AutoUpperLimit TRUE/FALSE Should the upper limit of recalRange be set automatically to the m/z where a set percentage of the cumulative intensity is reached? (default TRUE)
#' @param AutoUpperLimitTarget Target percentage of cumulative intensity starting from the lowest detected m/z (default 0.8)
#' @param Nitrogen_recal Number of N atoms allowed in inital sum formula assignment (default 0)
#' @param max_loop Maximum number of attempts at fixing recalibration errors (default 3)
#' @param parallel Allow parallel execution? (default TRUE)
#' @return A list containing for each sample the recalibrated m/z values, the noise level, an indicator of recailbration success or failure, and the sample name.
#' @export


recalilbrate <- function(import,
                     ionizationMode,  #'neg' or 'pos'
                     LowMassLimit,  #What is the lower limit of the mass range in your scans?
                     HighMassLimit,  #What is the upper limit of the mass range in your scans?
                     CHO_assign_error = 3, #allowed error range (in ppm) for CHO sum formula assignment. Needs to be kind of high, because the data is not recalibrated at this point.
                     SN_recal = 10,  #Signal to noise ratio during import.  Low intensities commonly lead to high recalibration errors(>1ppm). Lower this value if you have trouble finding recalibrants.
                     RecalRangeCheck = TRUE, #require that the recalibration covers a specified mass range (set below)
                     recalRange = c(150,500), #Sample must be recalibrated within this entire range, else the sample is reported as failed recalibration.
                     AutoUpperLimit = TRUE, #should the upper limit of recalRange be set automatically to the m/z where a set percentage of the cumulative intensity is reached?
                     AutoUpperLimitTarget = 0.80,
                     Nitrogen_recal = 0, #Number of nitrogen atoms allowed during homologous series search. Set to zero for CHO only.
                     max_loop = 3, # maximum number of retries during recalibration
                     parallel = TRUE){
  `%>%` <- magrittr::`%>%`
  time1 <- Sys.time()
  if (parallel) {
    require(future)
    future::plan(multisession) #use plan(multicore) under linux
  } else {
    future::plan(sequential)
  }
  suppressWarnings({
  progressr::handlers(progressr::handler_progress(format="[:bar] :percent :eta :message")) #formatting of the progress bar
  progressr::with_progress({
    p <- progressr::progressor(along = import)
    recal_result <- future.apply::future_lapply(import, FUN = function(x) {
      p("recalibrating...")
      # Assign CHO sum formulas. This is needed as an
      # input for the identification of the recalibration
      # series.
      
      assignCHO <- try(MFAssignR::MFAssign_RMD(peaks = x[[1]]$Mono, 
                                    isopeaks = x[[1]]$Iso, ionMode = ionizationMode, 
                                    lowMW = LowMassLimit, highMW = HighMassLimit, 
                                    ppm_err = CHO_assign_error,
                                    SN = SN_recal * x[[2]],
                                    Nx  = Nitrogen_recal
      )$Unambig,
      silent = TRUE)  
      df_recal <- try(MFAssignR::RecalList(assignCHO),
                      silent = TRUE)
      
      #check if any homologous series were found
      if (class(df_recal) != "try-error" &&
          nrow(df_recal) > 0) {
        
        #extract their mass ranges (used later for ordering the series)
        
        df_recal2 <- df_recal %>%
          dplyr::mutate(min_mass = as.numeric(stringr::str_extract(df_recal$`Mass Range`, pattern = "^[[:digit:]]+\\.[[:digit:]]+")),
                 max_mass = as.numeric(stringr::str_extract(df_recal$`Mass Range`, pattern = "[[:digit:]]+\\.[[:digit:]]+$"))) %>%
          dplyr::arrange(min_mass)
        
        recal_series_initial <- c() #initialize
        start_mass_series_from <- c(0,25,50,100,150,200,250,300,400,500)
        start_mass_series_to <- c(25,50,100,150,200,250,300,400,500,600)
        for (i in 1:10){    #select series section-wise over the mass range, ca 50Da per section
          recal_series_initial[i] <- ifelse(any(dplyr::between(df_recal2$min_mass, LowMassLimit + start_mass_series_from[i], LowMassLimit + start_mass_series_to[i])),
                                            df_recal2 %>%
                                              dplyr::filter(dplyr::between(df_recal2$`Tall Peak`, LowMassLimit + start_mass_series_from[i], LowMassLimit + start_mass_series_to[i])) %>%
                                              dplyr::slice_min(order_by = `Series Index`, n = 1) %>%
                                              dplyr::select("Series") %>%
                                              unlist(),
                                            NA)
        }
        recal_series <- recal_series_initial[!is.na(recal_series_initial)] 
        while (length(recal_series) < 10 &
               length(recal_series) < nrow({df_recal2 %>% dplyr::filter(min_mass < 400)})){
          recal_series[length(recal_series) + 1] <- df_recal2 %>%
            dplyr::filter(!Series %in% recal_series) %>%
            dplyr::slice_min(order_by = `Series Index`, n = 1) %>%
            dplyr::pull(Series)
        }
        
        
        #initialize an empty error object for the recalibration
        recal_output <- "placeholder for object initialization"
        class(recal_output) <- "try-error"
        
        #HERE THE ACTUAL RECALIBRATION STARTS
        #we go through all combinations of the top10 homologous series from df_recal
        #we start with all ten, then reduce to 9, 8, 7...
        #on each level we go through all possible combinations starting with low-mass combinations first
        
        
        #find first stable recalibration solution
        i <- length(recal_series)
        while (i >= 1){
          recal_series_combinations <- utils::combn(recal_series, m = i)
          k <- 1
          while (k <= ncol(recal_series_combinations)) {
            recal_output <- try(
              MFAssignR::Recal(
                df = assignCHO,
                peaks = x[[1]]$Mono,
                isopeaks = x[[1]]$Iso,
                mode = ionizationMode,
                series1 = ifelse(length(recal_series_combinations[,1])>=1,recal_series_combinations[1,k],NA),
                series2 = ifelse(length(recal_series_combinations[,1])>=2,recal_series_combinations[2,k],NA),
                series3 = ifelse(length(recal_series_combinations[,1])>=3,recal_series_combinations[3,k],NA),
                series4 = ifelse(length(recal_series_combinations[,1])>=4,recal_series_combinations[4,k],NA),
                series5 = ifelse(length(recal_series_combinations[,1])>=5,recal_series_combinations[5,k],NA),
                series6 = ifelse(length(recal_series_combinations[,1])>=6,recal_series_combinations[6,k],NA),
                series7 = ifelse(length(recal_series_combinations[,1])>=7,recal_series_combinations[7,k],NA),
                series8 = ifelse(length(recal_series_combinations[,1])>=8,recal_series_combinations[8,k],NA),
                series9 = ifelse(length(recal_series_combinations[,1])>=9,recal_series_combinations[9,k],NA),
                series10 = ifelse(length(recal_series_combinations[,1])>=10,recal_series_combinations[10,k],NA),
                SN = SN_recal * x[[2]]
              )[2:4],
              silent = TRUE)
            if (class(recal_output) == "list"){ if(stringr::str_detect(x[[3]],"[[B,b]]lank")) break
              if (RecalRangeCheck & (min(recal_output$RecalList$exp_mass) <= {recalRange[1]}) & (max(recal_output$RecalList$exp_mass) >= {
                ifelse(AutoUpperLimit,
                       x[[1]]$Mono %>%
                       dplyr::mutate(cumsum = cumsum(abundance),
                              sum = sum(abundance)) %>%
                       dplyr::filter(cumsum <= AutoUpperLimitTarget * sum(abundance)) %>%
                       dplyr::pull(exp_mass) %>%
                       max(),
                       recalRange[2])
              })) break else {
                k <- k+1
                class(recal_output) <- "try-error"
              }
            } else {
              k <- k+1
            }
            if (k >= max_loop) break #maximum number of retries
          }
          if (class(recal_output) == "list") break
          if (i == 1) break
          i <- i-1
        }
        
        if (class(recal_output) == "list") {
          howdiditgo <- paste("Successfully recalibrated.")
        }
        
        if (class(recal_output) == "try-error") {
          howdiditgo <- paste("No stable recalibration solution in range.")
        }
      }
      
      if (class(df_recal) == "try-error") {
        howdiditgo <- paste("The recalibration step was skipped because no homologous series were found.")
        recal_output <- list(Mono = NA, Iso = NA, RecalList = NA)
      }
      
      if (class(recal_output) == "try-error"){
        recal_output <- "try-error" }
      
      tmp <- list(recal_output,
                  x[[2]],
                  howdiditgo,
                  x[[3]])
      return(tmp)
    }#,
    #future.chunk.size = 5  #during parallel processing, this controls how many tasks(samples) are sent to each worker at a time
    )
  })
  return(recal_result)
  time2 <- Sys.time()
  round(difftime(time2,time1),1)
  })
}
