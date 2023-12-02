#' Recalibration ranges
#'
#' Displays the m/z of the lightest and heaviest masses used as recalibrants.
#' @param recal_result The result of the recalibrate function
#' @return A dataframe.
#' @export

recalRange <- function(recal_result) {
  `%>%` <- magrittr::`%>%`
minmz <- data.frame(sampleName = character(),
                    sampleNr = numeric(),
                    min_recal_mz = numeric(),
                    max_recal_mz = numeric())

for (i in 1:length(recal_result)){
  #print(paste0("SampleNr", i,": recalibrated from ",recal_result[[i]][[1]]$RecalList$exp_mass %>% min() %>% round()," to ",
  #             recal_result[[i]][[1]]$RecalList$exp_mass %>% max() %>% round()," m/z"))
  minmz[i,1] <- recal_result[[i]][[4]]
  minmz[i,2] <- i
  minmz[i,3] <- recal_result[[i]][[1]]$RecalList$exp_mass %>% min()
  minmz[i,4] <- recal_result[[i]][[1]]$RecalList$exp_mass %>% max()}

return(minmz)
}