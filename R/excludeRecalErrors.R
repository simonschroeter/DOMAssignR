#' Exclude recalibration errors
#'
#' Removes spectra that couldn't be recalibrated from the dataset
#' @param recal_result The result of the recalibrate function
#' @return A list containing for each sample the recalibrated m/z values, the noise level, an indicator of recailbration success or failure, and the sample name.
#' @export

excludeRecalErrors <- function(recal_result) {
  recal_result <- recal_result[which(sapply(recal_result, function(x){x[[3]]}) %in% c("Successfully recalibrated.","Successfully recalibrated with only one homologous series. Proceed with caution."))]
  return(recal_result)
}
