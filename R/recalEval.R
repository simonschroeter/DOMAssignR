#' Evaluate recalibrations
#'
#' Prints the recalibration outcome for all samples.
#' @param import The result of the import Mzml function
#' @param recal_result The result of the recalibrate function
#' @return Printout of an overview of the recalibration success or failure for each sample.
#' @export

recalEval <- function(import,
                       recal_result){
  `%>%` <- magrittr::`%>%`
  #WERE ANY UNSUCCESSFUL?
  print(ifelse(any(!sapply(recal_result, function(x){x[[3]]}) %in% c("Successfully recalibrated.","Successfully recalibrated with only one homologous series. Proceed with caution.")),
         "Some samples contain recalibration errors, please check and remove before continuing.",
         "All samples okay to conitnue."))
  #CHECK IF ALL SAMPLES COULD BE RECALIBRATED
  #for (i in 1:length(recal_result)){
  #  print(paste0("Sample",i, ": ", recal_result[[i]][3]))
  #}
  #
  if(any(!sapply(recal_result, function(x){x[[3]]}) %in% c("Successfully recalibrated.","Successfully recalibrated with only one homologous series. Proceed with caution."))){
  for (i in 1:length(recal_result)){
    if(i==1){print("The following samples could not be recalibrated:")}
    if(recal_result[[i]][[3]] %in% c("Successfully recalibrated.","Successfully recalibrated with only one homologous series. Proceed with caution.")){
      print(paste0("Sample",i, ": ", recal_result[[i]][4]))
    }}
  if(!any(!sapply(recal_result, function(x){x[[3]]}) %in% c("Successfully recalibrated.","Successfully recalibrated with only one homologous series. Proceed with caution."))) {print("none")}
  }
  print("_____________________________")
  print("Spectra info:")
  for (i in 1:length(import)){
    print(paste("Sample",i, import[[i]][[3]],"contains", import[[i]][[1]]$Mono %>% nrow(), "monoisotopic masses above noiselevel =", round(import[[i]][[2]])))
    }
    }


