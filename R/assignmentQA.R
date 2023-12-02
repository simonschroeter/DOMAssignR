#' Sum formula assignment QA
#'
#' Computes the ratio of the summed intensity of peaks with sum formula assignment to the raw data
#' @param assignment_output The result of the assignSumFormulas function
#' @param junction_output The result of the junction function
#' @return A data frame containing the intensity ratios of sum formulas relative to the raw data.
#' @export


assignmentQA <- function(assignment_output,
                         junction_output) {
  `%>%` <- magrittr::`%>%`

  QC_sampleNames <- colnames(assignment_output)[stringr::str_detect(colnames(assignment_output),"^Sample") & !stringr::str_detect(colnames(assignment_output),"[[B,b]]lank")] #select a sample name

  intensityRatios <- vector(length = length(QC_sampleNames))

for (sample in 1:length(QC_sampleNames)){
  intensityRatios[sample] <-  sum(assignment_output %>% dplyr::pull(QC_sampleNames[sample]), na.rm = T) / sum(junction_output %>% dplyr::pull(QC_sampleNames[sample]), na.rm = T) * 100
}
intensityRatios_df <- data.frame(sample = QC_sampleNames,
                                 intensityRatio = intensityRatios)

print(paste0("Average post-assignment intensity ratio: ", round(mean(intensityRatios)), "±", round(stats::sd(intensityRatios)), "%"))

return(intensityRatios_df)
}