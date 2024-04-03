#' Junction of m/z tables
#'
#' Joins together the m/z tables of all recalibrated mass spectra.
#' @param recal_result The result of the recalibrate function
#' @param tolerance Tolerance in joining close m/z values (in ppm, default 0.5)
#' @return A dataframe containing one column of m/z values followed by their respective intensities in all spectra.
#' @export


junction <- function(recal_result,
                     tolerance = 0.5) {

`%>%` <- magrittr::`%>%`

junction_input <-  dplyr::bind_rows(
  lapply(recal_result, function(x){
    colnames(x[[1]]$Mono)[1:2] <- c("mz", "I")
    colnames(x[[1]]$Iso)[1:2] <- c("mz", "I")
    rbind(
      x[[1]]$Mono[1:2],
      x[[1]]$Iso[1:2]) %>%
      as.data.frame() %>%
      dplyr::arrange(`mz`) %>%
      dplyr::mutate(`MDL_2.5` = 10,
             `ResPow` = 480000,
             index = fs::path_ext_remove(basename(x[[4]])))}))

mdlname<-names(junction_input)[3]
names(junction_input)[3]<-"MDL"

junction_output <- dplyr::arrange(junction_input,mz) %>%
  # Start with the smallest mass and check if the next mass is the tolerance range defined by 'tolp'. If Yes group togetherinto a mass cluster, if not start a new group (mass cluster). 
  dplyr::mutate(mz1=mz,group = cumsum( abs(((mz - dplyr::lag(mz, default = mz[1], order_by = mz) )/dplyr::lag(mz, default = mz[1], order_by = mz))*10^6) > tolerance)) %>%
  # Group by found mass cluster and check for duplicate samples within a mass cluster
  dplyr::group_by(group) %>%
  dplyr::mutate(max=abs(mz[which.max(I)]-mz)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(group, index) %>% dplyr::arrange(max)%>% dplyr::mutate(gjk=as.character(ifelse(dplyr::row_number()==1,paste0(group,"x",dplyr::row_number()),paste0(group,"x","x",dplyr::row_number()))))%>%            #dplyr::summarize(mz=(mz[which.min(max)]),I=(I[which.min(max)]),MDL=(MDL[which.min(max)]),ResPow=(ResPow[which.min(max)]),mz1=mz1[which.min(max)],refe=refe[which.min(max)],m1=m1[which.min(max)])%>%
  dplyr::ungroup() %>% 
  dplyr::select(c(-max,-group))%>%
  dplyr::group_by(gjk) %>%
  dplyr::mutate(
    # calculate weighted average mass within each mass group or use arithmethic mean with mean()
    mz = weighted.mean(mz, sqrt(I), na.rm=T),MDL=mean(MDL), ResPow=mean(ResPow),SE=sd(mz1)/mz*1e6,mz1=mean(mz1)) %>%  
  # Take out the  group assignment and columns not needed
  dplyr::ungroup() %>%
  dplyr::select(c(-mz1,-gjk)) %>%
  # Add the "B" leader to the index
  dplyr::mutate(index = paste0("Sample", index)) %>%
  #intensity filter
  #dplyr::filter(I > (SN_import * noise_abs_intensity)) %>%
  #remove singlets
  dplyr::group_by(mz) %>%
  dplyr::mutate(mzCount = dplyr::n()) %>%
  dplyr::filter(mzCount > 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(-mzCount) %>%
  # Put intensity values of each sample into separate columns
  tidyr::pivot_wider(names_from = index, values_from = I, values_fn = max) %>%
  #remove duplicate masses (weird error)
  dplyr::filter(!duplicated(mz)) %>%
  dplyr::arrange(mz)

duplicates5digits <- junction_output$mz[duplicated(round(junction_output$mz,5))]

junction_output <- junction_output %>%
  dplyr::filter(!mz %in% duplicates5digits) %>%
  dplyr::mutate(mz = round(mz,5))

return(junction_output)
}
