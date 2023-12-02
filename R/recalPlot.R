#' Recalibration plots
#'
#' Plots a recalibrated mass spectrum, highlighting the peaks used as recalibrants and their respective mass error in ppm.
#' @param recal_result The result of the recalibrate function
#' @param sample_nr The index of the spectrum to be plotted.
#' @return A plot of a mass spectrum and recalibration errors.
#' @export

recalPlot <- function(recal_result, sample_nr) {

#auto upper limit recal target
#recal_target_mz <- ifelse(AutoUpperLimit,
#                          as.numeric(recal_result[[sample_nr]][[1]]$Mono$exp_mass[min(which(cumsum(recal_result[[sample_nr]][[1]]$Mono$abundance) / sum(recal_result[[sample_nr]][[1]]$Mono$abundance) > AutoUpperLimitTarget))]),
#                          recalRange[2])

#OPTIONAL MANUAL CHECKS
#DETAILED QUALITY ASSESSMENT WITH RECALIBRATION MASS ERROR PLOTS
#resulting mass_error should be less than 1ppm

ggplot2::ggplot(data.frame(mass = recal_result[[sample_nr]][[1]]$RecalList$exp_mass,
                  intensity = recal_result[[sample_nr]][[1]]$RecalList$abundance), ggplot2::aes(x=mass, y = intensity))+
  
  ggplot2::geom_col(data = data.frame(mass = recal_result[[sample_nr]][[1]]$Mono$exp_mass,
                             intensity = recal_result[[sample_nr]][[1]]$Mono$abundance), ggplot2::aes(x=mass, y = intensity),
           color = "red")+
  ggplot2::geom_col(color="black")+
  #ggplot2::geom_vline(xintercept = recal_target_mz, color = "blue")+
  ggplot2::ggtitle(recal_result[[sample_nr]][[4]])+
  ggplot2::scale_x_continuous(limits = c(100,1000),
                     breaks = seq(100,1000, by = 100))+
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank())+
  ggplot2::ggplot(data.frame(mass_error = recal_result[[sample_nr]][[1]]$RecalList$recal_err,
                    mass = recal_result[[sample_nr]][[1]]$RecalList$exp_mass), ggplot2::aes(x=mass, y = mass_error))+
  ggplot2::geom_point(shape = 1)+
  ggplot2::geom_smooth(color="red")+
  ggplot2::scale_x_continuous(limits = c(100,1000),
                     breaks = seq(100,1000, by = 100))+
  ggplot2::scale_y_continuous(name="m/z error [ppm]")+
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())+
  patchwork::plot_layout(ncol = 1)

}