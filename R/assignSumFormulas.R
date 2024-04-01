#' Sum Formula assignment
#'
#' Assigns sum formulas to m/z values using MFAssignR. Refer to MFAssignR::MFAssign for detailed instructions and help.
#' @param junction_output The result of the junction function
#' @param ionizationMode "neg" or "pos"
#' @param LowMassLimit Lower end of the measured mass range in Da (default 100)
#' @param HighMassLimit Upper end of the measured mass range in Da (default 1000)
#' @param assign_error Error tolerance for sum formula assignment in ppm (default 0.7)
#' @param set_Nx Number of Nitrogen atoms allowed (default 4)
#' @param set_Sx Number of Sulphur atoms allowed (default 1)
#' @param set_Px Number of Phosphorous atoms allowed (default 1)
#' @param set_S34x Number of Sulphur-34 atoms allowed (default 0)
#' @param set_N15x Number of Nitrogen-15 atoms allowed (default 0)
#' @param set_Dx Number of Deuterium atoms allowed (default 0)
#' @param set_Ex Number of Carbon-13 atoms allowed (default 0)
#' @param set_Clx Number of Chlorine atoms allowed (default 0)
#' @param set_Cl37x Number of Chlorine-37 atoms allowed (default 0)
#' @param set_Brx Number of Bromine atoms allowed (default 0)
#' @param set_Br81x Number of Bromine-81 atoms allowed (default 0)
#' @param set_Ix Number of Iodine atoms allowed (default 0)
#' @param set_Mx Number of Sodium adducts allowed (default 0)
#' @param set_NH4x Number of Ammonium adducts allowed (default 0)
#' @param set_Zx Amount of charge allowed (default 1)
#' @return A data frame containing the assigned sum formulas and their intensities per sample
#' @export


assignSumFormulas <- function(junction_output,
                               ionizationMode,  #'neg' or 'pos'
                               LowMassLimit = 100,  
                               HighMassLimit = 1000,
                               assign_error = 0.7,
                              
                              #ELEMENT SETTINGS
                               set_Nx = 4,
                               set_Sx = 1,
                               set_Px = 1,
                               set_S34x = 0,
                               set_N15x = 0,
                               set_Dx = 0,
                               set_Ex = 0,
                               set_Clx = 0,
                               set_Cl37x = 0, 
                               set_Fx = 0,
                               set_Brx = 0, 
                               set_Br81x = 0, 
                               set_Ix = 0,
                               set_Mx = 0,
                               set_NH4x = 0,
                               set_Zx = 1 ,
                               set_Sval = 2,
                               set_Nval = 3,
                               set_S34val = 2,
                               set_N15val = 3,
                               set_Pval = 5,
                               set_Ox = 30,
                               set_O_Cmin = 0,
                               set_O_Cmax = 2.5,
                               set_H_Cmin = 0.3,
                               set_H_Cmax = 3,
                               set_DBEOmin = -13,
                               set_DBEOmax = 13,
                               set_Omin = 0,
                               
                               # DETAIL SETTINGS####
                               set_POEx = 0,
                               set_NOEx = 0,
                               set_max_def = 0.9,
                               set_min_def = 0.5,
                               set_HetCut = "off",
                               set_NMScut = "on",
                               set_DeNovo = 300,
                               set_nLoop = 5,
                               set_SulfCheck = "on",
                               set_Ambig = "off",
                               set_MSMS = "off",
                               set_S34_abund = 30,
                               set_C13_abund = 60,
                               set_N3corr = "on"){
  
  `%>%` <- magrittr::`%>%`
  
assignMonoIso <- MFAssignR::IsoFiltR(peaks = data.frame(mz=junction_output$mz,
                                                        abundance = junction_output %>% dplyr::select(dplyr::starts_with(c("mz","Sample"))) %>% tidyr::pivot_longer(cols = dplyr::starts_with("Sample"), values_drop_na = TRUE) %>% dplyr::group_by(mz) %>% dplyr::summarize(abundance_max = max(value)) %>% dplyr::pull(abundance_max)))

assigned <- MFAssignR::MFAssign(peaks = assignMonoIso$Mono, isopeaks = assignMonoIso$Iso, 
                     ionMode = ionizationMode, lowMW = LowMassLimit, 
                     highMW = HighMassLimit, 
                     ppm_err = assign_error,
                     POEx = set_POEx,
                     NOEx = set_NOEx,
                     Nx = set_Nx,
                     Sx = set_Sx,
                     Px = set_Px,
                     S34x = set_S34x,
                     N15x = set_N15x,
                     Dx = set_Dx,
                     Ex = set_Ex,
                     Clx = set_Clx,
                     Cl37x = set_Cl37x,
                     Fx = set_Fx,
                     Brx = set_Brx,
                     Br81x = set_Br81x,
                     Ix = set_Ix,
                     Mx = set_Mx,
                     NH4x = set_NH4x,
                     Zx = set_Zx,
                     Sval = set_Sval,
                     Nval = set_Nval,
                     S34val = set_S34val,
                     N15val = set_N15val,
                     Pval = set_Pval,
                     Ox = set_Ox,
                     O_Cmin = set_O_Cmin,
                     O_Cmax = set_O_Cmax,
                     H_Cmin = set_H_Cmin,
                     H_Cmax = set_H_Cmax,
                     DBEOmin = set_DBEOmin,
                     DBEOmax = set_DBEOmax,
                     Omin = set_Omin,
                     max_def = set_max_def,
                     min_def = set_min_def,
                     HetCut = set_HetCut,
                     NMScut = set_NMScut,
                     DeNovo = set_DeNovo,
                     nLoop = set_nLoop,
                     SulfCheck = set_SulfCheck,
                     Ambig = set_Ambig,
                     MSMS = set_MSMS,
                     S34_abund = set_S34_abund,
                     C13_abund = set_C13_abund,
                     N3corr = set_N3corr)


if(any(duplicated(join$formula))){
  join <-
    join %>%
    dplyr::filter(!formula %in% unique(join$formula[duplicated(join$formula)]))
}

join <- dplyr::inner_join(x=assigned$Unambig, y=junction_output, by = c("exp_mass" = "mz")) %>%
  dplyr::select(-abundance)


return(join)
}
