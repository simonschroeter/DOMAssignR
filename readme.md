# DOMAssignR
DOMAssignR is a wrapper around the MFAssignR package.

The purpose of this package is to provide a standardized and semi-automated way of recalibrating any number of MS1 mass spectra of Dissolved Organic Matter (DOM) locally on your device. This package allows for parallelized use of functions provided in package MFAssignR. Requires averaged and centroided MS1 spectra in .mzml format. Example data is available from: [https://doi.org/10.17617/3.OZVHVP](https://doi.org/10.17617/3.OZVHVP).

## Installation
You can install DOMAssignR from Github with:
```R
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("simonschroeter/DOMAssignR")
```
DOMAssignR requires packages MFAssignR, MALDIquant, future, future.apply, progressr, and the tidyverse.

## Usage
For testing purposes, download some example data from [https://doi.org/10.17617/3.OZVHVP](https://doi.org/10.17617/3.OZVHVP). Import the spectra with importMzml(). Set the signal to noise ratio (SN) to a value of your choice depending on the noisiness of your spectra. Peaks below the SN are not imported. Therefore, the choice strongly affects expected runtimes. Values around SN=5 give reasonable runtimes on the example data set.Recalibrate the imported spectra using recal(). Specify if positive or negative ionization mode was used and fix the m/z range to be evaluated, default is 100 Da to 1000 Da. The recal() function automates the recalibartion workflow of MFAssignR. Briefly, sum formulas containing only elements C, H, and O are assigned. Based on these, homologous series of molecules are extracted. These are used as recalibrants. recal() attempts to find up to ten homologous series, which together cover most of the m/z range of the respective samples (option RecalRangeCheck). The expected m/z range to be covered by the recalibrants can either be defined manually or based on a cumulative intensity threshold of e.g. 80%, summing from low to high m/z (option AutoUpperLimit). The success of the recal() function for each mass spectrum can be checked with recalEval().

The recalibration result file needs to be cleaned with excludeRecalErrors() in case any samples could not be successfully reaclibrated.

Further options for inspecting the recalibration quality are recalRange() to show the m/z ranges covered by the recalibrants and recalPlot(), which plots the mass spectrum, its recalibrants, and their respective mass error. In general, most of the m/z range should display errors < 1 ppm. In low intensity regions, such as high m/z ranges, this criterion is often not fulfilled. Careful inspection of the recalibration quality is advised and finetuning of the recalibration parameters is often required. Detailed information on all recalibration parameters is provided in the documentation of MFAssignR.

Before proceeding to sum formula assignment, all recalibrated spectra are joined with junction(). The default tolerance for m/z values of different spectra to be considered equal is 0.5 ppm.

Sum formula assignment with assignSumFormulas() works equal to MFAssign(). The proportion of intensity remaining among the assigned sum formulas compared to the raw data can be assessed with intensityRatio() and commonly ranges between 60 and 80%.

## Acknowledgments
This package was developed at the Max Planck Institute for Biogeochemistry, Jena, Germany within the Collaborative Research Centre AquaDiva of the Friedrich Schiller University Jena, funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation)—SFB 1076—Project Number 218627073.

## Disclaimer
Functionality on different machines is still being evaluated. Any feedback is highly appreciated.
