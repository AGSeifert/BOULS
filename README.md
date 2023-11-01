# SOULS
## Installation
```{r setup}
devtools::install_github("AGSeifert/RFSurrogates")
library(SOULS)
```

This is the introduction to the LCMS data processing approach SOULS (Segmentation of untargeted LCMS spectra). It is based on the xcms package [1-3] and particularly suitable for the development of large machine learning models over time. In contrast to the xcms workflow, it allows separate processing, as a unique peak list independent of the processing batch is achieved by summing up the signal intensities within defined segments in the LCMS spectra. 


## Data
To demonstrate the functionality of this approach, [example data is provided here](https://www.fdr.uni-hamburg.de/record/13535).
Please download the and unzip the file. The phenodata folder contains a csv file with information about the samples, which can be customized individually (e.g. name, variety, geogr. origin, instrument, harvesting year etc). The other folders provide mzML files (we recommend converting the vendor specific files to mzML format with MSConvert [4]). The example files have been shortened to one minute to enable processing on a normal laptop. For real data sets, we recommend using a workstation. The development of large machine learning models over time requires the definition of one sample as a reference sample for the retention time alignment. For example, this could be the sample with the most detected peaks. This reference sample is added to each new data set to be processed, making the retention time alignment independent of the processing batch. The "Samples" folder contains the samples to be processed in the respective batch.

```{r}
# Please insert your paths to the 'Sample', the 'Reference' and the 'phenodata' folders.
path_mzMLs <- '/home/hansen/example_data/Samples'
path_Ref <- '/home/hansen/example_data/Reference'
path_csv <- '/home/hansen/example_data/phenodata'

```

## Data processing using SOULS
This package provides two functions. To facilitate a first try-out, the first function `process_souls()` includes all steps from data import to R to the segmentation of the spectra. In this function, the settings for the xcms functions are predefined. To adjust the xcms settings, the second function `souls()` can be included in the [xcms workflow](https://bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms.html), replacing the correspondence step. 

### General processing using the SOULS approach
Here, the example data is processed using the predefined xcms settings. Two CPUs are used in parallel (num_workers parameter). In this case, a retention time range of 300 s to 360 s and a mass range of 250 Da to 750 Da is processed. The size of the segments is 10 s in retention time dimension and 5 Da in mass dimension. The result is a matrix with summed intensities for the respective segments. The segments are named (rownames). The first value corresponds to the beginning of the segment in retention time dimension and the second value to the beginning of the segment in mass dimension. For example, the segment 300-310 s and 250-255 Da would have the name "300.250". 

```{r , warning=FALSE}
seg.result <- process_souls(path_mzMLs = path_mzMLs, 
                            path_Ref = path_Ref, 
                            path_csv = path_csv,
                            num_workers = 2,
                            RT_range = c(300, 360),
                            size_RT_segs = 10,
                            mz_range = c(250, 750),
                            size_mz_segs = 5)
head(seg.result)
```

### Individual processing using xcms and SOULS approach
For individual settings in peak picking and retention time alignment, the `souls()` function can replace the [correspondence step](https://bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms.html#25_Correspondence) in the xcms workflow. The XCMSnExp object resulting from the retention time alignment should be used as object. This is an example, how you can process your LCMS data using xcms and SOULS.

```{r , warning=FALSE}
## load mzMLs
  mzMLs <- list.files(path_mzMLs, full.names = TRUE, recursive = TRUE)

  ## load mzML file of reference sample for retention time alignment
  mzMLs_Ref <- list.files(path_Ref, full.names = TRUE, recursive = TRUE)

  ## assembling of sample and reference mzMLs (last sample)
  mzMLs <- c(mzMLs, mzMLs_Ref[1])

  ## read phenodata from csv files
  pd <- read.csv(paste0(path_csv, "/phenodata.csv"))
  ## add reference sample to the pd (last row)
  pd <- rbind(pd, c("Ref.mzML", rep("Ref",times = ncol(pd)-1)))

  # read raw data
  raw_data <- MSnbase::readMSData(
    files = mzMLs,
    pdata = new("NAnnotatedDataFrame", pd),
    mode = "onDisk"
  )

  # Peak picking settings
  cwp <- xcms::CentWaveParam(
    peakwidth = c(5, 20),
    noise = 100,
    ppm = 5,
    snthresh = 5
  )

  # Peak detection
  xdata <- xcms::findChromPeaks(raw_data, param = cwp)

  # Retention time alignment with reference sample as center sample
  index.ref <- which(pd$Name == "Ref.mzML")

  xdata_adj <- xcms::adjustRtime(xdata,
                    param = xcms::ObiwarpParam(
                      binSize = 0.1,
                      centerSample = index.ref,
                      localAlignment = TRUE,
                      distFun = "cor_opt"
                    )
  )
  
  seg.result <- souls(object = xdata_adj,
                      num_workers = 2,
                      RT_range = c(300, 360),
                      size_RT_segs = 10,
                      mz_range = c(250, 750),
                      size_mz_segs = 5)
head(seg.result)
```

[1] Smith, C.A., Want, E.J., O'Maille, G., Abagyan,R., Siuzdak, G. (2006). “XCMS: Processing mass spectrometry data for metabolite profiling using nonlinear peak alignment, matching and identification.” Analytical Chemistry, 78, 779–787.

[2] Tautenhahn R, Boettcher C, Neumann S (2008). “Highly sensitive feature detection for high resolution LC/MS.” BMC Bioinformatics, 9, 504.

[3] Benton HP, Want EJ, Ebbels TMD (2010). “Correction of mass calibration gaps in liquid chromatography-mass spectrometry metabolomics data.” BIOINFORMATICS, 26, 2488.

[4] Chambers, M.C.; Maclean, B.; Burke, R.; Amodei, D.; Ruderman, D.L.; Neumann, S.; Gatto, L.; Fischer, B.; Pratt, B.; Egertson, J.; et al. A Cross-Platform Toolkit for Mass Spectrometry and Proteomics. Nat. Biotechnol. 2012, 30, 918–920, doi:10.1038/nbt.2377.
