#' Processing of raw LCMS data using BOULS (Bucketing of untargeted LCMS spectra)
#'
#' @param path_mzMLs Path to the folder where the mzML files of the samples are located.
#' @param path_Ref Path to the folder containing the reference sample.
#' @param index_ref Index of the reference sample in the folder; Default: 1
#' @param path_csv Path to the csv file containing the pheno-data such as name, instrument, origin, variety...
#' @param num_workers Number of cores used for parallel processing. Recommendation: 6/96 or 2/8 Cores; Default: 2
#' @param RT_range Retention time range of the LCMS spectra; Default: c(300, 360)
#' @param size_RT_bins Size of the buckets at retention time level; Default: 10
#' @param mz_range m/z range of the LCMS spectra; Default: c(250, 750)
#' @param size_mz_bins Size of the buckets at m/z level; Default: 5
#' @return Matrix showing summed intensities in each bucket.
#' @export
#'
#' @examples
#' bin.result <- process_bouls(path_mzMLs = "path/to/folder/Samples",
#' path_Ref = "path/to/folder/Ref",
#' path_csv = "path/to/folder/phenodata")



process_bouls <- function(path_mzMLs,
                  path_Ref,
                  index_ref = 1,
                  path_csv,
                  num_workers,
                  RT_range,
                  size_RT_bins,
                  mz_range,
                  size_mz_bins
                  ){

  ## load mzMLs
  mzMLs <- list.files(path_mzMLs, full.names = TRUE, recursive = TRUE)

  ## load mzML file of reference sample for retention time alignment
  mzMLs_Ref <- list.files(path_Ref, full.names = TRUE, recursive = TRUE)

  ## assembling of sample and reference mzMLs (last sample)
  mzMLs <- c(mzMLs, mzMLs_Ref[index_ref])

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


  # Set up parallel processing/amount of cores used
  options(MulticoreParam = BiocParallel::MulticoreParam(workers = num_workers))

  bin.rt <- diff(RT_range) / size_RT_bins # bin.rt = number of RT bins
  bin.rt.list <- as.list(seq(from = RT_range[1], to = (RT_range[2] - size_RT_bins), by = size_RT_bins))
  names(bin.rt.list) <- seq(from = RT_range[1], to = (RT_range[2] - size_RT_bins), by = size_RT_bins)

  bin.mz <- diff(mz_range) / size_mz_bins # bin.mz = number of mz bins
  bin.mz.list <- as.list(seq(from = mz_range[1], to = (mz_range[2] - size_mz_bins), by = size_mz_bins))
  names(bin.mz.list) <- seq(from = mz_range[1], to = (mz_range[2] - size_mz_bins), by = size_mz_bins)
  bin.mz.list <- lapply(bin.mz.list, function(mz) {
    mz <- c(mz, mz + size_mz_bins)
  })


  bin.result <- bplapply(bin.mz.list, function(mz.bin) {
    TIC <- xcms::chromatogram(xdata_adj,
                              mz = as.vector(mz.bin),
                              aggregationFun = "sum"
    )

    lapply(bin.rt.list, function(rt.bin) {
      apply(TIC, 2, function(TIC_x) {
        sum(TIC_x[[1]]@intensity[which(dplyr::between(TIC_x[[1]]@rtime, rt.bin, (rt.bin + size_RT_bins)) == TRUE)], na.rm = TRUE)
      })
    }) # End lapply bin.rt.list
  }) # End of bplapply mz


  # Creating a matrix with the results
  result.unlist <- unlist(bin.result)
  number.bins <- (diff(RT_range) / size_RT_bins) * ((mz_range[2] - mz_range[1]) / size_mz_bins)
  bindata <- matrix(result.unlist, nrow = number.bins, ncol = nrow(pd), byrow = TRUE)
  colnames(bindata) <- basename(fileNames(xdata_adj))

  # Rownames

  prefix <- rep(names(bin.rt.list), times = length(bin.mz.list))
  suffix <- rep(names(bin.mz.list), each = length(bin.rt.list))

  rownames(bindata) <- paste0(prefix, sep = ".", suffix)

  bindata

}
