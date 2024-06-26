
#' Bucketing of untargeted LCMS data
#'
#' @param object XCMSnExp object from xcms workflow (https://bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms.html) after retention time alignment (the BOULS step replaces the correspondence step in the xcms workflow).
#' @param RT_range Retention time range of the spectra to be processed.
#' @param size_RT_bins Size of the buckets in retention time dimension (RT_range should be divisible by size_RT_bins).
#' @param mz_range Mass range of the spectra to be processed.
#' @param size_mz_bins Size of the buckets in mass dimension (mz_range should be divisible by size_mz_bins)
#'
#' @return Matrix showing summed intensities in each bucket.
#' @export
#'
#' @examples
#' bin.result <- bouls(object = xdata_adj,
#' RT_range = c(300, 360),
#' size_RT_bins = 10,
#' mz_range = c(250, 750),
#' size_mz_bins = 5)



bouls <- function(object,
                  num_workers,
                  RT_range,
                  size_RT_bins,
                  mz_range,
                  size_mz_bins){

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
      TIC <- xcms::chromatogram(object,
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
    bindata <- matrix(result.unlist, nrow = number.bins, ncol = nrow(object), byrow = TRUE)
    colnames(bindata) <- basename(fileNames(object))

    # Rownames

    prefix <- rep(names(bin.rt.list), times = length(bin.mz.list))
    suffix <- rep(names(bin.mz.list), each = length(bin.rt.list))

    rownames(bindata) <- paste0(prefix, sep = ".", suffix)

    bindata
  }
