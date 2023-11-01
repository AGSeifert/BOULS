
#' Segmentation of untargeted LCMS data
#'
#' @param object XCMSnExp object from xcms workflow (https://bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms.html) after retention time alignment (the SOULS step replaces the correspondence step in the xcms workflow).
#' @param RT_range Retention time range of the spectra to be processed.
#' @param size_RT_segs Size of the segments in retention time dimension (RT_range should be divisible by size_RT_segs).
#' @param mz_range Mass range of the spectra to be processed.
#' @param size_mz_segs Size of the segments in mass dimension (mz_range should be divisible by size_mz_segs)
#'
#' @return Matrix showing summed intensities in each segment.
#' @export
#'
#' @examples seg.result <- souls(object = xdata_adj, RT_range = c(300, 360), size_RT_segs = 10, mz_range = mz_range = c(250, 750), size_mz_segs = 5)


souls <- function(object,
                  num_workers,
                  RT_range,
                  size_RT_segs,
                  mz_range,
                  size_mz_segs){

  options(MulticoreParam = BiocParallel::MulticoreParam(workers = num_workers))

    seg.rt <- diff(RT_range) / size_RT_segs # seg.rt = number of RT segs
    seg.rt.list <- as.list(seq(from = RT_range[1], to = (RT_range[2] - size_RT_segs), by = size_RT_segs))
    names(seg.rt.list) <- seq(from = RT_range[1], to = (RT_range[2] - size_RT_segs), by = size_RT_segs)

    seg.mz <- diff(mz_range) / size_mz_segs # seg.mz = number of mz segs
    seg.mz.list <- as.list(seq(from = mz_range[1], to = (mz_range[2] - size_mz_segs), by = size_mz_segs))
    names(seg.mz.list) <- seq(from = mz_range[1], to = (mz_range[2] - size_mz_segs), by = size_mz_segs)
    seg.mz.list <- lapply(seg.mz.list, function(mz) {
      mz <- c(mz, mz + size_mz_segs)
    })


    seg.result <- bplapply(seg.mz.list, function(mz.seg) {
      TIC <- xcms::chromatogram(object,
                                mz = as.vector(mz.seg),
                                aggregationFun = "sum"
      )

      lapply(seg.rt.list, function(rt.seg) {
        apply(TIC, 2, function(TIC_x) {
          sum(TIC_x[[1]]@intensity[which(dplyr::between(TIC_x[[1]]@rtime, rt.seg, (rt.seg + size_RT_segs)) == TRUE)], na.rm = TRUE)
        })
      }) # End lapply seg.rt.list
    }) # End of bplapply mz


    # Creating a matrix with the results
    result.unlist <- unlist(seg.result)
    number.segs <- (diff(RT_range) / size_RT_segs) * ((mz_range[2] - mz_range[1]) / size_mz_segs)
    segdata <- matrix(result.unlist, nrow = number.segs, ncol = nrow(object), byrow = TRUE)
    colnames(segdata) <- basename(fileNames(object))

    # Rownames

    prefix <- rep(names(seg.rt.list), times = length(seg.mz.list))
    suffix <- rep(names(seg.mz.list), each = length(seg.rt.list))

    rownames(segdata) <- paste0(prefix, sep = ".", suffix)

    segdata
  }
