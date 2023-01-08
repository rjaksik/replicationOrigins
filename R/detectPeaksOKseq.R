#' Detection of replication origins based on OK-seq data
#'
#' @description
#' Detection of replication origins based on OK-seq data using RFD profile changes. For detailes on how to calculate RFD profiles see: Petryk, N., Kahli, M., d'Aubenton-Carafa, Y. et al. Replication landscape of the human genome. Nat Commun 7, 10208 (2016). https://doi.org/10.1038/ncomms10208
#'
#' @param Data Data frame containing position-specific RFD values. Four columns are expected: chr, start, end, RFD. Column names can be different but the order needs to preserved.
#' @param IntSize Size of the window used to calculate linear regression.
#' @param StepSize Distance between intervals used to calculate linear regression.
#' @param MovingAvgPts Number of neighboring points used to calculate the moving average for the purpose of linear regression slope smoothing.
#' @param Verbose Show chromosome progress details.
#'
#' @return Table with coordinates of the RFD shifts associated with DNA replication origins
#' @export
#'
#' @examples
#' rfd = readRDS(system.file("testdata", "rfd_profile.RDS", package = "replicationOrigins"))
#' rfd_ori_pos = detectPeaksOKseq(rfd, Verbose = FALSE)
detectPeaksOKseq = function(Data,
                            IntSize = 200000,
                            StepSize = 1000,
                            MovingAvgPts = 200,
                            Verbose = TRUE) {


      #extract all chromosomes
      chrs= unique(Data[,1])

      #Calculate the lm slopes for each chromosome
      FitTable = data.frame()
      for (i in 1:length(chrs)) {
        if (Verbose) print(paste0("Calculating slopes: ",chrs[i]))
        tData = Data[Data[,1]==chrs[i],]
        Nint = round(max(tData[,2])/StepSize)+1
        tFitTable = matrix(NA,ncol=2,nrow=Nint)
        for (k in 1:Nint) {
          start = (k-1)*StepSize+1-IntSize/2
          end = k*StepSize+IntSize/2
          ptData = tData[tData[,3]<end & tData[,2]>start & !is.na(tData[,4]), ]
          tFitTable[k,1] = round((start+end)/2)
          if (dim(ptData)[1]>0) {
            model <- stats::lm(1:dim(ptData)[1] ~ poly(ptData[,4], 1, raw=TRUE))
            tFitTable[k,2] = model$coefficients[2][[1]]
          } else {
            tFitTable[k,2] = NA
          }
        }
        tFitTable = data.frame(tFitTable)
        colnames(tFitTable) = c('pos','slope')
        tFitTable$chr = chrs[i]
        FitTable = rbind(FitTable,tFitTable)
      }
      FitTable = FitTable[!is.na(FitTable$slope),]


      #moving average
      ma <- function(x, n = 5){stats::filter(x, rep(1 / n, n), sides = 2)}

      #peak detection
      peaks = rep(FALSE,dim(FitTable)[1])
      for (chr in unique(FitTable$chr)) {
        if (Verbose) print(paste0("Detecting peaks: ",chr))
        chridx = FitTable$chr == chr
        FitTable$slopeSmooth[chridx] = as.numeric(ma(FitTable$slope[chridx], n=MovingAvgPts))
        peaks[chridx] <- peakPick::peakpick(matrix(FitTable$slopeSmooth[chridx], ncol=1), neighlim=100, peak.npos=100, deriv.lim = 1) &  FitTable$slopeSmooth[chridx]>0
      }
      FitTablePeak = FitTable
      FitTablePeak$peak=0
      FitTablePeak$peak[peaks]=FitTablePeak$pos[peaks]

      return(FitTablePeak[FitTablePeak$peak>0,])
}
