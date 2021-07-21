#' Detection of RFD profile changes associated with DNA replication origins, based on OK-seq data
#'
#' @param Data data frame containing position-specific RFD values. The format is as following, 4 columns chr, start, end, RFD
#' @param IntSize size of the window used to calculate linear regression (default:200000)
#' @param StepSize distance between intervals used to calculate linear regression (default:1000)
#'
#' @return Table with coordinates of the RFD shifts associated with DNA replication origins
#' @export
detectPeaksOKseq = function(Data,IntSize=200000,StepSize=1000) {

      require('peakPick')

      #extract all chromosomes
      chrs= unique(Data[,1])

      #exclude chr1
      chrs = chrs[chrs != "chrM"]

      #Calculate the lm slopes for each chromosome
      FitTable = data.frame()
      for (i in 1:length(chrs)) {
        print(paste0("Calculating slopes: ",chrs[i]))
        tData = Data[Data[,1]==chrs[i],]
        Nint = round(max(tData[,2])/StepSize)+1
        tFitTable = matrix(NA,ncol=2,nrow=Nint)
        for (k in 1:Nint) {
          start = (k-1)*StepSize+1-IntSize/2
          end = k*StepSize+IntSize/2
          ptData = tData[tData[,3]<end & tData[,2]>start & !is.na(tData$RFD), ]
          tFitTable[k,1] = round((start+end)/2)
          if (dim(ptData)[1]>0) {
            model <- lm(1:dim(ptData)[1] ~ poly(ptData[,4], 1, raw=TRUE))
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


      #peak detection
      peaks = rep(FALSE,dim(FitTable)[1])
      ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}
      dist = 500
      for (chr in unique(FitTable$chr)) {
        print(paste0("Detecting peaks: ",chr))
        chridx = FitTable$chr == chr
        FitTable$slopeSmooth[chridx] = as.numeric(ma(FitTable$slope[chridx],n=ceiling(100000/dist)))  ### verify for other dist values
        peaks[chridx] <- peakpick(matrix(FitTable$slopeSmooth[chridx], ncol=1), neighlim=100, peak.npos=100,deriv.lim = 1) &  FitTable$slopeSmooth[chridx]>0  ##neglim should be calculated
      }
      FitTablePeak = FitTable
      FitTablePeak$peak=0
      FitTablePeak$peak[peaks]=FitTablePeak$pos[peaks]
      dim(FitTablePeak[FitTablePeak$peak>0,])


      return(FitTablePeak[FitTablePeak$peak>0,])
}
