#' Detection of replication origins based on mutation patterns associated with POLE-exo mutants
#'
#' @param MUTGR GenomicRanges object with somatic SNV positions, the object should additionally contain an "mutation" column specifying mutation type including its context e.g. TCT->TAT (this format is not obligatory and can be different as long as its the same for sel_mutations_* variables)
#' @param sel_mutations_norm vector of mutations found upstream of the replication origins (default: TCT->TAT, TCG->TTG)
#' @param sel_mutations_revcomp vector of mutations found downstream of the replication origins (default: AGA->ATA, CGA->CAA)
#' @param win window size used to calculate the PMA score (default: 200000)
#' @param dist distance between regions for which to calculate the PMA score (default: 1000)
#' @param PeakCut minimum PMA score required to consider a peak an origin (default: 0.1)
#' @param PvalCut adjusted p-value cutoff used to filter the ORI set using Fisher's exact test (default: 0.01)
#' @param useChrs perform the algorithm only for a specific set of chromosomes (default: NULL - use all basic chromosomes only)
#' @param refGenome version of the reference genome supported values are hg19 and hg38 (default: hg19)
#'
#' @return Returns a list of 2 variables PeaksGR and PMAscores. PeaksGR is a GRanges object with positions of all identified peaks and corresponding PMA score and adjusted p-value of the Fisher's exact test. PMAscores is a data.frame object with detailed statistics obtained for each of the genomic regions, including the number of mutations from each group used to calculate PMA score, smoothed PMA score (used in peak detection), positions of peaks and adjusted p-values of the Fisher's exact test.
#' @export
#'
pmaORIdetection = function(MUTGR, sel_mutations_norm = c('TCT->TAT','TCG->TTG'), sel_mutations_revcomp = c('AGA->ATA','CGA->CAA'),
                         win = 200000, dist = 1000, PeakCut = 0.1, PvalCut = 0.01, useChrs=NULL, refGenome = 'hg19') {

  require('peakPick')

  sel_mutations = c(sel_mutations_norm,sel_mutations_revcomp)
  win2 = floor(win/2)

  #check if the selected reference genome is supported
  if(!GenomeVer %in% c("hg19","hg38")) {
    stop('Only hg19 and hg38 reference genomes are supported')
  }

  #get the gap lengths
  con <- gzcon(url(paste0('http://hgdownload.cse.ucsc.edu/goldenpath/',GenomeVer,'/database/gap.txt.gz')))
  txt <- readLines(con)
  Gaps <- read.table(textConnection(txt),header=F)
  GapsGR = GRanges(Gaps[,2],IRanges(Gaps[,3]+1,Gaps[,4]),'+')  ## +1 - bed format

  #get the chromosome lengths
  con <- gzcon(url(paste0('http://hgdownload.cse.ucsc.edu/goldenpath/',GenomeVer,'/database/chromInfo.txt.gz')))
  txt <- readLines(con)
  tChrLengths <- read.table(textConnection(txt),header=F)
  #create a matrix object for safe data extraction
  ChrLengths = matrix(tChrLengths$V2,ncol=1,dimnames=list(tChrLengths$V1,'ChrLen'))

  #extract chromosome IDs to be used
  if (is.null(useChrs)) {
    chrs <- c(paste0("chr", 1:22),"chrX","chrY")
  } else {
    chrs = useChrs
  }
  missingChr = chrs[!chrs %in% rownames(ChrLengths)]

  #check if the reference contains all selected chromosomes and extract their lengths
  if (length(missingChr)>0) {
    cMissingChr = paste0(missingChr,collapse = ', ')
    stop(paste0('The following  chromosomes are missing in the selected ',GenomeVer,' reference genome: ',cMissingChr))
  }

  #calculate the PMA score for each region in specific chromosome
  ResMutFreq=data.frame()
  for (chr in chrs) {

    #determine the positions to be checked
    start=1
    end = ChrLengths[chr,'ChrLen']
    NclPositions = seq(start+win2,end,dist)
    N = length(NclPositions)
    MUTGR_chr = MUTGR[seqnames(MUTGR)==chr & MUTGR$mutation %in% sel_mutations]

    #create genomic ranges objects out of selected positions
    NclPositionsGR_up   = GRanges(chr,IRanges(NclPositions-win2 ,NclPositions-1),strand='+')
    NclPositionsGR_down = GRanges(chr,IRanges(NclPositions ,NclPositions+win2-1),strand='+')

    #check the number of specific mutations in the selected posisitons
    MutUp=MutDown=data.frame(pos=NclPositions)
    for (muttype in sel_mutations) {
      OverUp   = countOverlaps(NclPositionsGR_up,MUTGR_chr[MUTGR_chr$mutation==muttype])
      MutUp[,muttype]=OverUp

      OverDown = countOverlaps(NclPositionsGR_down,MUTGR_chr[MUTGR_chr$mutation==muttype])
      MutDown[,muttype]=OverDown
    }
    head(MutUp)

    #normalize the counts by their total number
    Tsum = rowSums(cbind(MutDown[-1],MutUp[-1]))
    TsumN=Tsum
    TsumN[TsumN==0]=1
    MutDownNorm=MutDown[-1]/TsumN
    MutUpNorm=MutUp[-1]/TsumN

    #calculate the frequency statistic
    kLeadUp = rowSums(MutUp[sel_mutations_norm])
    kLeadDown = rowSums(MutDown[sel_mutations_norm])
    kLagUp = rowSums(MutUp[sel_mutations_revcomp])
    kLagDown = rowSums(MutDown[sel_mutations_revcomp])
    PMA = 4/(TsumN^2) * (kLeadUp * kLagDown - kLagUp * kLeadDown)

    #save the results to a data.frame
    ResMutFreq = rbind(ResMutFreq, data.frame(chr=chr, start=NclPositions-win2,
                                              end=NclPositions+win2-1,
                                              mid=NclPositions,
                                              PMA=PMA,
                                              MutNr=Tsum,
                                              kLeadUp=kLeadUp,kLeadDown=kLeadDown,kLagUp=kLagUp,kLagDown=kLagDown))
  }
  ResMutFreq = ResMutFreq[!is.na(ResMutFreq$PMA) &  !is.na(ResMutFreq$MutNr),]

  #peak detection using peakPick
  ResMutFreq$peak=0
  ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}
  peaks = rep(FALSE,dim(ResMutFreq)[1])
  #for each chromosome separetly
  for (chr in chrs) {
    chridx = ResMutFreq$chr == chr
    ResMutFreq$PMASmooth[chridx] = as.numeric(ma(ResMutFreq$PMA[chridx],n=ceiling(100000/dist)))  ### verify for other dist values
    peaks[chridx] <- peakpick(matrix(ResMutFreq$PMASmooth[chridx], ncol=1), neighlim=100, peak.npos=50) & ResMutFreq$PMA[chridx]>PeakCut  ##neglim should be calculated
  }
  ResMutFreq$peak[peaks]=ResMutFreq$mid[peaks]

  #remove peaks located in genome gaps
  PeakCandidatesGR = GRanges(ResMutFreq$chr,IRanges(ResMutFreq$start,ResMutFreq$end))
  GapOverlaps = data.frame(findOverlaps(PeakCandidatesGR,GapsGR,ignore.strand=T))
  ResMutFreq$peak[unique(GapOverlaps$queryHits)]=0

  #extract the peak positions
  Peaks = ResMutFreq[ResMutFreq$peak>0,]
  for (i in 1:dim(Peaks)[1]) {
    tab <- matrix(c(Peaks$kLeadUp[i],Peaks$kLeadDown[i],Peaks$kLagUp[i],Peaks$kLagDown[i]),nrow = 2,dimnames = list(Lead = c("Up", "Down"),Lag = c("Up", "Down")))
    Peaks$pval[i] = fisher.test(tab, alternative = "greater")$p.value
  }
  Peaks$padj = p.adjust(Peaks$pval,method="BH")


  #adjusted p-value filtering and GRanges opject creation
  PeaksSig = Peaks[Peaks$padj<PvalCut,]
  if(dim(PeaksSig)[1]>1) {
    PeaksGR = GRanges(PeaksSig$chr,IRanges(PeaksSig$peak-dist,PeaksSig$peak+dist),strand="+",PMA=PeaksSig$PMA,Padj=PeaksSig$padj)
  } else {
    PeaksGR = GRanges()
  }

  #add adjusted p-values to the details table
  ResMutFreq$padj = 1
  ResMutFreq[ResMutFreq$peak>0,'padj'] = Peaks$padj

  Result = list(PeaksGR=PeaksGR,PMAscores=ResMutFreq)
  return(Result)
}

