#' Compares two numeric vectors using Nedleman-Wunsh alignment algorithm
#'
#' @param vec1 first numerical vector (e.g. genomic coordinates)
#' @param vec2 second numerical vector
#' @param gap_penalty penalty for gap insertion (default:1000000)
#' @param sort if True the vectors will be position sorted before conducting the alignemnt (in this case they cannot contain NA values) (default: False)
#'
#' @return Returns a list of  variables: Score, AlignTab and Cons. Score is the genreral similarity score between the sequences. AlignTab is a table providing paired coordinates from both vectors with NA indicating gaps. Cons is the consensus vector obtained by averaging the paired coordinates.
#' @export
#'
#' @examples
#' pos_sets = readRDS(system.file("testdata", "coordinate_sets.RDS", package = "replicationOrigins"))
#' pos_sets_num = lapply(pos_sets, function(x) x@ranges@start)
#' nva_res = numericVectorAlign(pos_sets_num$set1,pos_sets_num$set2, 100000)
numericVectorAlign = function(vec1,
                              vec2,
                              gap_penalty=1000000,
                              sort=FALSE) {

  #length of both vectors
  Nrow=length(vec1)
  Ncol=length(vec2)

  #sort the vectors if required
  if(sort) {
    NAelem = is.na(c(vec1,vec1))
    if (sum(NAelem)>0) {
      warning("Cannot sort the vectors since they contain NA values.")
    } else {
      vec1 = sort(vec1)
      vec2 = sort(vec2)
    }
  }

  #declare the score matrix (+1 to both dimentions)
  ScoreMat = matrix(NA,nrow=Nrow+1,ncol=Ncol+1)
  ScoreMat[1,1] = 0

  #first row and column of the score matrix is based only on the gap_penalty
  ScoreMat[2:(Nrow+1),1] = gap_penalty*(1:Nrow)
  ScoreMat[1,2:(Ncol+1)] = gap_penalty*(1:Ncol)

  #declare the trace table to track the movement in the score table, there are up to 3 possible movement directions if the scores are equal
  #following values are used: 1-diagonal; 2-up; 3-left
  # 1 2
  # 3
  TraceTab = array(NA,dim=c(Nrow+1,Ncol+1,3))
  TraceTab[1,1,1] = 1
  #first row tells to move to the left, first column always to the top
  TraceTab[2:(Nrow+1),1,1] = 2
  TraceTab[1,2:(Ncol+1),1] = 3

  #Step1: calculate the score matrix while tracing the movement in it
  for(i in 2:dim(ScoreMat)[1]) {
    for(j in 2:dim(ScoreMat)[2]) {
      #difference between the vector values
      dif = abs(vec1[i-1]-vec2[j-1])
      if(is.na(dif)) dif=0
      #vector of scores for diagonal, top and left neighbouring fields
      vec = c(ScoreMat[i-1,j-1]+dif,ScoreMat[i-1,j]+gap_penalty,ScoreMat[i,j-1]+gap_penalty)
      #get the smallest value and save it to the table
      minScore = min(vec)
      ScoreMat[i,j]=minScore
      #save the information on which field had the smallest value
      elem = (1:3)[vec == minScore]
      TraceTab[i,j,1:length(elem)] = elem
    }
  }

  #Step2: based on the trace table construct the alignment info
  i=Nrow+1
  j=Ncol+1
  #declare the final trace table
  AlignTab = matrix(NA,nrow=max(c(i,j))*2,ncol=2)
  idx = 0
  while (i>1 | j>1) {
    idx=idx+1

    #depending on the value in the trace table move in it in specific direction (diagonal, up down)
    if (TraceTab[i,j,1]==1) {  ## for diagonal (match)
      AlignTab[idx,1]=vec1[i-1]
      AlignTab[idx,2]=vec2[j-1]
      j=j-1
      i=i-1
    } else if (TraceTab[i,j,1]==2) { ## for up (gap)
      AlignTab[idx,1]=vec1[i-1]
      AlignTab[idx,2]=NA
      i=i-1
    } else if (TraceTab[i,j,1]==3) { ## for down (gap)
      AlignTab[idx,1]=NA
      AlignTab[idx,2]=vec2[j-1]
      j=j-1
    }
  }
  #trim and reverse the alignment table
  AlignTab = AlignTab[rev(1:idx),]
  if(idx==1) {
    AlignTab = matrix(AlignTab,ncol=2)
  }

  #calculate the alignment score
  tScore = abs(AlignTab[,1]-AlignTab[,2])
  tScore[is.na(tScore)]=gap_penalty
  Score = sum(tScore)

  #create the consensus
  Cons = rowMeans(AlignTab)
  Cons[is.na(Cons)] = AlignTab[is.na(Cons),1]
  Cons[is.na(Cons)] = AlignTab[is.na(Cons),2]

  return(list(Score=Score,AlignTab=AlignTab,Cons=Cons))
}
