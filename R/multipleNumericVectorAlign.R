#' Performs a comparison of multiple numerical vectors (e.g. genomic coordinates) using an algorithm similar to ClustalW
#'
#' @param vector_list List of numeric vectors to be compared.
#' @param gap_penalty Penalty for gap insertion. Use lower values to match more similar positions, at the cost of increased number of gaps.
#' @param gap_symbol Value used to represent gaps
#'
#' @return Matrix with matched coordinates from all vectors with NA (or other chosen symbol) indicating gaps.
#' @export
#'
#' @examples
#' pos_sets = readRDS(system.file("testdata", "coordinate_sets.RDS", package = "replicationOrigins"))
#' pos_sets_num = lapply(pos_sets, function(x) x@ranges@start)
#' mnva_res = multipleNumericVectorAlign(pos_sets_num, 100000, '-')
multipleNumericVectorAlign = function(vector_list,
                                      gap_penalty=1000000,
                                      gap_symbol = NA) {


  #1. Calculate all possible pairwise alignments, record the score for each pair.
  Nseq = length(vector_list)
  ScoreMatrix = matrix(NA,ncol=Nseq,nrow=Nseq)
  for (i in 1:Nseq) {
    for (j in 1:Nseq) {
      if (i>=j) {
        ScoreMatrix[i,j] = ScoreMatrix[j,i] = numericVectorAlign(vector_list[[i]],vector_list[[j]],gap_penalty)$Score
      }
    }
  }
  dimnames(ScoreMatrix) <- list(1:Nseq, 1:Nseq)


  #2. Calculate a guide tree based on the pairwise distances (algorithm: Neighbor Joining).
  tree = ape::nj(ScoreMatrix)

  #create alignment order matrix
  Nnode=tree$Nnode
  Edges = tree$edge
  ResIdx = Nseq+Nnode+1
  AlnMatrix = matrix(NA,nrow=Nnode+1,ncol=3)
  for(i in 1:(Nnode)) {
    Nidx = Nseq+i
    tEdges = Edges[Edges[,1]==Nidx,]
    tEdges = tEdges[order(tEdges[,2]),]
    if(dim(tEdges)[1]==3) {
      AlnMatrix[Nnode+1,1:2] = sort(tEdges[3,])
      AlnMatrix[Nnode+1,3] = ResIdx
    }
    AlnMatrix[i,1:2] = sort(tEdges[1:2,2])
    AlnMatrix[i,3] = Nidx
  }
  AlnMatrix = AlnMatrix[order(rowSums(AlnMatrix[,1:2])),] #### Warning the order might be wrong!!!!


  #sort the order matrix (bubble sort)
  changed = T
  Noper = dim(AlnMatrix)[1]
  while (changed==T) {
    changed=F

    for (i in 1:(Noper-1)) {
      replace=F
      if (AlnMatrix[i,2]>Nseq) {
        res2idx = which(AlnMatrix[i,2]==AlnMatrix[,3])
        if(res2idx>i) {
          replace=T
        }
      }
      if (AlnMatrix[i,1]>Nseq) {
        res1idx = which(AlnMatrix[i,1]==AlnMatrix[,3])
        if(res1idx>i) {
          replace=T
        }
      }

      if(replace) {
        tmp = AlnMatrix[i,]
        AlnMatrix[i,] = AlnMatrix[i+1,]
        AlnMatrix[i+1,] = tmp
        changed=T
      }
    }
  }



  #check order
  for (i in 1:Noper) {
    if (AlnMatrix[i,2]>Nseq & !AlnMatrix[i,2] %in% AlnMatrix[1:i,3]) {
      print(i)
    }
    if (AlnMatrix[i,1]>Nseq & !AlnMatrix[i,1] %in% AlnMatrix[1:i,3]) {
      print(i)
    }
  }


  #3. Align the sequences by progressive method
  ext_vector_list = c(vector_list,as.list(rep(NA,Nnode+1)))
  for (i in 1:dim(AlnMatrix)[1]) {
    idx1=AlnMatrix[i,1]
    idx2=AlnMatrix[i,2]
    idx3=AlnMatrix[i,3]
    Aln = numericVectorAlign(ext_vector_list[[idx1]],ext_vector_list[[idx2]],gap_penalty)
    ext_vector_list[[idx1]]=Aln$AlignTab[,1]
    ext_vector_list[[idx2]]=Aln$AlignTab[,2]
    ext_vector_list[[idx3]]=Aln$Cons
  }


  #4. Expand the consensus sequences with the (gapped) original sequences
  ConLen = length(ext_vector_list[[ResIdx]])
  MSAres=matrix(NA,nrow=Nseq,ncol=ConLen)
  for (i in 1:Nseq) {
    AlignTab = numericVectorAlign(ext_vector_list[[ResIdx]],ext_vector_list[[i]],gap_penalty^2)$AlignTab
    MSAres[i,1:ConLen] = AlignTab[,2]
  }


  #5. Report the multiple sequence alignment
  MSAres[is.na(MSAres)] = gap_symbol
  return(data.frame(t(MSAres),check.names = F))
}
