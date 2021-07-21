# replicationOrigins
Functions used for replication origin detection and evaluation

**Instalation**:

>install.packages("devtools")  
>devtools::install_github("rjaksik/replicationOrigins")

**Available functions:**

_pmaORIdetection_ - Detection of replication origins based on mutation patterns associated with POLE-exo mutants  
_detectPeaksOKseq_ - Detection of RFD profile changes associated with DNA replication origins, based on OK-seq data  
_numericVectorAlign_ - Compares two numeric vectors using Nedleman-Wunsh alignment algorithm  
_multipleNumericVectorAlign_ - Performs a comparison of multiple numerical vectors (e.g. genomic coordinates) using an algorithm similar to ClustalW  
