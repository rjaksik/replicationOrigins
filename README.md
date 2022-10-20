
<!-- README.md is generated from README.Rmd. Please edit that file -->

# replicationOrigins

<!-- badges: start -->
<!-- badges: end -->

Functions used for replication origin detection and evaluation

## Instalation:

``` r
install.packages("devtools")  
devtools::install_github("cran/peakPick")
devtools::install_github("rjaksik/replicationOrigins")
```

## pmaORIdetection

<table style="width: 100%;">
<tr>
<td>
pmaORIdetection
</td>
<td style="text-align: right;">
R Documentation
</td>
</tr>
</table>
<h4>
Detection of replication origins based on mutation patterns associated
with POLE-exo mutants
</h4>
<h4>
Description
</h4>
<p>
Detection of replication origins based on mutation patterns associated
with POLE-exo mutants
</p>
<h4>
Usage
</h4>
<pre><code class='language-R'>pmaORIdetection(
  MUTGR,
  sel_mutations_norm = c("TCT-&gt;TAT", "TCG-&gt;TTG"),
  sel_mutations_revcomp = c("AGA-&gt;ATA", "CGA-&gt;CAA"),
  win = 2e+05,
  dist = 1000,
  PeakCut = 0.1,
  PvalCut = 0.01,
  useChrs = NULL,
  refGenome = "hg19"
)
</code></pre>
<h4>
Arguments
</h4>
<table>
<tr style="vertical-align: top;">
<td>
<code>MUTGR</code>
</td>
<td>
<p>
GenomicRanges object with somatic SNV positions, the object should
additionally contain an "mutation" column specifying mutation type
including its context e.g. TCT-\>TAT (this format is not obligatory and
can be different as long as its the same for sel_mutations\_\*
variables)
</p>
</td>
</tr>
<tr style="vertical-align: top;">
<td>
<code>sel_mutations_norm</code>
</td>
<td>
<p>
vector of mutations found upstream of the replication origins (default:
TCT-\>TAT, TCG-\>TTG)
</p>
</td>
</tr>
<tr style="vertical-align: top;">
<td>
<code>sel_mutations_revcomp</code>
</td>
<td>
<p>
vector of mutations found downstream of the replication origins
(default: AGA-\>ATA, CGA-\>CAA)
</p>
</td>
</tr>
<tr style="vertical-align: top;">
<td>
<code>win</code>
</td>
<td>
<p>
window size used to calculate the PMA score (default: 200000)
</p>
</td>
</tr>
<tr style="vertical-align: top;">
<td>
<code>dist</code>
</td>
<td>
<p>
distance between regions for which to calculate the PMA score (default:
1000)
</p>
</td>
</tr>
<tr style="vertical-align: top;">
<td>
<code>PeakCut</code>
</td>
<td>
<p>
minimum PMA score required to consider a peak an origin (default: 0.1)
</p>
</td>
</tr>
<tr style="vertical-align: top;">
<td>
<code>PvalCut</code>
</td>
<td>
<p>
adjusted p-value cutoff used to filter the ORI set using Fisher’s exact
test (default: 0.01)
</p>
</td>
</tr>
<tr style="vertical-align: top;">
<td>
<code>useChrs</code>
</td>
<td>
<p>
perform the algorithm only for a specific set of chromosomes (default:
NULL - use all basic chromosomes only)
</p>
</td>
</tr>
<tr style="vertical-align: top;">
<td>
<code>refGenome</code>
</td>
<td>
<p>
version of the reference genome supported values are hg19 and hg38
(default: hg19)
</p>
</td>
</tr>
</table>
<h4>
Value
</h4>
<p>
Returns a list of 2 variables PeaksGR and PMAscores. PeaksGR is a
GRanges object with positions of all identified peaks and corresponding
PMA score and adjusted p-value of the Fisher’s exact test. PMAscores is
a data.frame object with detailed statistics obtained for each of the
genomic regions, including the number of mutations from each group used
to calculate PMA score, smoothed PMA score (used in peak detection),
positions of peaks and adjusted p-values of the Fisher’s exact test.
</p>
</div>

## detectPeaksOKseq

<table style="width: 100%;">
<tr>
<td>
detectPeaksOKseq
</td>
<td style="text-align: right;">
R Documentation
</td>
</tr>
</table>
<h4>
Detection of RFD profile changes associated with DNA replication
origins, based on OK-seq data
</h4>
<h4>
Description
</h4>
<p>
Detection of RFD profile changes associated with DNA replication
origins, based on OK-seq data
</p>
<h4>
Usage
</h4>
<pre><code class='language-R'>detectPeaksOKseq(Data, IntSize = 2e+05, StepSize = 1000)
</code></pre>
<h4>
Arguments
</h4>
<table>
<tr style="vertical-align: top;">
<td>
<code>Data</code>
</td>
<td>
<p>
data frame containing position-specific RFD values. The format is as
following, 4 columns chr, start, end, RFD
</p>
</td>
</tr>
<tr style="vertical-align: top;">
<td>
<code>IntSize</code>
</td>
<td>
<p>
size of the window used to calculate linear regression (default:200000)
</p>
</td>
</tr>
<tr style="vertical-align: top;">
<td>
<code>StepSize</code>
</td>
<td>
<p>
distance between intervals used to calculate linear regression
(default:1000)
</p>
</td>
</tr>
</table>
<h4>
Value
</h4>
<p>
Table with coordinates of the RFD shifts associated with DNA replication
origins
</p>
</div>

## multipleNumericVectorAlign

<table style="width: 100%;">
<tr>
<td>
multipleNumericVectorAlign
</td>
<td style="text-align: right;">
R Documentation
</td>
</tr>
</table>
<h4>
Performs a comparison of multiple numerical vectors (e.g. genomic
coordinates) using an algorithm similar to ClustalW
</h4>
<h4>
Description
</h4>
<p>
Performs a comparison of multiple numerical vectors (e.g. genomic
coordinates) using an algorithm similar to ClustalW
</p>
<h4>
Usage
</h4>
<pre><code class='language-R'>multipleNumericVectorAlign(vector_list, gap_penalty = 1e+06)
</code></pre>
<h4>
Arguments
</h4>
<table>
<tr style="vertical-align: top;">
<td>
<code>vector_list</code>
</td>
<td>
<p>
list of numeric vectors to be compared
</p>
</td>
</tr>
<tr style="vertical-align: top;">
<td>
<code>gap_penalty</code>
</td>
<td>
<p>
penalty for gap insertion (default:1000000)
</p>
</td>
</tr>
</table>
<h4>
Value
</h4>
<p>
Matrix with matched coordinates from all vectors with NA indicating
gaps.
</p>
</div>

## numericVectorAlign

<table style="width: 100%;">
<tr>
<td>
numericVectorAlign
</td>
<td style="text-align: right;">
R Documentation
</td>
</tr>
</table>
<h4>
Compares two numeric vectors using Nedleman-Wunsh alignment algorithm
</h4>
<h4>
Description
</h4>
<p>
Compares two numeric vectors using Nedleman-Wunsh alignment algorithm
</p>
<h4>
Usage
</h4>
<pre><code class='language-R'>numericVectorAlign(vec1, vec2, gap_penalty = 1e+06, sort = FALSE)
</code></pre>
<h4>
Arguments
</h4>
<table>
<tr style="vertical-align: top;">
<td>
<code>vec1</code>
</td>
<td>
<p>
first numerical vector (e.g. genomic coordinates)
</p>
</td>
</tr>
<tr style="vertical-align: top;">
<td>
<code>vec2</code>
</td>
<td>
<p>
second numerical vector
</p>
</td>
</tr>
<tr style="vertical-align: top;">
<td>
<code>gap_penalty</code>
</td>
<td>
<p>
penalty for gap insertion (default:1000000)
</p>
</td>
</tr>
<tr style="vertical-align: top;">
<td>
<code>sort</code>
</td>
<td>
<p>
if True the vectors will be position sorted before conducting the
alignemnt (in this case they cannot contain NA values) (default: False)
</p>
</td>
</tr>
</table>
<h4>
Value
</h4>
<p>
Returns a list of variables: Score, AlignTab and Cons. Score is the
genreral similarity score between the sequences. AlignTab is a table
providing paired coordinates from both vectors with NA indicating gaps.
Cons is the consensus vector obtained by averaging the paired
coordinates.
</p>
</div>
