#########################################################################/**
# @RdocFunction adjustBinCount
#
# @alias adjustBinCount,QDNAseqReadCounts-method
#
# @title "Adjust binned read counts for bins that have no counts."
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
#\arguments{
#  \item{object}{An @see "QDNAseqReadCounts" object with \code{counts} data.}
#  \item{method}{A character(1) string speficying the correction method.
#    \code{median} calculates the median of \code{counts} and devides this 
#    by used.reads/1000, which is then added to \code{counts}.
#    \code{mean} calculates the mean of \code{counts} and devides this 
#    by used.reads/1000, which is then added to \code{counts}.
#    \code{map} as \code{median(fit) - fit(gc,map)}, which is added to
#    \code{counts}. Method \code{none} leaves \code{counts} untouched.}
#  \item{...}{Additional arguments passed to @see "estimateCorrection".}
#}
#
# \value{
#     Returns a @see "QDNAseqReadCounts" object with assay data element
#     \code{counts}.
# }
#
# \examples{
# data(LGG150)
# readCounts <- LGG150
# 
# readCounts<- adjustBinCount(readCounts)
# readCountsFiltered <- applyFilters(readCounts)
# readCountsFiltered <- estimateCorrection(readCountsFiltered)
# copyNumbers <- correctBins(readCountsFiltered)
# 
# or if adjustment wants to be done after applyFilters
#
# readCountsFiltered <- applyFilters(readCounts)
# readCountsFiltered <- adjustBinCount(readCountsFiltered)
# readCountsFiltered <- estimateCorrection(readCountsFiltered)
# copyNumbers <- correctBins(readCountsFiltered
# }
#
# @author "HFvE"
#
# \seealso{
#     Discussion, should the bins with 0 counts be included or excluded 
#     in the calculation of the mean/median?
# }
# @keyword manip
# 
#*/#########################################################################
setMethod("adjustBinCount", signature=c(object="QDNAseqReadCounts"),
          definition=function(object, 
                              method=c("median", "mean"), ...) {
            counts <- assayDataElement(object, "counts")
            used.reads <- object@phenoData@data$used.reads
            method <- match.arg(method)
            for (s in sampleNames(object)) {
              binSelection <- object@featureData@data$use == TRUE & assayDataElement(object, "counts")[,s] > 0
              if(method == "median") { 
                centralTendency <- median(counts[binSelection,s], na.rm = TRUE)
              } else { 
                centralTendency <- mean(counts[binSelection,s], na.rm = TRUE)
              }
              counts[,s] <- counts[,s] + (centralTendency / (used.reads / 10))
            }
            new("QDNAseqReadCounts", bins=featureData(object), counts=counts,
                phenodata=phenoData(object))
          })

# EOF