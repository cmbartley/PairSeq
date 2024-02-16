# PairSeq
One Line Explanation:
PairSeq is an R library that extracts and ranks mutually enriched, overlapping linear peptides from phage display sequencing data.

Goal:
To develop an R library that allows non-computational scientists to select candidate antigens from phage display sequencing data for downstream research assays. 

Logic of Approach:
PairSeq operationalizes two logical assumptions. First, that antibodies that enrich linear phage display peptides primarily recognize contiguous epitopes. Second, that in large phage libraries, it is unlikely that epitope-sharing phage will be stochastically co-enriched above the reference background. 

Abstract:
Peptide phage display is a powerful affinity-based method that has been used for linear antigen discovery. However, when screening IgG-containing biospecimens, candidate antigen lists often include many off targets. Although selecting antigens shared across a related sample cohort can overcome this problem, antigen identification in individual biospecimens remains a significant challenge. Common methods for generating candidate antigen lists utilize a single parameter (e.g. fold change), however in our experience univariate thresholds include many off targets and can exclude valid antigens. To address this issue, we developed PairSeq, a tool implemented in R that generates lists of candidate antigens from peptide phage display sequencing data. In contrast to existing methods, PairSeq extracts and ranks mutually enriched peptides that share epitopes. In replicated experiments using commercial polyclonal antibodies to glial acidic fibrillary protein (GFAP) and gephyrin (GPHN) PairSeq significantly improved the proportion of true positives compared to univariate thresholding methods. Finally, we back tested PairSeq against paraneoplastic biospecimens and confirm its ability to identify disease-relevant autoantibodies in human biospecimens.
