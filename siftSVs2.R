#siftSVs2.R
source("SVSiftFun.R")

MutationMath(Strip_to_Integers("Full_CDS_sets.bed","overlapping_SVs_and_GOIsandCs_CDS.bed"))
SVs_intersect_RE_FUNC("SVs_upstream_of_genes.bed","SVs_only_intersect_fives.bed")
Generate_genelistsFUNC("SVs_intersect_CDS.txt","SVs_intersect_RE.txt")
