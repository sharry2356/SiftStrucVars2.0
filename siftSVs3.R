#siftSVs3.R
source("SVSiftFun.R")

Relevent_Expression_SVs_intersect_CDS<-Relevent_Expression(Expression_Analysis_Func(ReadFiles(Meristems_file = "Meristems_CDS.txt", Anthesis_file = "Anthesis_CDS.txt")
  ,zero.rm = T),15,1.25,10,1.25,1.25,0.5)
Relevent_Expression_SVs_intersect_RE<-Relevent_Expression(Expression_Analysis_Func(ReadFiles(Meristems_file = "Meristems_RE.txt", Anthesis_file = "Anthesis_RE.txt")
  ,zero.rm = T),15,1.25,10,1.25,1.25,0.5)

write.table(Relevent_Expression_SVs_intersect_CDS, "Relevent_Expression_SVs_intersect_CDS.txt",
            sep = "\t")
write.table(Relevent_Expression_SVs_intersect_RE, "Relevent_Expression_SVs_intersect_RE.txt",
            sep = "\t")

SVs_intersect_CDS<-AddReleventExpressiontoSVs("Relevent_Expression_SVs_intersect_CDS.txt","SVs_intersect_CDS.txt")

SVs_intersect_RE<-AddReleventExpressiontoSVs("Relevent_Expression_SVs_intersect_RE.txt","SVs_intersect_RE.txt")

write.table(SVs_intersect_CDS,"SVs_intersect_CDS.txt", sep = "\t")
write.table(SVs_intersect_RE,"SVs_intersect_RE.txt", sep = "\t")
