library(ggplot2)
library(reshape)
library(data.table)

# #SURVIVOR merge if start and end are within 1000bp
# survivor1k<-read.delim("/Users/jzook/Documents/AJTrio/SVs/triounion_171212/SURVIVOR_dist1k_min20bp_suppvec_quotehead.csv", stringsAsFactors=FALSE, header = FALSE, sep = ",", na.strings = c("NA", NA, "NaN"))
# survivor1kvcf<-read.delim("/Users/jzook/Documents/AJTrio/SVs/triounion_171212/SURVIVOR_dist1k_min20bp.tsv", stringsAsFactors=FALSE, header = TRUE, sep = "\t", na.strings = c("NA", NA, "NaN"))
# colnames(survivor1k)<-survivor1k[1,]
# survivor1k<-survivor1k[2:nrow(survivor1k),1:118]
# indiv<-sapply(strsplit(as.character(sapply(strsplit(as.character(colnames(survivor1k)),'_'), "[", 1)),'_'), "[", 1)
# tech<-sub("10X","TenX",sapply(strsplit(as.character(sapply(strsplit(as.character(colnames(survivor1k)),'_'), "[", 2)),'_'), "[", 1))
# caller<-sub(".*_.*_(.*)","\\1",colnames(survivor1k))
# techs<-unique(tech)
# indivs<-unique(indiv)
# for (i in techs) {
#   survivor1k[,paste0(i,"calls")]<-rowSums(survivor1k[,tech==i]==1)
#   tech=c(tech,"") #these lines keep tech and indiv vectors the same length as the columns in suvivor1k
#   indiv=c(indiv,"")
#   survivor1kvcf$INFO<-paste0(survivor1kvcf$INFO,";",paste0(i,"calls"),"=",survivor1k[,paste0(i,"calls")])
# }
# 
# for (i in indivs) {
#   survivor1k[,paste0(i,"count")]<-rowSums(survivor1k[,indiv==i]==1)
#   indiv=c(indiv,"")
#   survivor1kvcf$INFO<-paste0(survivor1kvcf$INFO,";",paste0(i,"count"),"=",survivor1k[,paste0(i,"count")])
# }
# 
# survivor1k$NumTechs<-(survivor1k$TenXcalls>0)+(survivor1k$Bionanocalls>0)+(survivor1k$CGcalls>0)+(survivor1k$Illcalls>0)+(survivor1k$PBcalls>0)
# survivor1k$MultiTech<-(survivor1k$NumTechs>1)
# 
# survivor1kvcf$INFO<-paste0(survivor1kvcf$INFO,";NumTechs=",survivor1k$NumTechs,";MultiTech=",survivor1k$MultiTech,";ClusterIDs=")
# for (i in 1:118) {
#   survivor1kvcf[survivor1k[,i]==1,]$INFO<-paste0(survivor1kvcf[survivor1k[,i]==1,]$INFO,colnames(survivor1k)[i],",")
# }
# survivor1kvcf$INFO<-sub(",$","",survivor1kvcf$INFO)
# survivor1kvcf$INFO<-sub("SUPP=","NumClusterSVs=",survivor1kvcf$INFO)
# survivor1kvcf[1:5,]$INFO
# 
# survivor1kvcf$FILTER<-ifelse(survivor1k$MultiTech,"PASS","NOT2TECH")
# survivor1kvcf$QUAL<-ifelse(survivor1k$MultiTech,20,10)
# 
# contingencytable <- xtabs(~(TenXcalls>0)+(Bionanocalls>0)+(CGcalls>0)+(Illcalls>0)+(PBcalls>0), data=survivor1k)
# ftable(contingencytable,row.vars=c(2,5))
# 
# hist(survivor1k$NumTechs)
# sum(survivor1k$NumTechs>1)
# 
# survivor1kSVTYPE<-read.delim("/Users/jzook/Documents/AJTrio/SVs/triounion_171212/SURVIVOR_dist1k_min20bp_SVTYPE.csv", stringsAsFactors=FALSE, header = TRUE, sep = ",")
# survivor1k$SVTYPE <- survivor1kSVTYPE$ALT
# sum(survivor1k$NumTechs>1 & survivor1k$SVTYPE=="DEL")
# sum(survivor1k$NumTechs>1 & survivor1k$SVTYPE=="INS")
# sum(survivor1k$NumTechs>1 & survivor1k$SVTYPE=="DUP")
# sum(survivor1k$NumTechs>1 & survivor1k$SVTYPE=="INV")
# sum(survivor1k$NumTechs>1 & survivor1k$SVTYPE=="OTHER")
# 
# survivor1kvcf$FORMAT<-"."
# survivor1kvcf$AJTRIO<-"."
# #write.table(survivor1kvcf[,c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","AJTRIO")],sep="\t",row.names=FALSE,col.names=FALSE,file="/Users/jzook/Documents/AJTrio/SVs/triounion_171212/SURVIVOR_dist1k_min20bp.techcounts.vcf",quote=FALSE)
# 
# # sed 's/^X/23/;s/^Y/24/' SURVIVOR_dist1k_min20bp.techcounts.vcf | sort -k1,1n -k2,2n | sed 's/^23/X/;s/^24/Y/' > SURVIVOR_dist1k_min20bp.techcounts.sort.vcf
# # /Applications/bioinfo/tabix/tabix -f SURVIVOR_dist1k_min20bp.techcounts.vcf.gz
# # rm SURVIVOR_dist1k_min20bp.techcounts.TR*
# # /Applications/bioinfo/rtg-tools-3.7.1/rtg vcfannotate --bed-info AllRepeats_gt95percidentity_slop5.annotate.bed.gz -i SURVIVOR_dist1k_min20bp.techcounts.vcf.gz -o SURVIVOR_dist1k_min20bp.techcounts.TR.vcf.gz
# # /Applications/bioinfo/rtg-tools-3.7.1/rtg vcfannotate --bed-info segdupall.annotate.bed.gz  -i SURVIVOR_dist1k_min20bp.techcounts.TR.vcf.gz -o SURVIVOR_union_171212_v0.3.0a.vcf.gz


##SVcomp analysis

#vcf clustering output for all distance metrics within 0.2
# zgrep -v ^# union_171212_refalt.2.2.2.clustered.vcf.gz | awk 'BEGIN {FS=OFS="\t"} {if(!($1 ~ /^#/)) $7="PASS"; print}' | sed 's/PASS.*ClusterIDs=\(.*\);NumClusterSVs=\(.*\);ExactMatchIDs=\(.*\);NumExactMatchSVs=\(.*\);ClusterMaxShiftDist=\(.*\);ClusterMaxSizeDiff=\(.*\);ClusterMaxEditDist=\(.*\)    .*      .*/PASS \1      \2      \3      \4      \5      \6      \7/;' > union_171212_refalt.2.2.2.clustered.simpleINFO.tsv
# cut -f8 union_171212_refalt.2.2.2.clustered.simpleINFO.tsv | awk -F '_PB_' '{print NF-1}' > union_171212_refalt.2.2.2.clustered.simpleINFO.PBcount.tsv
# cut -f8 union_171212_refalt.2.2.2.clustered.simpleINFO.tsv | awk -F '_Ill' '{print NF-1}' > union_171212_refalt.2.2.2.clustered.simpleINFO.Illcount.tsv
# cut -f8 union_171212_refalt.2.2.2.clustered.simpleINFO.tsv | awk -F '_10X_' '{print NF-1}' > union_171212_refalt.2.2.2.clustered.simpleINFO.10Xcount.tsv
# cut -f8 union_171212_refalt.2.2.2.clustered.simpleINFO.tsv | awk -F '_CG_' '{print NF-1}' > union_171212_refalt.2.2.2.clustered.simpleINFO.CGcount.tsv
# cut -f10 union_171212_refalt.2.2.2.clustered.simpleINFO.tsv | awk -F '_PB_' '{print NF-1}' > union_171212_refalt.2.2.2.clustered.simpleINFO.PBexactcount.tsv
# cut -f10 union_171212_refalt.2.2.2.clustered.simpleINFO.tsv | awk -F '_Ill' '{print NF-1}' > union_171212_refalt.2.2.2.clustered.simpleINFO.Illexactcount.tsv
# cut -f10 union_171212_refalt.2.2.2.clustered.simpleINFO.tsv | awk -F '_10X_' '{print NF-1}' > union_171212_refalt.2.2.2.clustered.simpleINFO.10Xexactcount.tsv
# cut -f10 union_171212_refalt.2.2.2.clustered.simpleINFO.tsv | awk -F '_CG_' '{print NF-1}' > union_171212_refalt.2.2.2.clustered.simpleINFO.CGexactcount.tsv
# cut -f8 union_171212_refalt.2.2.2.clustered.simpleINFO.tsv | awk -F 'HG2_' '{print NF-1}' > union_171212_refalt.2.2.2.clustered.simpleINFO.HG2count.tsv
# cut -f8 union_171212_refalt.2.2.2.clustered.simpleINFO.tsv | awk -F 'HG3_' '{print NF-1}' > union_171212_refalt.2.2.2.clustered.simpleINFO.HG3count.tsv
# cut -f8 union_171212_refalt.2.2.2.clustered.simpleINFO.tsv | awk -F 'HG4_' '{print NF-1}' > union_171212_refalt.2.2.2.clustered.simpleINFO.HG4count.tsv
# paste union_171212_refalt.2.2.2.clustered.simpleINFO.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.PBcount.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.Illcount.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.10Xcount.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.CGcount.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.PBexactcount.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.Illexactcount.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.10Xexactcount.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.CGexactcount.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.HG2count.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.HG3count.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.HG4count.tsv | awk '$9>1' > union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.tsv
#note that awk removes sites from one callset that are not within 20% of any other callset
# paste union_171212_refalt.2.2.2.clustered.simpleINFO.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.PBcount.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.Illcount.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.10Xcount.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.CGcount.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.PBexactcount.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.Illexactcount.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.10Xexactcount.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.CGexactcount.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.HG2count.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.HG3count.tsv union_171212_refalt.2.2.2.clustered.simpleINFO.HG4count.tsv  > union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.tsv
#note that this command does not remove sites from one callset that are not within 20% of any other callset


svcompthresh2<-fread("/Users/jzook/Documents/AJTrio/SVs/triounion_171212/union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.tsv", header = FALSE, sep = "\t")
colnames(svcompthresh2)<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","ClusterIDs","NumClusterSVs","ExactMatchIDs","NumExactMatchSVs","ClusterMaxShiftDist","ClusterMaxSizeDiff","ClusterMaxEditDist","PBcalls","Illcalls","TenXcalls","CGcalls","PBexactcalls","Illexactcalls","TenXexactcalls","CGexactcalls","HG2count","HG3count","HG4count")
svcompthresh2$NumTechs<-(svcompthresh2$PBcalls>0)+(svcompthresh2$Illcalls>0)+(svcompthresh2$TenXcalls>0)+(svcompthresh2$CGcalls>0)
svcompthresh2$NumTechsExact<-(svcompthresh2$PBexactcalls>0)+(svcompthresh2$Illexactcalls>0)+(svcompthresh2$TenXexactcalls>0)+(svcompthresh2$CGexactcalls>0)
svcompthresh2$SVLEN<-(nchar(svcompthresh2$ALT)-nchar(svcompthresh2$REF))
svcompthresh2$DistBack<-c(99999,svcompthresh2[2:nrow(svcompthresh2),]$POS-svcompthresh2[1:(nrow(svcompthresh2)-1),]$POS-nchar(svcompthresh2[1:(nrow(svcompthresh2)-1),]$REF))
svcompthresh2$DistForward<-c(svcompthresh2[2:nrow(svcompthresh2),]$POS-svcompthresh2[1:(nrow(svcompthresh2)-1),]$POS-nchar(svcompthresh2[1:(nrow(svcompthresh2)-1),]$REF),99999)
svcompthresh2$DistMin<-pmin(svcompthresh2$DistBack,svcompthresh2$DistForward)
svcompthresh2$DistMinlt1000<-(svcompthresh2$DistMin<1000)
svcompthresh2$MultiTech<-(svcompthresh2$NumTechs>1)
svcompthresh2$MultiTechExact<-(svcompthresh2$NumTechsExact>1)
svcompthresh2$SVTYPE<-ifelse(nchar(svcompthresh2$ALT)>19 & nchar(svcompthresh2$REF)>19,"COMPLEX",ifelse(nchar(svcompthresh2$ALT)>19,"INS","DEL"))
svcompthresh2$END<-svcompthresh2$POS + nchar(svcompthresh2$REF) - 1
svcompthresh2$FILTER<-ifelse(svcompthresh2$MultiTech,"PASS","NOT2TECH")
svcompthresh2$QUAL<-ifelse(svcompthresh2$MultiTech,20,10)
svcompthresh2$sizecat <- ifelse(abs(svcompthresh2$SVLEN)>=1000,"gt1000",ifelse(abs(svcompthresh2$SVLEN)>=300,"300to999",ifelse(abs(svcompthresh2$SVLEN)>=100,"100to299",ifelse(abs(svcompthresh2$SVLEN)>=50,"50to99","20to49"))))

#svcompthresh2[svcompthresh2$ID=="HG4_Ill_FB_8673","ALT"]

hist(svcompthresh2$NumTechs)
hist(svcompthresh2$NumTechsExact)
sum(svcompthresh2$NumTechs>1)
sum(svcompthresh2$NumTechsExact>1)
sum(svcompthresh2[svcompthresh2$SVLEN>0,]$NumTechs>1)
sum(svcompthresh2[svcompthresh2$SVLEN<0,]$NumTechs>1)
sum(svcompthresh2[svcompthresh2$SVLEN>0,]$NumTechsExact>1)
sum(svcompthresh2[svcompthresh2$SVLEN<0,]$NumTechsExact>1)
hist(svcompthresh2[svcompthresh2$NumTechs>1,]$NumClusterSVs,breaks=(-1:400),xlim=c(0,60))
sum(svcompthresh2$SVLEN<(-299))

sum(svcompthresh2$SVLEN>0 & (svcompthresh2$NumTechs>1 ))
sum(svcompthresh2$SVLEN<0 & (svcompthresh2$NumTechs>1 ))
sum(svcompthresh2$SVLEN>0 & (svcompthresh2$NumTechs>1 | svcompthresh2$HG2count>=5 | svcompthresh2$HG3count>=5 | svcompthresh2$HG4count>=5))
sum(svcompthresh2$SVLEN<0 & (svcompthresh2$NumTechs>1 | svcompthresh2$HG2count>=5 | svcompthresh2$HG3count>=5 | svcompthresh2$HG4count>=5))
sum(svcompthresh2$SVLEN>0 & !(svcompthresh2$NumTechs>1) & (svcompthresh2$HG2count>=5 | svcompthresh2$HG3count>=5 | svcompthresh2$HG4count>=5))
sum(svcompthresh2$SVLEN<0 & !(svcompthresh2$NumTechs>1) & (svcompthresh2$HG2count>=5 | svcompthresh2$HG3count>=5 | svcompthresh2$HG4count>=5))
sum(svcompthresh2$SVLEN>300 & !(svcompthresh2$NumTechs>1) & (svcompthresh2$HG2count>=5 | svcompthresh2$HG3count>=5 | svcompthresh2$HG4count>=5))
sum(svcompthresh2$SVLEN< -300 & !(svcompthresh2$NumTechs>1) & (svcompthresh2$HG2count>=5 | svcompthresh2$HG3count>=5 | svcompthresh2$HG4count>=5))

sum(svcompthresh2$SVLEN>0 & !(svcompthresh2$NumTechs>1) & (svcompthresh2$NumClusterSVs>=4))
sum(svcompthresh2$SVLEN<0 & !(svcompthresh2$NumTechs>1) & (svcompthresh2$NumClusterSVs>=4))
sum(svcompthresh2$SVLEN>0 & (svcompthresh2$NumTechs>1 | svcompthresh2$NumClusterSVs>=4))
sum(svcompthresh2$SVLEN<0 & (svcompthresh2$NumTechs>1 | svcompthresh2$NumClusterSVs>=4))
sum(svcompthresh2$SVLEN>5000 & (svcompthresh2$NumTechs>1 | svcompthresh2$NumClusterSVs>=4))
sum(svcompthresh2$SVLEN< -5000 & (svcompthresh2$NumTechs>1 | svcompthresh2$NumClusterSVs>=4))

#how many bps?
sum(svcompthresh2[svcompthresh2$SVLEN>0 & svcompthresh2$NumTechs>1,]$SVLEN)
sum(svcompthresh2[svcompthresh2$SVLEN<0 & svcompthresh2$NumTechs>1,]$SVLEN)
sum(svcompthresh2[svcompthresh2$SVLEN>50 & svcompthresh2$NumTechs>1,]$SVLEN)
sum(svcompthresh2[svcompthresh2$SVLEN< -50 & svcompthresh2$NumTechs>1,]$SVLEN)

sum(svcompthresh2[svcompthresh2$SVLEN>0 & (svcompthresh2$NumTechs>1 | svcompthresh2$NumClusterSVs>=4),]$SVLEN)
sum(svcompthresh2[svcompthresh2$SVLEN<0 & (svcompthresh2$NumTechs>1 | svcompthresh2$NumClusterSVs>=4),]$SVLEN)
sum(svcompthresh2[svcompthresh2$SVLEN>50 & (svcompthresh2$NumTechs>1 | svcompthresh2$NumClusterSVs>=4),]$SVLEN)
sum(svcompthresh2[svcompthresh2$SVLEN< -50 & (svcompthresh2$NumTechs>1 | svcompthresh2$NumClusterSVs>=4),]$SVLEN)

hist(log10(svcompthresh2$DistMin),breaks=(-1:200)/20,xlim=c(0,6))
hist(log10(svcompthresh2[svcompthresh2$SVLEN>0,]$DistMin),breaks=(-1:200)/20,xlim=c(0,6))
sum(svcompthresh2$DistMin>1000)
sum(svcompthresh2$DistMin>1000 & svcompthresh2$NumTechs>1)
sum(svcompthresh2$DistMin>1000 & svcompthresh2$NumTechsExact>1)
sum(svcompthresh2$DistMin>1000 & svcompthresh2$NumTechsExact>1 & svcompthresh2$SVLEN>0)
sum(svcompthresh2$DistMin>1000 & svcompthresh2$NumTechsExact>1 & svcompthresh2$SVLEN<0)

hist(svcompthresh2[svcompthresh2$NumTechs>1,]$ClusterMaxShiftDist,breaks=(-1:400)/200,xlim=c(0,.5))
hist(svcompthresh2[svcompthresh2$NumTechsExact>1,]$ClusterMaxShiftDist,breaks=(-1:400)/200,xlim=c(0,.5))

hist(log10(svcompthresh2[svcompthresh2$SVLEN>0,]$SVLEN),breaks=(-1:200)/20,xlim=c(1,5),main="All candidate insertions",xlab="log10(size)")
hist(log10(-svcompthresh2[svcompthresh2$SVLEN<0,]$SVLEN),breaks=(-1:200)/20,xlim=c(1,5),main="All candidate deletions",xlab="log10(size)")
hist(log10(svcompthresh2[svcompthresh2$SVLEN>0 & svcompthresh2$NumTechs>1,]$SVLEN),breaks=(-1:200)/20,xlim=c(1,5),main="Insertions supported by 2+ techs",xlab="log10(size)")
hist(log10(-svcompthresh2[svcompthresh2$SVLEN<0 & svcompthresh2$NumTechs>1,]$SVLEN),breaks=(-1:200)/20,xlim=c(1,5),main="Deletions supported by 2+ techs",xlab="log10(size)")
hist(log10(svcompthresh2[svcompthresh2$SVLEN>0 & svcompthresh2$NumTechsExact>1,]$SVLEN),breaks=(-1:200)/20,xlim=c(1,5),main="Insertions supported by 2+ techs",xlab="log10(size)")
hist(log10(-svcompthresh2[svcompthresh2$SVLEN<0 & svcompthresh2$NumTechsExact>1,]$SVLEN),breaks=(-1:200)/20,xlim=c(1,5),main="Deletions supported by 2+ techs",xlab="log10(size)")
hist(log10(svcompthresh2[svcompthresh2$SVLEN>0 & svcompthresh2$DistMin>1000 & svcompthresh2$NumTechsExact>1,]$SVLEN),breaks=(-1:200)/20,xlim=c(1,5),main="Insertions supported by 2+ techs",xlab="log10(size)")
hist(log10(-svcompthresh2[svcompthresh2$SVLEN<0 & svcompthresh2$DistMin>1000 & svcompthresh2$NumTechsExact>1,]$SVLEN),breaks=(-1:200)/20,xlim=c(1,5),main="Deletions supported by 2+ techs",xlab="log10(size)")

#display contingency table for support from each tech, with short reads as rows and long reads/linked reads as columns, only for high-confidence calls
contingencytable <- xtabs(~(svcompthresh2$PBcalls>0)+(svcompthresh2$Illcalls>0)+(svcompthresh2$TenXcalls>0)+(svcompthresh2$CGcalls>0), data=svcompthresh2)
ftable(contingencytable,row.vars=c(1,3))

svcompthresh2<-data.frame(svcompthresh2)
svcompthresh2$INFO<-""
colno<-ncol(svcompthresh2)
for (i in 8:(colno-1)) {
  svcompthresh2$INFO<-paste0(svcompthresh2$INFO,colnames(svcompthresh2)[i],"=",svcompthresh2[,i],";")
}
svcompthresh2$INFO<-sub(";;$","",svcompthresh2$INFO)
svcompthresh2$FORMAT<-"."
svcompthresh2$AJTRIO<-"."
svcompthresh2<-data.table(svcompthresh2)


#output new vcf
write.table(svcompthresh2[,c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","AJTRIO")],sep="\t",row.names=FALSE,col.names=FALSE,file="/Users/jzook/Documents/AJTrio/SVs/triounion_171212/union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.vcf",quote=FALSE)
# cat header.txt union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.vcf | /Applications/bioinfo/htslib-1.3/bgzip -c > union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.vcf.gz
# /Applications/bioinfo/htslib-1.3/tabix -f union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.vcf.gz
# /Applications/bioinfo/bcftools/bcftools norm -D -c x -f /Users/jzook/Documents/references/human_g1k_v37.fasta union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.vcf.gz | /Applications/bioinfo/htslib-1.3/bgzip -c > union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.norm.vcf.gz
# /Applications/bioinfo/htslib-1.3/tabix -f union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.norm.vcf.gz
# rm union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.norm.TR*
# /Applications/bioinfo/rtg-tools-3.7.1/rtg vcfannotate --bed-info AllRepeats_gt95percidentity_slop5.annotate.bed.gz -i union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.norm.vcf.gz -o union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.norm.TR.vcf.gz
# /Applications/bioinfo/rtg-tools-3.7.1/rtg vcfannotate --bed-info segdupall.annotate.bed.gz  -i union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.norm.TR.vcf.gz -o svanalyzer_union_171212_v0.5.0.vcf.gz

# awk '$3-$2>10000' /Users/jzook/Downloads/segdups_selfchain_merged50.bed  | awk '{FS=OFS="\t"} { print $1, $2-50, $3+50}' | awk '{FS=OFS="\t"} { if($2<0) $2=0; print}' > /Users/jzook/Downloads/segdups_selfchain_merged50_slop50_gt10k.bed
# awk '$3-$2>10000' /Users/jzook/Downloads/repeats_trfSimplerepLowcomplex_merged50.bed  | awk '{FS=OFS="\t"} { print $1, $2-50, $3+50}' | awk '{FS=OFS="\t"} { if($2<0) $2=0; print}' > /Users/jzook/Downloads/repeats_trfSimplerepLowcomplex_merged50_slop50_gt10k.bed
# awk '$3-$2>100' /Users/jzook/Downloads/repeats_trfSimplerepLowcomplex_merged50.bed  | awk '{FS=OFS="\t"} { print $1, $2-50, $3+50}' | awk '{FS=OFS="\t"} { if($2<0) $2=0; print}' > /Users/jzook/Downloads/repeats_trfSimplerepLowcomplex_merged50_slop50_gt100.bed
# awk '{FS=OFS="\t"} { print $1, $2-5, $3+5}'  /Users/jzook/Downloads/repeats_trfSimplerepLowcomplex_merged50.bed | awk '{FS=OFS="\t"} { if($2<0) $2=0; print}' >  /Users/jzook/Downloads/repeats_trfSimplerepLowcomplex_merged50_slop5.bed
# awk '{FS="\t";OFS="\t"} {print $1,$2,$2+length($4),$3}' union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.vcf > svanalyzer_union_171212_v0.5.0.bed
# /Applications/bioinfo/bedtools2.26.0/bin/annotateBed -i svanalyzer_union_171212_v0.5.0.bed -both -files /Users/jzook/Downloads/repeats_trfSimplerepLowcomplex_merged50_slop5.bed /Users/jzook/Downloads/repeats_trfSimplerepLowcomplex_merged50_slop50_gt100.bed /Users/jzook/Downloads/repeats_trfSimplerepLowcomplex_merged50_slop50_gt10k.bed /Users/jzook/Downloads/segdups_selfchain_merged50_slop50_gt10k.bed /Applications/bioinfo/nist-integration-v3.2.2/resources/example_of_no_ref_regions_input_file_b37.bed -names  TRall TRgt100 TRgt10k segdup refNs > svanalyzer_union_171212_v0.5.0_TRall_TRgt100_TRgt10k_segdup_refN.bed
annbeds <- fread("/Users/jzook/Documents/AJTrio/SVs/triounion_171212/svanalyzer_union_171212_v0.5.0_TRall_TRgt100_TRgt10k_segdup_refN.bed", stringsAsFactors=FALSE, header = FALSE, sep = "\t")
colnames(annbeds) <- c("1","2","3","ID","TRall_cnt","TRall_pct","TRgt100_cnt","TRgt100_pct","TRgt10k_cnt","TRgt10k_pct","segdup_cnt",	"segdup_pct",	"refNs_cnt",	"refNs_pct")
setkey(annbeds,"ID")
setkey(svcompthresh2,"ID")
svcompthresh2annbeds <- merge(svcompthresh2,annbeds,by="ID",all.x=TRUE,sort=FALSE)
svcompthresh2annbeds$TRall <- (svcompthresh2annbeds$TRall_pct>0.2)
svcompthresh2annbeds$TRgt100 <- (svcompthresh2annbeds$TRgt100_pct>0.2)
svcompthresh2annbeds$TRgt10k <- (svcompthresh2annbeds$TRgt10k_pct>0.2)
svcompthresh2annbeds$segdup <- (svcompthresh2annbeds$segdup_pct>0.2)
sum(svcompthresh2annbeds$TRall,na.rm=TRUE)
sum(svcompthresh2annbeds$TRgt100,na.rm=TRUE)
sum(svcompthresh2annbeds$TRgt10k,na.rm=TRUE)
sum(svcompthresh2annbeds$segdup,na.rm=TRUE)

#output vcf with only 2+ techs or 4+ callsets for svviz eval
write.table(svcompthresh2[(svcompthresh2$NumTechs>1 | svcompthresh2$NumClusterSVs>=4),c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","AJTRIO")][order(CHROM,POS)],sep="\t",row.names=FALSE,col.names=FALSE,file="/Users/jzook/Documents/AJTrio/SVs/triounion_171212/union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.2techor4callset.vcf",quote=FALSE)
sum((svcompthresh2$NumTechs>1 | svcompthresh2$NumClusterSVs>=4))

# cat header.txt union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.2techor5caller.vcf | /Applications/bioinfo/htslib-1.3/bgzip -c > union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.2techor5caller.vcf.gz
 # /Applications/bioinfo/htslib-1.3/tabix -f union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.2techor5caller.vcf.gz


#add bionano comparison performed by Alex Hastie and Joyce Lee
#/Applications/bioinfo/bedtools2.26.0/bin/intersectBed -wa -wb -a BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/Seq_based/GM24385_HG2_seq_overlap_deletion.bed -b BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_overlap_deletion.bed > BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/GM24385_HG2_seq_BNG_overlap_deletion.bed
#/Applications/bioinfo/bedtools2.26.0/bin/intersectBed -wa -wb -a BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/Seq_based/GM24385_HG2_seq_overlap_insertion.bed -b BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_overlap_insertion.bed > BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/GM24385_HG2_seq_BNG_overlap_insertion.bed

# cat BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_overlap_deletion.bed BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_uniq_deletion.bed | sed 's/^X/23/;s/^Y/24/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^23/X/;s/^24/Y/' > BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_deletion.bed
#sed 's/_/***tab****/g' BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_deletion.bed > BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_deletion_sep.bed
#awk '$4>0.9' BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_deletion_sep.bed | /Applications/bioinfo/bedtools2.26.0/bin/mergeBed -c 4,5,6,7,8 -o max,collapse,collapse,distinct,max -i stdin > BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_deletion_sep_mergedmaxsize.bed
#/Applications/bioinfo/bedtools2.26.0/bin/intersectBed -wa -wb -a svanalyzer_union_171212_v0.5.0.bed -b ../triounion_170313/BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_deletion_sep_mergedmaxsize.bed > svanalyzer_union_171212_v0.5.0_HG2_BNG_overlap_deletion.bed
# cat BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_overlap_insertion.bed BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_uniq_insertion.bed | sed 's/^X/23/;s/^Y/24/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^23/X/;s/^24/Y/' > BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_insertion.bed
#sed 's/_/***tab****/g' BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_insertion.bed > BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_insertion_sep.bed
#awk '$4>0.9' BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_insertion_sep.bed | /Applications/bioinfo/bedtools2.26.0/bin/mergeBed -c 4,5,6,7,8 -o max,collapse,collapse,distinct,max -i stdin > BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_insertion_sep_mergedmaxsize.bed
#/Applications/bioinfo/bedtools2.26.0/bin/intersectBed -wa -wb -a svanalyzer_union_171212_v0.5.0.bed -b ../triounion_170313/BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_insertion_sep_mergedmaxsize.bed > svanalyzer_union_171212_v0.5.0_HG2_BNG_overlap_insertion.bed

# cat BNG_overlap_NIST_UnionsSV_v0.3.0/GM24149_HG3/BNG_based/GM24149_HG3_BNG_overlap_deletion.bed BNG_overlap_NIST_UnionsSV_v0.3.0/GM24149_HG3/BNG_based/GM24149_HG3_BNG_uniq_deletion.bed | sed 's/^X/23/;s/^Y/24/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^23/X/;s/^24/Y/' > BNG_overlap_NIST_UnionsSV_v0.3.0/GM24149_HG3/BNG_based/GM24149_HG3_BNG_deletion.bed
#sed 's/_/***tab****/g' BNG_overlap_NIST_UnionsSV_v0.3.0/GM24149_HG3/BNG_based/GM24149_HG3_BNG_deletion.bed > BNG_overlap_NIST_UnionsSV_v0.3.0/GM24149_HG3/BNG_based/GM24149_HG3_BNG_deletion_sep.bed
#awk '$4>0.9' BNG_overlap_NIST_UnionsSV_v0.3.0/GM24149_HG3/BNG_based/GM24149_HG3_BNG_deletion_sep.bed | /Applications/bioinfo/bedtools2.26.0/bin/mergeBed -c 4,5,6,7,8 -o max,collapse,collapse,distinct,max -i stdin > BNG_overlap_NIST_UnionsSV_v0.3.0/GM24149_HG3/BNG_based/GM24149_HG3_BNG_deletion_sep_mergedmaxsize.bed
#/Applications/bioinfo/bedtools2.26.0/bin/intersectBed -wa -wb -a svanalyzer_union_171212_v0.5.0.bed -b ../triounion_170313/BNG_overlap_NIST_UnionsSV_v0.3.0/GM24149_HG3/BNG_based/GM24149_HG3_BNG_deletion_sep_mergedmaxsize.bed > svanalyzer_union_171212_v0.5.0_HG3_BNG_overlap_deletion.bed
# cat BNG_overlap_NIST_UnionsSV_v0.3.0/GM24149_HG3/BNG_based/GM24149_HG3_BNG_overlap_insertion.bed BNG_overlap_NIST_UnionsSV_v0.3.0/GM24149_HG3/BNG_based/GM24149_HG3_BNG_uniq_insertion.bed | sed 's/^X/23/;s/^Y/24/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^23/X/;s/^24/Y/' > BNG_overlap_NIST_UnionsSV_v0.3.0/GM24149_HG3/BNG_based/GM24149_HG3_BNG_insertion.bed
#sed 's/_/***tab****/g' BNG_overlap_NIST_UnionsSV_v0.3.0/GM24149_HG3/BNG_based/GM24149_HG3_BNG_insertion.bed > BNG_overlap_NIST_UnionsSV_v0.3.0/GM24149_HG3/BNG_based/GM24149_HG3_BNG_insertion_sep.bed
# awk '$4>0.9' BNG_overlap_NIST_UnionsSV_v0.3.0/GM24149_HG3/BNG_based/GM24149_HG3_BNG_insertion_sep.bed | /Applications/bioinfo/bedtools2.26.0/bin/mergeBed -c 4,5,6,7,8 -o max,collapse,collapse,distinct,max -i stdin > BNG_overlap_NIST_UnionsSV_v0.3.0/GM24149_HG3/BNG_based/GM24149_HG3_BNG_insertion_sep_mergedmaxsize.bed
# /Applications/bioinfo/bedtools2.26.0/bin/intersectBed -wa -wb -a svanalyzer_union_171212_v0.5.0.bed -b ../triounion_170313/BNG_overlap_NIST_UnionsSV_v0.3.0/GM24149_HG3/BNG_based/GM24149_HG3_BNG_insertion_sep_mergedmaxsize.bed > svanalyzer_union_171212_v0.5.0_HG3_BNG_overlap_insertion.bed

# cat BNG_overlap_NIST_UnionsSV_v0.3.0/GM24143_HG4/BNG_based/GM24143_HG4_BNG_overlap_deletion.bed BNG_overlap_NIST_UnionsSV_v0.3.0/GM24143_HG4/BNG_based/GM24143_HG4_BNG_uniq_deletion.bed | sed 's/^X/23/;s/^Y/24/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^23/X/;s/^24/Y/' > BNG_overlap_NIST_UnionsSV_v0.3.0/GM24143_HG4/BNG_based/GM24143_HG4_BNG_deletion.bed
# sed 's/_/***tab****/g' BNG_overlap_NIST_UnionsSV_v0.3.0/GM24143_HG4/BNG_based/GM24143_HG4_BNG_deletion.bed > BNG_overlap_NIST_UnionsSV_v0.3.0/GM24143_HG4/BNG_based/GM24143_HG4_BNG_deletion_sep.bed
# awk '$4>0.9' BNG_overlap_NIST_UnionsSV_v0.3.0/GM24143_HG4/BNG_based/GM24143_HG4_BNG_deletion_sep.bed | /Applications/bioinfo/bedtools2.26.0/bin/mergeBed -c 4,5,6,7,8 -o max,collapse,collapse,distinct,max -i stdin > BNG_overlap_NIST_UnionsSV_v0.3.0/GM24143_HG4/BNG_based/GM24143_HG4_BNG_deletion_sep_mergedmaxsize.bed
# /Applications/bioinfo/bedtools2.26.0/bin/intersectBed -wa -wb -a svanalyzer_union_171212_v0.5.0.bed -b ../triounion_170313/BNG_overlap_NIST_UnionsSV_v0.3.0/GM24143_HG4/BNG_based/GM24143_HG4_BNG_deletion_sep_mergedmaxsize.bed > svanalyzer_union_171212_v0.5.0_HG4_BNG_overlap_deletion.bed
# cat BNG_overlap_NIST_UnionsSV_v0.3.0/GM24143_HG4/BNG_based/GM24143_HG4_BNG_overlap_insertion.bed BNG_overlap_NIST_UnionsSV_v0.3.0/GM24143_HG4/BNG_based/GM24143_HG4_BNG_uniq_insertion.bed | sed 's/^X/23/;s/^Y/24/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^23/X/;s/^24/Y/' > BNG_overlap_NIST_UnionsSV_v0.3.0/GM24143_HG4/BNG_based/GM24143_HG4_BNG_insertion.bed
# sed 's/_/***tab****/g' BNG_overlap_NIST_UnionsSV_v0.3.0/GM24143_HG4/BNG_based/GM24143_HG4_BNG_insertion.bed > BNG_overlap_NIST_UnionsSV_v0.3.0/GM24143_HG4/BNG_based/GM24143_HG4_BNG_insertion_sep.bed
# awk '$4>0.9' BNG_overlap_NIST_UnionsSV_v0.3.0/GM24143_HG4/BNG_based/GM24143_HG4_BNG_insertion_sep.bed | /Applications/bioinfo/bedtools2.26.0/bin/mergeBed -c 4,5,6,7,8 -o max,collapse,collapse,distinct,max -i stdin > BNG_overlap_NIST_UnionsSV_v0.3.0/GM24143_HG4/BNG_based/GM24143_HG4_BNG_insertion_sep_mergedmaxsize.bed
# /Applications/bioinfo/bedtools2.26.0/bin/intersectBed -wa -wb -a svanalyzer_union_171212_v0.5.0.bed -b ../triounion_170313/BNG_overlap_NIST_UnionsSV_v0.3.0/GM24143_HG4/BNG_based/GM24143_HG4_BNG_insertion_sep_mergedmaxsize.bed > svanalyzer_union_171212_v0.5.0_HG4_BNG_overlap_insertion.bed

#make bed for svrefine
# cat BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_uniq_deletion.bed BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_overlap_deletion.bed BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_uniq_insertion.bed BNG_overlap_NIST_UnionsSV_v0.3.0/GM24385_HG2/BNG_based/GM24385_HG2_BNG_overlap_insertion.bed BNG_overlap_NIST_UnionsSV_v0.3.0/GM24149_HG3/BNG_based/GM24149_HG3_BNG_uniq_deletion.bed BNG_overlap_NIST_UnionsSV_v0.3.0/GM24149_HG3/BNG_based/GM24149_HG3_BNG_overlap_deletion.bed BNG_overlap_NIST_UnionsSV_v0.3.0/GM24149_HG3/BNG_based/GM24149_HG3_BNG_uniq_insertion.bed BNG_overlap_NIST_UnionsSV_v0.3.0/GM24149_HG3/BNG_based/GM24149_HG3_BNG_overlap_insertion.bed  BNG_overlap_NIST_UnionsSV_v0.3.0/GM24143_HG4/BNG_based/GM24143_HG4_BNG_uniq_deletion.bed BNG_overlap_NIST_UnionsSV_v0.3.0/GM24143_HG4/BNG_based/GM24143_HG4_BNG_overlap_deletion.bed BNG_overlap_NIST_UnionsSV_v0.3.0/GM24143_HG4/BNG_based/GM24143_HG4_BNG_uniq_insertion.bed BNG_overlap_NIST_UnionsSV_v0.3.0/GM24143_HG4/BNG_based/GM24143_HG4_BNG_overlap_insertion.bed | sed 's/^X/23/;s/^Y/24/' | sort -k1,1n -k2,2n -k3,3n | sed 's/^23/X/;s/^24/Y/' > BNG_overlap_NIST_UnionsSV_v0.3.0/AJTrio_BNGbased_union.bed

#bionanoHG2del<-fread("/Users/jzook/Documents/AJTrio/SVs/triounion_171212/BNG_overlap_NIST_UnionsSV_v0.3.0/HG2BNG_overlap_all_seq/GM24385_overlap_all_seq_del.bed", stringsAsFactors=FALSE, header = FALSE, sep = "\t")
#bionanoHG2ins<-fread("/Users/jzook/Documents/AJTrio/SVs/triounion_171212/BNG_overlap_NIST_UnionsSV_v0.3.0/HG2BNG_overlap_all_seq/GM24385_overlap_all_seq_ins.bed", stringsAsFactors=FALSE, header = FALSE, sep = "\t")
bionanoHG2del<-fread("/Users/jzook/Documents/AJTrio/SVs/triounion_171212/svanalyzer_union_171212_v0.5.0_HG2_BNG_overlap_deletion.bed", stringsAsFactors=FALSE, header = FALSE, sep = "\t")
bionanoHG2ins<-fread("/Users/jzook/Documents/AJTrio/SVs/triounion_171212/svanalyzer_union_171212_v0.5.0_HG2_BNG_overlap_insertion.bed", stringsAsFactors=FALSE, header = FALSE, sep = "\t")
colnames(bionanoHG2del)<-c("SEQ_CHROM","SEQ_START","SEQ_END","ID","BNG_CHROM_HG2DEL","BNG_START_HG2DEL","BNG_END_HG2DEL","BNG_QUAL_MAX_HG2DEL","BspQIIDs_HG2DEL","BssSIIDs_HG2DEL","BNG_SVTYPE_HG2DEL","BNG_LEN_HG2DEL")
colnames(bionanoHG2ins)<-c("SEQ_CHROM","SEQ_START","SEQ_END","ID","BNG_CHROM_HG2INS","BNG_START_HG2INS","BNG_END_HG2INS","BNG_QUAL_MAX_HG2INS","BspQIIDs_HG2INS","BssSIIDs_HG2INS","BNG_SVTYPE_HG2INS","BNG_LEN_HG2INS")
bionanoHG3del<-fread("/Users/jzook/Documents/AJTrio/SVs/triounion_171212/svanalyzer_union_171212_v0.5.0_HG3_BNG_overlap_deletion.bed", stringsAsFactors=FALSE, header = FALSE, sep = "\t")
bionanoHG3ins<-fread("/Users/jzook/Documents/AJTrio/SVs/triounion_171212/svanalyzer_union_171212_v0.5.0_HG3_BNG_overlap_insertion.bed", stringsAsFactors=FALSE, header = FALSE, sep = "\t")
colnames(bionanoHG3del)<-c("SEQ_CHROM","SEQ_START","SEQ_END","ID","BNG_CHROM_HG3DEL","BNG_START_HG3DEL","BNG_END_HG3DEL","BNG_QUAL_MAX_HG3DEL","BspQIIDs_HG3DEL","BssSIIDs_HG3DEL","BNG_SVTYPE_HG3DEL","BNG_LEN_HG3DEL")
colnames(bionanoHG3ins)<-c("SEQ_CHROM","SEQ_START","SEQ_END","ID","BNG_CHROM_HG3INS","BNG_START_HG3INS","BNG_END_HG3INS","BNG_QUAL_MAX_HG3INS","BspQIIDs_HG3INS","BssSIIDs_HG3INS","BNG_SVTYPE_HG3INS","BNG_LEN_HG3INS")
bionanoHG4del<-fread("/Users/jzook/Documents/AJTrio/SVs/triounion_171212/svanalyzer_union_171212_v0.5.0_HG4_BNG_overlap_deletion.bed", stringsAsFactors=FALSE, header = FALSE, sep = "\t")
bionanoHG4ins<-fread("/Users/jzook/Documents/AJTrio/SVs/triounion_171212/svanalyzer_union_171212_v0.5.0_HG4_BNG_overlap_insertion.bed", stringsAsFactors=FALSE, header = FALSE, sep = "\t")
colnames(bionanoHG4del)<-c("SEQ_CHROM","SEQ_START","SEQ_END","ID","BNG_CHROM_HG4DEL","BNG_START_HG4DEL","BNG_END_HG4DEL","BNG_QUAL_MAX_HG4DEL","BspQIIDs_HG4DEL","BssSIIDs_HG4DEL","BNG_SVTYPE_HG4DEL","BNG_LEN_HG4DEL")
colnames(bionanoHG4ins)<-c("SEQ_CHROM","SEQ_START","SEQ_END","ID","BNG_CHROM_HG4INS","BNG_START_HG4INS","BNG_END_HG4INS","BNG_QUAL_MAX_HG4INS","BspQIIDs_HG4INS","BssSIIDs_HG4INS","BNG_SVTYPE_HG4INS","BNG_LEN_HG4INS")
bionanomatch<-merge(merge(merge(merge(merge(bionanoHG2del,bionanoHG2ins,by=c("SEQ_CHROM","SEQ_START","SEQ_END","ID"),all=TRUE,sort=FALSE),bionanoHG3del,by=c("SEQ_CHROM","SEQ_START","SEQ_END","ID"),all=TRUE,sort=FALSE),bionanoHG3ins,by=c("SEQ_CHROM","SEQ_START","SEQ_END","ID"),all=TRUE,sort=FALSE),bionanoHG4del,by=c("SEQ_CHROM","SEQ_START","SEQ_END","ID"),all=TRUE,sort=FALSE),bionanoHG4ins,by=c("SEQ_CHROM","SEQ_START","SEQ_END","ID"),all=TRUE,sort=FALSE)
setkey(bionanomatch, "ID")
svcompthresh2bionano <- merge(svcompthresh2annbeds,bionanomatch,by="ID",all.x=TRUE,sort=FALSE)
svcompthresh2bionano$BNGLENHG2DELvsSEQLEN <- 1-(svcompthresh2bionano$BNG_LEN_HG2DEL/-svcompthresh2bionano$SVLEN)
svcompthresh2bionano$BNGLENHG2INSvsSEQLEN <- 1-(svcompthresh2bionano$BNG_LEN_HG2INS/svcompthresh2bionano$SVLEN)
svcompthresh2bionano$BNGLENHG2DELdiffSEQLEN <- log10(abs((svcompthresh2bionano$BNG_LEN_HG2DEL + svcompthresh2bionano$SVLEN)))
svcompthresh2bionano$BNGLENHG2INSdiffSEQLEN <- log10(abs((svcompthresh2bionano$BNG_LEN_HG2INS - svcompthresh2bionano$SVLEN)))
plot(log10(svcompthresh2bionano$BNG_LEN_HG2DEL),log10(-svcompthresh2bionano$SVLEN))
plot(log10(svcompthresh2bionano$BNG_LEN_HG2INS),log10(svcompthresh2bionano$SVLEN))
svcompthresh2bionano$BNGLENHG3DELvsSEQLEN <- 1-(svcompthresh2bionano$BNG_LEN_HG3DEL/-svcompthresh2bionano$SVLEN)
svcompthresh2bionano$BNGLENHG3INSvsSEQLEN <- 1-(svcompthresh2bionano$BNG_LEN_HG3INS/svcompthresh2bionano$SVLEN)
svcompthresh2bionano$BNGLENHG3DELdiffSEQLEN <- log10(abs((svcompthresh2bionano$BNG_LEN_HG3DEL + svcompthresh2bionano$SVLEN)))
svcompthresh2bionano$BNGLENHG3INSdiffSEQLEN <- log10(abs((svcompthresh2bionano$BNG_LEN_HG3INS - svcompthresh2bionano$SVLEN)))
plot(log10(svcompthresh2bionano$BNG_LEN_HG3DEL),log10(-svcompthresh2bionano$SVLEN))
plot(log10(svcompthresh2bionano$BNG_LEN_HG3INS),log10(svcompthresh2bionano$SVLEN))
svcompthresh2bionano$BNGLENHG4DELvsSEQLEN <- 1-(svcompthresh2bionano$BNG_LEN_HG4DEL/-svcompthresh2bionano$SVLEN)
svcompthresh2bionano$BNGLENHG4INSvsSEQLEN <- 1-(svcompthresh2bionano$BNG_LEN_HG4INS/svcompthresh2bionano$SVLEN)
svcompthresh2bionano$BNGLENHG4DELdiffSEQLEN <- log10(abs((svcompthresh2bionano$BNG_LEN_HG4DEL + svcompthresh2bionano$SVLEN)))
svcompthresh2bionano$BNGLENHG4INSdiffSEQLEN <- log10(abs((svcompthresh2bionano$BNG_LEN_HG4INS - svcompthresh2bionano$SVLEN)))
plot(log10(svcompthresh2bionano$BNG_LEN_HG4DEL),log10(-svcompthresh2bionano$SVLEN))
plot(log10(svcompthresh2bionano$BNG_LEN_HG4INS),log10(svcompthresh2bionano$SVLEN))

contingencytable <- xtabs(~(svcompthresh2[abs(svcompthresh2$SVLEN)>999,]$PBcalls>0)+(svcompthresh2[abs(svcompthresh2$SVLEN)>999,]$Illcalls>0)+(svcompthresh2[abs(svcompthresh2$SVLEN)>999,]$TenXcalls>0)+(svcompthresh2[abs(svcompthresh2$SVLEN)>999,]$CGcalls>0), data=svcompthresh2)
ftable(contingencytable,row.vars=c(1,3))

contingencytable <- xtabs(~(svcompthresh2bionano$PBcalls>0)+(svcompthresh2bionano$Illcalls>0)+(svcompthresh2bionano$TenXcalls>0)+(svcompthresh2bionano$CGcalls>0)+(svcompthresh2bionano$BNGLENHG2DELdiffSEQLEN<log10(300)), data=svcompthresh2bionano)
ftable(contingencytable,row.vars=c(1,3,5))

contingencytable <- xtabs(~(svcompthresh2bionano$PBcalls>0)+(svcompthresh2bionano$Illcalls>0)+(svcompthresh2bionano$TenXcalls>0)+(svcompthresh2bionano$CGcalls>0)+(svcompthresh2bionano$BNGLENHG2INSdiffSEQLEN<log10(300)), data=svcompthresh2bionano)
ftable(contingencytable,row.vars=c(1,3,5))

sum(svcompthresh2$NumTechs>1 & svcompthresh2$SVLEN>999)
sum(svcompthresh2$NumTechs>1 & svcompthresh2$SVLEN < -999)

sum(svcompthresh2bionano$NumTechs==1 & svcompthresh2bionano$NumClusterSVs > 2 & svcompthresh2bionano$SVLEN>999)
sum(svcompthresh2bionano$NumTechs==1 & svcompthresh2bionano$NumClusterSVs > 2 & svcompthresh2bionano$SVLEN < -999)

length(unique(svcompthresh2bionano[svcompthresh2bionano$NumTechs==1 & svcompthresh2bionano$NumClusterSVs > 2 & svcompthresh2bionano$SVLEN>999,]$ID))
length(unique(svcompthresh2bionano[svcompthresh2bionano$NumTechs==1 & svcompthresh2bionano$NumClusterSVs > 2 & svcompthresh2bionano$SVLEN>999,]$BNG_QUAL_BspQIID_BssSIID_SVTYPE))
#dedup BNG INS confirm
length(unique(svcompthresh2bionano[svcompthresh2bionano$NumTechs==1 & svcompthresh2bionano$NumClusterSVs > 2 & svcompthresh2bionano$SVLEN>999,]$BNG_QUAL_BspQIID_BssSIID_SVTYPE))-sum(svcompthresh2bionano$NumTechs==1 & svcompthresh2bionano$NumClusterSVs > 2 & svcompthresh2bionano$SVLEN>999)+length(unique(svcompthresh2bionano[svcompthresh2bionano$NumTechs==1 & svcompthresh2bionano$NumClusterSVs > 2 & svcompthresh2bionano$SVLEN>999,]$ID))

#dedup BNG DEL confirm
length(unique(svcompthresh2bionano[svcompthresh2bionano$NumTechs==1 & svcompthresh2bionano$NumClusterSVs > 2 & svcompthresh2bionano$SVLEN< -999,]$BNG_QUAL_BspQIID_BssSIID_SVTYPE))-sum(svcompthresh2bionano$NumTechs==1 & svcompthresh2bionano$NumClusterSVs > 2 & svcompthresh2bionano$SVLEN< -999)+length(unique(svcompthresh2bionano[svcompthresh2bionano$NumTechs==1 & svcompthresh2bionano$NumClusterSVs > 2 & svcompthresh2bionano$SVLEN< -999,]$ID))

hist(log10(-svcompthresh2bionano[svcompthresh2bionano$NumTechs==1 & svcompthresh2bionano$BNGLENdiffSEQLEN<log10(300) & svcompthresh2bionano$NumClusterSVs > 2,]$SVLEN))
hist(log10(svcompthresh2bionano[svcompthresh2bionano$NumTechs==1 & svcompthresh2bionano$BNGLENdiffSEQLEN<log10(300) & svcompthresh2bionano$NumClusterSVs > 2,]$SVLEN))
sum(svcompthresh2bionano$NumTechs==1 & svcompthresh2bionano$BNGLENdiffSEQLEN<log10(300) & svcompthresh2bionano$NumClusterSVs > 2 & svcompthresh2bionano$SVLEN>999,na.rm =TRUE)
sum(svcompthresh2bionano$NumTechs==1 & svcompthresh2bionano$BNGLENdiffSEQLEN<log10(300) & svcompthresh2bionano$NumClusterSVs > 2 & svcompthresh2bionano$SVLEN < -999,na.rm =TRUE)

length(unique(svcompthresh2bionano[svcompthresh2bionano$NumTechs>1 & svcompthresh2bionano$HG2count>0 & svcompthresh2bionano$BNGLENdiffSEQLEN<log10(300) & svcompthresh2bionano$SVLEN>999,]$ID))
length(unique(svcompthresh2bionano[svcompthresh2bionano$NumTechs>1 & svcompthresh2bionano$HG2count>0 & svcompthresh2bionano$BNGLENdiffSEQLEN<log10(300) & svcompthresh2bionano$SVLEN < -999,]$ID))

length(unique(svcompthresh2bionano[svcompthresh2bionano$NumTechs>1 & svcompthresh2bionano$HG2count>0 & svcompthresh2bionano$BNGLENdiffSEQLEN>=log10(300) & svcompthresh2bionano$SVLEN>999,]$ID))
length(unique(svcompthresh2bionano[svcompthresh2bionano$NumTechs>1 & svcompthresh2bionano$HG2count>0 & svcompthresh2bionano$BNGLENdiffSEQLEN>=log10(300) & svcompthresh2bionano$SVLEN < -999,]$ID))

sum(svcompthresh2bionano$NumTechs>1 & svcompthresh2bionano$HG2count>0 & is.na(svcompthresh2bionano$BNGLENdiffSEQLEN) & svcompthresh2bionano$SVLEN>999,na.rm =TRUE)
sum(svcompthresh2bionano$NumTechs>1 & svcompthresh2bionano$HG2count>0 & is.na(svcompthresh2bionano$BNGLENdiffSEQLEN) & svcompthresh2bionano$SVLEN < -999,na.rm =TRUE)

test<- svcompthresh2bionano[svcompthresh2bionano$NumTechs>1 & svcompthresh2bionano$HG2count>0 & is.na(svcompthresh2bionano$BNGLENHG2DELdiffSEQLEN) & svcompthresh2bionano$SVLEN< -999,]
test<- svcompthresh2bionano[svcompthresh2bionano$NumTechs>1 & svcompthresh2bionano$HG2count>0 & is.na(svcompthresh2bionano$BNGLENHG2INSdiffSEQLEN) & svcompthresh2bionano$SVLEN> 999,]

#add Nabsys HG2 deletion evaluations
nabsysHG2del<-fread("/Users/jzook/Documents/AJTrio/SVs/triounion_170313/nabsys/ALLTN_Results_sizing.bed", stringsAsFactors=FALSE, header = TRUE, sep = "\t")
setkey(nabsysHG2del, "CustomerSVID")
svcompthresh2bionanonabsys <- unique(merge(svcompthresh2bionano,nabsysHG2del,all.x=TRUE,sort=FALSE,by.x="ID",by.y="CustomerSVID"),by="ID")
svcompthresh2bionanonabsys$svm_score<- (as.numeric(as.character(svcompthresh2bionanonabsys$svm_score)))
svcompthresh2bionanonabsys$MeasDelSize<- (as.numeric(as.character(svcompthresh2bionanonabsys$MeasDelSize)))
svcompthresh2bionanonabsys$NABSYSLENdiffSEQLEN <- log10(abs((svcompthresh2bionanonabsys$MeasDelSize + svcompthresh2bionanonabsys$SVLEN)))
svcompthresh2bionanonabsys$noHG2 <- (svcompthresh2bionanonabsys$HG2count==0 & svcompthresh2bionanonabsys$DistMinlt1000==FALSE)
plot(log10(svcompthresh2bionanonabsys$MeasDelSize),log10(-(svcompthresh2bionanonabsys$SVLEN)))
ggplot(svcompthresh2bionanonabsys, aes(x=NABSYSLENdiffSEQLEN)) +  geom_histogram(binwidth=0.1,aes(colour = factor(MultiTech))) + facet_grid( noHG2 ~ MultiTech)

#nrow(duplicated(svcompthresh2bionanonabsys, by="ID"))
#nrow(unique(svcompthresh2bionanonabsys, by=c("ID")))

ggplot(svcompthresh2bionanonabsys, aes(x=BNGLENHG2DELdiffSEQLEN)) +  geom_histogram(binwidth=0.1,aes(colour = factor(MultiTech))) + facet_grid( noHG2 ~ MultiTech)
ggplot(svcompthresh2bionanonabsys[!is.na(svcompthresh2bionanonabsys$BNG_LEN_HG2DEL),], aes(x=log10(BNG_LEN_HG2DEL),y=log10(-SVLEN))) +  geom_point(aes(colour = factor(MultiTech)),alpha = 1/10) + facet_grid( HG2count ~ MultiTech)
ggplot(svcompthresh2bionanonabsys[!is.na(svcompthresh2bionanonabsys$BNG_LEN_HG2DEL) & svcompthresh2bionanonabsys$HG2count<10,], aes(x=log10(BNG_LEN_HG2DEL),y=log10(-SVLEN))) +  geom_point(aes(colour = factor(MultiTech)),alpha = 1/10) + facet_grid( HG2count ~ MultiTech)
ggplot(svcompthresh2bionanonabsys, aes(x=BNGLENHG2INSdiffSEQLEN)) +  geom_histogram(binwidth=0.1,aes(colour = factor(MultiTech))) + facet_grid( noHG2 ~ MultiTech)
ggplot(svcompthresh2bionanonabsys[!is.na(svcompthresh2bionanonabsys$BNG_LEN_HG2INS),], aes(x=log10(BNG_LEN_HG2INS),y=log10(SVLEN))) +  geom_point(aes(colour = factor(MultiTech)),alpha = 1/10) + facet_grid( HG2count ~ MultiTech)
ggplot(svcompthresh2bionanonabsys[!is.na(svcompthresh2bionanonabsys$BNG_LEN_HG2INS) & svcompthresh2bionanonabsys$HG2count<10,], aes(x=log10(BNG_LEN_HG2INS),y=log10(SVLEN))) +  geom_point(aes(colour = factor(MultiTech)),alpha = 1/20) + facet_grid( HG2count ~ MultiTech + TR)

ggplot(svcompthresh2bionanonabsys[!is.na(svcompthresh2bionanonabsys$BNG_LEN_HG2DEL),], aes(x=log10(BNG_LEN_HG2DEL),y=log10(abs(SVLEN)))) +  geom_point(aes(colour = factor(MultiTech)),alpha = 1/10) + facet_grid( SVTYPE ~ MultiTech)
ggplot(svcompthresh2bionanonabsys[!is.na(svcompthresh2bionanonabsys$BNG_LEN_HG2INS),], aes(x=log10(BNG_LEN_HG2INS),y=log10(abs(SVLEN)))) +  geom_point(aes(colour = factor(MultiTech)),alpha = 1/10) + facet_grid( SVTYPE ~ MultiTech)

ggplot(svcompthresh2bionanonabsys[(svcompthresh2bionano$BNGLENHG2INSdiffSEQLEN<log10(300) | svcompthresh2bionano$BNGLENHG2DELdiffSEQLEN<log10(300)) & abs(svcompthresh2bionanonabsys$SVLEN)>100,], aes(x=log10(abs(SVLEN)))) +  geom_histogram(binwidth=0.1,aes(colour = factor(MultiTech))) + facet_grid( SVTYPE ~ MultiTech + TR)

sum(svcompthresh2bionanonabsys$NumTechs>1 | svcompthresh2bionanonabsys$NumClusterSVs>=4)
sum(svcompthresh2bionanonabsys$NumTechs>1 | svcompthresh2bionanonabsys$HG2count>=5 | svcompthresh2bionanonabsys$HG3count>=5 | svcompthresh2bionanonabsys$HG4count>=5 | svcompthresh2bionano$BNGLENHG2INSdiffSEQLEN<log10(300) | svcompthresh2bionano$BNGLENHG2DELdiffSEQLEN<log10(300),na.rm=TRUE)
sum(svcompthresh2bionanonabsys$NumTechs>1 | svcompthresh2bionanonabsys$HG2count>=5 | svcompthresh2bionanonabsys$HG3count>=5 | svcompthresh2bionanonabsys$HG4count>=5 | svcompthresh2bionano$BNGLENHG2INSdiffSEQLEN<log10(300) | svcompthresh2bionano$BNGLENHG2DELdiffSEQLEN<log10(300) | svcompthresh2bionanonabsys$svm_score> 0.9,na.rm=TRUE)
test<-svcompthresh2bionanonabsys[!(svcompthresh2bionanonabsys$NumTechs>1 | svcompthresh2bionanonabsys$NumClusterSVs>=4)  & (svcompthresh2bionano$BNGLENHG2INSdiffSEQLEN<log10(500) | svcompthresh2bionano$BNGLENHG2DELdiffSEQLEN<log10(500) | svcompthresh2bionanonabsys$svm_score> 0.9),]
setorder(test,CHROM,POS)

sum(svcompthresh2bionanonabsys$NumTechs==1 & svcompthresh2bionanonabsys$NumClusterSVs > 2 & svcompthresh2bionanonabsys$svm_score> 0.9 & svcompthresh2bionanonabsys$SVLEN < -299,na.rm =TRUE)
sum(svcompthresh2bionanonabsys$NumTechs==1 & svcompthresh2bionano$BNGLENdiffSEQLEN<log10(300) & svcompthresh2bionanonabsys$NumClusterSVs > 2 & svcompthresh2bionanonabsys$svm_score> 0.9 & svcompthresh2bionanonabsys$SVLEN < -999,na.rm =TRUE)

sum(svcompthresh2bionanonabsys$NumTechs>1 & svcompthresh2bionanonabsys$SVLEN < -299,na.rm =TRUE)
sum(svcompthresh2bionanonabsys$NumTechs>1 & svcompthresh2bionanonabsys$svm_score> 0.9 & svcompthresh2bionanonabsys$SVLEN < -299,na.rm =TRUE)
sum(svcompthresh2bionanonabsys$NumTechs>1 & svcompthresh2bionanonabsys$svm_score<= 0.9 & svcompthresh2bionanonabsys$SVLEN < -299,na.rm =TRUE)
sum(svcompthresh2bionanonabsys$NumTechs>1 & svcompthresh2bionanonabsys$svm_score< 0.2 & svcompthresh2bionanonabsys$SVLEN < -299,na.rm =TRUE)
sum(svcompthresh2bionanonabsys$NumTechs>1 & svcompthresh2bionanonabsys$HG2count>0 & svcompthresh2bionanonabsys$svm_score< 0.2 & svcompthresh2bionanonabsys$SVLEN < -299,na.rm =TRUE)

length(unique(svcompthresh2bionanonabsys[svcompthresh2bionanonabsys$NumTechs==1 & svcompthresh2bionano$BNGLENdiffSEQLEN<log10(300) & svcompthresh2bionanonabsys$NumClusterSVs > 2 & svcompthresh2bionanonabsys$svm_score> 0.9 & svcompthresh2bionanonabsys$SVLEN < -999,]$BNG_QUAL_BspQIID_BssSIID_SVTYPE))-sum(svcompthresh2bionanonabsys$NumTechs==1 & svcompthresh2bionano$BNGLENdiffSEQLEN<log10(300) & svcompthresh2bionanonabsys$NumClusterSVs > 2 & svcompthresh2bionanonabsys$svm_score> 0.9 & svcompthresh2bionanonabsys$SVLEN < -999,na.rm =TRUE)+length(unique(svcompthresh2bionanonabsys[svcompthresh2bionanonabsys$NumTechs==1 & svcompthresh2bionano$BNGLENdiffSEQLEN<log10(300) & svcompthresh2bionanonabsys$NumClusterSVs > 2 & svcompthresh2bionanonabsys$svm_score> 0.9 & svcompthresh2bionanonabsys$SVLEN < -999,]$ID))
length(unique(svcompthresh2bionanonabsys[svcompthresh2bionanonabsys$NumTechs==1 & svcompthresh2bionanonabsys$NumClusterSVs > 2 & svcompthresh2bionanonabsys$svm_score> 0.9 & svcompthresh2bionanonabsys$SVLEN < -999,]$BNG_QUAL_BspQIID_BssSIID_SVTYPE))-sum(svcompthresh2bionanonabsys$NumTechs==1 & svcompthresh2bionanonabsys$NumClusterSVs > 2 & svcompthresh2bionanonabsys$svm_score> 0.9 & svcompthresh2bionanonabsys$SVLEN < -999,na.rm =TRUE)+length(unique(svcompthresh2bionanonabsys[svcompthresh2bionanonabsys$NumTechs==1 &  svcompthresh2bionanonabsys$NumClusterSVs > 2 & svcompthresh2bionanonabsys$svm_score> 0.9 & svcompthresh2bionanonabsys$SVLEN < -999,]$ID))
length(unique(svcompthresh2bionanonabsys[svcompthresh2bionanonabsys$NumTechs==1 & svcompthresh2bionanonabsys$NumClusterSVs > 2 & svcompthresh2bionanonabsys$svm_score> 0.9 & svcompthresh2bionanonabsys$SVLEN >= -999,]$ID))
length(unique(svcompthresh2bionanonabsys[svcompthresh2bionanonabsys$NumTechs==1 & svcompthresh2bionanonabsys$NumClusterSVs > 2 & svcompthresh2bionanonabsys$svm_score> 0.9 & svcompthresh2bionanonabsys$SVLEN < -999,]$ID))

ggplot(svcompthresh2bionanonabsys[svcompthresh2bionanonabsys$NumClusterSVs > 2 & svcompthresh2bionanonabsys$svm_score> 0.9,], aes(x=log10(-SVLEN))) +  geom_histogram(binwidth=0.1,colour="white") + facet_grid( ~ MultiTech)
ggplot(svcompthresh2bionanonabsys[svcompthresh2bionanonabsys$BNGLENdiffSEQLEN<log10(300) & svcompthresh2bionanonabsys$NumClusterSVs > 2,], aes(x=log10(abs(SVLEN)))) +  geom_histogram(binwidth=0.1,colour="white") + facet_grid( SVTYPE ~ MultiTech)

ggplot(svcompthresh2bionanonabsys, aes(x=svm_score)) +  geom_histogram(binwidth=0.05,colour="white") + facet_grid( noHG2 ~ MultiTech)
ggplot(svcompthresh2bionanonabsys, aes(x=BNGLENdiffSEQLEN)) +  geom_histogram(binwidth=0.2,colour="white") + facet_grid(  SVTYPE ~ MultiTech )

ggplot(svcompthresh2bionanonabsys, aes(x=svm_score, ..density..)) +  geom_histogram(binwidth=0.05,colour="white") + facet_grid( noHG2 ~ MultiTech)
ggplot(svcompthresh2bionanonabsys, aes(x=BNGLENdiffSEQLEN, ..density..)) +  geom_histogram(binwidth=0.2,colour="white") + facet_grid(  SVTYPE ~ MultiTech )

#output new vcfs
#non-PASS sites supported by BioNano or Nabsys but not called by 5 callers
sum(!(svcompthresh2bionanonabsys$NumTechs>1 | svcompthresh2bionanonabsys$NumClusterSVs>=4)  & (svcompthresh2bionano$BNGLENHG2INSdiffSEQLEN<log10(500) | svcompthresh2bionano$BNGLENHG2DELdiffSEQLEN<log10(500) | svcompthresh2bionanonabsys$svm_score> 0.9),na.rm=TRUE)
write.table(svcompthresh2bionanonabsys[!(svcompthresh2bionanonabsys$NumTechs>1 | svcompthresh2bionanonabsys$NumClusterSVs>=4)  & (svcompthresh2bionano$BNGLENHG2INSdiffSEQLEN<log10(500) | svcompthresh2bionano$BNGLENHG2DELdiffSEQLEN<log10(500) | svcompthresh2bionanonabsys$svm_score> 0.9),c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","AJTRIO")][order(CHROM,POS)],sep="\t",row.names=FALSE,col.names=FALSE,file="/Users/jzook/Documents/AJTrio/SVs/triounion_171212/union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.nabsysorbionanoandnot2techor4callset.vcf",quote=FALSE)
# sed 's/^X/23/;s/^Y/24/;s/^MT/25/' union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.nabsysorbionanoandnot2techor5caller.vcf  | sort -k1,1n -k2,2n -k4,4 -k5,5 | sed 's/^23/X/;s/^24/Y/;s/^25/MT/'  > union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.nabsysorbionanoandnot2techor5caller.sort.vcf
# cat header.txt union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.nabsysorbionanoandnot2techor5caller.sort.vcf | /Applications/bioinfo/htslib-1.3/bgzip -c > union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.nabsysorbionanoandnot2techor5caller.vcf.gz
# /Applications/bioinfo/htslib-1.3/tabix -f union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.nabsysorbionanoandnot2techor5caller.vcf.gz
#test<-svcompthresh2bionanonabsys[!(svcompthresh2bionanonabsys$NumTechs>1 | svcompthresh2bionanonabsys$HG2count>=5 | svcompthresh2bionanonabsys$HG3count>=5 | svcompthresh2bionanonabsys$HG4count>=5)  & (svcompthresh2bionano$BNGLENHG2INSdiffSEQLEN<log10(500) | svcompthresh2bionano$BNGLENHG2DELdiffSEQLEN<log10(500) | svcompthresh2bionanonabsys$svm_score> 0.9),]

ggplot(svcompthresh2bionanonabsys, aes(x=NumClusterSVs, ..density..)) +  geom_histogram(binwidth=1,colour="black") + facet_grid( SVTYPE ~ sizecat) + xlim(c(0,50))
ggplot(svcompthresh2bionanonabsys, aes(x=log10(abs(SVLEN)), ..density..)) +  geom_histogram(binwidth=0.03,colour="black") + facet_grid( SVTYPE ~ sizecat, scales = "free") 
ggplot(svcompthresh2bionanonabsys[abs(svcompthresh2bionanonabsys$SVLEN)>999], aes(x=(abs(SVLEN)))) +  geom_histogram(binwidth=100,colour="black") + facet_grid( SVTYPE ~ sizecat, scales = "free") + xlim(c(0,20000))
ggplot(svcompthresh2bionanonabsys[(svcompthresh2bionanonabsys$NumTechs>1 | svcompthresh2bionanonabsys$NumClusterSVs>=4) & abs(svcompthresh2bionanonabsys$SVLEN)>999,], aes(x=(abs(SVLEN)))) +  geom_histogram(binwidth=100,colour="black") + facet_grid( SVTYPE ~ sizecat, scales = "free") + xlim(c(0,20000))
ggplot(svcompthresh2bionanonabsys[(svcompthresh2bionanonabsys$NumTechs>1 | svcompthresh2bionanonabsys$NumClusterSVs>=4 | svcompthresh2bionano$BNGLENHG2INSdiffSEQLEN<log10(500) | svcompthresh2bionano$BNGLENHG2DELdiffSEQLEN<log10(500) | svcompthresh2bionanonabsys$svm_score> 0.9) & abs(svcompthresh2bionanonabsys$SVLEN)>999,], aes(x=(abs(SVLEN)))) +  geom_histogram(binwidth=100,colour="black") + facet_grid( SVTYPE ~ sizecat) + xlim(c(0,20000))


#PASS dels called as incorrect by nabsys (svmscore<0.2)
write.table(svcompthresh2bionanonabsys[svcompthresh2bionanonabsys$NumTechs>1 & !is.na(svcompthresh2bionanonabsys$svm_score) & svcompthresh2bionanonabsys$svm_score< 0.2 & svcompthresh2bionanonabsys$SVLEN < -299,c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","AJTRIO")],sep="\t",row.names=FALSE,col.names=FALSE,file="/Users/jzook/Documents/AJTrio/SVs/triounion_171212/union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts_HG2DELgt299.nabsyssvmlt0.2.vcf",quote=FALSE)
 # sed 's/^X/23/;s/^Y/24/;s/^MT/25/' union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts_HG2DELgt299.nabsyssvmlt0.2.vcf  | sort -k1,1n -k2,2n -k4,4 -k5,5 | sed 's/^23/X/;s/^24/Y/;s/^25/MT/'  > union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts_HG2DELgt299.nabsyssvmlt0.2.sort.vcf
 # cat header.txt union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts_HG2DELgt299.nabsyssvmlt0.2.sort.vcf | /Applications/bioinfo/htslib-1.3/bgzip -c > union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts_HG2DELgt299.nabsyssvmlt0.2.vcf.gz
 # /Applications/bioinfo/htslib-1.3/tabix -f union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts_HG2DELgt299.nabsyssvmlt0.2.vcf.gz
View(svcompthresh2bionanonabsys[svcompthresh2bionanonabsys$NumTechs>1 & !is.na(svcompthresh2bionanonabsys$svm_score) & svcompthresh2bionanonabsys$svm_score< 0.2 & svcompthresh2bionanonabsys$SVLEN < -599,c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","PBcalls","ClusterIDs","SVLEN")])

#PASS ins missing from bionano
write.table(svcompthresh2bionano[svcompthresh2bionano$NumTechs>1 & svcompthresh2bionano$HG2count>0 & !is.na(svcompthresh2bionano$BNGLENdiffSEQLEN) & svcompthresh2bionano$BNGLENdiffSEQLEN>=log10(300),c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","AJTRIO")],sep="\t",row.names=FALSE,col.names=FALSE,file="/Users/jzook/Documents/AJTrio/SVs/triounion_171212/union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.bionanodiffgt300.vcf",quote=FALSE)
 # sed 's/^X/23/;s/^Y/24/;s/^MT/25/' union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.bionanodiffgt300.vcf  | sort -k1,1n -k2,2n -k4,4 -k5,5 | sed 's/^23/X/;s/^24/Y/;s/^25/MT/'  > union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.bionanodiffgt300.sort.vcf
 # cat header.txt union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.bionanodiffgt300.sort.vcf | /Applications/bioinfo/htslib-1.3/bgzip -c > union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.bionanodiffgt300.vcf.gz
 # /Applications/bioinfo/htslib-1.3/tabix -f union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.bionanodiffgt300.vcf.gz
length(unique(svcompthresh2bionano[svcompthresh2bionano$NumTechs>1 & svcompthresh2bionano$HG2count>0 & svcompthresh2bionano$BNGLENdiffSEQLEN>=log10(300),]$ID))
length(unique(svcompthresh2bionano[svcompthresh2bionano$NumTechs>1 & svcompthresh2bionano$HG2count>0 & svcompthresh2bionano$BNGLENdiffSEQLEN>=log10(300) & svcompthresh2bionano$SVLEN < -999,]$ID))

length(unique(svcompthresh2bionano[svcompthresh2bionano$NumTechs>1 & svcompthresh2bionano$HG2count>0 & is.na(svcompthresh2bionano$BNGLENdiffSEQLEN) & svcompthresh2bionano$SVLEN < -999,]$ID))
length(unique(svcompthresh2bionano[svcompthresh2bionano$NumTechs>1 & svcompthresh2bionano$HG2count>0 & is.na(svcompthresh2bionano$BNGLENdiffSEQLEN) & svcompthresh2bionano$SVLEN < -999,]$ID))/length(unique(svcompthresh2bionano[svcompthresh2bionano$NumTechs>1 & svcompthresh2bionano$HG2count>0 & svcompthresh2bionano$SVLEN < -999,]$ID))
hist(log10(-svcompthresh2bionano[svcompthresh2bionano$NumTechs>1 & svcompthresh2bionano$HG2count>0 & is.na(svcompthresh2bionano$BNGLENdiffSEQLEN) & svcompthresh2bionano$SVLEN < -999,]$SVLEN))
hist(log10(-svcompthresh2bionano[svcompthresh2bionano$NumTechs>1 & svcompthresh2bionano$HG2count>0 & !is.na(svcompthresh2bionano$BNGLENdiffSEQLEN) & svcompthresh2bionano$SVLEN < -999,]$SVLEN))


sum(svcompthresh2bionanonabsys$NumTechs>1 & svcompthresh2bionanonabsys$PBcalls==0 & svcompthresh2bionanonabsys$SVLEN < -299,na.rm =TRUE)
sum(svcompthresh2bionanonabsys$NumTechs>1 & svcompthresh2bionanonabsys$PBcalls==0 & svcompthresh2bionanonabsys$SVLEN > 299,na.rm =TRUE)
sum(svcompthresh2bionanonabsys$NumTechs>1 & svcompthresh2bionanonabsys$PBcalls==0 & svcompthresh2bionanonabsys$SVLEN < -49,na.rm =TRUE)
sum(svcompthresh2bionanonabsys$NumTechs>1 & svcompthresh2bionanonabsys$PBcalls==0 & svcompthresh2bionanonabsys$SVLEN > 49,na.rm =TRUE)
sum(svcompthresh2bionanonabsys$NumTechs>1 & svcompthresh2bionanonabsys$PBcalls==0 & svcompthresh2bionanonabsys$SVLEN < 0,na.rm =TRUE)
sum(svcompthresh2bionanonabsys$NumTechs>1 & svcompthresh2bionanonabsys$PBcalls==0 & svcompthresh2bionanonabsys$SVLEN > 0,na.rm =TRUE)


#Output bed for merging calls not tested by genosv
write.table(svcompthresh2bionanonabsys[!((svcompthresh2bionanonabsys$NumTechs>1 | svcompthresh2bionanonabsys$NumClusterSVs>=4)  | ((!is.na(svcompthresh2bionanonabsys$BNGLENHG2INSdiffSEQLEN) & svcompthresh2bionano$BNGLENHG2INSdiffSEQLEN<log10(500)) | (!is.na(svcompthresh2bionanonabsys$BNGLENHG2DELdiffSEQLEN) & svcompthresh2bionano$BNGLENHG2DELdiffSEQLEN<log10(500)) | (!is.na(svcompthresh2bionanonabsys$svm_score) & svcompthresh2bionanonabsys$svm_score> 0.9))),c("CHROM","POS","END","SVLEN","HG2count","PBcalls","Illcalls","TenXcalls","CGcalls")][order(CHROM,POS,END)],sep="\t",row.names=FALSE,col.names=FALSE,file="/Users/jzook/Documents/AJTrio/SVs/triounion_171212/union_171212_refalt.2.2.2.clustered.simpleINFO.techcounts.notnabsysorbionanoor2techor4callset.bed",quote=FALSE)
