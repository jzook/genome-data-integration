# Compares summary of annotations for different VCF files
#  - uses SQLite database that contains:
#     - consensus genotype calls at each position (chrposall)
#     - annotations, including VQSR tranches, for each dataset at each postion (vcfannot)
#     - genotype calls using different datasets and variant callers (vcfhomhet)
#     - genotype calls from validation methods (e.g., OMNI microarray) - (validation)
#  - generates 3x4 matrix of histograms of annotations for each category of overlap between consensus and individual dataset's genotype calls
#  - generates tables comparing individual datasets to consensus genotypes
#  - generates tables comparing individual datasets and consensus genotypes to validation genotypes
# 
# Author: Justin Zook, NIST
###############################################################################
library(RSQLite)

drv <- dbDriver("SQLite")
#con <- dbConnect(drv,dbname="data/NA12878vcfCompare_120315/NA12878vcfs.db")
con <- dbConnect(drv,dbname="/Volumes/SSD960/workspace/data/NA12878vcfCompare_120315/NA12878vcfs.db")
#con <- dbConnect(drv,dbname="/Volumes/UHTS2/workspace/data/NA12878vcfCompare_120315/NA12878vcfs.db")
#con <- dbConnect(drv,dbname="/data/results/justin/SRA/NA12878/NA12878vcfs.db")
dbListTables(con)

#after running sqlupdateVcfUsingVQSR_v1.2_unix.pl, run these sql commands
#sqlite> update chrposall set genoSeqV=20000+genoAll,genoMapGoodSeqV=genoMapGood where genoSeqV=0;
#sqlite> update chrposall set genoSeqV=30000+genoSeq,genoMapGoodSeqV=genoMapGoodSeq where genoSeqV=0;
#sqlite> update chrposall set genoSeqV=30000+genoSeq,genoMapGoodSeqV=genoMapGoodSeq where genoSeqV=20000;
#sqlite> update chrposall set genoSeqV=0,genoMapGoodSeqV=0 where genoSeqV=50000 or (genoSeqV>20000 and genoSeqV<30000);


## Order of datasets in database
metaInf <- c("HSWG","HSWEx","ILLWG","ILLWEx","XIll","SOLWG","454WG","454WEx","CG","ILLCLIA","IllPCRFreeNoBQSR","XSolWG")
nf <- length(metaInf)



##plot 3x4 matrix of histograms of annotations for each category of overlap between consensus and individual dataset's genotype calls
	
	plotHistTruthTable <- function(dataset, vcf, annType, minxlimfactor, maxxlimfactor,addwhere="") {
		#select only rows with coverage >= cov
	
    nann <- length(annType)
    annlist <- paste("vcfannot.",annType#[1],sep="")
    for (i in 2:nann) {
      annlist <- paste(annlist,",vcfannot.",annType[i],sep="")
    } 
    #print(annlist)
    xv <- dbGetQuery(con,paste("SELECT cast(vcfhomhet.chrompos as real) as chrompos, vcfhomhet.geno genov FROM vcfhomhet WHERE vcfhomhet.vcf=",vcf," AND vcfhomhet.dataset=",dataset,sep=""))
		xall <- dbGetQuery(con,paste("SELECT cast(vcfannot.chrompos as real) as chrompos, ",annlist,", (chrposall.genoAll+chrposall.genoSeqV) genoc, (chrposall.genoMapGood+chrposall.genoMapGoodSeqV) genoMapGood, TrancheSSEmin2,TrancheABQDmin2,TrancheAlignmin2,TrancheMapmin2 FROM chrposall, vcfannot WHERE chrposall.chrompos=vcfannot.chrompos AND vcfannot.bam=",dataset,addwhere,sep=""))
		mm <- match(xall$chrompos,xv$chrompos)
		xall$genov <- xv[mm,"genov"]
	
		r=3
	#	if (nrow(subset(xall,xall$genov%%10==0))>0) {r=4} 
			
    for (i in 1:nann) {
        print(annType[i])
        if (annType[i]=="DP") {cov=0} else {cov=5}
          maxann<-max(xall[!is.na(xall[,annType[i]]),annType[i]])
          minann<-min(xall[!is.na(xall[,annType[i]]),annType[i]])
        if (minann==maxann) {next}
        breakslist <- ((minann/(maxann-minann)*100):(maxann/(maxann-minann)*100))*(maxann-minann)/100
        if (annType[i]=="AB" | annType[i]=="ABCI5" |annType[i]=="ABCI50" |annType[i]=="ABCI95" |annType[i]=="OND") {
          if (annType[i]=="OND") {xlims <- c(0,0.5)} else {xlims <- c(0,1)}
        } else {
          xlims <- c(minann*minxlimfactor[i],maxann*maxxlimfactor[i])
        }
        
        # generate figures for filtered consensus genotypes without mapping filter for Hom calls
        tiff(paste("plots/TruthTableHist_genoAllSeq_Sol_WEx_1.2V_NoHomMapfilt_dataset",dataset,"_vcf",vcf,"_Cov",cov,annType[i],".tif",sep=""),width=900,height=900*r/4)
        layout(matrix(1:(4*r), ncol=4, byrow=TRUE))
        op <- par(mar=rep(2,4))
        
        ##HomRef in VCF
        tempsub <- subset(xall,is.na(xall$genov) & (xall$genoc%%10==1 & xall$genoc%%1000>30 & xall$TrancheAlignmin2<99) & xall$DP>=cov) 
        hist(tempsub[,annType[i]], breaks=breakslist, xlim=xlims, main="Consensus HomRef and VCF HomRef")
        tempsub <- subset(xall,is.na(xall$genov) & (xall$genoc%%10==2 & xall$genoMapGood>=3 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95) & xall$DP>=cov) 
        hist(tempsub[,annType[i]], breaks=breakslist, xlim=xlims, main="Consensus Het and VCF HomRef")
        tempsub <- subset(xall,is.na(xall$genov) & (xall$genoc%%10==3 & xall$genoc%%1000>30 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99) & xall$DP>=cov) 
        hist(tempsub[,annType[i]], breaks=breakslist, xlim=xlims, main="Consensus HomVar and VCF HomRef")
        tempsub <- subset(xall,is.na(xall$genov) & (xall$genoc%%10==0 | (xall$genoc%%10==1 & (xall$TrancheAlignmin2>=99 | xall$genoc%%1000<=30)) | (xall$genoc%%10==2 & (xall$TrancheSSEmin2>=95 | xall$TrancheABQDmin2>=95 | xall$TrancheAlignmin2>=95 | xall$TrancheMapmin2>=95 | xall$genoMapGood<3)) | (xall$genoc%%10==3 & (xall$TrancheABQDmin2>=99 | xall$TrancheAlignmin2>=99 | xall$genoc%%1000<=30))) & xall$DP>=cov)
        hist(tempsub[,annType[i]], breaks=breakslist, xlim=xlims, main="Consensus Uncertain and VCF HomRef")
        
        ##Het in VCF
        tempsub <- subset(xall,xall$genov==2 & (xall$genoc%%10==1 & xall$genoMapGood>=3 & xall$TrancheAlignmin2<99 & xall$TrancheMapmin2<99) & xall$DP>=cov) 
        hist(tempsub[,annType[i]], breaks=breakslist, xlim=xlims, main="Consensus HomRef and VCF Het")
        tempsub <- subset(xall,xall$genov==2 & (xall$genoc%%10==2 & xall$genoMapGood>=3 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95) & xall$DP>=cov) 
        hist(tempsub[,annType[i]], breaks=breakslist, xlim=xlims, main="Consensus Het and VCF Het")
        tempsub <- subset(xall,xall$genov==2 & (xall$genoc%%10==3 & xall$genoc%%1000>30 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99) & xall$DP>=cov)
        hist(tempsub[,annType[i]], breaks=breakslist, xlim=xlims, main="Consensus HomVar and VCF Het")
        tempsub <- subset(xall,xall$genov==2 & (xall$genoc%%10==0 | (xall$genoc%%10==1 & (xall$TrancheAlignmin2>=99 | xall$genoc%%1000<=30)) | (xall$genoc%%10==2 & (xall$TrancheSSEmin2>=95 | xall$TrancheABQDmin2>=95 | xall$TrancheAlignmin2>=95 | xall$TrancheMapmin2>=95 | xall$genoMapGood<3)) | (xall$genoc%%10==3 & (xall$TrancheABQDmin2>=99 | xall$TrancheAlignmin2>=99 | xall$genoc%%1000<=30))) & xall$DP>=cov)
        hist(tempsub[,annType[i]], breaks=breakslist, xlim=xlims, main="Consensus Uncertain and VCF Het")
        
        ##HomVar in VCF
        tempsub <- subset(xall,xall$genov==3 & (xall$genoc%%10==1 & xall$genoc%%1000>30 & xall$TrancheAlignmin2<99) & xall$DP>=cov) 
        hist(tempsub[,annType[i]], breaks=breakslist, xlim=xlims, main="Consensus HomRef and VCF HomVar")
        tempsub <- subset(xall,xall$genov==3 & (xall$genoc%%10==2 & xall$genoMapGood>=3 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95) & xall$DP>=cov) 
        hist(tempsub[,annType[i]], breaks=breakslist, xlim=xlims, main="Consensus Het and VCF HomVar")
        tempsub <- subset(xall,xall$genov==3 & (xall$genoc%%10==3 & xall$genoc%%1000>30 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99) & xall$DP>=cov)
        hist(tempsub[,annType[i]], breaks=breakslist, xlim=xlims, main="Consensus HomVar and VCF HomVar")
        tempsub <- subset(xall,xall$genov==3 & (xall$genoc%%10==0 | (xall$genoc%%10==1 & (xall$TrancheAlignmin2>=99 | xall$genoc%%1000<=30)) | (xall$genoc%%10==2 & (xall$TrancheSSEmin2>=95 | xall$TrancheABQDmin2>=95 | xall$TrancheAlignmin2>=95 | xall$TrancheMapmin2>=95 | xall$genoMapGood<3)) | (xall$genoc%%10==3 & (xall$TrancheABQDmin2>=99 | xall$TrancheAlignmin2>=99 | xall$genoc%%1000<=30))) & xall$DP>=cov)
        hist(tempsub[,annType[i]], breaks=breakslist, xlim=xlims, main="Consensus Uncertain and VCF HomVar")
        
        ##Other in VCF
          tempsub <- subset(xall,xall$genov%%10==0 & (xall$genoc%%10==1 & xall$genoc%%1000>30 & xall$TrancheAlignmin2<99) & xall$DP>=cov) 
        if (r>3 & nrow(tempsub)>0) {
          hist(tempsub[,annType[i]], breaks=breakslist, xlim=xlims, main="Consensus HomRef and VCF Other")
          tempsub <- subset(xall,xall$genov%%10==0 & (xall$genoc%%10==2 & xall$genoMapGood>=3 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95) & xall$DP>=cov) 
          hist(tempsub[,annType[i]], breaks=breakslist, xlim=xlims, main="Consensus Het and VCF Other")
          tempsub <- subset(xall,xall$genov%%10==0 & (xall$genoc%%10==3 & xall$genoc%%1000>30 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99) & xall$DP>=cov)
          hist(tempsub[,annType[i]], breaks=breakslist, xlim=xlims, main="Consensus HomVar and VCF Other")
          tempsub <- subset(xall,xall$genov%%10==0 & (xall$genoc%%10==0 | (xall$genoc%%10==1 & (xall$TrancheAlignmin2>=99 | xall$genoc%%1000<=30)) | (xall$genoc%%10==2 & (xall$TrancheSSEmin2>=95 | xall$TrancheABQDmin2>=95 | xall$TrancheAlignmin2>=95 | xall$TrancheMapmin2>=95 | xall$genoMapGood<3)) | (xall$genoc%%10==3 & (xall$TrancheABQDmin2>=99 | xall$TrancheAlignmin2>=99 | xall$genoc%%1000<=30))) & xall$DP>=cov)
          hist(tempsub[,annType[i]], breaks=breakslist, xlim=xlims, main="Consensus Uncertain and VCF Other")
          
        }
        
        par(op)
        dev.off()
        gc()
        
        
        print(i)
    }
	}
dataset=1
vcf=2
#plotHistTruthTable(dataset,vcf,c("DP","AB","ReadPosRankSum","ReadPosEndDist","ReadMeanPos","ReadMeanLen","QD","HaplotypeScore","MapQRankSum","MQ","FS","ABCI50","ABCI5","ABCI95","OND","MQ0frac","VQSLODSSEHet","VQSLODABQDHet","VQSLODAlignHet","VQSLODMapHet"),c(0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,1,1,1,1),c(1,1,1,1,1,1,1,0.1,1,1,0.3,1,1,1,1,1,1,1,1,1))
plotHistTruthTable(dataset,vcf,c("DP","AB","ReadPosRankSum","ReadPosEndDist","ReadMeanPos","ReadMeanLen","QD","HaplotypeScore","MapQRankSum","MQ","FS","ABCI50","ABCI5","ABCI95","OND","MQ0frac"),c(0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0),c(1,1,1,1,1,1,1,0.1,1,1,0.3,1,1,1,1,1))
dataset=9
vcf=4
plotHistTruthTable(dataset,vcf,c("DP","AB","ReadPosRankSum","ReadPosEndDist","ReadMeanPos","ReadMeanLen","QD","HaplotypeScore","MapQRankSum","MQ","FS","ABCI50","ABCI5","ABCI95","OND","MQ0frac"),c(0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0),c(1,1,1,1,1,1,1,0.1,1,1,0.3,1,1,1,1,1))
dataset=10
vcf=6
plotHistTruthTable(dataset,vcf,c("DP","AB","ReadPosRankSum","ReadPosEndDist","ReadMeanPos","ReadMeanLen","QD","HaplotypeScore","MapQRankSum","MQ","FS","ABCI50","ABCI5","ABCI95","OND","MQ0frac"),c(0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0),c(1,1,1,1,1,1,1,0.1,1,1,0.3,1,1,1,1,1))
dataset=7
vcf=5
plotHistTruthTable(dataset,vcf,c("DP","AB","ReadPosRankSum","ReadPosEndDist","ReadMeanPos","ReadMeanLen","QD","HaplotypeScore","MapQRankSum","MQ","FS","ABCI50","ABCI5","ABCI95","OND","MQ0frac"),c(0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0),c(1,1,1,1,1,1,1,0.1,1,1,0.3,1,1,1,1,1))
dataset=1
vcf=5
plotHistTruthTable(dataset,vcf,c("DP","AB","ReadPosRankSum","ReadPosEndDist","ReadMeanPos","ReadMeanLen","QD","HaplotypeScore","MapQRankSum","MQ","FS","ABCI50","ABCI5","ABCI95","OND","MQ0frac"),c(0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0),c(1,1,1,1,1,1,1,0.1,1,1,0.3,1,1,1,1,1))
dataset=12
vcf=1
plotHistTruthTable(dataset,vcf,c("DP","AB","ReadPosRankSum","ReadPosEndDist","ReadMeanPos","ReadMeanLen","QD","HaplotypeScore","MapQRankSum","MQ","FS","ABCI50","ABCI5","ABCI95","OND","MQ0frac"),c(0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0),c(1,1,1,1,1,1,1,0.1,1,1,0.3,1,1,1,1,1))




##Count variants in each truth table category comparing vcfs to consensus calls, not filtering by mapping bias for hom calls
datavcf <- dbGetQuery(con,paste("SELECT DISTINCT vcfhomhet.dataset, vcfhomhet.vcf FROM vcfhomhet",sep=""))
ndv <- nrow(datavcf)
writeTruthTableSum <- matrix(0,ndv,16)
xall <- dbGetQuery(con,paste("SELECT ref,alt,cast(chrposall.chrompos as real) as chrompos, (chrposall.genoAll+chrposall.genoSeqV) genoc, (chrposall.genoMapGood+chrposall.genoMapGoodSeqV) genoMapGood, TrancheSSEmin2,TrancheABQDmin2,TrancheAlignmin2,TrancheMapmin2 FROM chrposall",sep=""))
for (i in 1:(ndv)) {
  xv <- dbGetQuery(con,paste("SELECT cast(vcfhomhet.chrompos as real) as chrompos, vcfhomhet.geno genov FROM vcfhomhet WHERE vcfhomhet.vcf=",datavcf$vcf[i]," AND vcfhomhet.dataset=",datavcf$dataset[i],sep=""))
  mm <- match(xall$chrompos,xv$chrompos)
  xall$genov <- 0
  xall$genov <- xv[mm,"genov"]
  

  genos=c("0/0","0/1","1/1")
  genotxt=c("HomRef","Het","HomVar")
  k=0
  for (j in c(1,2,3)) {
    tmpsel <- subset(xall[,c("chrompos","genoc","ref","alt")],((j==1 & is.na(xall$genov)) | (j!=1 & xall$genov==j)) & (xall$genoc%%10==1 & xall$TrancheAlignmin2<99 & xall$genoc%%1000>30))
    exportvcf(tmpsel,paste("/Volumes/SSD960/workspace/plots/TruthTableSum_genoAllSeq_Sol_WEx_v1.2V_NoHomMapfilt_consHomRef_vcf",i,genotxt[k+1],".vcf",sep=""),genos[k+1])
    writeTruthTableSum[i,k*4+1] <- nrow(tmpsel) 
    k=k+1
  }
  
  k=0
  for (j in c(1,2,3)) {
    tmpsel <- subset(xall[,c("chrompos","genoc","ref","alt")],((j==1 & is.na(xall$genov)) | (j!=1 & xall$genov==j)) & (xall$genoc%%10==2 & xall$genoMapGood>=3 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95))
    exportvcf(tmpsel,paste("/Volumes/SSD960/workspace/plots/TruthTableSum_genoAllSeq_Sol_WEx_v1.2V_NoHomMapfilt_consHet_vcf",i,genotxt[k+1],".vcf",sep=""),genos[k+1])
    writeTruthTableSum[i,k*4+2] <- nrow(tmpsel) 
    k=k+1
  }
  
  k=0
  for (j in c(1,2,3)) {
    tmpsel <- subset(xall[,c("chrompos","genoc","ref","alt")],((j==1 & is.na(xall$genov)) | (j!=1 & xall$genov==j)) & (xall$genoc%%10==3 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$genoc%%1000>30))
    exportvcf(tmpsel,paste("/Volumes/SSD960/workspace/plots/TruthTableSum_genoAllSeq_Sol_WEx_v1.2V_NoHomMapfilt_consHomVar_vcf",i,genotxt[k+1],".vcf",sep=""),genos[k+1])
    writeTruthTableSum[i,k*4+3] <- nrow(tmpsel) 
    k=k+1
  }
  
  k=0
  for (j in c(1,2,3)) {
    tmpsel <- subset(xall[,c("chrompos","genoc","ref","alt")],((j==1 & is.na(xall$genov)) | (j!=1 & xall$genov==j)) & (xall$genoc%%10==0 | (xall$genoc%%10==1 & (xall$TrancheAlignmin2>=99 | xall$genoc%%1000<=30)) | (xall$genoc%%10==2 & (xall$TrancheSSEmin2>=95 | xall$TrancheABQDmin2>=95 | xall$TrancheAlignmin2>=95 | xall$TrancheMapmin2>=95 | xall$genoMapGood<3)) | (xall$genoc%%10==3 & (xall$TrancheABQDmin2>=99 | xall$TrancheAlignmin2>=99 | xall$genoc%%1000<=30))))
    exportvcf(tmpsel,paste("/Volumes/SSD960/workspace/plots/TruthTableSum_genoAllSeq_Sol_WEx_v1.2V_NoHomMapfilt_consUncert_vcf",i,genotxt[k+1],".vcf",sep=""),genos[k+1])
    writeTruthTableSum[i,k*4+4] <- nrow(tmpsel) 
    k=k+1
  }
  
  ##HomRef in VCF
  writeTruthTableSum[i,1] <- nrow(subset(xall,is.na(xall$genov) & (xall$genoc%%10==1 & xall$TrancheAlignmin2<99 & xall$genoc%%1000>30))) #& xall$genoMapGood>=3 & xall$TrancheAlignmin2<99 & xall$TrancheMapmin2<99))) 
  writeTruthTableSum[i,2] <- nrow(subset(xall,is.na(xall$genov) & (xall$genoc%%10==2 & xall$genoMapGood>=3 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95))) 
  writeTruthTableSum[i,3] <- nrow(subset(xall,is.na(xall$genov) & (xall$genoc%%10==3 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$genoc%%1000>30)) ) #& xall$genoMapGood>=3 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$TrancheMapmin2<99)) )
  #	writeTruthTableSum[i,4] <- nrow(subset(xall,is.na(xall$genov) & (xall$genoc%%10==0 | xall$genoMapGood<3 | (xall$genoc%%10==1 & (xall$TrancheAlignmin2>=99 | xall$TrancheMapmin2>=99)) | (xall$genoc%%10==2 & (xall$TrancheSSEmin2>=95 | xall$TrancheABQDmin2>=95 | xall$TrancheAlignmin2>=95 | xall$TrancheMapmin2>=95)) | (xall$genoc%%10==3 & (xall$TrancheABQDmin2>=99 | xall$TrancheAlignmin2>=99 | xall$TrancheMapmin2>=99)))) )
  writeTruthTableSum[i,4] <- nrow(subset(xall,is.na(xall$genov) & (xall$genoc%%10==0 | (xall$genoc%%10==1 & (xall$TrancheAlignmin2>=99 | xall$genoc%%1000<=30)) | (xall$genoc%%10==2 & (xall$TrancheSSEmin2>=95 | xall$TrancheABQDmin2>=95 | xall$TrancheAlignmin2>=95 | xall$TrancheMapmin2>=95 | xall$genoMapGood<3)) | (xall$genoc%%10==3 & (xall$TrancheABQDmin2>=99 | xall$TrancheAlignmin2>=99 | xall$genoc%%1000<=30)))) )
 # tmpsel <- subset(xall[,c("chrompos","genoc","ref","alt")],is.na(xall$genov) & (xall$genoc%%10==2 & xall$genoMapGood>=3 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95),paste("/Volumes/SSD960/workspace/plots/TruthTableSum_genoAllSeq_Sol_WEx_v1.2V_NoHomMapfilt_consHet_vcf",i,"HomRef.vcf",sep=""),"0/1")
  
  ##Het in VCF
  writeTruthTableSum[i,5] <- nrow(subset(xall,xall$genov==2 & (xall$genoc%%10==1 & xall$TrancheAlignmin2<99 & xall$genoc%%1000>30))) #& xall$genoMapGood>=3 & xall$TrancheAlignmin2<99 & xall$TrancheMapmin2<99)))
  writeTruthTableSum[i,6] <- nrow(subset(xall,xall$genov==2 & (xall$genoc%%10==2 & xall$genoMapGood>=3 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95)) )
  writeTruthTableSum[i,7] <- nrow(subset(xall,xall$genov==2 & (xall$genoc%%10==3 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$genoc%%1000>30)) ) #& xall$genoMapGood>=3 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$TrancheMapmin2<99)) )
  #	writeTruthTableSum[i,8] <- nrow(subset(xall,xall$genov==2 & (xall$genoc%%10==0 | xall$genoMapGood<3 | (xall$genoc%%10==1 & (xall$TrancheAlignmin2>=99 | xall$TrancheMapmin2>=99)) | (xall$genoc%%10==2 & (xall$TrancheSSEmin2>=95 | xall$TrancheABQDmin2>=95 | xall$TrancheAlignmin2>=95 | xall$TrancheMapmin2>=95)) | (xall$genoc%%10==3 & (xall$TrancheABQDmin2>=99 | xall$TrancheAlignmin2>=99 | xall$TrancheMapmin2>=99)))) )
  writeTruthTableSum[i,8] <- nrow(subset(xall,xall$genov==2 & (xall$genoc%%10==0 | (xall$genoc%%10==1 & (xall$TrancheAlignmin2>=99 | xall$genoc%%1000<=30)) | (xall$genoc%%10==2 & (xall$TrancheSSEmin2>=95 | xall$TrancheABQDmin2>=95 | xall$TrancheAlignmin2>=95 | xall$TrancheMapmin2>=95 | xall$genoMapGood<3)) | (xall$genoc%%10==3 & (xall$TrancheABQDmin2>=99 | xall$TrancheAlignmin2>=99 | xall$genoc%%1000<=30)))) )
  
  ##HomVar in VCF
  writeTruthTableSum[i,9] <- nrow(subset(xall,xall$genov==3 & (xall$genoc%%10==1 & xall$TrancheAlignmin2<99 & xall$genoc%%1000>30))) #& xall$genoMapGood>=3 & xall$TrancheAlignmin2<99 & xall$TrancheMapmin2<99)))
  writeTruthTableSum[i,10] <- nrow(subset(xall,xall$genov==3 & (xall$genoc%%10==2 & xall$genoMapGood>=3 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95)) )
  writeTruthTableSum[i,11] <- nrow(subset(xall,xall$genov==3 & (xall$genoc%%10==3 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$genoc%%1000>30)) ) #& xall$genoMapGood>=3 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$TrancheMapmin2<99)) )
  #	writeTruthTableSum[i,12] <- nrow(subset(xall,xall$genov==3 & (xall$genoc%%10==0 | xall$genoMapGood<3 | (xall$genoc%%10==1 & (xall$TrancheAlignmin2>=99 | xall$TrancheMapmin2>=99)) | (xall$genoc%%10==2 & (xall$TrancheSSEmin2>=95 | xall$TrancheABQDmin2>=95 | xall$TrancheAlignmin2>=95 | xall$TrancheMapmin2>=95)) | (xall$genoc%%10==3 & (xall$TrancheABQDmin2>=99 | xall$TrancheAlignmin2>=99 | xall$TrancheMapmin2>=99)))) )
  writeTruthTableSum[i,12] <- nrow(subset(xall,xall$genov==3 & (xall$genoc%%10==0 | (xall$genoc%%10==1 & (xall$TrancheAlignmin2>=99 | xall$genoc%%1000<=30)) | (xall$genoc%%10==2 & (xall$TrancheSSEmin2>=95 | xall$TrancheABQDmin2>=95 | xall$TrancheAlignmin2>=95 | xall$TrancheMapmin2>=95 | xall$genoMapGood<3)) | (xall$genoc%%10==3 & (xall$TrancheABQDmin2>=99 | xall$TrancheAlignmin2>=99 | xall$genoc%%1000<=30)))) )
  
  ##Other in VCF
  writeTruthTableSum[i,13] <- nrow(subset(xall,xall$genov%%10==0 & (xall$genoc%%10==1 & xall$TrancheAlignmin2<99 & xall$genoc%%1000>30))) #& xall$genoMapGood>=3 & xall$TrancheAlignmin2<99 & xall$TrancheMapmin2<99)))
  writeTruthTableSum[i,14] <- nrow(subset(xall,xall$genov%%10==0 & (xall$genoc%%10==2 & xall$genoMapGood>=3 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95)) )
  writeTruthTableSum[i,15] <- nrow(subset(xall,xall$genov%%10==0 & (xall$genoc%%10==3 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$genoc%%1000>30)) ) #& xall$genoMapGood>=3 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$TrancheMapmin2<99)) )
  #		writeTruthTableSum[i,16] <- nrow(subset(xall,xall$genov%%10==0 & (xall$genoc%%10==0 | xall$genoMapGood<3 | (xall$genoc%%10==1 & (xall$TrancheAlignmin2>=99 | xall$TrancheMapmin2>=99)) | (xall$genoc%%10==2 & (xall$TrancheSSEmin2>=95 | xall$TrancheABQDmin2>=95 | xall$TrancheAlignmin2>=95 | xall$TrancheMapmin2>=95)) | (xall$genoc%%10==3 & (xall$TrancheABQDmin2>=99 | xall$TrancheAlignmin2>=99 | xall$TrancheMapmin2>=99)))) )
  writeTruthTableSum[i,16] <- nrow(subset(xall,xall$genov%%10==0 & (xall$genoc%%10==0 | (xall$genoc%%10==1 & (xall$TrancheAlignmin2>=99 | xall$genoc%%1000<=30)) | (xall$genoc%%10==2 & (xall$TrancheSSEmin2>=95 | xall$TrancheABQDmin2>=95 | xall$TrancheAlignmin2>=95 | xall$TrancheMapmin2>=95 | xall$genoMapGood<3)) | (xall$genoc%%10==3 & (xall$TrancheABQDmin2>=99 | xall$TrancheAlignmin2>=99 | xall$genoc%%1000<=30)))) )
  
  
  print(i)
}
write.csv(writeTruthTableSum,file="plots/TruthTableSum_genoAllSeq_Sol_WEx_v1.2V_NoHomMapfilt.csv")

exportvcf <- function(tmpsel,fname,geno) {
  tmpsel$alt1=tmpsel$alt
  tmpsel$alt=tmpsel$ref
  tmpsel$ref="."
  tmpsel$qual=tmpsel$genoc
  tmpsel$genoc=tmpsel$chrompos%%1000000000
  tmpsel$chrompos=floor(tmpsel$chrompos/1000000000)
  tmpsel$filt="PASS"
  tmpsel$info="AC=2"
  tmpsel$format="GT"
  tmpsel$na12878=geno
  
  system(paste("cp /Volumes/SSD960/workspace/data/vcfannot/vcf.head",fname))
  write.table(tmpsel,file=fname,sep="\t",append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
}



##Count variants in each truth table category comparing vcfs to OMNI calls
datavcf <- dbGetQuery(con,paste("SELECT DISTINCT vcfhomhet.dataset, vcfhomhet.vcf FROM vcfhomhet",sep=""))
ndv <- nrow(datavcf)
writeTruthTableSum <- matrix(0,ndv+1,28)
valdataset = 13
xall <- dbGetQuery(con,paste("SELECT cast(validation.chrompos as real) as chrompos, validation.geno genoc,ref,alt FROM validation where dataset=",valdataset,sep=""))
for (i in (ndv+1):(ndv+1)) {
  if (i==ndv+1) {
    xv <- dbGetQuery(con,paste("SELECT cast(chrposall.chrompos as real) as chrompos, (chrposall.genoAll+chrposall.genoSeqV) genov, (chrposall.genoMapGood+chrposall.genoMapGoodSeqV) genoMapGood,TrancheSSEmin2,TrancheABQDmin2,TrancheAlignmin2,TrancheMapmin2 FROM chrposall",sep=""))
  } else {
    xv <- dbGetQuery(con,paste("SELECT cast(vcfhomhet.chrompos as real) as chrompos, vcfhomhet.geno genov FROM vcfhomhet WHERE vcfhomhet.vcf=",datavcf$vcf[i]," AND vcfhomhet.dataset=",datavcf$dataset[i],sep=""))
  }
  mm <- match(xall$chrompos,xv$chrompos)
  xall$genov <- NA
  xall$genov <- xv[mm,"genov"]
  if (i==ndv+1) { 
    xall$genoMapGood <- xv[mm,"genoMapGood"] 
    xall$TrancheSSEmin2 <- xv[mm,"TrancheSSEmin2"] 
    xall$TrancheABQDmin2 <- xv[mm,"TrancheABQDmin2"] 
    xall$TrancheAlignmin2 <- xv[mm,"TrancheAlignmin2"] 
    xall$TrancheMapmin2 <- xv[mm,"TrancheMapmin2"] 
  } else {
    xall$genoMapGood <- 3
    xall$TrancheSSEmin2 <- 0
    xall$TrancheABQDmin2 <- 0
    xall$TrancheAlignmin2 <- 0
    xall$TrancheMapmin2 <- 0
  }
  
  ##HomRef in VCF
  genos=c("0/0","0/1","1/1")
  genotxt=c("HomRef","Het","HomVar")
  for (j in 1:3) {
    tmpsel <- subset(xall[,c("chrompos","genoc","ref","alt")],((xall$genov%%10==1 & (xall$genov%%1000>30 | xall$genov==1)) | is.na(xall$genov)) & ((xall$TrancheAlignmin2<99) | is.na(xall$genoMapGood)) & (xall$genoc%%10==j))
    exportvcf(tmpsel,paste("/Volumes/SSD960/workspace/plots/TruthTableSumVal",valdataset,"_GenoAllSeq_Sol_WEx_v1.2V_NoHomMapfilt_val",genotxt[j],"_vcf",i,"HomRef.vcf",sep=""),genos[j])
    writeTruthTableSum[i,j] <- nrow(tmpsel) 
  }
  writeTruthTableSum[i,4] <- nrow(subset(xall,((xall$genov%%10==1 & (xall$genov%%1000>30 | xall$genov==1)) | is.na(xall$genov)) & ((xall$TrancheAlignmin2<99) | is.na(xall$genoMapGood)) & xall$genoc==4) )
  writeTruthTableSum[i,5] <- nrow(subset(xall,((xall$genov%%10==1 & (xall$genov%%1000>30 | xall$genov==1)) | is.na(xall$genov)) & ((xall$TrancheAlignmin2<99) | is.na(xall$genoMapGood)) & xall$genoc==5) )
  writeTruthTableSum[i,6] <- nrow(subset(xall,((xall$genov%%10==1 & (xall$genov%%1000>30 | xall$genov==1)) | is.na(xall$genov)) & ((xall$TrancheAlignmin2<99) | is.na(xall$genoMapGood)) & xall$genoc==6) )
  writeTruthTableSum[i,7] <- nrow(subset(xall,((xall$genov%%10==1 & (xall$genov%%1000>30 | xall$genov==1)) | is.na(xall$genov)) & ((xall$TrancheAlignmin2<99) | is.na(xall$genoMapGood)) & xall$genoc%%10==0) )
  
  ##Het in VCF
  for (j in 1:3) {
    tmpsel <- subset(xall[,c("chrompos","genoc","ref","alt")],xall$genov%%10==2 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95 & xall$genoMapGood>=3 & (xall$genoc%%10==j))
    exportvcf(tmpsel,paste("/Volumes/SSD960/workspace/plots/TruthTableSumVal",valdataset,"_GenoAllSeq_Sol_WEx_v1.2V_NoHomMapfilt_val",genotxt[j],"_vcf",i,"Het.vcf",sep=""),genos[j])
    writeTruthTableSum[i,j+7] <- nrow(tmpsel) 
  }
  writeTruthTableSum[i,11] <- nrow(subset(xall,xall$genov%%10==2 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95 & xall$genoMapGood>=3 & xall$genoc==4) )
  writeTruthTableSum[i,12] <- nrow(subset(xall,xall$genov%%10==2 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95 & xall$genoMapGood>=3 & xall$genoc==5) )
  writeTruthTableSum[i,13] <- nrow(subset(xall,xall$genov%%10==2 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95 & xall$genoMapGood>=3 & xall$genoc==6) )
  writeTruthTableSum[i,14] <- nrow(subset(xall,xall$genov%%10==2 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95 & xall$genoMapGood>=3 & xall$genoc%%10==0) )
  
  ##HomVar in VCF
  for (j in 1:3) {
    tmpsel <- subset(xall[,c("chrompos","genoc","ref","alt")],xall$genov%%10==3 & (xall$genov%%1000>30 | xall$genov==3) & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & (xall$genoc%%10==j))
    exportvcf(tmpsel,paste("/Volumes/SSD960/workspace/plots/TruthTableSumVal",valdataset,"_GenoAllSeq_Sol_WEx_v1.2V_NoHomMapfilt_val",genotxt[j],"_vcf",i,"HomVar.vcf",sep=""),genos[j])
    writeTruthTableSum[i,j+14] <- nrow(tmpsel) 
  }
  writeTruthTableSum[i,18] <- nrow(subset(xall,xall$genov%%10==3 & (xall$genov%%1000>30 | xall$genov==3) & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$genoc==4) )
  writeTruthTableSum[i,19] <- nrow(subset(xall,xall$genov%%10==3 & (xall$genov%%1000>30 | xall$genov==3) & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$genoc==5) )
  writeTruthTableSum[i,20] <- nrow(subset(xall,xall$genov%%10==3 & (xall$genov%%1000>30 | xall$genov==3) & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$genoc==6) )
  writeTruthTableSum[i,21] <- nrow(subset(xall,xall$genov%%10==3 & (xall$genov%%1000>30 | xall$genov==3) & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$genoc%%10==0) )
  
  ##Other in VCF
  for (j in 1:3) {
    tmpsel <- subset(xall[,c("chrompos","genoc","ref","alt")],(xall$genov%%10==0 | (xall$genov%%10==1 & (xall$TrancheAlignmin2>=99 | (xall$genov%%1000<=30 & xall$genov%%1000>10))) | (xall$genov%%10==2 & (xall$TrancheABQDmin2>=95 | xall$TrancheAlignmin2>=95 | xall$TrancheMapmin2>=95 | xall$genoMapGood<3)) | (xall$genov%%10==3 & (xall$TrancheABQDmin2>=99 | xall$TrancheAlignmin2>=99 | (xall$genov%%1000<=30 & xall$genov%%1000>10)))) & (xall$genoc%%10==j))
    exportvcf(tmpsel,paste("/Volumes/SSD960/workspace/plots/TruthTableSumVal",valdataset,"_GenoAllSeq_Sol_WEx_v1.2V_NoHomMapfilt_val",genotxt[j],"_vcf",i,"Other.vcf",sep=""),genos[j])
    writeTruthTableSum[i,j+21] <- nrow(tmpsel) 
  }
  writeTruthTableSum[i,25] <- nrow(subset(xall,(xall$genov%%10==0 | (xall$genov%%10==1 & (xall$TrancheAlignmin2>=99 | (xall$genov%%1000<=30 & xall$genov%%1000>10))) | (xall$genov%%10==2 & (xall$TrancheABQDmin2>=95 | xall$TrancheAlignmin2>=95 | xall$TrancheMapmin2>=95 | xall$genoMapGood<3)) | (xall$genov%%10==3 & (xall$TrancheABQDmin2>=99 | xall$TrancheAlignmin2>=99 | (xall$genov%%1000<=30 & xall$genov%%1000>10)))) & xall$genoc==4) )
  writeTruthTableSum[i,26] <- nrow(subset(xall,(xall$genov%%10==0 | (xall$genov%%10==1 & (xall$TrancheAlignmin2>=99 | (xall$genov%%1000<=30 & xall$genov%%1000>10))) | (xall$genov%%10==2 & (xall$TrancheABQDmin2>=95 | xall$TrancheAlignmin2>=95 | xall$TrancheMapmin2>=95 | xall$genoMapGood<3)) | (xall$genov%%10==3 & (xall$TrancheABQDmin2>=99 | xall$TrancheAlignmin2>=99 | (xall$genov%%1000<=30 & xall$genov%%1000>10)))) & xall$genoc==5) )
  writeTruthTableSum[i,27] <- nrow(subset(xall,(xall$genov%%10==0 | (xall$genov%%10==1 & (xall$TrancheAlignmin2>=99 | (xall$genov%%1000<=30 & xall$genov%%1000>10))) | (xall$genov%%10==2 & (xall$TrancheABQDmin2>=95 | xall$TrancheAlignmin2>=95 | xall$TrancheMapmin2>=95 | xall$genoMapGood<3)) | (xall$genov%%10==3 & (xall$TrancheABQDmin2>=99 | xall$TrancheAlignmin2>=99 | (xall$genov%%1000<=30 & xall$genov%%1000>10)))) & xall$genoc==6) )
  writeTruthTableSum[i,28] <- nrow(subset(xall,(xall$genov%%10==0 | (xall$genov%%10==1 & (xall$TrancheAlignmin2>=99 | (xall$genov%%1000<=30 & xall$genov%%1000>10))) | (xall$genov%%10==2 & (xall$TrancheABQDmin2>=95 | xall$TrancheAlignmin2>=95 | xall$TrancheMapmin2>=95 | xall$genoMapGood<3)) | (xall$genov%%10==3 & (xall$TrancheABQDmin2>=99 | xall$TrancheAlignmin2>=99 | (xall$genov%%1000<=30 & xall$genov%%1000>10)))) & xall$genoc%%10==0) )
  
  
  print(i)
}
write.csv(writeTruthTableSum,file=paste("/Volumes/SSD960/workspace/plots/TruthTableSumVal",valdataset,"_GenoAllSeq_Sol_WEx_v1.2V_NoHomMapfilt.csv",sep=""))
write.csv(datavcf,file=paste("plots/TruthTableSum_Datasetnums.csv",sep=""))


datavcfannot <- dbGetQuery(con,paste("SELECT DISTINCT vcfannot.bam FROM vcfannot",sep=""))
nf <- nrow(datavcfannot)

#calculate totals before and after arbitration
xall <- dbGetQuery(con,paste("SELECT chrposall.genoAll,chrposall.genoSeqV, (chrposall.genoMapGood+chrposall.genoMapGoodSeqV) genoMapGood, TrancheSSEmin2,TrancheABQDmin2,TrancheAlignmin2,TrancheMapmin2 FROM chrposall",sep=""))
nrow(subset(xall,(xall$genoSeqV%%10==1 & xall$genoSeqV%%1000>30 & xall$TrancheAlignmin2<99)))
#[1] 2819418
nrow(subset(xall,(xall$genoAll%%10==1 & xall$genoAll%%1000>30 & xall$TrancheAlignmin2<99)))
#[1] 629827
nrow(subset(xall,(xall$genoSeqV%%10==1 & xall$genoSeqV>15000 & xall$genoSeqV%%1000>30 & xall$TrancheAlignmin2<99)))
#[1] 2077317
2077317+629827
#[1] 2707144
nrow(subset(xall,(xall$genoSeqV%%10==2 & xall$genoSeqV>15000 & xall$genoMapGood>=3 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95)))
#[1] 93401
nrow(subset(xall,(xall$genoAll%%10==2 & xall$genoMapGood>=3 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95)))
#[1] 1979968
1979968+93401
#[1] 2073369
> 
  nrow(subset(xall,(xall$genoSeqV%%10==3 & xall$genoSeqV>15000 & xall$genoMapGood>=3 & xall$TrancheABQDmin2<99 &xall$TrancheAlignmin2<99 & xall$TrancheMapmin2<99 & xall$TrancheMapmin2<99)))
#[1] 144737
nrow(subset(xall,(xall$genoAll%%10==3 & xall$genoMapGood>=3 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$TrancheMapmin2<99 & xall$TrancheMapmin2<99)))
#[1] 1031937
1031937+144737
#[1] 1176674
521001+1979968+1031937
1971491+93401+144737
xall$genoc%%10==1 & xall$TrancheAlignmin2<99 & xall$genoc%%1000>30

#count variants arbitrated in different ways
xall <- dbGetQuery(con,paste("SELECT cast(chrposall.chrompos as real) as chrompos, chrposall.genoAll, (chrposall.genoAll+chrposall.genoSeqV) genoc, (chrposall.genoMapGood+chrposall.genoMapGoodSeqV) genoMapGood, TrancheSSEmin2,TrancheABQDmin2,TrancheAlignmin2,TrancheMapmin2 FROM chrposall",sep=""))
nrow(subset(xall,xall$genoc<1000 & xall$genoAll==0 & ((xall$genoc%%10==1 & xall$TrancheAlignmin2<99 & xall$genoc%%1000>30) | (xall$genoc%%10==2 & xall$genoMapGood>=3 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95) | (xall$genoc%%10==3 & xall$genoMapGood>=3 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$TrancheMapmin2<99))))
#[1] 809384
nrow(subset(xall,xall$genoc>1000 & xall$genoc<2000 &xall$genoAll==0 & ((xall$genoc%%10==1 & xall$TrancheAlignmin2<99 & xall$genoc%%1000>30) | (xall$genoc%%10==2 & xall$genoMapGood>=3 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95) | (xall$genoc%%10==3 & xall$genoMapGood>=3 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$TrancheMapmin2<99))))
#[1] 113888
nrow(subset(xall,xall$genoc>2000 & xall$genoc<3000 &xall$genoAll==0 & ((xall$genoc%%10==1 & xall$TrancheAlignmin2<99 & xall$genoc%%1000>30) | (xall$genoc%%10==2 & xall$genoMapGood>=3 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95) | (xall$genoc%%10==3 & xall$genoMapGood>=3 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$TrancheMapmin2<99))))
#[1] 113950
nrow(subset(xall,xall$genoc>10000 &xall$genoc<11000 & xall$genoAll==0 & ((xall$genoc%%10==1 & xall$TrancheAlignmin2<99 & xall$genoc%%1000>30) | (xall$genoc%%10==2 & xall$genoMapGood>=3 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95) | (xall$genoc%%10==3 & xall$genoMapGood>=3 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$TrancheMapmin2<99))))
#[1] 68408
nrow(subset(xall,xall$genoc>11000 & xall$genoc<12000 &xall$genoAll==0 & ((xall$genoc%%10==1 & xall$TrancheAlignmin2<99 & xall$genoc%%1000>30) | (xall$genoc%%10==2 & xall$genoMapGood>=3 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95) | (xall$genoc%%10==3 & xall$genoMapGood>=3 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$TrancheMapmin2<99))))
#[1] 14181
nrow(subset(xall,xall$genoc>12000 & xall$genoc<13000 &xall$genoAll==0 & ((xall$genoc%%10==1 & xall$TrancheAlignmin2<99 & xall$genoc%%1000>30) | (xall$genoc%%10==2 & xall$genoMapGood>=3 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95) | (xall$genoc%%10==3 & xall$genoMapGood>=3 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$TrancheMapmin2<99))))
#[1] 61639
#>  809384+68408
#[1] 877792
#>  113888+14181
#[1] 128069
#>  113950+61639
#[1] 175589
nrow(subset(xall,xall$genoAll>0 & ((xall$genoc%%10==1 & xall$TrancheAlignmin2<99 & xall$genoc%%1000>30) | (xall$genoc%%10==2 & xall$genoMapGood>=3 & xall$TrancheSSEmin2<95 & xall$TrancheABQDmin2<95 & xall$TrancheAlignmin2<95 & xall$TrancheMapmin2<95) | (xall$genoc%%10==3 & xall$genoMapGood>=3 & xall$TrancheABQDmin2<99 & xall$TrancheAlignmin2<99 & xall$TrancheMapmin2<99))))

877792+128069+175589

dbDisconnect(con)
