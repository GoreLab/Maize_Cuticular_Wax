###############################################################
###############################################################
# HapMap3 GWAS
###############################################################
################################################################
chr=2 # chr=(2:12) chr 1 is split into 11 and 12

#trait<-c(1,2:13) #cbsulm16 & 12 PCs
#trait<-c(1,14:25) #cbsulm17
#trait<-c(1,26:37) #cbsulm18
#trait<-c(1,38:49) #cbsulm20
#trait<-c(1,50,59:65) #gore02
#trait<-c(1,11:12)# for model testing with FT as covariate: FFA_C20, FFA_C22


#trait<-c(1,33,34,46,49:54) # traits that need to be updated after imputation

setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT")
#setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT")
#setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT")

#pheno.all<-read.table("Transformed_BLUP_wax_imputed.txt",header=T,sep="\t")
pheno.all<-read.table("Transformed_BLUP_wax_imputed_OLRM.txt",header=T,sep="\t")
#pheno.all<-read.table("PCs_from_53waxes.txt",header=T,sep="\t")

olrm<-which(colnames(pheno.all) %in% c("FFA_C36_tr","HC_C23_tr","HC_C31_tr","HC_C37_tr","WE_C50_tr","AC_Unk4_tr"))
trait<-c(1,olrm)# outlier removed traits: 


## if need FT as covariate, use the following line
#FT<-read.table("FT_as_covariate.txt",header=T,sep="\t")

#####################################################
## no pruning, 450 imputation, 310 filter
geno.all<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/Hapmap3/UpliftTo_AGPv4/Imputed_GBS/forGAPIT/MLC_GBSSNP310_468K_v4_chr",chr,"_filter.hmp.txt",sep=""),
                    header=F,sep="\t",comment.char="")

## no pruning, 450 imputation, 450 filter
#geno.all<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/Hapmap3/UpliftTo_AGPv4/Imputed_GBS/forGAPIT/MLC_GBSSNP450_468K_v4_chr",chr,"_filter.hmp.txt",sep=""),
#                     header=F,sep="\t",comment.char="")


myG<-geno.all

#trait<-c(2,12:20)# 09242018 CE rate, just AZ, SD and all4
#trait<-c(2,12:19) # 12:15 for single env, 18:19 for 16 combined and 17 combined

myY<-pheno.all[,trait]
##### if log CE rate ###
#trait<-c(21:23)
#myY<-cbind(pheno.all[,1],log(pheno.all[,trait]))
#colnames(myY)[1]<-"Taxa"
#######################

#myKI<-read.table("centeredIBS_HMP3_AGPv4_mac_LD02.txt")
#myKI<-read.table("centeredIBS_HMP3_AGPv4_mac_LD02_320lines.txt")

#setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT")
#myKI<-read.table("centeredIBS_HMP3_AGPv4_LD02_310lines.txt")
#myKI<-read.table("centeredIBS_HMP3_AGPv4_LD01_310lines.txt")
#myKI<-read.table("centeredIBS_HMP3_AGPv4_LD02_310lines_1kNE.txt")
#myKI<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT/kin_VR/prun02/GAPIT.Kin.VanRaden.csv",sep=",")

## 450 imputation, 310 filtering
myKI<-read.table("centeredIBS_HMP3_AGPv4_LD02_450to310.txt")
## 450 imputation, 450 filtering
#myKI<-read.table("centeredIBS_HMP3_AGPv4_LD02_450.txt")

# source("http://www.bioconductor.org/biocLite.R") 
# biocLite("multtest")
# install.packages("gplots") 
# install.packages("LDheatmap") 
# install.packages("genetics")
# install.packages("EMMREML") 
# install.packages("scatterplot3d") #The downloaded link at: http://cran.r-project.org/package=scatterplot3d

library(multtest) 
library(gplots) 
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler) #this library is already installed in R 
library("scatterplot3d")
#source("http://zzlab.net/GAPIT/gapit_functions.txt")
#source("/Users/Meng/Google Drive/MLC_AZ_2017/GAPIT/gapit_functions.R")
source("gapit_functions.R") # 2017 GAPIT
#source('http://www.zzlab.net/GAPIT/previous/gapit_functions20190714.txt') ## 2018 GAPIT, from Di
source("http://zzlab.net/GAPIT/emma.txt")

################################
# re-direct to folders of chr
###############################
#mainDir<-"/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/imputed_wax"
#mainDir<-"/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/" # cbsu servers
#mainDir<-"/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/PC12"
#mainDir<-"/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/adjusted_models"
mainDir<-"/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/OLRM_traits"

subDir<-paste("chr",chr,sep="")
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))

log_file <- paste('SAM_log_GAPIT_chr',chr,'.txt',  sep = '')
sink(log_file)
#### without FT ######################################
nm_ind<-nrow(myKI)

myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  #CV=myCV, ## to be ditermined
  KI=myKI,
  #kinship.algorithm="Zhang",
  group.from=nm_ind, # 353 for (lsa adjusted ce;RNA), 359 for (lsa adjusted ce;GBS)
  group.to=nm_ind,
  group.by=1,
  Major.allele.zero=T,
  Model.selection=FALSE
) #automatically include best number of K groups

##############################


#########################################
#FT[is.na(FT[,2]),2]<-mean(FT[,2],na.rm=T)
#FT[is.na(FT[,3]),3]<-mean(FT[,3],na.rm=T)

##### with FT ######################################

myCV=FT[,c(1,4)]

nm_ind<-nrow(myKI)

myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  CV=myCV, ## to be ditermined
  KI=myKI,
  #kinship.algorithm="Zhang",
  group.from=nm_ind, # 353 for (lsa adjusted ce;RNA), 359 for (lsa adjusted ce;GBS)
  group.to=nm_ind,
  group.by=1,
  Major.allele.zero=T,
  Model.selection=FALSE
) #automatically include best number of K groups

####################################
# stacking result and Manhattan plot
####################################
library(qqman)
library(qvalue)
library(data.table)

## waxes
#Files<-read.table("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/pheno/wax_names.txt",sep="\t",stringsAsFactors=F)
#Files<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno/wax_names_subDir.txt",header=T,sep="\t",stringsAsFactors=F)
#Files<-t(Files[1,]) #64 wax names

Files<-read.table("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/pheno/wax_names_subDir.txt",header=T,sep="\t",stringsAsFactors=F)
Files<-Files[!is.na(Files[,2]),]


#### ****only for six traits with outlier removed**** ######
Files<-Files[which(Files$subDir=="OLRM_traits"),]
#################################################### 
  

## PCs
#Files<-paste("PC",1:12,sep="")

###
path="/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/"

#for(j in 1:12){ ## 12 PCs
#for (j in c(32,33,45,49:50,52,53,56,57)){ # imputed wax components
for (j in 1:nrow(Files)){
#for (j in c(25:48)){
#for (j in c(1:49,58:64)){
  f=Files[j,1]
  print(f)
  # if (j %in% 1:12){
  #   server="cbsulm16/"
  #   
  # }else if (j %in% 13:24){
  #   server="cbsulm17/"
  # }else if (j %in% 25:36){
  #   server="cbsulm18/"
  # }else if (j %in% 37:48){
  #   server="cbsulm20/"
  # }else if (j %in% c(49,58:64)){
  #   server="GORE/"
  # }
  #server="imputed_wax/"
  #server="PC12/"
  server=as.character(Files[j,2])
  
  results<-matrix(nrow=0,ncol=10)
  
  for (i in 2:12){
    
    setwd(paste(path,server,"/chr",i,sep=""))
  
    #file<-paste("GAPIT.MLM.ce_",f,"_untr.GWAS.Results.csv",sep="")
    #file<-paste("GAPIT.CMLM.ce_",f,"_untr.GWAS.Results.csv",sep="")
    file<-paste("GAPIT.MLM.",f,"_tr.GWAS.Results.csv",sep="")
    #gemma<-read.csv(file,header=T)
    #file<-paste("GAPIT.MLM.",f,".GWAS.Results.csv",sep="")
    gemma<-fread(file,data.table=F)
    
    if(i ==11 | i==12){
      gemma$Chromosome<-1
    }else{
      gemma$Chromosome<-i
    }
    
    results<-rbind(results,gemma)
  }
  results<-results[order(results$Chromosome,results$Position),]
  
  setwd(paste(path,server,sep=""))
  
  qobj <- qvalue(p = results$P.value)
  qvalues <- qobj$qvalues
  results$qvalues<-qvalues
  #results<-results[order(results$qvalues),]
 
   ## pruning
  GI.MP=results
  GI.MP=GI.MP[order(GI.MP$P.value),]
  topSNP<-GI.MP[1:5000,] # top 5000, ~0.1%
  topSNP.2<-GI.MP[1:ceiling(nrow(GI.MP)*0.1),] # top 10% 
  topSNP.3<-GI.MP[1:ceiling(nrow(GI.MP)*0.0001),] # top 0.01%
  
  p0005<--log10(GI.MP$P.value[ceiling(nrow(GI.MP)*0.00005)])
  p0002<--log10(GI.MP$P.value[ceiling(nrow(GI.MP)*0.00002)])
  p0001<--log10(GI.MP$P.value[ceiling(nrow(GI.MP)*0.00001)])
  
  #restSNP<-GI.MP[-(1:5000),]
  #set.seed(89898)
  #keep<-sample(1:nrow(restSNP),round(nrow(restSNP)*0.1,0),replace=F)
  #restSNP<-restSNP[keep,]
  #GI.MP.pruning<-rbind(topSNP,restSNP)
  GI.MP.pruning<-GI.MP[which(GI.MP$P.value<0.01),]
    
  # output for top SNPs
  #write.table(topSNP,paste("CE_",f,"_HapMap3_topSNPs.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F) # for making tables
  #write.table(topSNP.2,paste("CE_",f,"_HapMap3_topSNPs01.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F) # for finsher's combined test

  #write.table(topSNP,paste(f,"_tr_HapMap3_topSNPs.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F) # for making tables
  #write.table(topSNP.2,paste(f,"_tr_HapMap3_topSNPs01.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F) # for finsher's combined test
  #write.table(topSNP.3,paste(f,"_tr_HapMap3_topSNPs_forCand.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F) # for making tables

  write.table(topSNP,paste(f,"_tr_HapMap3_topSNPs.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F) # for making tables
  write.table(topSNP.2,paste(f,"_tr_HapMap3_topSNPs01.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F) # for finsher's combined test
  write.table(topSNP.3,paste(f,"_tr_HapMap3_topSNPs_forCand.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F) # for making tables
  
  #}
  ### plot
  if (length(results$P.value[which(results$qvalues<0.1)])==0){
    #pdf(paste("CE_",f,"_HapMap3_man_01.pdf",sep=""),width=8,height=4)  # "_1": 1% SNPs are highlighted; "_01": 5000 SNPs are highlighted
    #pdf(paste(f,"_tr_man_SNPforCand.pdf",sep=""),width=8,height=4)
    png(paste(f,"_tr_man_SNPforCand.png",sep=""), width=12,height=5,units="in",res=300,family="serif")
    manhattan(GI.MP.pruning, chr = "Chromosome", bp = "Position", p = "P.value", snp = "SNP",
            col = c("navy", "darkorange1"), chrlabs = NULL,
            suggestiveline = FALSE,genomewideline =FALSE, 
            highlight = as.character(topSNP.3$SNP),
            main=paste("GWAS ",f," HMP3",sep=""))
    dev.off()
  }else {
    threshold1<-(max(GI.MP$P.value[which(GI.MP$qvalues<0.1)])+min(GI.MP$P.value[which(GI.MP$qvalues>0.1)]))/2
    try(threshold2<-(max(GI.MP$P.value[which(GI.MP$qvalues<0.05)])+min(GI.MP$P.value[which(GI.MP$qvalues>0.05)]))/2)
  
    #pdf(paste("CE_",f,"_HapMap3_man_01.pdf",sep=""),width=8,height=4) # "_1": 1% SNPs are highlighted; "_01": 5000 SNPs are highlighted
    #pdf(paste(f,"_tr_man_SNPforCand.pdf",sep=""),width=8,height=4)
    png(paste(f,"_tr_man_SNPforCand.png",sep=""), width=12,height=5,units="in",res=300,family="serif")
    manhattan(GI.MP.pruning, chr = "Chromosome", bp = "Position", p = "P.value", snp = "SNP",
            col = c("navy", "darkorange1"), chrlabs = NULL,
            suggestiveline = -log10(threshold2), 
            genomewideline = -log10(threshold1),
            highlight = as.character(topSNP.3$SNP),
            main=paste("GWAS ",f," HMP3",sep=""))
    dev.off()
  }
  #jpeg(paste("CE_",f,"_HapMap3_qq.jpg",sep=""))
  jpeg(paste(f,"_tr_HapMap3_qq.jpg",sep=""))
  qq(results$P.value)
  dev.off()

}


