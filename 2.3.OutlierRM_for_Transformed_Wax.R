###############################################
# for transformed waxes
###############################################
setwd("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/Wax_TWAS/GAPIT")
#setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT")

wax<-read.table("Transformed_BLUP_wax_imputed.txt",header=T,row.names=1,sep="\t")
wax<-cbind.data.frame(rownames(wax),wax)
colnames(wax)[1]<-"MLC_STANDARD"

taxa<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/taxa310_pheno_geno_rna.txt")
#taxa<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/taxa310_pheno_geno_rna.txt")
taxa<-as.character(taxa[,1])
wax<-wax[which(wax$MLC_STANDARD %in% taxa),]
################################################
# outlier removal
################################################
source("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/script/outlier_removal_function_Meng.R")

nbvariables=1
nbtraits=ncol(wax)-1
random=c()
fixed=c()
BLUx="BLUE" # because the model used for outlier removal only contains grand mean, and there is no random factors in the model

clean.dataset <- wax[,c(1:nbvariables)]
alltraitnames <- vector()
#for (i in 1:10){
for (i in 1:nbtraits){
  curr.trait <- colnames(wax[ (nbvariables + i) ])
  alltraitnames <- c(alltraitnames, curr.trait)
  
  transfpheno=wax[,c(1:nbvariables,i+nbvariables)]
  initial.outliers(transfpheno, curr.trait, random, fixed, nbvariables, BLUx) -> out
  as.matrix(out) -> out
  colnames(out) <- curr.trait
  
  if(i == 1) {
    cleanedpheno <- cbind(clean.dataset, out)
  } else {
    cleanedpheno <- cbind(cleanedpheno, out)
  }
  
}
cleanedpheno<-as.data.frame(cleanedpheno)
cleanedpheno[,1]<-wax[,1]
colnames(cleanedpheno)[1]<-"MLC_STANDARD"
write.table(cleanedpheno,"Transformed_BLUP_wax_imputed_OLRM.txt",col.names=T,row.names=F,sep="\t",quote=F)
write.table(wax,"Transformed_BLUP_wax_imputed_310.txt",col.names=T,row.names=F,sep="\t",quote=F) ## for supplemental table


####### check traits that were changed
for(i in 2:ncol(cleanedpheno)){
  temp<-cleanedpheno[,i]
  wax_name<-colnames(cleanedpheno)[i]
  nm_na<-length(temp[is.na(temp)])
  if(nm_na>0){
    print(paste(wax_name,": ",nm_na,sep=""))
  }
  
}

# FFA_C36_tr: 1
# HC_C23_tr: 1
# HC_C31_tr: 2
# HC_C37_tr: 1
# WE_C50_tr: 1
# AC_Unk4_tr: 2
########################################
