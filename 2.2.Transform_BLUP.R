#### transform, record lamda

## format wax BLUPs from wrapper

file<-read.table("/Users/menglin/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno/wrapper/result_wax/BLUPS_all_traits_wax_imputed_BLUPinput_02192020.txt",header=T,sep="\t")
#file<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno/wrapper/result_wax/BLUPS_all_traits_wax_imputed_BLUPinput_02192020.txt",header=T,sep="\t")
file$predict<-as.character(file$predict)

index<-grep("MLC_STANDARD_[a-zA-Z0-9_.-]*:IS_EXPERIMENTAL$",file$predict)
blup<-file[index,]
taxa<-strsplit(blup$predict,split=":",fixed = T)
taxa<-matrix(unlist(taxa),ncol=2,byrow=T)
taxa<-substr(taxa[,1],14,nchar(taxa[,1]))
blup<-cbind.data.frame(taxa,blup[,-1])
colnames(blup)[1]<-"MLC_STANDARD"

for (i in 2:ncol(blup)){
  if(sum(blup[,i]<0)){
    print(i)
  }
}
# col 50, 51 and 57 still have negative values
blup[,49]<-blup[,49]+2*abs(min(blup[,49],na.rm=T))
blup[,50]<-blup[,50]+2*abs(min(blup[,50],na.rm=T))
blup[,53]<-blup[,53]+2*abs(min(blup[,53],na.rm=T))
## double check
for (i in 2:ncol(blup)){
  if(sum(blup[,i]<0)){
    print(i)
  }
}


blup<-blup[-which(blup[,1]=="MO17"),] ## remove blup

### box-cox transfromation
library(MASS)

trans_info<-data.frame(matrix(nrow=(ncol(blup)-1),ncol=2),stringsAsFactors=FALSE)
colnames(trans_info)<-c("expt","lambda")
blup_tr<-blup$MLC_STANDARD

for (i in 2:ncol(blup)){
  trans_info$expt[i-1]=colnames(blup)[i]
  
  temp<-blup[,i]
  trans <- boxcox(temp ~ 1)
  
  lambda <- trans$x[which.max(trans$y)]
  easy<-seq(-3,3, by=0.5)
  lambda.1<-easy[which.min(abs(easy - lambda))]
  trans_info$lambda[(i-1)]=lambda.1
  if (lambda.1!=0){
    temp_tr <-(temp^lambda.1-1)/lambda.1
    
  }else {
    temp_tr <- log(temp)
  }
  blup_tr<-cbind.data.frame(blup_tr,temp_tr) 
  colnames(blup_tr)[i]<-paste(colnames(blup)[i],'tr',sep="_")
}
colnames(blup_tr)[1]<-"MLC_STANDARD"
setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT")
write.table(trans_info,"Transform_info_lamda.txt",col.names=T,row.names=F,sep="\t",quote=F)
write.table(blup_tr,"Transformed_BLUP_wax_imputed.txt",col.names=T,row.names=F,sep="\t",quote=F)

### plot transformed wax blups
setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/plot")
pdf("Dist_wax_BLUPs_trans.pdf",height=3,width=3)
for (i in 2:ncol(blup_tr)){
  hist(blup_tr[,i],main=colnames(blup_tr)[i])
}
dev.off()




