#### check GC-FID design
GCdesign<-read.table("/Users/menglin/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/raw_data/design_machine_column.txt",
                     header=T,sep="\t",comment.char="",stringsAsFactors=F)
# GCdesign<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/raw_data/design_machine_column.txt",
#                      header=T,sep="\t",comment.char="",stringsAsFactors=F)

## remove "Fresh", keep Frozen only
GCdesign<-GCdesign[-grep("Fresh",GCdesign$Barcode_temp),]
## plot 228 keep the "rerun" sample in batch 56?57?
GCdesign<-GCdesign[-intersect(grep("228",GCdesign$Barcode_temp),which(GCdesign$batch=="52")),]

col_nm<-unique(GCdesign$column_ser_No)
for (c in col_nm){
  batch<-unique(GCdesign$batch[which(GCdesign$column_ser_No==c)])
  print(paste(c,": ",paste(batch,collapse = ","),sep=""))
}
# [1] "SN 1579334: 1,2,3,4,5,6,7,9,11,13,15,19,21"
# [1] "SN 1579331: 8,10,12,14,16,17,18,20"
# [1] "SN 1579333: 22,24,25,26,27,28,29,31,32,37,38,41,44,45,48,49,50,51,52"
# [1] "SN 1579332: 23"
# [1] "SN 1626961: 30,33,34,35,36,37,38,39,40,42,43,46,47,64,65,66,67,75,76"
# [1] "S/N 1647227: 53,54,55,56,57,58,59,60,61,62,63,68,69,70,71,72,73,74,75"
# [1] "SN 1644772: 77,79,81,83,85,87,89,91,93,95,97"
# [1] "SN 1644771: 78,80,82,84,86,88,90,92,94,96,98,99,100,101,102"

Barcode<-strsplit(GCdesign$Barcode_temp,split=".",fixed=T)
Barcode<-as.data.frame(matrix(unlist(Barcode),ncol=2,byrow=T))
Barcode[,1]<-as.character(Barcode[,1])
for (i in 1:nrow(Barcode)){
  if(nchar(Barcode[i,1])==3){
    Barcode[i,1]<-paste("18SD",Barcode[i,1],sep="")
  }else{
    Barcode[i,1]<-paste("18SD",paste(rep("0",3-nchar(Barcode[i,1])),collapse=""),Barcode[i,1],sep="")
  }
  
}
GCdesign$Barcode<-Barcode[,1]
for (bar in unique(GCdesign$Barcode)){
  sub<-GCdesign[which(GCdesign$Barcode==bar),]
  if (length(unique(sub$batch))>1 | length(unique(sub$run_date))>1|length(unique(sub$machine))>1|length(unique(sub$column_ser_No))>1){
    print(bar)
  }
}
# [1] "18SD396" # duplicated samples
# [1] "18SD110" # overnight run, okay
# [1] "18SD513" # overnight run, okay

GCdesign<-GCdesign[-intersect(which(GCdesign$run_date=="27-Jun-19"),which(GCdesign$Barcode_temp=="110.12 from seq 53")),]
GCdesign<-GCdesign[-intersect(which(GCdesign$batch==67),grep("396",GCdesign$Barcode_temp,fixed=T)),] # correspondingly, only keep plot 396 in row 605 for wax mean (remove row 382)
GCdesign<-GCdesign[-intersect(which(GCdesign$batch==58),grep("513",GCdesign$Barcode_temp,fixed=T)),]

##
GCdesign<-GCdesign[,-(1:2)]
uniq_barcode<-unique(GCdesign$Barcode)
GC_design<-GCdesign[match(uniq_barcode,GCdesign$Barcode),]

setwd("/Users/menglin/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/raw_data")
#setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/raw_data")
write.table(GC_design,"GC_design_machine_column_cleaned.txt",col.names=T,row.names=F,sep="\t",quote=F)


## check col-No within rep: SN 1579333 and SN 1626961 were used across rep!! not desirable...
unique(GCdesign$column_ser_No[which(GCdesign$Barcode %in% paste("18SD",1:342,sep=""))])
#"SN 1579334"  "SN 1579331"  "SN 1579333"  "SN 1579332"  "SN 1626961"  "S/N 1647227"
unique(GCdesign$column_ser_No[which(GCdesign$Barcode %in% paste("18SD",343:684,sep=""))])
#"SN 1579332"  "SN 1579333"  "S/N 1647227" "SN 1626961"  "SN 1644772"  "SN 1644771" 

##################################
# incorperate field design and wax
###################################
setwd("/Users/menglin/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/raw_data")
# setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/raw_data")
# setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/raw_data")
# setwd("/Users/zhenghao/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/raw_data")

## plot mean from Dropbox, chemicals with only 0's have been removed
wax<-read.table("panel_wax_line_mean.txt",header=T,sep="\t")
wax<-wax[-382,] ## remove one duplicated sample for plot 396
wax$Plot<-as.character(wax$Plot)

Barcode<-vector()
for (i in 1:nrow(wax)){
  if(nchar(wax$Plot[i])==3){
    barcode<-paste("18SD",wax$Plot[i],sep="")
  }else{
    barcode<-paste("18SD",paste(rep("0",3-nchar(wax$Plot[i])),collapse=""),wax$Plot[i],sep="")
  }
  Barcode<-c(Barcode,barcode)
}
wax<-cbind(Barcode,wax)
wax<-wax[,-2]
wax$Barcode[duplicated(wax$Barcode)] #18SD133 18SD588!!! dropbox summary table: second sheet, two errors, corrected according to the first sheet.

## calculate category total
FA<-wax[,substring(colnames(wax),1,3)=="FA_"]
FA$Total<-apply(FA,1,sum)

FFA<-wax[,substring(colnames(wax),1,4)=="FFA_"]
FFA$Total<-apply(FFA,1,sum)

HC<-wax[,substring(colnames(wax),1,3)=="HC_"]
HC$Total<-apply(HC,1,sum)

AD<-wax[,substring(colnames(wax),1,3)=="AD_"]
AD$Total<-apply(AD,1,sum)

WE<-wax[,substring(colnames(wax),1,3)=="WE_"]
WE$Total<-apply(WE,1,sum)

AC<-wax[,substring(colnames(wax),1,3)=="AC_"]
AC$Total<-apply(AC,1,sum)

wax_total<-apply(wax[,-1],1,sum)

## add the total values to wax dataset
wax$FA_total<-FA$Total
wax$FFA_total<-FFA$Total
wax$HC_total<-HC$Total
wax$AD_total<-AD$Total
wax$WE_total<-WE$Total
wax$AC_total<-AC$Total
wax$Total<-wax_total


##### GC design ########

GC<-read.table("GC_design_machine_column_cleaned.txt",header=T,sep="\t",comment.char="")
GC$order_within_batch<-as.character(GC$order_within_batch)
GC$order_within_batch<-substr(GC$order_within_batch,2,nchar(GC$order_within_batch))

##### field design ########
#FieldDesign<-read.table("/Users/menglin/Desktop/Labserver/OfficeCmp/GoogleDrive/TWAS_2018/ExpDesign/SD18_Design_Chk_Barcode_forBLUP.txt",header=T,sep="\t")
FieldDesign<-read.table("/home/ml2498/Desktop/Labserver/OfficeCmp/GoogleDrive/TWAS_2018/ExpDesign/SD18_Design_Chk_Barcode_forBLUP.txt",header=T,sep="\t")
FieldDesign<-FieldDesign[,c(2,3,5:7)]


Field_GC<-merge(FieldDesign,GC,by="Barcode",all.y=T)

colnm<-which(colnames(Field_GC)=="MLC_STANDARD")
for (i in 1:dim(Field_GC)[1]){
  if (Field_GC[i,colnm]=="MO17") {Field_GC$CHECK[i]<-1}
  else Field_GC$CHECK[i]<-99
}
for (i in 1:dim(Field_GC)[1]){
  if (Field_GC[i,colnm]=="MO17") {Field_GC$IS_EXEXPERIMENTAL[i]<-0
  } else {Field_GC$IS_EXEXPERIMENTAL[i]<-1
  }
}

Field_GC$Year<-18
Field_GC$LOC<-"SD"

colnames(Field_GC)<-c("Barcode","COL","BLOCK","ENV","MLC_STANDARD","Batch","OrderWithinBatch","RunDate","Machine",
                       "ColSerNo","CHECK","IS_EXPERIMENTAL","Year","LOC")
Field_GC$COL1<-Field_GC$COL
## SD18 only
Field_GC$COL1[which(Field_GC$COL>9)]<-19-Field_GC$COL[which(Field_GC$COL>9)]

Filed_GC_wax<-merge(Field_GC,wax,by="Barcode",all=T)

setwd("/Users/menglin/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno")
# setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno")
write.table(Filed_GC_wax,"wax_BLUPinput_02102020.txt",col.names=T,row.names=F,sep="\t",quote=FALSE)
#############################################


###### impute zero values (within ENV) ###
rawwax<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno/wax_BLUPinput_02102020.txt",header=T,sep="\t")
# rawwax<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno/wax_BLUPinput_02102020.txt",header=T,sep="\t")

rawwax$WE_C56[which(rawwax$WE_C56==0)]<-NA # special process for WE_C56 for now
rawwax$WE_C40[which(rawwax$WE_C40==0)]<-NA
rawwax$WE_C41[which(rawwax$WE_C41==0)]<-NA

## remove componds with number of zero's > 40% of the population
rm_wax<-c(63) # col number for Campesterol, not stable, Isabel suggested to remove it
for (j in 2:ncol(rawwax)){
  nm_zero<-length(rawwax[which(rawwax[,j]==0),j])
  if(nm_zero>nrow(rawwax)*0.4){
    rm_wax<-c(rm_wax,j)
  }
}
rawwax<-rawwax[,-rm_wax]

imputed<-vector()
for (i in 1:2){
  raw_sub<-rawwax[which(rawwax$ENV==i),]
  
  for (w in 16:68){ # the remaining AC columns, because they contains zero's
    if (min(raw_sub[,w],na.rm=T)==0){
      
      wax_temp<-raw_sub[,w]
      wax_temp[which(wax_temp==0)]<-NA
      
      print(colnames(raw_sub)[w])
      print(length(wax_temp[is.na(wax_temp)]))
      
      ceil<-min(wax_temp,na.rm=T)
      print(ceil)
      raw_sub[is.na(wax_temp),w]<-runif(length(wax_temp[is.na(wax_temp)]), min = 0, max = ceil)
    }
  }
  imputed<-rbind(imputed,raw_sub)
}
setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno")
# setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno")
write.table(imputed,"wax_imputed_BLUPinput_02192020.txt",col.names=T,row.names=F,sep="\t",quote=FALSE)

#for(i in 61:68){
#  hist(rawwax[,i],breaks=50,main=colnames(rawwax)[i])
#}

## WE_C41&40&56, AC_F, AC_unk need GWAS





# vidualize raw wax data
Filed_GC_wax<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno/wax_BLUPinput_02102020.txt",header=T,sep="\t")
# Filed_GC_wax<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno/wax_BLUPinput_02102020.txt",header=T,sep="\t")
setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/plot")
pdf("Dist_raw_wax.pdf")
for (i in 16:ncol(Filed_GC_wax)){
  hist(Filed_GC_wax[,i],main=colnames(Filed_GC_wax)[i])
}
dev.off()

#### check distribution of zero's for each wax
sub1<-Filed_GC_wax[which(Filed_GC_wax$ENV==1),];length(unique(sub1$MLC_STANDARD)) # 306
sub2<-Filed_GC_wax[which(Filed_GC_wax$ENV==2),] ;length(unique(sub2$MLC_STANDARD))# 316
common<-intersect(as.character(sub1$MLC_STANDARD),as.character(sub2$MLC_STANDARD)) # 298 in both env

len_zero1<-vector()
for (i in 16:79){
  zero<-length(which(sub1[,i]==0))
  len_zero1<-c(len_zero1,zero)
}
len_zero1<-as.data.frame(len_zero1)
rownames(len_zero1)<-colnames(sub1)[16:79]

len_zero2<-vector()
for (i in 16:79){
  zero<-length(which(sub2[,i]==0))
  len_zero2<-c(len_zero2,zero)
}
len_zero2<-as.data.frame(len_zero2)
rownames(len_zero2)<-colnames(sub2)[16:79]




##### check machine effect
Filed_GC_wax<-Filed_GC_wax[order(Filed_GC_wax$ENV,Filed_GC_wax$Machine,Filed_GC_wax$ColSerNo),]
Filed_GC_wax$Barcode <- factor(Filed_GC_wax$Barcode, levels = Filed_GC_wax$Barcode)
library(ggplot2)

pdf("col_machine_effect_on_wax.pdf",height=4,width=7)
for(i in 16:ncol(Filed_GC_wax)){
  print(ggplot(Filed_GC_wax, aes(x = Barcode, y = Filed_GC_wax[,i], color=ColSerNo,shape=Machine,group=interaction(ColSerNo,Machine)))+
          geom_point()+
          theme(axis.text.x = element_blank())+
          labs(y = colnames(Filed_GC_wax)[i])
  )
  
}
dev.off()

pdf("col_machine_effect_on_wax_boxplot.pdf",height=4,width=7)
for(i in 16:ncol(Filed_GC_wax)){
  print(ggplot(Filed_GC_wax, aes(x = ColSerNo, y = Filed_GC_wax[,i], fill=ColSerNo))+
          geom_boxplot()+
          #theme(axis.text.x = element_blank())+
          labs(y = colnames(Filed_GC_wax)[i])
  )
  
}
dev.off()





