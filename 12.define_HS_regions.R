library(data.table)
setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/GWAS_hotspot")

file<-"GWAS_peak_counts_byCategory_WS2e+05_SS40000_04042022";ws<-200000 ## use this one to define hotspots
REGION<-fread(paste(file,".txt",sep=""),data.table=F)

#check<-REGION[which(REGION$AC>0),]
check<-REGION[which(REGION$FA>0),]


i=4
HS<-c()
for(i in 4:10){
  cate<-colnames(REGION)[i]
  
  if(i==10){
    threshold<-4 ## across category
  }else{
    threshold<-3  ## within category
  }
  
  j=1
  while (j<nrow(REGION)){
    if(REGION[j,i]>=threshold){
      hs_st<-REGION[j,1]
      hs_chr<-REGION[j,3]
      nm_peak<-REGION[j,i]
      
      while (REGION[j,i]>=threshold){
        hs_end<-REGION[j,2]
        nm_peak2<-REGION[j,i]
        if(nm_peak2>nm_peak){ ## find the maximum number of peak SNPs in any 200Kb window for a hotspot
          nm_peak<-nm_peak2
        }
        j<-j+1
      }
      print(j)
      hs<-c(hs_st,hs_end,hs_chr,cate,nm_peak)
      HS<-rbind.data.frame(HS,hs)
    }else{
      j<-j+1
    }
  }
}
colnames(HS)<-c("HS_ST","HS_END","HS_Chr","HS_Category","HS_max_number_peakSNP_200Kb")
setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/GWAS_hotspot")
#write.table(HS,"GWAS_hotspot_genomic_regions_Cate_and_All.txt",row.names=F,sep="\t",quote=F)
write.table(HS,"GWAS_hotspot_genomic_regions_Cate_and_All_04042022.txt",row.names=F,sep="\t",quote=F)
#######################

## all GWAS, TWAS and FCT candidates in hotspots

library(data.table)
setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/GWAS_hotspot")
#HS<-fread("GWAS_hotspot_genomic_regions_Cate_and_All.txt",data.table=F)
#HS<-fread("GWAS_hotspot_genomic_regions_Cate_and_All_04042022.txt",data.table=F)
HS<-fread("GWAS_hotspot_genomic_regions_Cate_and_All_04042022_noDup.txt",data.table=F)
## ***manually adjusted two intervals and removed 9 intervalls in ALL that were duplicated for within class HS

setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
G_genes<-fread("GWAS_wax_Top0002SNPs_all_cand_genes_Dist.txt",data.table=F)
T_GENEs<-fread("TWAS_0.0025_cand_for_SuppTable.txt",data.table=F)
F_GENEs<-fread("FCT_top0.0025_for_SuppTable.txt",data.table=F)


############# overlap between hotspot and TWAS/FCT (by category) #########
TWAS_hot<-c()
FCT_hot<-c()
GWAS_hot<-c()
i=1
for (i in 1:nrow(HS)){
  chr<-HS$HS_Chr[i]
  start<-HS$HS_ST[i]
  end<-HS$HS_END[i]
  cate<-HS$HS_Category[i]
  
  
  twas<-T_GENEs[which(T_GENEs$chr==chr &((T_GENEs$start<end & T_GENEs$start>start)|(T_GENEs$end<end & T_GENEs$end>start))),]
  fct<-F_GENEs[which(F_GENEs$chr==chr & ((F_GENEs$start<end & F_GENEs$start>start)|(F_GENEs$end<end & F_GENEs$end>start))),]
  gwas<-G_genes[which(G_genes$Chromosome==chr & ((G_genes$start<end & G_genes$start>start)|(G_genes$end<end & G_genes$end>start))),]
  if(nrow(twas)>0){
    #print(twas)
    twas$region_l<-start
    twas$region_r<-end
    twas$region_ch<-chr
    twas$nm_peakSNP<-HS[i,5]
    twas$hotspot_cate<-cate
    TWAS_hot<-rbind.data.frame(TWAS_hot,twas)
  }
  if(nrow(fct)>0){
    #print(fct)
    fct$region_l<-start
    fct$region_r<-end
    fct$region_ch<-chr
    fct$nm_peakSNP<-HS[i,5]
    fct$hotspot_cate<-cate
    FCT_hot<-rbind.data.frame(FCT_hot,fct)
    
  }
  if(nrow(gwas)>0){
    #print(gwas)
    gwas$region_l<-start
    gwas$region_r<-end
    gwas$region_ch<-chr
    gwas$nm_peakSNP<-HS[i,5]
    gwas$hotspot_cate<-cate
    GWAS_hot<-rbind.data.frame(GWAS_hot,gwas)
    
  }
}

TWAS_hot_uniq<-unique(TWAS_hot$gene_name)
TWAS_hot.1<-c()
for (tgene in TWAS_hot_uniq){
  temp<-TWAS_hot[which(TWAS_hot$gene_name==tgene),]
  all_HS_cate<-unique(temp$hotspot_cate)
  temp$all_HS_cate<-paste(all_HS_cate,collapse = ";")
  TWAS_hot.1<-rbind.data.frame(TWAS_hot.1,temp)
}

FCT_hot_uniq<-unique(FCT_hot$gene_name)
FCT_hot.1<-c()
for (tgene in FCT_hot_uniq){
  temp<-FCT_hot[which(FCT_hot$gene_name==tgene),]
  all_HS_cate<-unique(temp$hotspot_cate)
  temp$all_HS_cate<-paste(all_HS_cate,collapse = ";")
  FCT_hot.1<-rbind.data.frame(FCT_hot.1,temp)
}

GWAS_hot_uniq<-unique(GWAS_hot$gene_name)
GWAS_hot.1<-c()
for (tgene in GWAS_hot_uniq){
  temp<-GWAS_hot[which(GWAS_hot$gene_name==tgene),]
  all_HS_cate<-unique(temp$hotspot_cate)
  temp$all_HS_cate<-paste(all_HS_cate,collapse = ";")
  GWAS_hot.1<-rbind.data.frame(GWAS_hot.1,temp)
}

setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/GWAS_hotspot")
#write.table(TWAS_hot.1,"TWAS_cand_genes_inHotspot.txt",row.names=F,sep="\t",quote=F)
#write.table(FCT_hot.1,"FCT_cand_genes_inHotspot.txt",row.names=F,sep="\t",quote=F)
#write.table(GWAS_hot.1,"GWAS_cand_genes_inHotspot.txt",row.names=F,sep="\t",quote=F)

write.table(TWAS_hot.1,"TWAS_cand_genes_inHotspot_noDup.txt",row.names=F,sep="\t",quote=F)
write.table(FCT_hot.1,"FCT_cand_genes_inHotspot_noDup.txt",row.names=F,sep="\t",quote=F)
write.table(GWAS_hot.1,"GWAS_cand_genes_inHotspot_noDup.txt",row.names=F,sep="\t",quote=F)

# add annotation
annotation<-read.delim("/Users/ml2498/Desktop/Labserver2/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)
#annotation<-read.delim("/Users/Meng/Desktop/Labserver/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)
#annotation<-read.delim("/workdir/ml2498/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)

colnames(annotation)[1]<-"gene_name"
## need import from csv
TWAS_hot.1$gene_name<-as.character(TWAS_hot.1$gene_name)
FCT_hot.1$gene_name<-as.character(FCT_hot.1$gene_name)
GWAS_hot.1$gene_name<-as.character(GWAS_hot.1$gene_name)

TWAS_hot_ANNO<-merge(TWAS_hot.1,annotation,by="gene_name",all.y=F,all.x=T)
FCT_hot_ANNO<-merge(FCT_hot.1,annotation,by="gene_name",all.y=F,all.x=T)
GWAS_hot_ANNO<-merge(GWAS_hot.1,annotation,by="gene_name",all.y=F,all.x=T)

#setwd("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/GWAS_hotspot")
setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/GWAS_hotspot")

write.table(TWAS_hot_ANNO,paste("TWAS_cand_genes_inHotspot_noDup_wANNO.txt",sep=""),row.names=F,col.names=T,sep="\t",quote=F)
write.table(FCT_hot_ANNO,paste("FCT_cand_genes_inHotspot_noDup_wANNO.txt",sep=""),row.names=F,col.names=T,sep="\t",quote=F)
write.table(GWAS_hot_ANNO,paste("GWAS_cand_genes_inHotspot_noDup_wANNO.txt",sep=""),row.names=F,col.names=T,sep="\t",quote=F)




## all peak SNPs in hotspots

library(data.table)
setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/GWAS_hotspot")
#HS<-fread("GWAS_hotspot_genomic_regions_Cate_and_All.txt",data.table=F)
#HS<-fread("GWAS_hotspot_genomic_regions_Cate_and_All_04042022.txt",data.table=F)
HS<-fread("GWAS_hotspot_genomic_regions_Cate_and_All_04042022_noDup.txt",data.table=F)

cut="0002"
setwd(paste("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/peak_SNP_",cut,"_final/2.manual_peaks/",sep=""))

#gwas_peak<-fread("nonHeli_peaks_all_waxes_wGWAS_R2.txt",data.table=F)
gwas_peak<-fread("nonHeli_peaks_all_waxes_manual.txt",data.table=F)
## generate category
gwas_peak$trait[which(gwas_peak$trait=="Total")]<-"Total_total"
temp_name<-gwas_peak$trait
temp_name<-strsplit(temp_name,split="_")
temp_name<-matrix(unlist(temp_name),ncol=2,byrow=T)
gwas_peak$category<-temp_name[,1]

i=1
SNPs<-c()
for(i in 1:nrow(HS)){
  cate<-HS$HS_Category[i]
  chr<-HS$HS_Chr[i]
  st<-HS$HS_ST[i]
  end<-HS$HS_END[i]
  
  if(cate=="ALL"){
    snps<-gwas_peak[which(gwas_peak$Chromosome==chr & gwas_peak$Position>st & gwas_peak$Position<end),c(1:4,6,7,15)]
  }else{
    snps<-gwas_peak[which(gwas_peak$category==cate & gwas_peak$Chromosome==chr & gwas_peak$Position>st & gwas_peak$Position<end),c(1:4,6,7,15)]
  }
  
  if (nrow(snps)!=HS$HS_max_number_peakSNP_200Kb[i]){
    print(paste("warning for HS row ",i,sep=""))
  }
  
  temp<-HS[rep(i, nrow(snps)),]
  hs_snps<-cbind.data.frame(temp,snps)
  SNPs<-rbind.data.frame(SNPs,hs_snps)
}
setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/GWAS_hotspot")

#write.table(SNPs,"GWAS_peak_SNPs_inHotspot.txt",row.names=F,sep="\t",quote=F)
write.table(SNPs,"GWAS_peak_SNPs_inHotspot_noDup.txt",row.names=F,sep="\t",quote=F)

