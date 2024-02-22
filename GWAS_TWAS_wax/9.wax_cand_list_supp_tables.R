library(data.table)

############# GWAS SNPs ########
###### need to recalculate distance between peaks and other SNPs
setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/")
Files<-read.table("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/pheno/wax_names_subDir.txt",header=T,sep="\t",stringsAsFactors=F)
#Files<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno/wax_names_subDir.txt",header=T,sep="\t",stringsAsFactors=F)
#Files<-t(Files[1,]) #64 wax names

main<-"/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/"

gwas_res0=fread('/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/peak_SNP_0002_final/wax_GWAS_SD18_both_HapMap3_for_peakSNPs_final.txt',data.table = F)
traits<-unique(gwas_res0$trait)
gwas_res<-c()
for (t in traits){
  temp<-gwas_res0[which(gwas_res0$trait==t),]
  temp<-temp[order(temp$P.value),]
  temp<-temp[1:ceiling(971.5*2/10),]
  gwas_res<-rbind.data.frame(gwas_res,temp)
}


#gwas_snp<-fread("ce_SD18_both_HapMap3_topSNPs_forCand.txt",data.table=F)
t=traits[1]
GWAS_SNP<-c()
for (t in traits){
  print(t)
  gwas_snp<-gwas_res[which(gwas_res$trait==t),]
  gwas_snp<-gwas_snp[order(gwas_snp$P.value),]
  gwas_snp$SNP_rank<-c(1:nrow(gwas_snp))
  gwas_snp$SNP_ranking_perc<-gwas_snp$SNP_rank*100/9745050 # generate percentile
  
  ### allelic effect ###
  Effect<-c()
  folder<-Files$subDir[which(Files$Wax_trait==t)]
  CHR<-unique(gwas_snp$Chromosome)
  for (i in CHR){
    sub<-paste(folder,"/chr",i,sep="")
    setwd(paste(main,sub,sep=""))
    effect<-fread(paste("GAPIT.MLM.",t,"_tr.Allelic_Effect_Estimates.csv",sep=""),data.table=F)
    effect_sub<-effect[which(effect$SNP %in% gwas_snp$SNP),c(1,4)] #SNP ID and allelic effect estimate
    Effect<-rbind.data.frame(Effect,effect_sub)
    
  }
  gwas_snp<-merge(gwas_snp,Effect,by="SNP",all.x=T)
  GWAS_SNP<-rbind.data.frame(GWAS_SNP,gwas_snp)
}

## sort by chr and position in each trait
GWAS_SNP1<-c()
for(t in traits){
  print(t)
  gwas_temp<-GWAS_SNP[which(GWAS_SNP$trait==t),]
  #gwas_temp<-gwas_temp[order(gwas_temp$Chromosome.x, gwas_temp$Position.x),]
  gwas_temp<-gwas_temp[order(gwas_temp$Chromosome, gwas_temp$Position),]
  GWAS_SNP1<-rbind.data.frame(GWAS_SNP1,gwas_temp)
}

## generate category
GWAS_SNP1$trait[which(GWAS_SNP1$trait=="Total")]<-"Total_total"
temp_name<-GWAS_SNP1$trait
temp_name<-strsplit(temp_name,split="_")
temp_name<-matrix(unlist(temp_name),ncol=2,byrow=T)
GWAS_SNP1$category<-temp_name[,1]

setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
write.table(GWAS_SNP1,"Supp_Tablex_wax_GWAS_top0002_allSNPs.txt",row.names=F,sep="\t",quote=F)




############# GWAS candidate genes ########
library(data.table)

### GWAS statistics
table<-fread("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/Cand_Selection/Supp_Tablex_wax_GWAS_top0002_allSNPs.txt",data.table=F)
table<-table[,-(14:15)]
colnames(table)[2:3]<-c("Chromosome","Position")
table$trait[which(table$trait=="Total_total")]<-"Total"
#check<-table[which(table$trait=="Total_total"),]


#v4_gene_gtf<-read.table("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/TWAS/v4.37_gene_pos_formatted.txt",header=T,sep="\t")
v4_gene_gtf<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS/v4.37_gene_pos_formatted.txt",header=T,sep="\t")

setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/peak_SNP_0002_final/2.manual_peaks")
Files<-read.table("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/pheno/wax_names_subDir.txt",header=T,sep="\t",stringsAsFactors=F)
Files<-Files[!is.na(Files[,2]),]

#all_files<-list.files(path=".",pattern="peak_snps_v3_")
f=Files[1,1]
ALL_genes<-c()

for(f in Files[,1]){
  print(f)
  GWAS_Peak<-fread(paste("peak_snps_v3_",f,".csv",sep=""),data.table = F)
  if(f %in% c("AC_b.Amyrin","AC_a.Amyrin","AC_Friedelin","AC_total","AC_Unk1","AC_Unk7","FFA_C30")){
    GWAS_Peak<-GWAS_Peak
  }else{
    GWAS_Peak<-GWAS_Peak[which(GWAS_Peak$comments==""),] ## exclude helicopters
  }
  
  GWAS_Peak$SNP<-as.character(GWAS_Peak$SNP)
  
  table_peak<-table[intersect(which(table$SNP %in% GWAS_Peak$SNP),which(table$trait==f)),]
  table_peak$Position<-as.numeric(as.character(table_peak$Position))
  table_peak<-table_peak[order(table_peak$Chromosome,table_peak$Position),]
  table_peak$SNP<-as.character(table_peak$SNP)  # v4 results
  table_peak$left250<-table_peak$Position-200000
  table_peak$right250<-table_peak$Position+200000
  
  gwas_hit<-table_peak
  all_genes<-matrix(ncol=15,nrow=0)
  for (i in 1:nrow(gwas_hit)){
    one_hit<-gwas_hit[i,]
    sub_genes<-v4_gene_gtf[which(v4_gene_gtf$chr==one_hit$Chromosome & v4_gene_gtf$start<one_hit$right250 & v4_gene_gtf$end>one_hit$left250),]
    nm_subgenes<-nrow(sub_genes)
    hit_info<-one_hit[rep(1, each=nm_subgenes),]
    sub_genes<-cbind(hit_info,sub_genes)
    all_genes<-rbind(all_genes,sub_genes)
  }
  
  ## remove duplicated genes
  print(length(unique(all_genes$gene_name)))
  keeplist2<-unique(all_genes$gene_name)
  all_genes<-all_genes[match(keeplist2,all_genes$gene_name),]
  
  ## Distance
  all_genes$Dist<-NA
  for (i in 1:nrow(all_genes)){
    if (all_genes$Position[i]>all_genes$start[i]&all_genes$Position[i]<all_genes$end[i]){
      all_genes$Dist[i]<-0
    }else {
      all_genes$Dist[i]<-min(c(abs(all_genes$Position[i]-all_genes$start[i]),abs(all_genes$Position[i]-all_genes$end[i])))
    }
  }
  
  ALL_genes<-rbind.data.frame(ALL_genes,all_genes)
}

setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
write.table(ALL_genes,"GWAS_wax_Top0002SNPs_all_cand_genes_Dist.txt",row.names=F,col.names=T,sep="\t",quote=F)

##### Add annotation ######

#annotation<-read.delim("/Users/Meng/Desktop/Labserver/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)
annotation<-read.delim("/workdir/ml2498/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)

colnames(annotation)[1]<-"gene_name"

keeplist3<-unique(annotation$gene_name)
annotation<-annotation[match(keeplist3,annotation$gene_name),]


ALL_genes$gene_name<-as.character(ALL_genes$gene_name)

ALL_genes_v4ANNO<-merge(ALL_genes,annotation,by="gene_name",all.y=F,all.x=T)

#setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/Cand_selection")
setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
write.table(ALL_genes_v4ANNO,"GWAS_wax_Top0002SNPs_all_cand_genes_Dist_wANNO.txt",row.names=F,col.names=T,sep="\t",quote=F)

# add "TWAS P-value,	 TWAS rank, 	TWAS ranking %" using vlookup ****

######## TWAS cand list #########
## total genes 20,013 tested in TWAS
library(data.table)
v4_gene_gtf<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS/v4.37_gene_pos_formatted.txt",header=T,sep="\t")

## add gwas statistics if common genes were found
gwas_cand<-fread("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/Cand_Selection/GWAS_wax_Top0002SNPs_all_cand_genes_Dist.txt",
                 data.table=F)

#### All TWAS candidates ####
setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/TWAS/Kin_P3D")
file<-"TWAS_cand_001_wax_nonP3D_K=20_rlog_OLRM"

TWAS0<-read.table(paste(file,".txt",sep=""),header=T,sep="\t")
TWAS0$wax<-as.character(TWAS0$wax)
TWAS0$wax<-gsub("_tr","",TWAS0$wax)

Traits<-unique(TWAS0$wax)
TWAS1<-c()
for(t in Traits){
  TWAS<-TWAS0[which(TWAS0$wax==t),]
  TWAS<-merge(TWAS[,c(1,4,5)],v4_gene_gtf,by="gene_name",all.x=T)
  TWAS<-TWAS[order(TWAS$pvalue),]
  TWAS<-TWAS[1:ceiling(20013*0.0025),]
  TWAS$TWAS_rank<-c(1:nrow(TWAS))
  TWAS$TWAS_rank_perc<-TWAS$TWAS_rank/20013
  
  gwas_sub<-gwas_cand[which(gwas_cand$trait==t),]
  gwas_sub<-gwas_sub[,c(1,12:13,24:25)]
  
  TWAS<-merge(TWAS,gwas_sub,by="gene_name",all.x=T)
  TWAS1<-rbind.data.frame(TWAS1,TWAS) ## 1994 unique gene, 3060 in 60 traits
}

## generate category
TWAS1$wax[which(TWAS1$wax=="Total")]<-"Total_wax"
temp_name<-TWAS1$wax
temp_name<-strsplit(temp_name,split="_")
temp_name<-matrix(unlist(temp_name),ncol=2,byrow=T)
TWAS1$category<-temp_name[,1]


##### Add annotation ######

#annotation<-read.delim("/Users/Meng/Desktop/Labserver/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)
annotation<-read.delim("/workdir/ml2498/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)

colnames(annotation)[1]<-"gene_name"

keeplist3<-unique(annotation$gene_name)
annotation<-annotation[match(keeplist3,annotation$gene_name),]

TWAS1$gene_name<-as.character(TWAS1$gene_name)

TWAS1_v4ANNO<-merge(TWAS1,annotation,by="gene_name",all.y=F,all.x=T)

setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
write.table(TWAS1_v4ANNO,"TWAS_0.0025_cand_for_SuppTable.txt",row.names=F,col.names=T,sep="\t",quote=F)

####################################
## FCT candidate gene list table ####
######################################
library(data.table)
setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")

##
cd /workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/TWAS/Fisher
wc -l FisherRes_*_nonP3D_TWAS.*_K20_HMP3_rlog_0715pair_short.txt > Number_genes_in_FCT
wc -l FisherRes_FA_total_nonP3D_TWAS.*_K20_HMP3_rlog_0715pair_short.txt
##

FCT0<-fread("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/Cand_Selection/Cand_FisherRes_allWax_kin_TWAS_K20_rlog_cut0.0025.txt",data.table=F)
FCT0<-fread("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/Cand_Selection/Cand_FisherRes_allWax_kin_TWAS_K20_rlog_cut0.0025.txt",data.table=F)
FCT_gene_nm<-fread("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/TWAS/Fisher/Number_genes_in_FCT.txt",data.table=F)

v4_gene_gtf<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS/v4.37_gene_pos_formatted.txt",header=T,sep="\t")
v4_gene_gtf<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS/v4.37_gene_pos_formatted.txt",header=T,sep="\t")
FCT0<-merge(v4_gene_gtf[,c(4:7)],FCT0[,-8],by="gene_name",all.y=T)

###### relative TWAS results

TWAS_all<-fread("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/TWAS/Kin_P3D/PEER_K20_GWAS_kin_wax_SD18_TWAS_combScores.txt",data.table=F)
TWAS_all<-TWAS_all[which(TWAS_all$chr %in% c(1:10)),]

for(i in 4:ncol(TWAS_all)){
  TWAS_all[,i]<-10^((-1)*TWAS_all[,i])
}
colnames(TWAS_all)<-gsub("_tr","",colnames(TWAS_all))


## gene to SNP distance
FCT0$gene_SNP_dist<-NA

for (i in 1:nrow(FCT0)){
  pos<-FCT0$position[i]
  start<-FCT0$start[i]
  end<-FCT0$end[i]
  
  if (pos>=start&pos<=end){
    FCT0$gene_SNP_dist[i]<-0
  }else {
    if(FCT0$strand[i]=="+"){
      if(pos<start){
        FCT0$gene_SNP_dist[i]<-pos-start
      }else if(pos>end){
        FCT0$gene_SNP_dist[i]<-pos-end
      }
      
    }else if(FCT0$strand[i]=="-"){
      if(pos<start){
        FCT0$gene_SNP_dist[i]<-start-pos ## SNP at 3' of a gene
      }else if(pos>end){
        FCT0$gene_SNP_dist[i]<-end-pos  ## SNP at 5' of a gene
      }
    }
  }
}

####### Distance between gene and closest peak SNP (could be no relative peak SNPs)

FCT0$peak_snp<-NA
FCT0$peak_snp_pos<-NA
FCT0$peak_snp_dist<-NA

traits<-unique(FCT0$trait)

FCT_ALL<-c()
for(t in traits){
  FCT<-FCT0[which(FCT0$trait==t),]
  GWAS_Peak<-fread(paste("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/peak_SNP_0002_final/2.manual_peaks/peak_snps_v3_",t,".csv",sep=""),data.table=F) # modified based on "ce_SD18_Peak_SNPs.txt" results
  
  if(t %in% c("AC_b.Amyrin","AC_a.Amyrin","AC_Friedelin","AC_total","AC_Unk1","AC_Unk7","FFA_C30")){
    GWAS_Peak<-GWAS_Peak[,1:3]
  }else{
    GWAS_Peak<-GWAS_Peak[which(GWAS_Peak$comments==""),1:3] ## exclude helicopters
  }
  
  ## identify closest peak SNP
  for (i in 1:nrow(FCT)){
    chr<-FCT$chr[i]
    start<-FCT$start[i]
    end<-FCT$end[i]
    
    sub_peak<-GWAS_Peak[which(GWAS_Peak$Chromosome==chr),]
    if(nrow(sub_peak)>0){
      DIST<-abs(sub_peak$Position-start)
      close_peak<-which.min(DIST)
      FCT$peak_snp_pos[i]<-sub_peak$Position[close_peak]
      #FCT$peak_snp_nm[i]<-sub_peak$peak_num_man[close_peak]
      FCT$peak_snp[i]<-as.character(sub_peak$SNP[close_peak])
    }
  }
  
  ## distance to closest Peak SNP
  
  for (j in 1:nrow(FCT)){
    if(!is.na(FCT$peak_snp[j])){
      pos<-FCT$peak_snp_pos[j]
      start<-FCT$start[j]
      end<-FCT$end[j]
      
      if (pos>=start&pos<=end){
        FCT$peak_snp_dist[j]<-0
      }else {
        if(FCT$strand[j]=="+"){
          if(pos<start){
            FCT$peak_snp_dist[j]<-pos-start
          }else if(pos>end){
            FCT$peak_snp_dist[j]<-pos-end
          }
          
        }else if(FCT$strand[j]=="-"){
          if(pos<start){
            FCT$peak_snp_dist[j]<-start-pos ## SNP at 3' of a gene
          }else if(pos>end){
            FCT$peak_snp_dist[j]<-end-pos  ## SNP at 5' of a gene
          }
        }
      }
    }
  }
  
  ###### relative TWAS results
  
  TWAS<-TWAS_all[,c(1,which(colnames(TWAS_all)==t))]
  colnames(TWAS)[2]<-"Pvalue_t"
  TWAS<-TWAS[order(TWAS$Pvalue_t),]
  TWAS$rank<-c(1:nrow(TWAS))
  TWAS$rank_perc<-TWAS$rank*100/20013
  
  FCT_final<-merge(FCT,TWAS,by="gene_name",all.x=T)
  
  FCT_final<-FCT_final[order(FCT_final$p),]
  FCT_final$FCT_rank<-c(1:nrow(FCT_final))
  fct_nm<-FCT_gene_nm[which(FCT_gene_nm[,2]==t),1]
  FCT_final$FCT_rank_perc<-FCT_final$FCT_rank*100/fct_nm
  
  FCT_ALL<-rbind(FCT_ALL,FCT_final)
}

## generate category
FCT_ALL$trait[which(FCT_ALL$trait=="Total")]<-"Total_wax"
temp_name<-FCT_ALL$trait
temp_name<-strsplit(temp_name,split="_")
temp_name<-matrix(unlist(temp_name),ncol=2,byrow=T)
FCT_ALL$category<-temp_name[,1]


##### Add annotation ######

annotation<-read.delim("/Users/Meng/Desktop/Labserver/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)
annotation<-read.delim("/workdir/ml2498/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)
colnames(annotation)[1]<-"gene_name"

keeplist3<-unique(annotation$gene_name)
annotation<-annotation[match(keeplist3,annotation$gene_name),]

FCT_ALL$gene_name<-as.character(FCT_ALL$gene_name)
FCT_ALL_v4ANNO<-merge(FCT_ALL,annotation,by="gene_name",all.y=F,all.x=T)

setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
write.table(FCT_ALL_v4ANNO,"FCT_top0.0025_for_SuppTable.txt",row.names=F,sep="\t",quote=F)
