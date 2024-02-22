library(VennDiagram)

cut_g="0002"
#cut_g="0001"
#cut2="FDR10"
#setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/")
#setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/")
#cand_hybrid<-read.table(paste("GWAS_hybrid_cand_",cut_g,"_",cut2,".txt",sep=""),header=T,sep="\t")
#cand_hybrid$gene_name<-as.character(cand_hybrid$gene_name)

#setwd(paste("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/peak_SNP_",cut_g,sep=""))
#setwd(paste("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/peak_SNP_",cut_g,sep=""))
#setwd(paste("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/peak_SNP_",cut_g,sep=""))

# final GWAS candidates
setwd(paste("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/peak_SNP_",cut_g,"_final/2.manual_peaks/",sep=""))
#G_genes<-read.table(paste("all_waxes_",cut_g,"SNPs_200K_cand_genes.txt",sep=""),header=T,sep="\t")
G_genes<-read.table(paste("all_waxes_",cut_g,"SNPs_200K_cand_genes_final.txt",sep=""),header=T,sep="\t")
G_genes$gene_name<-as.character(G_genes$gene_name)

setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
#setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
#setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
#setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
cut=0.0025

F_GENEs<-read.table(paste("Cand_FisherRes_allWax_kin_TWAS_K20_rlog_cut",cut,".txt",sep=""),header=T,sep="\t") # previous results were replaced
F_GENEs$gene_name<-as.character(F_GENEs$gene_name)
#T_GENEs<-read.table(paste("TWAS_cand_",cut,"_wax_nonP3D_K=20_rlog.txt",sep=""),header=T,sep='\t')
T_GENEs<-read.table(paste("TWAS_cand_",cut,"_wax_nonP3D_K=20_rlog_OLRM.txt",sep=""),header=T,sep='\t')
T_GENEs$gene_name<-as.character(T_GENEs$gene_name)

GT_genes<-intersect(G_genes$gene_name,T_GENEs$gene_name)
GF_genes<-intersect(G_genes$gene_name,F_GENEs$gene_name)
TF_genes<-intersect(T_GENEs$gene_name,F_GENEs$gene_name)
GTF_genes<-intersect(F_GENEs$gene_name,intersect(G_genes$gene_name,T_GENEs$gene_name))

setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
#setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
#write.table(GT_genes,paste("GWAS_hybrid_TWAS",cut,"_intersect_AgriGo.txt",sep=""),row.names=F,col.names=F,sep="\t",quote=F)
#write.table(GT_genes,paste("GWAS_0.002_TWAS",cut,"_intersect_AgriGo.txt",sep=""),row.names=F,col.names=F,sep="\t",quote=F)
write.table(GT_genes,paste("GWAS_0.002_TWAS",cut,"_intersect_AgriGo_final.txt",sep=""),row.names=F,col.names=F,sep="\t",quote=F)
################################
#### GTF genes with annotation

G_genes_sub<-G_genes[which(G_genes$gene_name %in% GTF_genes),]
GWAS_trait<-c()
for (g in GTF_genes){
  temp<-G_genes_sub[which(G_genes_sub$gene_name==g),]
  G_trait<-paste(temp$trait,collapse=";")
  GWAS_trait<-c(GWAS_trait,G_trait)
}
GWAS_res<-cbind.data.frame(GTF_genes,GWAS_trait)
  
#check<-T_GENEs[which(T_GENEs$gene_name=="Zm00001d046454"),] # fat1: HC_C23 rank 39; AD_total rank 36
T_GENEs_sub<-T_GENEs[which(T_GENEs$gene_name %in% GTF_genes),]
TWAS_trait<-c()
for (g in GTF_genes){
  temp<-T_GENEs_sub[which(T_GENEs_sub$gene_name==g),]
  T_trait<-paste(paste(temp$wax,"(",temp$rank,")",sep=""),collapse=";")
  TWAS_trait<-c(TWAS_trait,T_trait)
}
TWAS_res<-cbind.data.frame(GTF_genes,TWAS_trait)

F_GENEs_sub<-F_GENEs[which(F_GENEs$gene_name %in% GTF_genes),] ### need add rank for FCT results
FCT_trait<-c()
for (g in GTF_genes){
  temp<-F_GENEs_sub[which(F_GENEs_sub$gene_name==g),]
  F_trait<-paste(paste(temp$trait,"(",temp$rank,")",sep=""),collapse=";")
  FCT_trait<-c(FCT_trait,F_trait)
}
FCT_res<-cbind.data.frame(GTF_genes,FCT_trait)

all_res<-merge(GWAS_res,TWAS_res,by="GTF_genes",all=T)
all_res<-merge(all_res,FCT_res,by="GTF_genes",all=T)

#annotation<-read.delim("/workdir/ml2498/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)
annotation<-read.delim("/home/ml2498/Desktop/Labserver/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)
colnames(annotation)[1]<-"gene_name"
## need import from csv

all_genes_v4ANNO<-merge(all_res,annotation,by.y="gene_name",by.x="GTF_genes",all.y=F,all.x=T)

#setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
#setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
#write.table(all_genes_v4ANNO,paste("GTF_results_intersect_wAnno.txt",sep=""),row.names=F,col.names=T,sep="\t",quote=F)
write.table(all_genes_v4ANNO,paste("GTF_results_intersect_wAnno_final.txt",sep=""),row.names=F,col.names=T,sep="\t",quote=F)

######### grouping the core list into 4 priory classes (according to Laurie):
G_genes_sub<-G_genes[which(G_genes$gene_name %in% GTF_genes),]
T_GENEs_sub<-T_GENEs[which(T_GENEs$gene_name %in% GTF_genes),]
F_GENEs_sub<-F_GENEs[which(F_GENEs$gene_name %in% GTF_genes),] ### need add rank for FCT results

## generate categories
G_genes_sub$trait[which(G_genes_sub$trait=="Total")]<-"Total_total"
temp_name<-G_genes_sub$trait
temp_name<-strsplit(temp_name,split="_")
temp_name<-matrix(unlist(temp_name),ncol=2,byrow=T)
G_genes_sub$category<-temp_name[,1]

F_GENEs_sub$trait[which(F_GENEs_sub$trait=="Total")]<-"Total_total"
temp_name<-F_GENEs_sub$trait
temp_name<-strsplit(temp_name,split="_")
temp_name<-matrix(unlist(temp_name),ncol=2,byrow=T)
F_GENEs_sub$category<-temp_name[,1]

##### start grouping ###
g=GTF_genes[2]
which(GTF_genes=="Zm00001d044467")
g=GTF_genes[226]

group<-c()
for (g in GTF_genes){
  temp.G<-G_genes_sub[which(G_genes_sub$gene_name==g),]
  temp.T<-T_GENEs_sub[which(T_GENEs_sub$gene_name==g),]
  temp.F<-F_GENEs_sub[which(F_GENEs_sub$gene_name==g),]
  
  (trait.G<-unique(temp.G$trait))
  (trait.T<-unique(temp.T$wax))
  (trait.F<-unique(temp.F$trait))
  
  (cate.G<-unique(temp.G$category))
  (cate.T<-unique(temp.T$category))
  (cate.F<-unique(temp.F$category))
  
  if(length(intersect(trait.G, intersect(trait.T,trait.F)))==0){ 
    
    if(length(intersect(cate.G, intersect(cate.T,cate.F)))==0){
      
      res.G<-duplicated(temp.G$category)
      res.T<-duplicated(temp.T$category)
      res.F<-duplicated(temp.F$category)
      all_cate<-unique(c(cate.G,cate.T,cate.F))
      
      if(any(c(res.G,res.T,res.F))){ ## at least two traits in one category in at least one method
        group<-c(group,"Group3")
      }else{
        group<-c(group,"Group4") ## none of the other situations
      }
      
    }else{
      group<-c(group,"Group2") ## same category in three methods
    }
    
  }else{
    group<-c(group,"Group1") ## same trait in three methods
  }
  
}
group<-cbind.data.frame(GTF_genes,group)
setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
write.table(group,"GTF_intersect_group_info.txt",row.names=F,quote=F,sep="\t")

#### Supplemental Table: GFT list ###
library(data.table)
setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
all_genes_v4ANNO<-fread("GTF_results_intersect_wAnno_final.txt",data.table=F)
uniq_gene<-unique(all_genes_v4ANNO$GTF_genes)
all_genes<-all_genes_v4ANNO[match(uniq_gene,all_genes_v4ANNO$GTF_genes),1:5]

all_genes$nm_Gtrait<-1
all_genes$nm_Ttrait<-1
all_genes$nm_Ftrait<-1
all_genes$nm_ALLtrait<-1

for(i in 1:nrow(all_genes)){
  Gtrait<-all_genes$GWAS_trait[i]
  Ttrait<-all_genes$TWAS_trait[i]
  Ftrait<-all_genes$FCT_trait[i]
  
  ## count G_traits
  gtrait<-strsplit(Gtrait,split=";")
  gtrait<-data.frame(matrix(unlist(gtrait),ncol=1))
  all_genes$nm_Gtrait[i]<-nrow(gtrait)
  
  ## count T_traits
  ttrait<-strsplit(Ttrait,split=";")
  ttrait<-data.frame(matrix(unlist(ttrait),ncol=1))
  all_genes$nm_Ttrait[i]<-nrow(ttrait)
  ttrait.1<-ttrait[,1]
  
  ttrait.2<-gsub("\\s*\\([^\\)]+\\)", "", ttrait.1)
  
  ## count F_traits
  ftrait<-strsplit(Ftrait,split=";")
  ftrait<-data.frame(matrix(unlist(ftrait),ncol=1))
  all_genes$nm_Ftrait[i]<-nrow(ftrait)
  ftrait.1<-ftrait[,1]
  
  ftrait.2<-gsub("\\s*\\([^\\)]+\\)", "", ftrait.1)
  
  ## count All_traits
  all_traits<-unique(c(gtrait[,1],ttrait.2,ftrait.2))
  all_genes$nm_ALLtrait[i]<-length(all_traits)
}

group<-fread("GTF_intersect_group_info.txt",data.table=F)
all_genes_group<-merge(all_genes,group,by="GTF_genes",all=T)
write.table(all_genes_group,"SuppTable_GTF_intersect_group_counts.txt",row.names=F,sep="\t",quote=F)

#### update the table with final wax traits 
library(data.table)
setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
all_genes_group<-fread("SuppTable_GTF_intersect_group_counts.txt",data.table=F)
LCWE<-c("WE_C47","WE_C48","WE_C49","WE_C50","WE_C52","WE_C54","WE_C56")

#### how many genes were associated with longer chain WEs
geneID<-c()
for(i in 1:nrow(all_genes_group)){
  for(j in 2:4){ # trait columns
    temp_tr<-all_genes_group[i,j]
    for(w in LCWE){
      if(length(grep(w,temp_tr))>0){
        geneID<-c(geneID,all_genes_group$GTF_genes[i])
      }
    }
  }
}
geneID.1<-unique(geneID)
###### update trait names:***


