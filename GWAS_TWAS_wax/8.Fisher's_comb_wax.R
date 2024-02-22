library(reshape2)
library(data.table)
reads_cutoff=2

rep<-"BLUP" 

Fisher <- function(x){res <- sumlog(x); return(res$p)}

############# TWAS results ############

twas0<-read.table("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/TWAS/Kin_P3D/PEER_K20_GWAS_kin_wax_SD18_TWAS_combScores.txt",header=T,sep="\t")
#twas0<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/TWAS/Kin_P3D/PEER_K20_GWAS_kin_wax_SD18_TWAS.txt",header=T,sep="\t")
#twas0<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/TWAS/Kin_P3D/PEER_K20_GWAS_kin_wax_SD18_TWAS.txt",header=T,sep="\t")
#Note: K10 is the best PEERres for cpm, K20 is the best PEERres for rlog (use rlog)

for(i in 4:ncol(twas0)){
  twas0[,i]<-10^((-1)*twas0[,i])
}

#########################################

#Files<-read.table("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/Transformed_BLUP_wax_imputed.txt",header=T,sep="\t",stringsAsFactors=F)
#Files<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/Transformed_BLUP_wax_imputed.txt",header=T,sep="\t",stringsAsFactors=F)
#Files<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/Transformed_BLUP_wax_imputed.txt",header=T,sep="\t",stringsAsFactors=F)
#Files<-colnames(Files)[-1]


Files<-read.table("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/pheno/wax_names_subDir.txt",header=T,sep="\t",stringsAsFactors=F)
#Files<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno/wax_names_subDir.txt",header=T,sep="\t",stringsAsFactors=F)


path="/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/"
#path="/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/"
OLRM_rows<-which(Files$subDir=="OLRM_traits")
j=45

#for (j in 4:60){
#for (j in 21:40){
for (i in OLRM_rows){
#for (i in 1:nrow(Files)){
  if(!is.na(Files[i,2])){
    f=Files[i,1]
    print(f)
    sub_dir<-Files[i,2]
    
    ## gwas results
    gwas<-fread(paste(path,sub_dir,"/",f,"_tr_HapMap3_topSNPs01.txt",sep=""),data.table=F)
    #gwas<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT/HMP3_v4_NOPR_450to310_LD02/ce_SD18_both_HapMap3_topSNPs01.txt",header=T,sep="\t")
    
    gwas$SNP<-as.character(gwas$SNP)
    ### combine gwas results with SNP-gene pairs
    gwas_new<-vector()
    for (chr in 1:10){
      gwas_sub<-gwas[which(gwas$Chromosome==chr),]
      #pair<-read.table(paste("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/Cand_selection/imputedGBS_450to310_filtered_chr",chr,"_pos_ID_geneAssignment.txt",sep=""),header=T,sep="\t")
      #pair<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/Cand_selection/imputedGBS_450to310_filtered_chr",chr,"_pos_ID_geneAssignment.txt",sep=""),header=T,sep="\t")
      pair<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/Cand_selection/imputedGBS_450to310_filtered_chr",chr,"_pos_ID_SNP-Gene_merged1_0715.txt",sep=""),header=T,sep="\t")
      
      #pair.1<-melt(pair,id.vars=c("ID"),measure.vars = c("gene", "gene2"))
      gwas_sub<-merge(gwas_sub,pair,by.x="SNP",by.y="ID",all.x=T)
      gwas_new<-rbind.data.frame(gwas_new,gwas_sub)
    }
    
    GWAS_new<-gwas_new[,c(1:5,15:ncol(gwas_new))]
    #GWAS_new<-gwas_new[which(gwas_new$variable=="gene"),]
    
    
    wax_col<-which(colnames(twas0)==paste(f,"_tr",sep=""))
    twas<-twas0[,c(1:3,wax_col)]
    colnames(twas)[4]<-"pvalue"
    T_GWAS<-merge(twas,GWAS_new,by="gene_name",all.y=T)
    T_GWAS<-as.data.frame(T_GWAS)
    T_GWAS$pvalue[is.na(T_GWAS$pvalue)]<-1
    
    ###############################################################
    # Fisher's combined test
    ###############################################################
    #install.packages("metap")
    library(metap)
    
    df <- cbind(T_GWAS$P.value,T_GWAS$pvalue)
    a=apply(df, 1, Fisher)

    fisher<-cbind.data.frame(T_GWAS$gene_name,T_GWAS$SNP,T_GWAS$Chromosome,T_GWAS$Position,a,T_GWAS$P.value,T_GWAS$pvalue)
    colnames(fisher)<-c("gene_name","SNP","chr","position","p","pG","pT")
    fisher$DF<-4
    setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/TWAS/Fisher")
    #write.table(fisher,paste("FisherRes_",f,"_nonP3D_TWAS.",rep,"_K20_HMP3.450to310_rlog.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    write.table(fisher,paste("FisherRes_",f,"_nonP3D_TWAS.",rep,"_K20_HMP3.450to310_rlog_0715pair.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
  }
}
################################################
# Manhattan plot
################################################
library(qqman)
library(qvalue)
library(data.table)
K=20
rep="BLUP"

Files<-read.table("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/pheno/wax_names_subDir.txt",header=T,sep="\t",stringsAsFactors=F)
#Files<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno/wax_names_subDir.txt",header=T,sep="\t",stringsAsFactors=F)


path="/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/"
#path="/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/"
OLRM_rows<-which(Files$subDir=="OLRM_traits")


# Files<-read.table("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/Transformed_BLUP_wax_imputed.txt",header=T,sep="\t",stringsAsFactors=F)
# #Files<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/Transformed_BLUP_wax_imputed.txt",header=T,sep="\t",stringsAsFactors=F)
# #Files<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/Transformed_BLUP_wax_imputed.txt",header=T,sep="\t",stringsAsFactors=F)
# Files<-colnames(Files)[-1]
#Files<-Files[c(1:15,21:35,41:55)]


setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/TWAS/Fisher")
setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/TWAS/Fisher")

for (r in OLRM_rows){
    t<-Files$Wax_trait[r]
    filename<-paste("FisherRes_",t,"_nonP3D_TWAS.BLUP_K20_HMP3.450to310_rlog_0715pair",sep="")
    #twas<-read.table(paste(filename,".txt",sep=""),header=T,sep="\t")
    twas<-fread(paste(filename,".txt",sep=""),data.table=F)

    #pdf(paste(filename,"_man.pdf",sep=""),width=8,height=4)
    png(paste(filename,"_man.png",sep=""), width=12,height=5,units="in",res=300,family="serif")
    manhattan(twas, chr = "chr", bp = "position", p = "p", snp = "SNP",
              col = c("navy", "darkorange1"), chrlabs = NULL,
              suggestiveline = FALSE,genomewideline =FALSE, 
              main=paste(filename))
    dev.off()
    
    jpeg(paste(filename,"_qq.jpg",sep=""))
    qq(twas$p)
    dev.off()
}
#######################################
#####################################
# pull out important genes
#####################################
library(pkgcond)
library(data.table)
align<-"hisat"
#act_P=2
K=20
rep="BLUP" # 1, 2,"b","PEERresBLUE"
transf="rlog"

Files<-read.table("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/pheno/wax_names_subDir.txt",header=T,sep="\t",stringsAsFactors=F)
#Files<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno/wax_names_subDir.txt",header=T,sep="\t",stringsAsFactors=F)
traits<-Files[!is.na(Files$subDir),1]
#path="/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/"
#path="/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/"
#OLRM_rows<-which(Files$subDir=="OLRM_traits")


# traits<-read.table("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/Transformed_BLUP_wax_imputed.txt",header=T,sep="\t",stringsAsFactors=F)
# traits<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/Transformed_BLUP_wax_imputed.txt",header=T,sep="\t",stringsAsFactors=F)
# #traits<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/Transformed_BLUP_wax_imputed.txt",header=T,sep="\t",stringsAsFactors=F)
# traits<-colnames(traits)[-1]
# traits<-gsub("_tr","",traits)
#annotation<-read.delim("/home/ml2498/Desktop/Labserver/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)
#annotation<-read.delim("/Users/Meng/Desktop/Labserver/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)
annotation<-read.delim("/workdir/ml2498/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)

colnames(annotation)[1]<-"gene_name"
Models<-"Kin"
Cutoff<-c(0.005,0.004,0.0035,0.003,0.0025,0.002)


for(cut in Cutoff[5]){
  print(cut)
  F_GENEs<-c()
  Result_SUM<-c()
  for (t in traits){
    #for (model in Models){
    model=Models[1]
    print(t)
      #setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/TWAS/Fisher")
      setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/TWAS/Fisher")
      
      #setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/TWAS/Fisher")
      #filename<-paste("FisherRes_",t,"_",model,"_TWAS.",rep,"_HMP3.450to310",sep="")
      filename<-paste("FisherRes_",t,"_nonP3D_TWAS.",rep,"_K20_HMP3.450to310_rlog_0715pair",sep="")
  
      twas0<-fread(paste(filename,".txt",sep=""),header=T,sep="\t")
      
      twas0$gene_name<-as.character(twas0$gene_name)
      twas0<-twas0[order(twas0$gene_name,twas0$p),]
      uniq_gene<-unique(twas0$gene_name)
      twas1<-twas0[match(uniq_gene,twas0$gene_name),]
      twas1<-twas1[order(twas1$p),]
      if(cut==0.005){
        write.table(twas1,paste("FisherRes_",t,"_nonP3D_TWAS.",rep,"_K20_HMP3_rlog_0715pair_short.txt",sep=""),row.names=F,sep="\t",quote=F)
      }
      
      nm_cand<-ceiling(length(uniq_gene)*cut)
     
      F_genes<-twas1[1:nm_cand,]
      F_genes$rank<-1:nrow(F_genes)
      F_genes$trait<-t
      F_GENEs<-rbind.data.frame(F_GENEs,F_genes)
      #which(F_genes$gene_name != F_genes.2$gene_name)
      #which(F_genes$p != F_genes.2$p)
      result_sum<-merge(F_genes,annotation,by="gene_name",all.x=T)
      #result_sum$trait<-t
      
      Result_SUM<-rbind.data.frame(Result_SUM,result_sum)
    #}
  }
  setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
  #setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
  
    write.table(F_GENEs,paste("Cand_FisherRes_allWax_kin_TWAS_K20_rlog_cut",cut,".txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
  write.table(Result_SUM,paste("Cand_FisherRes_allWax_kin_TWAS_K20_rlog_cut",cut,"_anno.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
}

