library(qvalue)
library(rrBLUP)
library(matrixcalc)
#install.packages("Matrix")
library(Matrix)
############################
align<-"hisat"
rep<-"BLUP" #
transf<-"rlog" # or "untr"
K=20

## kinship
setwd("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/")
#setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS_v4/")
prun<-"LD02" # "LD01", "LD02"

geno_LD02<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT/centeredIBS_HMP3_AGPv4_LD02_450to310.txt",header=F,sep="\t")
#geno_LD02<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT/centeredIBS_HMP3_AGPv4_LD02_450to310.txt",header=F,sep="\t")#
#geno_LD02<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT/centeredIBS_HMP3_AGPv4_LD02_450to310.txt",header=F,sep="\t")
rownames(geno_LD02)<-geno_LD02[,1]
geno_LD02<-geno_LD02[,-1]
colnames(geno_LD02)<-rownames(geno_LD02)

### gene position
v4_gene_gff<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS/v4.37_gene_pos_formatted.txt",header=T,sep="\t") # the complete .gtf file was shorterned by picking lines for genes (server)
#v4_gene_gff<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS/v4.37_gene_pos_formatted.txt",header=T,sep="\t") # the complete .gtf file was shorterned by picking lines for genes (server)
#v4_gene_gff<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS/v4.37_gene_pos_formatted.txt",header=T,sep="\t") # the complete .gtf file was shorterned by picking lines for genes (server)

v4_gene_gff<-v4_gene_gff[which(v4_gene_gff$chr %in% c(1:10)),]
v4_gene_gff$chr<-as.numeric(as.character(v4_gene_gff$chr))
v4_gene_gff$start<-as.numeric(as.character(v4_gene_gff$start))
v4_gene_gff<-v4_gene_gff[,c(1,4,7)]

setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/TWAS/Kin_P3D")
#setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/TWAS/Kin_P3D")
#taxa<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/taxa310_pheno_geno_rna.txt")
taxa<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/taxa310_pheno_geno_rna.txt")
#taxa<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/taxa310_pheno_geno_rna.txt")
taxa<-as.character(taxa[,1])

kin<-geno_LD02[which(rownames(geno_LD02) %in% taxa),which(colnames(geno_LD02) %in% taxa)]
kin<-kin[order(rownames(kin)),order(colnames(kin))]
kin<-as.matrix(kin)

posdefmat <- function(mat) {
  if (is.positive.definite(round(mat, 18))) {
    g = mat
  }
  else {
    g <-nearPD(mat)$mat
    warning("The matrix was adjusted for the nearest positive definite matrix")
  }
  return(g)
}

kin.adj<-posdefmat(kin)
kin.test<-as.matrix(kin.adj)

### phenotype
#pheno.all<-read.table("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/Transformed_BLUP_wax_imputed.txt",header=T,sep="\t")
#pheno.all<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/Transformed_BLUP_wax_imputed.txt",header=T,sep="\t")
#pheno.all<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/Transformed_BLUP_wax_imputed_OLRM.txt",header=T,sep="\t")

pheno.all<-read.table("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/Transformed_BLUP_wax_imputed_OLRM.txt",header=T,sep="\t")


pheno.all<-pheno.all[which(pheno.all$MLC_STANDARD %in% taxa),] # only test SD18_both for CE and FT
pheno.all<-pheno.all[order(pheno.all$MLC_STANDARD),]

#######################################################

#pheno.all[is.na(pheno.all[,3]),3]<-mean(pheno.all[,3],na.rm=T)

counts.log<-read.table(paste("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/TWAS_v4/PEERfactorSel/",
                             transf,"_",align,"_PEER",K,"_PEERres_rep",rep,"_OLRM_final.txt",sep=""),header=T,sep="\t",row.names=1)
# counts.log<-read.table(paste("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/TWAS_v4/PEERfactorSel/",
#                              transf,"_",align,"_PEER",K,"_PEERres_rep",rep,"_OLRM_final.txt",sep=""),header=T,sep="\t",row.names=1)

#check<-as.data.frame(counts.log[,which(colnames(counts.log) %in% c("Zm00001d013798", "Zm00001d008671" ,"Zm00001d0086712", "Zm00001d008674"))])
#colnames(counts.log)[which(colnames(counts.log) %in% c("Zm00001d013798", "Zm00001d008671" ,"Zm00001d0086712", "Zm00001d008674"))]
# Zm00001d013798

counts.log<-t(counts.log) # row: marker; col: observation
counts.log<-counts.log[,order(colnames(counts.log))]
which(colnames(counts.log)!=pheno.all$MLC_STANDARD)
which(colnames(counts.log)!=colnames(kin))
counts.log<-cbind.data.frame(rownames(counts.log),counts.log)
colnames(counts.log)[1]<-"gene_name"
counts.log<-merge(v4_gene_gff,counts.log,by="gene_name",all=F)

#######################################
#
#pheno<-pheno.all[,1:3]
scores.noP3D<-GWAS(pheno.all, counts.log, fixed=NULL, K=kin.test, n.PC=0,min.MAF=-Inf, n.core=10, P3D=FALSE, plot=FALSE)
#write.table(scores.noP3D,"PEER_K10_GWAS_kin_wax_SD18_TWAS.txt",col.names=T,row.names=F,sep="\t",quote=F)
write.table(scores.noP3D,"PEER_K20_GWAS_kin_wax_SD18_TWAS.txt",col.names=T,row.names=F,sep="\t",quote=F)

#### ****only for six traits with outlier removed**** ######
scores.noP3D<-GWAS(pheno.all, counts.log, fixed=NULL, K=kin.test, n.PC=0,min.MAF=-Inf, n.core=10, P3D=FALSE, plot=FALSE)
#write.table(scores.noP3D,"PEER_K10_GWAS_kin_wax_SD18_TWAS.txt",col.names=T,row.names=F,sep="\t",quote=F)
write.table(scores.noP3D,"PEER_K20_GWAS_kin_wax_SD18_TWAS_6OLRM.txt",col.names=T,row.names=F,sep="\t",quote=F)

###### combine TWAS results #####
library(data.table)
setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/TWAS/Kin_P3D")

old_scores<-fread("PEER_K20_GWAS_kin_wax_SD18_TWAS.txt",data.table=F)
new_scores<-fread("PEER_K20_GWAS_kin_wax_SD18_TWAS_6OLRM.txt",data.table=F)
comb_scores<-old_scores
for(i in 3:ncol(new_scores)){
  trait<-colnames(new_scores)[i]
  comb_scores[,which(colnames(comb_scores)==trait)]<-new_scores[,i]
}
write.table(comb_scores,"PEER_K20_GWAS_kin_wax_SD18_TWAS_combScores.txt",col.names=T,row.names=F,sep="\t",quote=F)

#######################################
# Manhattan plot
######################################
library(qqman)
library(qvalue)

setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/TWAS/Kin_P3D")
#setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/TWAS/Kin_P3D")

#twas0<-read.table("PEER_K10_P3D_kin_wax_SD18_TWAS.txt",header=T,sep="\t")
#twas0<-read.table("PEER_K20_GWAS_kin_wax_SD18_TWAS.txt",header=T,sep="\t")
twas0<-read.table("PEER_K20_GWAS_kin_wax_SD18_TWAS_combScores.txt",header=T,sep="\t")

for(i in 4:ncol(twas0)){
  twas0[,i]<-10^((-1)*twas0[,i])
}

TWAS_Cand<-vector()

pdf(paste("Wax_TWAS_nonP3D_PEER20_man_rlog_combScores.pdf",sep=""),width=8,height=4)
#pdf(paste("Wax_TWAS_nonP3D_PEER20_man_rlog.pdf",sep=""),width=8,height=4)
#pdf(paste("Wax_TWAS_P3D_PEER10_man.pdf",sep=""),width=8,height=4)
for (i in 4:ncol(twas0)){
  twas<-twas0[,c(1:3,i)]
  twas<-twas[order(twas[,4]),]
  twas.cand<-twas[1:ceiling(nrow(twas)*0.01),]
  
  wax<-colnames(twas)[4]
  
  twas.cand<-cbind(twas.cand,rep(wax,nrow(twas.cand)))
  colnames(twas.cand)[4:5]<-c("pvalue","wax")
  TWAS_Cand<-rbind(TWAS_Cand,twas.cand)
  
  qobj <- qvalue(p = twas[,4])
  qvalues <- qobj$qvalues
  twas$qvalues<-qvalues
  colnames(twas)[4]<-"pvalue"
  
  ## start of if else
  if (length(twas[which(twas$qvalues<0.05),4])==0){
    
    manhattan(twas, chr = "chr", bp = "start", p = "pvalue", snp = "gene_name",
              col = c("navy", "darkorange1"), chrlabs = NULL,
              suggestiveline = FALSE,genomewideline =FALSE, 
              main=paste("TWAS kinship model, PEER K=10: ",wax,sep=""))
    
  }else {
    threshold<-(max(twas$pvalue[which(twas$qvalues<0.05)])+min(twas$pvalue[which(twas$qvalues>0.05)]))/2
    manhattan(twas, chr = "chr", bp = "start", p = "pvalue", snp = "gene_name",
              col = c("navy", "darkorange1"), chrlabs = NULL,
              suggestiveline = FALSE, genomewideline = -log10(threshold),
              main=paste("TWAS kinship model, PEER K=10: ",wax,sep=""))
  }
  ## end of if else
}
dev.off()
#write.table(TWAS_Cand,paste("TWAS_cand_001_wax_P3D_K=10.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
#write.table(TWAS_Cand,paste("TWAS_cand_001_wax_nonP3D_K=20_rlog.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
write.table(TWAS_Cand,paste("TWAS_cand_001_wax_nonP3D_K=20_rlog_OLRM.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)

pdf(paste("Wax_TWAS_nonP3D_PEER20_rlog_qq_combScores.pdf",sep=""),width=4,height=4)
#pdf(paste("Wax_TWAS_nonP3D_PEER20_rlog_qq.pdf",sep=""),width=4,height=4)
#pdf(paste("Wax_TWAS_P3D_PEER10_qq.pdf",sep=""),width=4,height=4)
for (i in 4:ncol(twas0)){
  twas<-twas0[,c(1:3,i)]
  wax<-colnames(twas)[4]
  colnames(twas)[4]<-"pvalue"
  
  qq(as.numeric(twas$pvalue),main=wax)
}
dev.off()
################################################
# add annotation
################################################
setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/TWAS/Kin_P3D")
#setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/TWAS/Kin_P3D")
#file<-"TWAS_cand_001_wax_nonP3D_K=20_rlog"
file<-"TWAS_cand_001_wax_nonP3D_K=20_rlog_OLRM"

TWAS<-read.table(paste(file,".txt",sep=""),header=T,sep="\t")

#annotation<-read.delim("/home/ml2498/Desktop/Labserver/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)
#annotation<-read.delim("/Users/Meng/Desktop/Labserver/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)
annotation<-read.delim("/workdir/ml2498/OfficeCmp/GoogleDrive/MLC_AZ_2017/gene_study/Susanne/Anno_database_maize_at_rice.txt",header=T,sep="\t",stringsAsFactors=FALSE)

colnames(annotation)[1]<-"gene_name"
## need import from csv
TWAS$gene_name<-as.character(TWAS$gene_name)

all_genes_v4ANNO<-merge(TWAS,annotation,by="gene_name",all.y=F,all.x=T)

#setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
#setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")
setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/Cand_Selection")

write.table(all_genes_v4ANNO,paste(file,"_wAnno.txt",sep=""),row.names=F,col.names=T,sep="\t",quote=F)





