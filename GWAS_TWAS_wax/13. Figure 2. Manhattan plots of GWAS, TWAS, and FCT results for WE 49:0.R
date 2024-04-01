##13. Figure 2. Manhattan plots of GWAS, TWAS, and FCT results for WE 49:0.

library(RColorBrewer)
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(compiler)
library(scatterplot3d)
library(EMMREML)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")
library(qvalue)
library(data.table)
library(tidyverse)
#install.packages("svglite")
library(svglite)

setwd("C:/Users/harel/OneDrive/Documents/")
# Open a PDF graphics device
#pdf("GTF_plots2.pdf", width = 16, height = 9)  # Adjust size as needed
svglite("GTF_plots2.svg", width = 16, height = 9) 

# Set up layout for 3 rows, 1 column
par(mfrow = c(3, 1), mar = c(4, 10, 4, 10))

# ---------------------
# GWAS Plot
# ---------------------
meth="GWAS_"
Files<-c("WE_C49")

#for (f in Files[4]){
f=Files[1]
results<-matrix(nrow=0,ncol=10)

setwd("C:/Users/harel/OneDrive/intermediate/Wax_MS/GTF_WE49/GTF_data_plot/plots/gwas_we49")
mainDir="C:/Users/harel/OneDrive/intermediate/Wax_MS/GTF_WE49/GTF_data_plot/plots/gwas_we49"
for (i in 1:10){
  subDir<-paste("chr",i,sep="")
  setwd(file.path(mainDir, subDir))
  
  file<-paste("GAPIT.MLM.",f,".GWAS.Results.csv",sep="")
  gemma<-fread(file,data.table=F)
  gemma$Chromosome<-i
  results<-rbind(results,gemma)
}
results<-results[order(results$Chromosome,results$Position),]


results$logP<- (-1)*log10(results$P.value)
results<-results[,c(2,3,10,1)]
colnames(results)<-c("Chromosome","Position","logP","gene_name") ## gene_name is SNP ID for gwas
#}

GI.MP <- results
plot_name<-paste(meth,".Manhattan.Genomewise_",f,"_geneLab_hl001.png" ,sep = "")

highlight<-c("4-14843410","5-215544967")
geneName<-c("TOL","GELP2")


#########################
# Plot Manhattan
#########################
y.lim <- ceiling(max(GI.MP$logP))

threshold_TWAS <- -log10(0.0032194578567087) # highest P-value of the top 0.25% genes
threshold_FCT <- -log10(0.000138567677322815) # highest P-value of the top 0.25% genes

#Remove all SNPs that do not have a chromosome, bp position and p value(NA)
GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
GI.MP <- GI.MP[!is.na(GI.MP[,2]),]
GI.MP <- GI.MP[!is.na(GI.MP[,3]),]
#Retain SNPs that have -logP values >=0
GI.MP <- GI.MP[GI.MP[,3]>=0,]

#Remove chr 0 and 99
GI.MP <- GI.MP[GI.MP[,1] %in% c(1:10),]
GI.MP[,1]<-as.numeric(as.character(GI.MP[,1]))
numMarker=nrow(GI.MP)

chm.to.analyze <- unique(GI.MP[,1])
chm.to.analyze=chm.to.analyze[order(chm.to.analyze)]
numCHR= length(chm.to.analyze)
nchr=length(chm.to.analyze)
#
GI.MP <- GI.MP[order(GI.MP[,2]),]
GI.MP <- GI.MP[order(GI.MP[,1]),]

GI.MP$color<-"lightskyblue"# "darkorange1"
GI.MP$color[which(GI.MP$Chromosome %in% c(2,4,6,8,10))]<-"slateblue" #"royalblue3"
GI.MP$color[which(GI.MP$gene_name %in% highlight)]<-"red"

ticks=NULL
lastbase=0
lastbase_chr=c()
#change base position to accumulatives (ticks)
for (i in chm.to.analyze){
  index=(GI.MP[,1]==i)
  ticks <- c(ticks, lastbase+mean(GI.MP[index,2]))
  GI.MP[index,2]=GI.MP[index,2]+lastbase
  lastbase_chr=c(lastbase_chr,lastbase)
  lastbase=max(GI.MP[index,2])
}

x0 <- as.numeric(GI.MP[,2])
y0 <- as.numeric(GI.MP[,3])
z0 <- as.numeric(GI.MP[,1])

position=order(y0,decreasing = TRUE)
index0=GAPIT.Pruning(y0[position],DPP=500000) # pruning for plot
index=position[index0]
GI.MP_plot<-GI.MP[index,]
color<-GI.MP_plot$color


if(length(which(color=="red"))!=0){
  HL_points<-GI.MP_plot[which(GI.MP_plot$color=="red"),]
}


plot(NULL,xlab="",ylab="",xaxt = 'n', yaxt = 'n',ylim=c(0,y.lim+1),xlim=c(0,lastbase),axes=TRUE)

threshold <- -log10(0.0000528071068872993) #highest P-value of the top 0.002% snps

# for Fisher and TWAS only
abline(v=lastbase_chr[-1],lwd=0.5,col='darkgray',lty="dotted")

#Set axis
axis(1, at=c(0,ticks,lastbase),cex.axis=1,labels=c("",chm.to.analyze,""),tick=FALSE,las=1)
axis(2, at=c(0:20),cex.axis=1,labels=c(0:20),tick=T,las=1)
mtext(side = 1, text = "Chromosome", line = 3,cex=1.2)

#plot snps
points(GI.MP_plot[,3]~GI.MP_plot[,2],col=color,pch=20,cex=0.2)
mtext("GWAS", side = 4, line = 2, las = 0, cex = 1.2)

if(length(which(color == "red")) != 0) {
  HL_points <- GI.MP_plot[which(GI.MP_plot$color == "red"),]
  points(HL_points[,3]~HL_points[,2], col="red", pch=20, cex=1.4) # Increase size with cex=0.6
  
  some_offset_y <-0.7
  # Loop through each highlighted point to add gene names
  for (i in 1:nrow(HL_points)) {
    # Adjust text position as needed
    text_x <- HL_points[i, "Position"] 
    text_y <- HL_points[i, "logP"] + some_offset_y
    text_label <- geneName[which(highlight == HL_points[i, "gene_name"])]
    
    text(text_x, text_y, labels = text_label, cex = 1.4, col = "black")
  }
}

abline(h = threshold, col = "red", lty = "dashed")
mtext(side = 2, line = 2.5, text=expression(paste(-log[10](italic(p)))),cex=1)

gene_lab<-cbind.data.frame(highlight,geneName)
colnames(gene_lab)<-c("gene_name","gene_lab")

gene_lab<-merge(HL_points,gene_lab,by="gene_name",all=T)
gene_lab$color<-"maroon3" ## "CuticleBio"

for (l in 1:nrow(gene_lab)){
  if(gene_lab$gene_lab[l]=="Ypt/Rab-GAP"){
    adj=-0.5
    text(x=gene_lab$Position[l],y=8.5+adj,labels=gene_lab$gene_lab[l],col=gene_lab$color[l],cex=0.7) 
    
  }
 
  else{
    adj=0
    text(x=gene_lab$Position[l],y=8.5+adj,labels=gene_lab$gene_lab[l],col=gene_lab$color[l],cex=0.7) #mute this line if for FCPU, for publication, use cex=0.7
    
  }
  
  
}

# ---------------------
# TWAS Plot
# ---------------------
meth="TWAS_"
results<-read.table("C:/Users/harel/OneDrive/intermediate/Wax_MS/GTF_WE49/GTF_data_plot/plots/twas_we_49/PEER_K20_GWAS_kin_wax_SD18_TWAS_combScores.txt"
, header=T,sep="\t")
# Select only the specified columns
results<- results[, c("gene_name", "chr", "start", "WE_C49")]
results <- rename(results, logP = WE_C49)

Files<-c("WE_C49")
f<-Files[1]
if (f==Files[1]){
  highlight<-c("Zm00001d038404", "Zm00001d045660")
  
  geneName<-c("Ypt/Rab-GAP","KCS22")
  results<-results[,c(2:4,1)]
  
}else if(f==Files[1]){
  highlight<-c("Zm00001d038404","Zm00001d038405", "Zm00001d045660")
  results<-results.all[,c(2:4,1)]
}

colnames(results)<-c("Chromosome","Position","logP","gene_name")
GI.MP <- results

plot_name<-paste(meth,".Manhattan.Kin_",f,"_rlog_K20_HL.png",sep="")

#########################
# Plot Manhattan
#########################
y.lim <- ceiling(max(GI.MP$logP))

#Remove all SNPs that do not have a chromosome, bp position and p value(NA)
GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
GI.MP <- GI.MP[!is.na(GI.MP[,2]),]
GI.MP <- GI.MP[!is.na(GI.MP[,3]),]
#Retain SNPs that have -logP values >=0
GI.MP <- GI.MP[GI.MP[,3]>=0,]

#Remove chr 0 and 99
GI.MP <- GI.MP[GI.MP[,1] %in% c(1:10),]
GI.MP[,1]<-as.numeric(as.character(GI.MP[,1]))
numMarker=nrow(GI.MP)

chm.to.analyze <- unique(GI.MP[,1])
chm.to.analyze=chm.to.analyze[order(chm.to.analyze)]
numCHR= length(chm.to.analyze)
nchr=length(chm.to.analyze)
#
GI.MP <- GI.MP[order(GI.MP[,2]),]
GI.MP <- GI.MP[order(GI.MP[,1]),]

#Remove all SNPs that do not have a choromosome, bp position and p value(NA)
GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
GI.MP <- GI.MP[!is.na(GI.MP[,2]),]
GI.MP <- GI.MP[!is.na(GI.MP[,3]),]
#Retain SNPs that have -logP values >=0
GI.MP <- GI.MP[GI.MP[,3]>=0,]

#Remove chr 0 and 99
GI.MP <- GI.MP[GI.MP[,1] %in% c(1:10),]
GI.MP[,1]<-as.numeric(as.character(GI.MP[,1]))
numMarker=nrow(GI.MP)

chm.to.analyze <- unique(GI.MP[,1])
chm.to.analyze=chm.to.analyze[order(chm.to.analyze)]
numCHR= length(chm.to.analyze)
nchr=length(chm.to.analyze)
#
GI.MP <- GI.MP[order(GI.MP[,2]),]
GI.MP <- GI.MP[order(GI.MP[,1]),]

GI.MP$color<-"lightskyblue"
GI.MP$color[which(GI.MP$Chromosome %in% c(2,4,6,8,10))]<-"slateblue" 

GI.MP$color[which(GI.MP$gene_name %in% highlight)]<-"red"

ticks=NULL
lastbase=0
lastbase_chr=c()
#change base position to accumulative (ticks)
for (i in chm.to.analyze){
  index=(GI.MP[,1]==i)
  ticks <- c(ticks, lastbase+mean(GI.MP[index,2]))
  GI.MP[index,2]=GI.MP[index,2]+lastbase
  lastbase_chr=c(lastbase_chr,lastbase)
  lastbase=max(GI.MP[index,2])
}

x0 <- as.numeric(GI.MP[,2])
y0 <- as.numeric(GI.MP[,3])
z0 <- as.numeric(GI.MP[,1])

GI.MP_plot<-GI.MP
color<-GI.MP$color


if(length(which(color=="red"))!=0){
  HL_points<-GI.MP_plot[which(GI.MP_plot$color=="red"),]
}

plot(NULL,xlab="",ylab="",xaxt = 'n', yaxt = 'n',ylim=c(0,y.lim+1),xlim=c(0,lastbase),axes=TRUE)
threshold <- threshold_FCT
mtext("TWAS", side = 4, line = 2, las = 0, cex = 1.2)

abline(v=lastbase_chr[-1],lwd=0.5,col='darkgray',lty="dotted")

#set axises
axis(1, at=c(0,ticks,lastbase),cex.axis=1,labels=c("",chm.to.analyze,""),tick=FALSE,las=1)
axis(2, at=c(0:20),cex.axis=1,labels=c(0:20),tick=T,las=1)
mtext(side = 1, text = "Chromosome", line = 3,cex=1.2)

#plot snps
points(GI.MP_plot[,3]~GI.MP_plot[,2],col=color,pch=20,cex=0.2)

if(length(which(color == "red")) != 0) {
  HL_points <- GI.MP_plot[which(GI.MP_plot$color == "red"),]
  points(HL_points[,3]~HL_points[,2], col="red", pch=20, cex=1.4) # Increase size with cex=0.6
  
  some_offset_y <-0.7
  
  # Loop through each highlighted point to add gene names
  for (i in 1:nrow(HL_points)) {
    # Adjust text position as needed
    text_x <- HL_points[i, "Position"]
    text_y <- HL_points[i, "logP"] + some_offset_y
    text_label <- geneName[which(highlight == HL_points[i, "gene_name"])]
    
    text(text_x, text_y, labels = text_label, cex = 1.4, col = "black")
  }
}

abline(h = threshold_TWAS, col = "red", lty = "dashed")
mtext(side = 2, line = 2.5, text=expression(paste(-log[10](italic(p)))),cex=1)

gene_lab<-cbind.data.frame(highlight,geneName)
colnames(gene_lab)<-c("gene_name","gene_lab")

HL_points<-HL_points[order(HL_points$gene_name,HL_points$logP,decreasing=T),]
temp_lab<-unique(HL_points$gene_name)
HL_points_reduced<-HL_points[match(temp_lab,HL_points$gene_name),]

gene_lab<-merge(HL_points_reduced,gene_lab,by="gene_name",all=T)

# ---------------------
# Fisher Plot
# ---------------------
meth="Fisher_"
Files<-c("WE_C49")
f<-Files[1]

if(f==Files[1]){
  
  ## final list
  highlight<-c("Zm00001d045660")
  geneName<-c("KCS22")
  results<-read.table("C:/Users/harel/OneDrive/intermediate/Wax_MS/GTF_WE49/GTF_data_plot/plots/fct_we49/FisherRes_WE_C49_nonP3D_TWAS.BLUP_K20_HMP3_rlog_0715pair_short.txt", header=T,sep="\t")
  
}else if(f==Files[1]){
  highlight<-c("Zm00001d045660")
  results<-read.table("C:/Users/harel/Downloads/fct_we49/FisherRes_WE_C49_nonP3D_TWAS.BLUP_K20_HMP3_rlog_0715pair_short.txt", header=T,sep="\t")
  
}

results$logP<- (-1)*log10(results$p)
results<-results[,c(3:4,9,1)]
colnames(results)<-c("Chromosome","Position","logP","gene_name")
GI.MP <- results

plot_name<-paste(meth,".Manhattan.Kin_",f,"_rlog_K20_HL.png",sep="")

y.lim <- ceiling(max(GI.MP$logP))

threshold_TWAS <- -log10(0.0032194578567087) #highest P-value of the top 0.25% genes
threshold_FCT <- -log10(0.000138567677322815) #highest P-value of the top 0.25% genes

#Remove all SNPs that do not have a chromosome, bp position and p value(NA)
GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
GI.MP <- GI.MP[!is.na(GI.MP[,2]),]
GI.MP <- GI.MP[!is.na(GI.MP[,3]),]
#Retain SNPs that have -logP values >=0
GI.MP <- GI.MP[GI.MP[,3]>=0,]

#Remove chr 0 and 99
GI.MP <- GI.MP[GI.MP[,1] %in% c(1:10),]
GI.MP[,1]<-as.numeric(as.character(GI.MP[,1]))
numMarker=nrow(GI.MP)

chm.to.analyze <- unique(GI.MP[,1])
chm.to.analyze=chm.to.analyze[order(chm.to.analyze)]
numCHR= length(chm.to.analyze)
nchr=length(chm.to.analyze)
#
GI.MP <- GI.MP[order(GI.MP[,2]),]
GI.MP <- GI.MP[order(GI.MP[,1]),]

GI.MP$color<-"lightskyblue"
GI.MP$color[which(GI.MP$Chromosome %in% c(2,4,6,8,10))]<-"slateblue"

GI.MP$color[which(GI.MP$gene_name %in% highlight)]<-"red"
#}

ticks=NULL
lastbase=0
lastbase_chr=c()
#change base position to accumulative (ticks)
for (i in chm.to.analyze){
  index=(GI.MP[,1]==i)
  ticks <- c(ticks, lastbase+mean(GI.MP[index,2]))
  GI.MP[index,2]=GI.MP[index,2]+lastbase
  lastbase_chr=c(lastbase_chr,lastbase)
  lastbase=max(GI.MP[index,2])
}

x0 <- as.numeric(GI.MP[,2])
y0 <- as.numeric(GI.MP[,3])
z0 <- as.numeric(GI.MP[,1])

GI.MP_plot<-GI.MP
color<-GI.MP$color

if(length(which(color=="red"))!=0){
  HL_points<-GI.MP_plot[which(GI.MP_plot$color=="red"),]
}

plot(NULL,xlab="",ylab="",xaxt = 'n', yaxt = 'n',ylim=c(0,y.lim+1),xlim=c(0,lastbase),axes=TRUE)
threshold <- threshold_FCT
mtext("FCT", side = 4, line = 2, las = 0, cex = 1.2)

abline(v=lastbase_chr[-1],lwd=0.5,col='darkgray',lty="dotted")

#Set axis
axis(1, at=c(0,ticks,lastbase),cex.axis=1,labels=c("",chm.to.analyze,""),tick=FALSE,las=1)
axis(2, at=c(0:20),cex.axis=1,labels=c(0:20),tick=T,las=1)
mtext(side = 1, text = "Chromosome", line = 3,cex=1.2)

#plot snps
points(GI.MP_plot[,3]~GI.MP_plot[,2],col=color,pch=20,cex=0.2)

if(length(which(color == "red")) != 0) {
  HL_points <- GI.MP_plot[which(GI.MP_plot$color == "red"),]
  points(HL_points[,3]~HL_points[,2], col="red", pch=20, cex=1.4) 
  
  some_offset_y <-0.7
  # Loop through each highlighted point to add gene names
  for (i in 1:nrow(HL_points)) {
    # Adjust text position as needed
    text_x <- HL_points[i, "Position"] 
    text_y <- HL_points[i, "logP"] + some_offset_y
    text_label <- geneName[which(highlight == HL_points[i, "gene_name"])]
    
    text(text_x, text_y, labels = text_label, cex = 1.4, col = "black")
  }
}
# Add the horizontal dashed red line at the specified threshold
abline(h = threshold, col = "red", lty = "dashed")
mtext(side = 2, line = 2.5, text=expression(paste(-log[10](italic(p)))),cex=1)

gene_lab<-cbind.data.frame(highlight,geneName)
colnames(gene_lab)<-c("gene_name","gene_lab")

HL_points<-HL_points[order(HL_points$gene_name,HL_points$logP,decreasing=T),]
temp_lab<-unique(HL_points$gene_name)
HL_points_reduced<-HL_points[match(temp_lab,HL_points$gene_name),]

gene_lab<-merge(HL_points_reduced,gene_lab,by="gene_name",all=T)

dev.off()
