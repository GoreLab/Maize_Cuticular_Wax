pheno.all<-read.table("/Users/menglin/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT/CE_FT_alllines_03092019_wrapper.txt",header=T,sep="\t")
#pheno.all<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT/CE_FT_alllines_03092019_wrapper.txt",header=T,sep="\t")
#pheno.all<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/GAPIT/CE_FT_alllines_03092019_wrapper.txt",header=T,sep="\t")

file<-read.table("/Users/menglin/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno/wrapper/result_wax/BLUPS_all_traits_wax_imputed_BLUPinput_02192020.txt",header=T,sep="\t")
#file<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno/wrapper/result_wax/BLUPS_all_traits_wax_imputed_BLUPinput_02192020.txt",header=T,sep="\t")
file$predict<-as.character(file$predict)

index<-grep("MLC_STANDARD_[a-zA-Z0-9_.-]*:IS_EXPERIMENTAL$",file$predict)
blup<-file[index,]
taxa<-strsplit(blup$predict,split=":",fixed = T)
taxa<-matrix(unlist(taxa),ncol=2,byrow=T)
taxa<-substr(taxa[,1],14,nchar(taxa[,1]))
blup<-cbind(taxa,blup[,-1])
colnames(blup)[1]<-"MLC_STANDARD"
blup[which(blup[,1]=="MO17"),2:ncol(blup)]<-NA ## remove blup
pheno.all<-merge(pheno.all,blup,by="MLC_STANDARD",all=T)

setwd("/Users/menglin/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT")
#setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT")
write.table(pheno.all,"CE_FT_wax_imputed_alllines_02192020.txt",col.names=T,row.names=F,sep="\t",quote=F)
##############################

## check distribution of wax BLUPs
setwd("/Users/menglin/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT")
#setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT")
#setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT")

pheno.all<-read.table("CE_FT_wax_imputed_alllines_02192020.txt",header=T,sep="\t")

taxa<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/taxa310_pheno_geno_rna.txt")
#taxa<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/taxa310_pheno_geno_rna.txt")
#taxa<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/RNAseq/GeneExpression/v4_counts/taxa310_pheno_geno_rna.txt")
taxa<-as.character(taxa[,1])


wax_start<-which(colnames(pheno.all)=="FA_C24")
wax_end<-which(colnames(pheno.all)=="Total")

#setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/plot")
setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/plot")
pdf("Dist_wax_BLUPs.pdf",height=3,width=3)
for (i in c(wax_start:wax_end)){
  hist(pheno.all[,i],main=colnames(pheno.all)[i])
}
dev.off()

########### make tables for wax paper ########
wax_ce<-pheno.all[which(pheno.all$MLC_STANDARD %in% taxa),c(1,23,26,wax_start:wax_end)]
setwd("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/Wax_TWAS/GAPIT/")
write.table(wax_ce,"gc_ft_wax_310_untrans.txt",row.names=F,sep="\t",quote=F)


## boxplot to check relative abundance of 60 wax traits
install.packages("Rfast")
library(Rfast)
setwd("/Users/ml2498/Desktop/LabServer/MaizeLeafCuticle/Wax_TWAS/GAPIT/")
wax_ce<-read.table("gc_ft_wax_310_untrans.txt",header=T,sep="\t")
wax_forstat<-wax_ce[,c(4:56)]
wax_mean<-colMeans(wax_forstat)
wax_var<-colVars(as.matrix(wax_forstat))
wax_var<-colVars(as.matrix(wax_forstat[sapply(wax_forstat, is.numeric)]))
wax_stat<-cbind.data.frame(names(wax_mean),wax_mean,wax_var)
colnames(wax_stat)[1]<-"Wax"

categ<-strsplit(names(wax_mean),split="_")
categ<-matrix(unlist(categ),ncol=2,byrow=T)
colnames(categ)[1:2]<-c("Category","chain_length")
wax_stat<-cbind.data.frame(categ,wax_stat)
categ2<-cbind.data.frame(names(wax_mean),categ)
colnames(categ2)[1]<-"trait"

library(ggplot2)
p <- ggplot(wax_stat, aes(Wax, wax_mean,fill=Category)) + 
  geom_bar(stat='identity',aes(color=Category))+
  scale_fill_manual(values=c("#88CCEE", "#CC6677", "#DDCC77","#117733","#332288","#882255"))+
  geom_errorbar(aes(ymin=wax_mean, ymax=wax_mean+wax_var), width=.2,
                position=position_dodge(.9)) +
  labs(x="Cuticular waxes",
       y=expression("Abundance ("~mu~g~dm^-2~")"))+
  #geom_text_repel()
  theme_light()+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.y = element_text(size=10,colour = 'black'),
        strip.text.x = element_text(size=5,colour = 'black'),
        strip.placement = "outside",
        axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size=7),
        text = element_text(size=10),
        #legend.title=element_blank(),
        legend.text.align = 0,
        legend.text=element_text(size=rel(1)))
print(p)




Wax_plot1<-c()
for (i in 4:56){
  wax_temp<-wax_ce[,c(1,i)]
  colnames(wax_temp)[2]<-"abundance"
  wax_temp$trait<-colnames(wax_ce)[i]
  Wax_plot1<-rbind.data.frame(Wax_plot1,wax_temp)
}
Wax_plot1<-merge(categ2,Wax_plot1,all.y=T,by="trait")

## update wax names (3/16/2023)
library(data.table)
wax<-fread("/Users/ml2498/Library/CloudStorage/GoogleDrive-lemonsquirrel1987@gmail.com/My Drive/Wax_manuscript/Tables/draft/wax_name_Isabel_20230215.txt",data.table=F)

Wax_plot1$Category[which(Wax_plot1$Category=="FFA")]<-"DD"
Wax_plot1$Category[which(Wax_plot1$Category=="FA")]<-"PA"
Wax_plot1$Category[which(Wax_plot1$Category=="DD")]<-"FA"

Wax_plot1$trait<-gsub("FFA_","DD_",Wax_plot1$trait)

i=1
for(i in 1:nrow(Wax_plot1)){
  trait<-Wax_plot1$trait[i]
  Wax_plot1$trait[i]<-wax$`Cuticular wax`[which(wax$old_wax3==trait)]
}
unique(Wax_plot1$trait)
### end of update

ggplot(Wax_plot1,aes(x=trait,y=abundance))+
  geom_boxplot(outlier.shape = NA)

p <- ggplot(Wax_plot1, aes(trait, abundance)) + 
  geom_boxplot(outlier.shape = NA,aes(color=Category))+
  scale_color_manual(values=c("#88CCEE", "#CC6677", "#DDCC77","#117733","#332288","#882255"))+
  labs(x="Cuticular waxes",
       y=expression("Abundance ("~mu~g~dm^-2~")"),color="Wax class")+
  #geom_text_repel()
  scale_y_continuous(limits = c(0, 32))+
  theme_light()+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.y = element_text(size=10,colour = 'black'),
        strip.text.x = element_text(size=5,colour = 'black'),
        strip.placement = "outside",
        axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size=7),
        text = element_text(size=10),
        #legend.title=element_blank(),
        legend.text.align = 0,
        legend.text=element_text(size=rel(1)))
print(p)
setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/plot")
pdf("SupplFig_wax_boxplot1.pdf",height=4,width=8,family="sans")
#png("SupplFig_wax_boxplot1.png",height=400,width=800,family="sans")
print(p)
dev.off()



Wax_plot2<-c()
for (i in 57:62){
  wax_temp<-wax_ce[,c(1,i)]
  colnames(wax_temp)[2]<-"abundance"
  wax_temp$trait<-colnames(wax_ce)[i]
  Wax_plot2<-rbind.data.frame(Wax_plot2,wax_temp)
}
categ3<-strsplit(Wax_plot2$trait,split="_")
categ3<-matrix(unlist(categ3),ncol=2,byrow=T)
Wax_plot2<-cbind.data.frame(Wax_plot2,categ3[,1])
colnames(Wax_plot2)[4]<-"Category"
Wax_plot2$trait[which(Wax_plot2$trait=="FA_total")]<-"Total PA"
Wax_plot2$trait[which(Wax_plot2$trait=="FFA_total")]<-"Total FA"
Wax_plot2$trait[which(Wax_plot2$trait=="AD_total")]<-"Total AD"
Wax_plot2$trait[which(Wax_plot2$trait=="HC_total")]<-"Total HC"
Wax_plot2$trait[which(Wax_plot2$trait=="WE_total")]<-"Total WE"
Wax_plot2$trait[which(Wax_plot2$trait=="AC_total")]<-"Total AC"

Wax_plot2$Category[which(Wax_plot2$Category=="FFA")]<-"DD"
Wax_plot2$Category[which(Wax_plot2$Category=="FA")]<-"PA"
Wax_plot2$Category[which(Wax_plot2$Category=="DD")]<-"FA"



p <- ggplot(Wax_plot2, aes(trait, abundance)) + 
  geom_boxplot(outlier.shape = NA,aes(color=Category))+
  scale_color_manual(values=c("#88CCEE", "#CC6677", "#DDCC77","#117733","#332288","#882255"))+
  labs(x="Cuticular waxes",
       y=expression("Abundance ("~mu~g~dm^-2~")"),color="Wax class")+
  #geom_text_repel()
  #scale_y_continuous(limits = c(0, 32))+
  theme_light()+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text.y = element_text(size=10,colour = 'black'),
        strip.text.x = element_text(size=5,colour = 'black'),
        strip.placement = "outside",
        axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size=7),
        text = element_text(size=10),
        #legend.title=element_blank(),
        legend.text.align = 0,
        legend.text=element_text(size=rel(1)))
print(p)
setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/plot")
pdf("SupplFig_wax_boxplot2.pdf",height=4,width=3.5,family="sans")
print(p)
dev.off()


### check pairwise comparison wax-ce-ft
library(corrplot)
forplot<-pheno.all[,c(23,26,wax_start:wax_end)]
corrplot.mixed(cor(forplot,use = "pairwise.complete.obs"), number.cex = 1,tl.cex=0.8)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
cor_pheno<-cor(forplot,use = "pairwise.complete.obs")

pdf("cor_ce_ft_wax_upper.pdf",height=20,width=20)
corrplot(cor_pheno, method="color", col=col(200),
         type="upper",number.cex = 0.6,tl.cex=0.7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
)
dev.off()






