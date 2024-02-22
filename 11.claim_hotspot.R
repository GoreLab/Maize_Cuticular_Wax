library(data.table)
cut="0002"
#setwd(paste("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/peak_SNP_",cut,sep=""))
#setwd(paste("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/peak_SNP_",cut,sep=""))
#setwd(paste("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/peak_SNP_",cut,sep=""))

setwd(paste("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/peak_SNP_",cut,"_final/2.manual_peaks/",sep=""))
#setwd(paste("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/peak_SNP_",cut,"_final/2.manual_peaks/",sep=""))

#gwas_peak<-fread("nonHeli_peaks_all_waxes_wGWAS_R2.txt",data.table=F)
gwas_peak<-fread("nonHeli_peaks_all_waxes_manual.txt",data.table=F)
## generate category
gwas_peak$trait[which(gwas_peak$trait=="Total")]<-"Total_total"
temp_name<-gwas_peak$trait
temp_name<-strsplit(temp_name,split="_")
temp_name<-matrix(unlist(temp_name),ncol=2,byrow=T)
gwas_peak$category<-temp_name[,1]


str(gwas_peak)
gwas_peak.sub<-gwas_peak[intersect(which(gwas_peak$Chromosome==1), grep("WE",gwas_peak$trait)),]
gwas_peak.sub2<-gwas_peak[intersect(which(gwas_peak$Chromosome==1), grep("FFA",gwas_peak$trait)),]
gwas_peak.sub3<-gwas_peak[intersect(which(gwas_peak$Chromosome==1), grep("FA_",gwas_peak$trait)),]

last_base<-c(306970776,244418435,235631188,246966153,223706444,173376415,181719261,181046077, 159687194,150930636)

window_size<-c(10000,50000,200000)


Cate<-unique(gwas_peak$category)
Cate<-Cate[-which(Cate=="Total")]
Chr<-unique(gwas_peak$Chromosome)

ws=10000
cate="AC"

setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/GWAS_hotspot")


for(ws in window_size[3]){
  #step_size<-0.5*ws
  step_size<-0.2*ws
  REGION<-c()
  for (chr in 1:10){ # generate sliding windows
    print(chr)
    start<-0
    end<-start+ws
    ST=0
    END<-end
    #while (start<(100001-ws)){
    while (start<(last_base[chr]-ws)){
      start<-start+step_size
      end<-start+ws
      #print(end)
      #print(start)
      ST<-c(ST,start)
      END<-c(END,end)
    }
    region<-cbind.data.frame(ST,END)
    region$chr<-chr
    REGION<-rbind.data.frame(REGION,region)
  }
  
  REGION$AC<-0
  REGION$AD<-0
  REGION$FA<-0
  REGION$FFA<-0
  REGION$HC<-0
  REGION$WE<-0
  REGION$ALL<-0
  
  # count peak SNPs in each region (peak SNPs fall into which regions)
  for (cate in Cate){
    gwas_cate<-gwas_peak[which(gwas_peak$category==cate),]
    
    for(i in 1:nrow(gwas_cate)){
      ch<-gwas_cate$Chromosome[i]
      pos<-gwas_cate$Position[i]
      find_index<-which(REGION$chr==ch & REGION$ST<pos & REGION$END>pos)
      REGION[find_index,which(colnames(REGION)==cate)]<-REGION[find_index,which(colnames(REGION)==cate)]+1
    }
  }
  
  ## number of peaks SNPs in a 200Kb sliding window (including total as a trait)
  for(i in 1:nrow(gwas_peak)){
    ch<-gwas_peak$Chromosome[i]
    pos<-gwas_peak$Position[i]
    find_index<-which(REGION$chr==ch & REGION$ST<pos & REGION$END>pos)
    REGION[find_index,which(colnames(REGION)=="ALL")]<-REGION[find_index,which(colnames(REGION)=="ALL")]+1
  }
  
  #REGION$ALL<-rowSums(REGION[,4:9])

  #write.table(REGION,paste("GWAS_peak_counts_byCategory_WS",ws,"_SS",step_size,".txt",sep=""),row.names=F,sep="\t",quote=F)  
  write.table(REGION,paste("GWAS_peak_counts_byCategory_WS",ws,"_SS",step_size,"_04042022.txt",sep=""),row.names=F,sep="\t",quote=F)  
  
}

#REGION[rowSums(REGION[,4:9])>0,]

####### plot distribution of peak SNP counts #####
library(data.table)
setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/GWAS_hotspot")
#setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/HMP3_NOPR_450to310_LD02/GWAS_hotspot")

#file<-"GWAS_peak_counts_byCategory_WS2e+05_SS1e+05";ws<-200000
#file<-"GWAS_peak_counts_byCategory_WS10000_SS5000";ws<-10000
#file<-"GWAS_peak_counts_byCategory_WS50000_SS25000";ws<-50000

#file<-"GWAS_peak_counts_byCategory_WS2e+05_SS40000";ws<-200000 ## use this one to define hotspots
file<-"GWAS_peak_counts_byCategory_WS2e+05_SS40000_04042022";ws<-200000 ## use this one to define hotspots
#file<-"GWAS_peak_counts_byCategory_WS10000_SS2000";ws<-10000
#file<-"GWAS_peak_counts_byCategory_WS50000_SS10000";ws<-50000

REGION<-fread(paste(file,".txt",sep=""),data.table=F)
colnames(REGION)[4:10]<-c("alicyclic compounds", "aldehydes","fatty alcohols","free fatty acids",
                    "hydrocarbons","wax esters","all waxes combined")


pdf(paste(file,".pdf"),family="sans",height=5,width=6)
for(i in 4:ncol(REGION)){
  temp<-REGION[which(REGION[,i]>0),i]
  cate<-colnames(REGION)[i]
  #(thres<-quantile(temp, probs = c(0.05,0.25,0.75,0.95)))
  #thres<-quantile(temp, probs = c(0.75,0.90))
  thres<-quantile(temp, probs = c(0.90))
  hist(temp,main=paste("Peak SNP counts for ",cate," (window size ",ws,")",sep=""),xlab=paste("Peak SNP counts for ",cate,sep=""))
  #abline(v=thres,col=c("red","black"))
  abline(v=thres,col=c("red"))
}
hist(REGION[,ncol(REGION)],main=paste("Peak SNP counts for ",cate," (window size ",ws,")",sep=""),xlab=paste("Peak SNP counts for ",cate,sep=""))
thres<-quantile(temp, probs = c(0.90))
abline(v=thres,col=c("red"))
dev.off()



pdf(paste(file,"_0.95QT.pdf",sep=""),family="sans",height=5,width=6)
for(i in 4:ncol(REGION)){
  temp<-REGION[which(REGION[,i]>0),i]
  temp.1<-REGION[,i]
  cate<-colnames(REGION)[i]
  #(thres<-quantile(temp, probs = c(0.05,0.25,0.75,0.95)))
  #thres<-quantile(temp, probs = c(0.75,0.90))
  #thres<-quantile(temp, probs = c(0.90))
  thres<-quantile(temp, probs = c(0.95))
  thres.1<-quantile(temp.1, probs = c(0.95))
  hist(temp,main=paste("Peak SNP counts for ",cate," (window size ",ws,")",sep=""),xlab=paste("Peak SNP counts for ",cate,sep=""))
  #abline(v=thres,col=c("red","black"))
  abline(v=thres,col=c("red"))
  hist(temp.1,main=paste("Peak SNP counts for ",cate," (window size ",ws,")",sep=""),xlab=paste("Peak SNP counts for ",cate,sep=""))
  abline(v=thres.1,col=c("red"))
}
dev.off()


pdf(paste(file,"_0.95QT_paper.pdf",sep=""),family="sans",height=15,width=6)
#par(mfrow = c(7, 2))
par(mfrow = c(4, 2))
for(i in 4:ncol(REGION)){
  temp<-REGION[which(REGION[,i]>0),i]
  temp.1<-REGION[,i]
  cate<-colnames(REGION)[i]
  thres<-quantile(temp, probs = c(0.95))
  thres.1<-quantile(temp.1, probs = c(0.95))
  print(cate)
  print(thres)
  hist(temp,main=NULL,xlab=paste("Number of peak SNPs for ",cate,sep=""))
  #abline(v=thres,col=c("red","black"))
  abline(v=thres,col=c("red"))
  #hist(temp.1,main=NULL,xlab=paste("Number of peak SNPs for ",cate,sep=""))
  #abline(v=thres.1,col=c("red"))
}
dev.off()

# with window size from 10K, 50K and 200K and step size of 0.2 and 0.5 of the window size:
# Peak SNPs from all traits showed similar distribution; also for each category
# the 90% quantile of all peak SNPs indicate that >2 peaks SNPs in a window could be declared as a hotspot

check<-REGION[which(REGION$ALL>0),]
check$sum<-rowSums(check[,4:9])
check[which(check$ALL!=check$sum),]

###################### Go to "define_HS_region.R" for region bounderis and candidate genes ######
