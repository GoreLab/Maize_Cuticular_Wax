###########################

suppressMessages(library(rpart))
suppressMessages(library(ggplot2))
suppressMessages(library(mlbench))
suppressMessages(library(caret))
suppressMessages(library(party))
suppressMessages(library(pROC))

### Plot: simple regression

setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/RF")
#setwd("/Users/zhenghao/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/RF")
data0<-read.table("PopTag_gc_wax_untr_imputed_forRF.txt",header=T,sep="\t")

wax_name<-read.table("/Users/zhenghao/Desktop/Labserver/OfficeCmp/GoogleDrive/MLC_AZ_2017/uplift_to_v4/wax_GeneExp/raw_wax/wax_cutin_gradient.txt",header=T,sep="\t")
#wax_name<-read.table("/workdir/ml2498/OfficeCmp/GoogleDrive/MLC_AZ_2017/uplift_to_v4/wax_GeneExp/raw_wax/wax_cutin_gradient.txt",header=T,sep="\t")
#wax_name<-read.table("/home/ml2498/Desktop/Labserver/OfficeCmp/GoogleDrive/MLC_AZ_2017/uplift_to_v4/wax_GeneExp/raw_wax/wax_cutin_gradient.txt",header=T,sep="\t")
#wax_name$Chemicals
wax_name<-as.character(wax_name[c(1:18,20:22,24,26,28:35,38:51,55:56,60,61,63,64,67:75),1])
trait_names<-cbind.data.frame(wax_name,colnames(data0)[7:ncol(data0)])
colnames(trait_names)<-c("wax_name","short_name")

table(data0$Category)
##############################
# test parameter combinations for mtry and 
##############################

setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/RF")
data<-read.table("PopTag_gc_wax_untr_imputed_forRF.txt",header=T,sep="\t")

##### handle missing data
for (i in 1:nrow(data)){
  len<-length(data[i,is.na(data[i,7:66])])
  if (len>0){
    print(paste(i,": ",len,sep=""))
  }
}



iter=50
#### generate a column for grouping (by subpop):sampletool ####
data[,2]<-as.character(data[,2])
data<-data[order(data[,2]),]
size_subpop<-table(data[,2])
uniq_tag<-unique(as.character(data[,2]))
which(names(size_subpop)!=uniq_tag)

for (tag in uniq_tag){
  if(tag!="subpop1"){
    # when use fastStructure K=6
    one_four=round(size_subpop[which(names(size_subpop)==tag)]/5,0)
    five=size_subpop[which(names(size_subpop)==tag)]-4*one_four
    
    #eval(parse(text=paste(tag,"_v1=",one_four,sep="")))
    #eval(parse(text=paste(tag,"_v5=",size_subpop[which(names(size_subpop)==tag)]-4*one_four,sep="")))
    eval(parse(text=paste(tag,"=c(",paste(rep(1:4,each=one_four),collapse = ","),",",paste(rep(5,five),collapse = ","),")",sep="")))
  }else{ # subpop1 only has 8 lines, mannually assign group
    subpop1=c(1,2,3,3,4,4,5,5)
  }
}

#
#seed<-sample(1:100000,iter,replace=F) #generate ten seeds for ten iterations

#SAMPLE_tool<-c(1:length(sample_tool))
SAMPLE_tool<-c(1:310)
set.seed(999)
for (l in 1:50){
  # shuffle the sampling index within subpopulation
  for (tag2 in uniq_tag){
    eval(parse(text=paste(tag2,"= sample(",tag2,")",sep="")))
  }
  
  eval(parse(text=paste("sample_tool<-c(",paste(uniq_tag,collapse = ","),")",sep=""))) #stack the shuffled sampling index
  
  SAMPLE_tool<-cbind.data.frame(SAMPLE_tool,sample_tool)
}
colnames(SAMPLE_tool)<-c("index",paste("iteration",1:50,sep=""))

set.seed(813)
for(mtry in c(10,20,30,40,50)){
  for(ntree in c(500,1000,2000)){
    
    print(paste("mtry=",mtry," ntree=",ntree,sep=""))

    Fitted.TE<-vector()
    IMP<-colnames(data)[c(7:66)]
    
    for (j in 1:iter){
      #set.seed(seed[j])
      print(paste("iteration ",j,sep=""))
      
      
      ######
      sample_tool<-SAMPLE_tool[,(l+1)]
      data_wk<-as.data.frame(cbind(data,sample_tool))
      
      fitted.TE<-vector()
      Imp<-colnames(data_wk)[c(7:66)]
      
      ## start 5-fold cross validation
      for (fd in 1:5){
        #print(paste("fold ",fd,sep=""))
        data.te<-data_wk[which(data_wk$sample_tool==fd),]
        data.tr<-data_wk[-which(data_wk$sample_tool==fd),]
        
        fml<-paste("cforest(ce_SD18_both_untr~",paste(colnames(data.tr)[c(7:66)],collapse="+"),
                   ",data=data.tr,controls=cforest_unbiased(ntree=",ntree,",mtry=",mtry,"))",sep="")
        model<-eval(parse(text=fml))
        #model
        
        #fitted.tr <- predict(model, subset(data.tr,select=c(5:67,69:71)), OOB=TRUE, type= "response")
        fitted.te <- predict(model, newdata= subset(data.te,select=c(7:66)),type='response')## not correct
        
        if(j==1){
          combine<-cbind.data.frame(data.te[,c(1:6)],fitted.te)
          fitted.TE<-rbind.data.frame(fitted.TE,combine)
        }else{
          combine<-cbind.data.frame(data.te$MLC_STANDARD,fitted.te)
          colnames(combine)[1]<-"MLC_STANDARD"
          fitted.TE<-rbind.data.frame(fitted.TE,combine)
          
        }
        imp<-varimp(model)
        Imp<-cbind.data.frame(Imp,imp)
        
      }# loop of 5-fold cross validation
      
      ## predicted gc
      if(j==1){
        Fitted.TE<-fitted.TE
      }else{
        Fitted.TE<-merge(Fitted.TE,fitted.TE,by="MLC_STANDARD",all=T) # accumulate fitted values from all five folds
      }
      
      ## importance score
      for (m in 2:6){
        Imp[,m]<-as.numeric(as.character(Imp[,m]))
      }
      
      ## mean of importance score in 6-fold CV
      RF.imp.1<-apply(Imp[,-1],1,mean)
      IMP<-cbind.data.frame(IMP,RF.imp.1)
    }  
    
    colnames(Fitted.TE)[7:ncol(Fitted.TE)]<-paste("predicted_",1:iter,sep="")
    colnames(IMP)<-c("Variable",paste("iteration_",1:iter,sep=""))
    
    #setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/RF")
    setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/RF")
    write.table(Fitted.TE,paste("predicted_gc_5fCV_310lines_mtry",mtry,"ntree",ntree,"_cforest.txt",sep=""),row.names=F,sep="\t",quote=F)
    write.table(IMP,paste("importanceScore_gc_5fCV_310lines_mtry",mtry,"ntree",ntree,"_cforest.txt",sep=""),row.names=F,sep="\t",quote=F)
    
  }
}




##########################################
# calculate RMSE and R2 for each mtry
##########################################
library(data.table)
RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}
setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/RF")


iter=50
for(mtry in c(10,20,30,40,50)){
  for(ntree in c(500,1000,2000)){
    
    Fitted.TE<-fread(paste("predicted_gc_5fCV_310lines_mtry",mtry,"ntree",ntree,"_cforest.txt",sep=""),data.table=F)
    
    #RMSE.TE<-vector()
    #R.TE<-vector()
    
    RMSE.te<-vector()
    R.te<-vector()
    for(k in 7:ncol(Fitted.TE)){
      pred0<-Fitted.TE[,k]
      #obs0<-Fitted.TE[,s]
      obs0<-Fitted.TE[,6]
      rm1<-which(is.na(pred0))
      rm2<-which(is.na(obs0))
      rm<-unique(c(rm1,rm2))
      if(length(rm)>0){
        pred<-pred0[-rm]
        obs<-obs0[-rm]
      }else{
        pred<-pred0
        obs<-obs0
      }
      
      rmse.te<-RMSE(obs,pred)
      RMSE.te<-c(RMSE.te,rmse.te)
      
      cor.te<-cor(Fitted.TE[,k],Fitted.TE[,6],use="complete")
      R.te<-c(R.te,cor.te)
    }
    #RMSE.TE<-cbind(RMSE.TE,RMSE.te)
    #R.TE<-cbind(R.TE,R.te)
    
    accuracy<-cbind.data.frame(RMSE.te,R.te,c(paste("iteration_",1:iter,sep="")))
    #colnames(accuracy)<-c("RMSE_AZ","RMSE_SD","RMSE_All","RMSE_SD18","R2_AZ","R2_SD","R2_All","R2_SD18","Iteration")
    colnames(accuracy)<-c(paste("RMSE_mtry",mtry,"ntree",ntree,sep=""),paste("R2_mtry",mtry,"ntree",ntree,sep=""),"Iteration")
    
    if(mtry==10 & ntree==500){
      ACCU<-accuracy
    }else{
      ACCU<-cbind.data.frame(ACCU,accuracy)
    }
    
  }
}

setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/RF")
#write.table(accuracy,"accuracy_gc_5fCV_310lines_mtry12_cforest.txt",row.names=F,sep="\t",quote=F)
write.table(ACCU,"accuracy_gc_5fCV_310lines_ALLmtry_ntree_cforest.txt",row.names=F,sep="\t",quote=F)

for(i in 1:6){
  print(range(accuracy[,i]))
  print(mean(accuracy[,i]))
}

##############################

## Importance
library(reshape2)
library(ggplot2)
#mtry=12
mtry=10
ntree=2000

#setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno")
setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno")
wax_name<-read.table("wax_full_short_names_conversion.txt",header=T,sep="\t")

#setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/RF")
#wax_name<-read.table("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/TWAS_2018/stomata/short_vs_complete_wax_names.txt",header=T,sep="\t")

setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/RF")
#setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno")

IMP<-read.table(paste("importanceScore_gc_5fCV_310lines_mtry",mtry,"ntree",ntree,"_cforest.txt",sep=""),head=T,sep="\t")
#IMP[-which(IMP$Variable %in% wax_name$short_name),]

IMP<-merge(wax_name,IMP,by.x="short_name",by.y="Variable",all.y=T)
IMP<-IMP[,-1]
colnames(IMP)[1]<-"Chemicals"
IMP <- melt(IMP, id.vars=c("Chemicals"),value.name = "Importance")
IMP$Importance<-as.numeric(as.character(IMP$Importance))

## update wax names ##
library(data.table)
wax<-fread("/Users/ml2498/Library/CloudStorage/GoogleDrive-lemonsquirrel1987@gmail.com/My Drive/Wax_manuscript/Tables/draft/wax_name_Isabel_20230215.txt",data.table=F)
for(i in 1:nrow(IMP)){
  trait<-IMP$Chemicals[i]
  IMP$Chemicals[i]<-wax$`Cuticular wax`[which(wax$old_wax2==trait)]
}
### end of updating ###

IMP$Chemicals<-as.factor(IMP$Chemicals)


#imp<-aggregate(Importance~Chemicals,data=IMP,median)
imp<-aggregate(Importance~Chemicals,data=IMP,mean)
imp<-imp[order(imp$Importance),]
wax<-as.character(imp$Chemicals)
IMP$Chemicals<-factor(IMP$Chemicals,levels=c(paste(wax,sep=",")))

#pdf(paste("Importance_5foldCV_mtry",mtry,"_cforest.pdf",sep=""),height=5,width=10)
#pdf(paste("Importance_5foldCV_mtry",mtry,"_ntree",ntree,"_cforest.pdf",sep=""),height=5,width=10,family="sans")
#pdf(paste("Importance_5foldCV_mtry",mtry,"_ntree",ntree,"_cforest_updated.pdf",sep=""),height=5,width=10,family="sans")
pdf(paste("Importance_5foldCV_mtry",mtry,"_ntree",ntree,"_cforest_updated_rankMean.pdf",sep=""),height=5,width=10,family="sans")

print(ggplot(IMP,aes(Chemicals,Importance))+geom_boxplot()+
        theme_light()+
        labs(x="Adult maize leaf cuticular wax traits",y="Importance score")+
        #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
        theme_light()+
        theme(strip.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              strip.text.y = element_text(size=5,colour = 'black'),
              strip.text.x = element_text(size=5,colour = 'black'),
              strip.placement = "outside",
              axis.text.x = element_text(angle = 90, hjust = 1,vjust=0),
              text = element_text(size=8),
              #legend.title=element_blank(),
              legend.text.align = 0,
              legend.text=element_text(size=rel(0.8)))
      
)
dev.off()

##### generate table ###
mtry=10
ntree=2000

setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/RF")
IMP_0<-read.table(paste("importanceScore_gc_5fCV_310lines_mtry",mtry,"ntree",ntree,"_cforest.txt",sep=""),head=T,sep="\t")
setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno")
wax_name<-read.table("wax_full_short_names_conversion.txt",header=T,sep="\t")

IMP_0<-merge(wax_name,IMP_0,by.x="short_name",by.y="Variable",all.y=T)

IMP<-IMP_0
rownames(IMP)<-IMP[,2]
IMP<-IMP[,-c(1:2)]
stat_mean<-apply(IMP,1,mean)
stat_var<-apply(IMP,1,var)
IMP_sum<-cbind.data.frame(IMP_0[,2],stat_mean,stat_var)
colnames(IMP_sum)<-c("Cuticular features","Mean of importance score","Variance of importance score")

#pearson<-read.table("/Users/Meng/Desktop/LabServer/MaizeLeafCuticle/TWAS_2018/stomata/plot/wax_stomata_ce_cor_naive.txt",header=T,sep="\t")
#pearson<-cbind.data.frame(rownames(pearson),pearson[,c(1,2)])
#colnames(pearson)<-c("Cuticular features","Pearson's correlation","P-value of correlation")
#IMP_Pearson<-merge(pearson,IMP_sum,by="Cuticular features",all=T)

setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/RF")
write.table(IMP_sum,"TableSx_importance_score_RF.txt",row.names=F,sep="\t",quote=F)
###### Compare importance rank ########################
setwd("/Users/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/RF")
comp<-read.table("compare_importance_rank.txt",header=T,sep="\t")
str(comp)
cor.res <- cor.test(comp$Rank310,comp$Rank51,method='pearson',use="pairwise.complete.obs")
cor.res$p.value # 0.0001038084


###### Ends here ########################

