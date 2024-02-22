### model selection, plot BIC

library(MASS)
library(asreml)

setwd("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/pheno/ModelSel")

## kinship
geno_LD02<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT/centeredIBS_HMP3_AGPv4_LD02_450to310.txt",header=F,sep="\t")

rownames(geno_LD02)<-geno_LD02[,1]
geno_LD02<-geno_LD02[,-1]
colnames(geno_LD02)<-rownames(geno_LD02)

## PCs
P=10
PCs<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT/PCs20_HMP3_prun02_450to310.txt",header=T,row.names=1,sep="\t")
## wax
wax_tr<-read.table("/workdir/ml2498/MaizeLeafCuticle/Wax_TWAS/GAPIT/Transformed_BLUP_wax_imputed_OLRM.txt",header=T,sep="\t")
Taxa<-intersect(rownames(PCs),wax_tr[,1])
wax_tr<-wax_tr[which(wax_tr$MLC_STANDARD %in% Taxa),]
wax_tr$MLC_STANDARD<-as.factor(as.character(wax_tr$MLC_STANDARD))

## FT (310 lines)
FT<-read.table("/workdir/ml2498/MaizeLeafCuticle/TWAS_2018/GAPIT/FT_as_covariate.txt",header=T,sep="\t")
FT<-FT[,c(1,4)]

### function to parse kinship matrix
toSparse <- function(m) {
  comb <- data.frame(row = rep(1:nrow(m), each = nrow(m)),
                     column = rep.int(1:nrow(m), nrow(m)))
  x <- comb[comb$row >= comb$column, ]
  x$value <- m[cbind(x$row, x$column)]
  attr(x, 'rowNames') <- rownames(m)
  return(x)
}
##################################
calc <- function (asreml)
{
  summ <- summary(asreml)
  vc <- summ$varcomp
  DF <- nrow(vc)
  if ("Fixed" %in% levels(vc$constraint))
    DF <- DF - table(vc$constraint)["Fixed"]
  if ("Constrained" %in% levels(vc$constraint))
    DF <- DF - table(vc$constraint)["Constrained"]
  names(DF) <- ""
  logREML <- summ$loglik
  AIC <- -2 * logREML + 2 * DF
  BIC <- -2 * logREML + DF * log(summ$nedf)
  data.frame(DF, AIC, BIC, logREML)
}

####################################
# start model selection
###################################
for (t in 2:ncol(wax_tr)){
  pheno<-wax_tr[!is.na(wax_tr[,t]),c(1,t)]
  
  trait<-colnames(pheno)[2]
  cv<-colnames(FT)[2]
  
  
  taxa<-intersect(pheno[,1],Taxa)
  
  pheno<-pheno[which(pheno[,1] %in% taxa),]
  pheno$MLC_STANDARD<-as.character(pheno$MLC_STANDARD)
  pheno<-pheno[order(pheno$MLC_STANDARD),]
  
  ft<-FT[which(FT[,1] %in% taxa),]
  ft$MLC_STANDARD<-as.character(ft$MLC_STANDARD)
  ft<-ft[order(ft$MLC_STANDARD),]
  
  if(any(ft[,1]!=pheno[,1])){
    break
  }else{
    pheno<-merge(ft,pheno,by="MLC_STANDARD",all=F)
  }
  
  #########################################
  
  
  
  pcs<-PCs[which(rownames(PCs) %in% taxa),]
  pcs<-pcs[order(rownames(pcs)),]
  
  which(rownames(pcs)!=pheno[,1])
  
  counts_PC<-cbind(pheno,pcs)
  counts_PC$MLC_STANDARD<-as.factor(counts_PC$MLC_STANDARD)
  
  kin<-geno_LD02[which(rownames(geno_LD02) %in% taxa),which(colnames(geno_LD02) %in% taxa)]
  kin<-kin[order(rownames(kin)),order(colnames(kin))]
  kin<-as.matrix(kin)
  kin.inv<-ginv(kin)
  kin.ainv<-toSparse(kin.inv)
  attr(kin.ainv,'rowNames')<-rownames(kin)
  
  ###### model selection ##################
  
  modelfitting<-vector()
  
  for (ft in c("","ft")){
    
    if (ft ==""){
      
      for (p in (-1):P){
        if (p==(-1)){
          fit_test <- eval(parse(text=paste("asreml(fixed = ",trait,"~ 1,data=counts_PC,na.method.X='include')",sep="")))
        }
        else if (p==0){
          fit_test <- eval(parse(text=paste("asreml(fixed = ",trait,"~ 1,random= ~ped(MLC_STANDARD,var=T),ginverse=list(MLC_STANDARD=kin.ainv),data=counts_PC,na.method.X='include')",sep="")))
        }else {
          PCvar <- paste("PC", 1:p, sep="")
          #PEERvar<-paste("PEER", 1:j, sep="")
          fit_test <- eval(parse(text=paste("asreml(fixed = ",trait,"~1 +",paste(PCvar,sep="", collapse= "+"),",random= ~ped(MLC_STANDARD,var=T),ginverse=list(MLC_STANDARD=kin.ainv),data=counts_PC,na.method.X='include')",sep="")))
          #fit_test <- eval(parse(text=paste("asreml(fixed = ",trait,"~",GeneID[1],"+",paste(PCvar,sep="", collapse= "+"),",data=counts_PC,na.method.X='include')",sep="")))
        }
        modelfitting<-rbind(modelfitting,calc(fit_test))
        
      }
      
      for (p in 1:P){
        PCvar <- paste("PC", 1:p, sep="")
        #PEERvar<-paste("PEER", 1:j, sep="")
        #fit_test <- eval(parse(text=paste("asreml(fixed = ",trait,"~",GeneID[1],"+",paste(PCvar,sep="", collapse= "+"),",random= ~ped(MLC_STANDARD,var=T),ginverse=list(MLC_STANDARD=kin.ainv),data=counts_PC,na.method.X='include')",sep="")))
        fit_test <- eval(parse(text=paste("asreml(fixed = ",trait,"~ 1 +",paste(PCvar,sep="", collapse= "+"),",data=counts_PC,na.method.X='include')",sep="")))
        
        modelfitting<-rbind(modelfitting,calc(fit_test))
        
      }
      
    }else if (ft=="ft"){
      for (p in (-1):P){
        if (p==(-1)){
          fit_test <- eval(parse(text=paste("asreml(fixed = ",trait,"~ 1+",cv,",data=counts_PC,na.method.X='include')",sep="")))
        }
        else if (p==0){
          fit_test <- eval(parse(text=paste("asreml(fixed = ",trait,"~ 1+",cv,",random= ~ped(MLC_STANDARD,var=T),ginverse=list(MLC_STANDARD=kin.ainv),data=counts_PC,na.method.X='include')",sep="")))
        }else {
          PCvar <- paste("PC", 1:p, sep="")
          #PEERvar<-paste("PEER", 1:j, sep="")
          fit_test <- eval(parse(text=paste("asreml(fixed = ",trait,"~1 +",cv,"+",paste(PCvar,sep="", collapse= "+"),",random= ~ped(MLC_STANDARD,var=T),ginverse=list(MLC_STANDARD=kin.ainv),data=counts_PC,na.method.X='include')",sep="")))
          #fit_test <- eval(parse(text=paste("asreml(fixed = ",trait,"~",GeneID[1],"+",paste(PCvar,sep="", collapse= "+"),",data=counts_PC,na.method.X='include')",sep="")))
        }
        modelfitting<-rbind(modelfitting,calc(fit_test))
        
      }
      
      for (p in 1:P){
        PCvar <- paste("PC", 1:p, sep="")
        #PEERvar<-paste("PEER", 1:j, sep="")
        #fit_test <- eval(parse(text=paste("asreml(fixed = ",trait,"~",GeneID[1],"+",paste(PCvar,sep="", collapse= "+"),",random= ~ped(MLC_STANDARD,var=T),ginverse=list(MLC_STANDARD=kin.ainv),data=counts_PC,na.method.X='include')",sep="")))
        fit_test <- eval(parse(text=paste("asreml(fixed = ",trait,"~ 1 +",cv,"+",paste(PCvar,sep="", collapse= "+"),",data=counts_PC,na.method.X='include')",sep="")))
        
        modelfitting<-rbind(modelfitting,calc(fit_test))
      }
      
    }
  }
  rownames(modelfitting)<-c("simple","K",paste("K+P",1:P,sep=""),paste("P",1:P,sep=""),"FT","K+FT",paste("K+FT+P",1:P,sep=""),paste("FT+P",1:P,sep=""))
  #write.table(modelfitting,paste(align,"_",transf,"_",trait,"_modelfitting_rep",rep,".txt",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
  #write.table(modelfitting,paste(trait,"_modelfitting_310kin_",prun,"_w-wout_FT.txt",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
  write.table(modelfitting,paste(trait,"_modelfitting_310kin_LD02_PC02_w-wout_FT.txt",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
  
}  

#write.table(FT,"FT_as_covariate.txt",col.names=T,row.names=F,sep="\t",quote=F)
#########################
# Plot
#########################
library(ggplot2)
wax_tr<-read.table("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/GAPIT/Transformed_BLUP_wax_imputed.txt",header=T,sep="\t")
Traits<-colnames(wax_tr)[-1]
FA<-Traits[c(1:7,54)]
FFA<-Traits[c(8:18,55)]
HC<-Traits[c(19:27,56)]
AD<-Traits[c(28:31,57)]
WE<-Traits[c(32:45,58)]
AC<-Traits[c(46:53,59)]
Total<-Traits[c(60)]
P=10

setwd("/home/ml2498/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno/ModelSel")
#setwd("/Users/Meng/Desktop/Labserver/MaizeLeafCuticle/Wax_TWAS/pheno/ModelSel")

selectedModel<-vector()
tra<-vector()
for (group in c("FA","FFA","HC","AD","WE","AC","Total")){
  eval(parse(text=paste("traits=",group,sep="")))
  allBIC<-vector()
  
  for (t in traits){
    #BIC<-read.table(paste(t,"_modelfitting_310kin_",prun,"_w-wout_FT.txt",sep=""),header=T,sep="\t")
    BIC<-read.table(paste(t,"_modelfitting_310kin_LD02_PC02_w-wout_FT.txt",sep=""),header=T,sep="\t")
    BIC$sel_model<-"N"
    BIC$sel_model[which.min(BIC$BIC)]<-"Y" # selected model is indicated as triangle
    
    # find the best model and write to a form
    selectedModel<-c(selectedModel,rownames(BIC)[which.min(BIC$BIC)])
    tra<-c(tra,t)
    
    #format BIC for plot
    bic<-BIC[,c("BIC","sel_model")]
    trait<-rep(t,nrow(BIC))
    bic<-cbind(bic,trait)
    allBIC<-rbind(allBIC,bic)
  }
  
  allBIC<-cbind(rep(rownames(BIC),length(traits)),allBIC)
  allBIC<-as.data.frame(allBIC)
  allBIC$BIC<-as.numeric(as.character(allBIC$BIC))
  colnames(allBIC)<-c("model","BIC","sel_model","Trait")
  allBIC$model<-factor(allBIC$model,levels=c("simple","K",paste("K+P",1:P,sep=""),paste("P",1:P,sep=""),"FT","K+FT",paste("K+FT+P",1:P,sep=""),paste("FT+P",1:P,sep="")))
  allBIC$sel_model<-factor(allBIC$sel_model,levels=c("N","Y"))
  limY<-range(allBIC$BIC)
  
  p<-ggplot(allBIC, aes(x = model, y = BIC, group = Trait,color=Trait))+
    geom_point(aes(shape=sel_model,size = sel_model))+
    scale_shape_manual(values=c(1, 17))+
    scale_size_manual(values=c(1,2))+
    geom_line(size=0.2)+
    ylim(limY[1]-50,limY[2]+50)+
    labs(x = "Model",y ="BIC",title='',color="Trait")+
    scale_color_manual(values=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB",
                                "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA"))+
    theme_light()+
    theme(strip.background = element_blank(),
          strip.text.y = element_text(size=5,colour = 'black'),
          strip.text.x = element_text(size=5,colour = 'black'),
          strip.placement = "outside",
          axis.text.x = element_text(angle = 90, hjust = 0,vjust=0.5),
          text = element_text(size=8),
          #legend.title=element_blank(),
          legend.text.align = 0,
          legend.text=element_text(size=rel(0.8)))
  p
  
  #pdf(paste("modelfitting_BIC_310kin_",prun,".pdf",sep=""),height=4,width=7)
  pdf(paste("modelfitting_BIC_450to310kin_",group,".pdf",sep=""),height=4,width=7)
  print(p)
  dev.off()
}

FinalModels<-cbind(tra,selectedModel)
colnames(FinalModels)<-c("trait","BestModel")
write.table(FinalModels,"Final_models_for_waxes.txt",col.names=T,row.names=F,sep="\t",quote=F)

