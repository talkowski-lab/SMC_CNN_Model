library(ggplot2)
library(plotROC)
library(pROC)

#setwd("") #This needs to be set to $MY_FOLDER/Step02_Build_CNN_model

checkM=function(a,th){
    tmp=strsplit(a,"\\.")[[1]]
    nf=paste(tmp[1:(length(tmp)-2)],collapse=".")
    fa=scan(nf,what="character",sep="\n")
    fi=(1:(length(fa)/2))*2-1
    fh=fa[fi]
    rn=sapply(strsplit(fh,">"),"[[",2)
    res=read.table(a)
    rownames(res)=rn
    colnames(res)=c("I","E","S")
    if (th=="local"){exp=sapply(strsplit(rn,"::"),"[[",2)}
    else{exp=th}
    #pred=apply(res,1,FUN=function(x){if(max(x)==x[1]){return("Inclusion")}else if(max(x)==x[2]){return("Exclusion")}else{return("Stable")}})
    res$Exp=factor(exp)
    new=cbind(res[1:4], sapply(levels(res$Exp), function(x) as.integer(x == res$Exp)))
    #res$Pred=pred
    return(new)
}

pi=38
pfx=".4Junc"
m=checkM(paste0("../Output_data/3C.test0.fa.param",pi,".prob"), "local")

roc1 = roc(m$Inclusion, m$I, percent=T, show.thres=T, print.auc=T)
roc2 = roc(m$Exclusion, m$E, percent=T, show.thres=T, print.auc=T)
roc3 = roc(m$Stable, m$S, percent=T, show.thres=T, print.auc=T)

d1=data.frame(Specificity=roc1$specificities/100,Sensitivity=roc1$sensitivities/100,Threshold=roc1$thresholds)
d2=data.frame(Specificity=roc2$specificities/100,Sensitivity=roc2$sensitivities/100,Threshold=roc2$thresholds)
d3=data.frame(Specificity=roc3$specificities/100,Sensitivity=roc2$sensitivities/100,Threshold=roc3$thresholds)
plt1=ggplot(d1,aes(x=Specificity,y=Sensitivity))+geom_line(color="red")+geom_line(data=d2,aes(x=Specificity,y=Sensitivity),color="blue")+geom_line(data=d3,aes(x=Specificity,y=Sensitivity),color="black")+geom_abline(intercept = 1, slope = -1, color="grey")+geom_noBG()
plt2=plt1+geom_text(x=0.2,y=0.2,label=paste0("AUC for Inclusion Vote: ",round(roc1$auc/100,2),"\nAUC for Exclusion Vote: ",round(roc2$auc/100,2),"\nAUC for Unchanged Vote: ",round(roc2$auc/100,2)))
plot(plt2)
ggsave(paste0("../Output_figures/aucs_param",pi,".pdf"))

roc1_sp_ind=tail(which(round(roc1$specificities)==95),1)
roc2_sp_ind=tail(which(round(roc2$specificities)==95),1)
roc3_sp_ind=tail(which(round(roc3$specificities)==95),1)
c1=coords(roc1, roc1$specificities[roc1_sp_ind], input="specificity")["threshold"]
c2=coords(roc2, roc2$specificities[roc2_sp_ind], input="specificity")["threshold"]
c3=coords(roc3, roc3$specificities[roc3_sp_ind], input="specificity")["threshold"]

c1 #cutoff for inclusion
c2 #cutoff for exclusion
c3 #cutoff for unchanged