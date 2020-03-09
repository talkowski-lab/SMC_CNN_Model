library(ggplot2)
source("mantelhaen_test.R")

myT=function(x){
    ith=c(1,3,5,7,9,11)
    psi_l1=mean(x[ith]/(x[ith]+x[ith+12]))
    psi_l2=mean(x[ith+1]/(x[ith+1]+x[ith+1+12]))
    exp=c(x[ith]+x[ith+12],x[ith+1]+x[ith+1+12])
    psi_ar=array(c(x[1],x[2],x[13],x[14],
                   x[3],x[4],x[15],x[16],
                   x[5],x[6],x[17],x[18],
                   x[7],x[8],x[19],x[20],
                   x[9],x[10],x[21],x[22],
                   x[11],x[12],x[23],x[24]),dim=c(2,2,6))
    #print(unname(exp))
    if (min(exp) >=5){pv=mantelhaen_test(psi_ar)$p.value
                     ch=unname(mantelhaen_test(psi_ar)$statistic)
                     if (is.na(pv)){pv=1;ch=0}}else{pv=ch=NaN}
    return(c(min(exp),psi_l1,psi_l2,psi_l1-psi_l2,ch,pv))
} 


mySummary=function(x){
  ith=c(1,3,5,7,9,11)
  N1=length(ith)/2
  N2=length(ith)/2
  psi1=x[ith]/(x[ith]+x[ith+12])
  psi2=x[ith+1]/(x[ith+1]+x[ith+1+12])
  u1=mean(psi1,na.rm=T)
  u2=mean(psi2,na.rm=T)
  sd1=sd(psi1,na.rm=T)
  sd2=sd(psi2,na.rm=T)
  se1=sd1/sqrt(N1)
  se2=sd2/sqrt(N2)
  ciM1=qt(0.95/2 + 0.5, N1)
  ciM2=qt(0.95/2 + 0.5, N2)
  ci1=se1*ciM1
  ci2=se2*ciM2
  
  return(c(N1,u1,sd1,se1,ci1,N2,u2,sd2,se2,ci2))
} 


a=read.table("../Input_data/SEMatrix_15477onWT.txt",header=T,row.names=1)
a=as.matrix(a)

sm=t(apply(a,1,FUN=function(x) mySummary(x)))
colnames(sm)=c("sampleSize_15477","average_psi_15477","sd_psi_15477","se_psi_15477","ci_psi_15477","sampleSize_ctrl","average_psi_ctrl","sd_psi_ctrl","se_psi_ctrl","ci_psi_ctrl")
sm=as.data.frame(sm)
write.table(sm,file="../Output_data/psiSummary_15477onWT.txt",sep="\t",row.names=T,col.names=T,quote=F)

res=t(apply(a,1,FUN=function(x) myT(x)))
colnames(res)=c("min_splicing","average_psi_15477","average_psi_ctrl","average_psi_change","chi-sq","p.value")
res=as.data.frame(res)
write.table(res,file="../Output_data/psiChange_15477onWT.txt",sep="\t",row.names=T,col.names=T,quote=F)

res0=res
res=res[which(res$min_splicing>=5),]

ne=setdiff(rownames(res0),rownames(res))
res$fdr=p.adjust(res$p.value,method="BH")
res$color="black"
dn0   =rownames(res[which(res$average_psi_change<=-0.1 & res$min_splicing>=5 & res$p.value<0.05),])
dn1   =rownames(res[which(res$average_psi_change<=-0.1 & res$min_splicing>=5 & res$fdr<0.1),])
up0   =rownames(res[which(res$average_psi_change>= 0.1 & res$min_splicing>=5 & res$p.value<0.05),])
up1   =rownames(res[which(res$average_psi_change>= 0.1 & res$min_splicing>=5 & res$fdr<0.1),])

res[up1,"color"]="red"
res[dn1,"color"]="blue"

write.table(res,file="../Input_data/psiChange_15477onWT_exp.txt",sep="\t",row.names=T,col.names=T,quote=F)


ggplot(res,aes(x=average_psi_change,y=-log10(fdr),color=color))+geom_point()+geom_noBG()+scale_colour_manual(values=c("black","blue","red"))+geom_vline(xintercept=c(-0.1,0.1),linetype="dashed")+geom_hline(yintercept=-log10(0.1),linetype="dashed")
ggsave("../Output_figures/volcano_fdr.pdf")

nc_up=rownames(res[which(abs(res$average_psi_change)<0.01 & (res$average_psi_ctrl<0.9 | res$average_psi_15477<0.9)),])
nc_dn=rownames(res[which(abs(res$average_psi_change)<0.01 & (res$average_psi_ctrl>0.1 | res$average_psi_15477>0.1)),])
nc = intersect(nc_dn,nc_up)
rf = rownames(res[which(res$average_psi_ctrl>=0.99 & res$average_psi_15477>=0.99),])
fl = rownames(res[which(res$average_psi_ctrl<=0.01 & res$average_psi_15477<=0.01),])
gy = setdiff(rownames(res),c(rf,fl,nc,up1,dn1))
write.table(res[gy,],file="../Output_data/psiChange_15477onWT_exp_grey.txt",sep="\t",row.names=T,col.names=T,quote=F)

set.seed(122)
nc_up1=sample(nc_up,length(up1))
set.seed(122)
nc_up0=sample(nc_up,length(up0))
set.seed(122)
nc_dn1=sample(nc_dn,length(dn1))
set.seed(122)
nc_dn0=sample(nc_dn,length(dn0))

write.table(as.data.frame(c(up1)),file="../Output_data/up01tripletNames_pos_fdr01.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(as.data.frame(c(dn1)),file="../Output_data/dn01tripletNames_pos_fdr01.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(as.data.frame(c(up0)),file="../Output_data/up01tripletNames_pos_p005.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(as.data.frame(c(dn0)),file="../Output_data/dn01tripletNames_pos_p005.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(as.data.frame(c(nc_up)),file="../Output_data/up01tripletNames_neg_all.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(as.data.frame(c(nc_dn)),file="../Output_data/dn01tripletNames_neg_all.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(as.data.frame(c(nc_up1)),file="../Output_data/up01tripletNames_neg_254.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(as.data.frame(c(nc_dn1)),file="../Output_data/dn01tripletNames_neg_680.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(as.data.frame(c(nc_up0)),file="../Output_data/up01tripletNames_neg_377.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(as.data.frame(c(nc_dn0)),file="../Output_data/dn01tripletNames_neg_970.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(as.data.frame(c(nc)),file="../Output_data/nc01tripletNames_all_382.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(as.data.frame(c(rf)),file="../Output_data/rf01tripletNames_all_138816.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(as.data.frame(c(fl)),file="../Output_data/fl01tripletNames_all_3515.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(as.data.frame(c(ne)),file="../Output_data/ne01tripletNames_all_173550.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(as.data.frame(c(gy)),file="../Output_data/gy01tripletNames_all_17450.txt",sep="\t",row.names=F,col.names=F,quote=F)

#./basset_train.lua -job /data/talkowski/dg520/projects/fd_embryo_fibro/output/fibro/STARpass2nd/ml_new/up01_params.txt -save /data/talkowski/dg520/projects/fd_embryo_fibro/output/fibro/STARpass2nd/ml_new/up01_cnn /data/talkowski/dg520/projects/fd_embryo_fibro/output/fibro/STARpass2nd/ml_new/up01_train.h5
#./basset_train.lua -job /data/talkowski/dg520/projects/fd_embryo_fibro/output/fibro/STARpass2nd/ml_new/dn01_params.txt -save /data/talkowski/dg520/projects/fd_embryo_fibro/output/fibro/STARpass2nd/ml_new/dn01_cnn /data/talkowski/dg520/projects/fd_embryo_fibro/output/fibro/STARpass2nd/ml_new/dn01_train.h5

#./basset_test.lua /data/talkowski/dg520/projects/fd_embryo_fibro/output/fibro/STARpass2nd/ml_new/up01_cnn_best.th /data/talkowski/dg520/projects/fd_embryo_fibro/output/fibro/STARpass2nd/ml_new/up01_train.h5 /data/talkowski/dg520/projects/fd_embryo_fibro/output/fibro/STARpass2nd/ml_new/up01_test
#./basset_test.lua /data/talkowski/dg520/projects/fd_embryo_fibro/output/fibro/STARpass2nd/ml_new/dn01_cnn_best.th /data/talkowski/dg520/projects/fd_embryo_fibro/output/fibro/STARpass2nd/ml_new/dn01_train.h5 /data/talkowski/dg520/projects/fd_embryo_fibro/output/fibro/STARpass2nd/ml_new/dn01_test

#./basset_motifs_predict.lua /data/talkowski/dg520/projects/fd_embryo_fibro/output/fibro/STARpass2nd/ml_new/up01_cnn_best.th /data/talkowski/dg520/projects/fd_embryo_fibro/output/fibro/STARpass2nd/ml_new/up01_train.h5 /data/talkowski/dg520/projects/fd_embryo_fibro/output/fibro/STARpass2nd/ml_new/up01_motif_predicted.h5
#./basset_motifs.py -s 100 -t -d /data/talkowski/dg520/projects/fd_embryo_fibro/output/fibro/STARpass2nd/ml_new/up01_motif_predicted.h5 -o /data/talkowski/dg520/projects/fd_embryo_fibro/output/fibro/STARpass2nd/ml_new/up01_motif /data/talkowski/dg520/projects/fd_embryo_fibro/output/fibro/STARpass2nd/ml_new/up01_cnn_best.th /data/talkowski/dg520/projects/fd_embryo_fibro/output/fibro/STARpass2nd/ml_new/up01_train.h5
