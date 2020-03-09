#!/usr/bin/env python

from __future__ import division
import numpy.random as npr 
import shlex, subprocess
import commands


def extract_info(a):
    al=a.split("ALLELEID=")[1].split(";")[0]
    if "CLNDN=" in a:
        ds=a.split("CLNDN=")[1].split(";")[0]
    else:
        ds=a.split("CLNDNINCL=")[1].split(";")[0]
    if "CLNDISDB=" in a:
        db=a.split("CLNDISDB=")[1].split(";")[0]
    else:
        db=a.split("CLNDISDBINCL=")[1].split(";")[0]
    if "CLNSIG=" in a:
        patho=a.split("CLNSIG=")[1].split(";")[0]
    else:
        patho=a.split("CLNSIGINCL=")[1].split(";")[0]
    if "MC=" in a:
        mc=a.split("MC=")[1].split(";")[0]
    else:
        mc="NA"
    
    return [patho,mc,al,ds,db]
    
    
def pub(var_citations):
    liter={}
    with open(var_citations,"r") as f:
        d=f.readlines()
    for x in d[1:]:
        t=x.strip().split("\t")
        pid=t[-1]
        src=t[-2]
        ids=t[1]
        if src=="PubMed":
            if (ids in liter)==False:
                liter[ids]=[pid,1]
            else:
                liter[ids][0]+=","+pid
                liter[ids][1]+=1
        
    return liter
    

def cp(sad_inp):
    pp={}
    with open(sad_inp,"r") as ff:
        dd=ff.readlines()
    for rec in dd:
        tt=rec.strip().split("\t")    
        if ((tt[2],tt[5]) in pp) == False:
            pp[(tt[2],tt[5])]=[tt[-3],tt[-2],tt[-1]]
        else:
            print "Error: Record duplication"
            
    return pp

    
def incl_filter(a,piece,c0,c1,c2):
    outSAD=a[:-4]+"_inclusion.txt"
    dicM=pub("../Input_data/var_citations_20190325.txt")
    dicP=cp("../Output_data/clinvar_20190325_SNV_INDEL1000_mainChr_headerLength_spliceAIpredicted_ds02_pathogenic_junctions_mEX_4BassetSad.vcf")
    with open(a,"r") as f:
        d=f.readlines()
    f2=open(outSAD,"w")
    f2.write(d[0].strip()+"\tPathogenic\tMC\tAlleleID\tDisease\tDisease_DB\tNum_paper\tPubMedID\n")
    i=1
    while i<len(d):
        hit=0
        cid,tid,dis,ref,mut,target0,pred_ref0,pred_mut0,chg0=d[i].strip().split("\t")
        cid,tid,dis,ref,mut,target1,pred_ref1,pred_mut1,chg1=d[i+1].strip().split("\t")
        cid,tid,dis,ref,mut,target2,pred_ref2,pred_mut2,chg2=d[i+2].strip().split("\t")
        if (float(pred_mut0)>=c0 and float(pred_mut0)>float(pred_mut1) and float(pred_mut0)>float(pred_mut2)):
            hit+=1
        if "seqDiff" in dis:
            hit+=1
        if ("mEXdrop3" in dis) and ("nonsense" in dis)==False and ("frameshift" in dis)==False and ("junctGone" in dis)==False:
            hit+=1
        if cid in dicM:
            hit+=1
        if hit ==4:
            info=extract_info(dicP[(cid,tid)][0])
            f2.write(d[i].strip()+"\t"+"\t".join(info)+"\t"+dicM[cid][0]+"\t"+str(dicM[cid][1])+"\n")
            
        i+=3
    f2.close()
    
    return 1
    
    
def excl_filter(a,piece,c0,c1,c2):
    outSAD=a[:-4]+"_exclusion.txt"
    dicM=pub("../Input_data/var_citations_20190325.txt")
    dicP=cp("../Output_data/clinvar_20190325_SNV_INDEL1000_mainChr_headerLength_spliceAIpredicted_ds02_pathogenic_junctions_mEX_4BassetSad.vcf")
    with open(a,"r") as f:
        d=f.readlines()
    f2=open(outSAD,"w")
    f2.write(d[0].strip()+"\tPathogenic\tMC\tAlleleID\tDisease\tDisease_DB\tNum_paper\tPubMedID\n")
    i=1
    while i<len(d):
        hit=0
        cid,tid,dis,ref,mut,target0,pred_ref0,pred_mut0,chg0=d[i].strip().split("\t")
        cid,tid,dis,ref,mut,target1,pred_ref1,pred_mut1,chg1=d[i+1].strip().split("\t")
        cid,tid,dis,ref,mut,target2,pred_ref2,pred_mut2,chg2=d[i+2].strip().split("\t")
        seq0=dicP[(cid,tid)][-2]
        seq1=dicP[(cid,tid)][-1]
        start =int(tid.split("@")[2].split("^")[1].split("&")[0])
        end   =int(tid.split("@")[2].split("&")[1].split("^")[0])-1
        if (float(pred_mut0)<float(pred_mut1) and float(pred_mut1)>=c1 and float(pred_mut2)<float(pred_mut1)):
            hit+=1
        if "seqDiff" in dis:
            hit+=1
        if (("mEXmore" in dis) and ("nonsense" in dis)==False and ("frameshift" in dis)==False and ("junctGone" in dis)==False) or (((end-start)%3)==0 and ("nonsense_mEX" in dis)):
            hit+=1
        if cid in dicM:
            hit+=1
        if hit ==4:
            info=extract_info(dicP[(cid,tid)][0])
            f2.write(d[i+1].strip()+"\t"+"\t".join(info)+"\t"+dicM[cid][0]+"\t"+str(dicM[cid][1])+"\n")
            
        i+=3
    f2.close()
    
    return 1
    
def doit(inp,region,cut0,cut1,cut2):
    incl_filter(inp,region,cut0,cut1,cut2)
    excl_filter(inp,region,cut0,cut1,cut2)
    
    return 1


doit("../Output_data/sad_table.txt",
     4,
     0.225832171738145,
     0.803149670362475,
     0.36993284523487)
