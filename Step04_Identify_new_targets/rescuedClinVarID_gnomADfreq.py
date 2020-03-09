#!/usr/bin/env python

from __future__ import division
import gzip


def loc(b):
    with open(b,"r") as f:
        d=f.readlines()
    cid=[x.split("\t")[0] for x in d[1:]]
    
    dic={}
    with open("../Output_data/clinvar_20190325_SNV_INDEL1000_mainChr_headerLength_spliceAIpredicted_ds02_pathogenic_junctions_mEX_4BassetSad.vcf","r") as f:
        for x in f:
            ch1= x.split("\t")[0]
            pos1=x.split("\t")[1]
            cid1=x.split("\t")[2]
            if cid1 in cid:
                if ((ch1,pos1) in dic)==False:
                    dic[(ch1,pos1)]={"cid":cid1,"freq":""}
    
    #"gnomAD genome"
    with gzip.open("../Input_data/gnomad.genomes.r2.1.1.sites.vcf.bgz","r") as f:
        for x in f:
            if x[0]!="#":
                ch2= x.split("\t")[0]
                pos2=x.split("\t")[1]
                if ((ch2,pos2) in dic):
                    dic[(ch2,pos2)]["freq"]=x.strip().split("\t")[-1]
    
    #"gnomAD exome"
    with gzip.open("../Input_data/gnomad.exomes.r2.1.1.sites.vcf.bgz","r") as f:
        for x in f:
            if x[0]!="#":
                ch2= x.split("\t")[0]
                pos2=x.split("\t")[1]
                if ((ch2,pos2) in dic):
                    dic[(ch2,pos2)]["freq"]+=",,,"+x.strip().split("\t")[-1]
    
    f=open("../Output_data/sad_table_rescue_alleleFreq_raw.txt","w")
    f.write(d[0][:-1]+"\tFreq\n")
    for x in d[1:]:
        cid3=x.split("\t")[0]
        for y in dic:
            if cid3==dic[y]["cid"]:
                freq=dic[y]["freq"]            
        f.write(x[:-1]+"\t"+freq+"\n")

    f.close()        
    
    return 1
    

def freq(a):
    with open(a,"r") as f:
        dd=f.readlines()

    f=open(a[:-8] + ".txt","w")
    f.write(dd[0])
    for x in dd[1:]:
        t=x.split("\t")[-1]
        if t=="\n":
            freq="-1"
        else:
            l=t.split(",,,")
            if l[0]=="" or len(l)==1:
                freq=t.split("AF=")[1].split(";")[0]
            else:
                ac1=int(l[0].split("AC=")[1].split(";")[0])
                an1=int(l[0].split("AN=")[1].split(";")[0])
                ac2=int(l[1].split("AC=")[1].split(";")[0])
                an2=int(l[1].split("AN=")[1].split(";")[0])
                freq=str((ac1+ac2)/(an1+an2))
                print x.split("\t")[0]
                print ac1,an1,ac2,an2,freq
        content = "\t".join(x.split("\t")[:-1]+[freq])+"\n"
        f.write(content)
    
    f.close()
    
    return 1
    

loc("../Output_data/sad_table_rescue.txt")
freq("../Output_data/sad_table_rescue_alleleFreq_raw.txt")