#!/usr/bin/env python

from __future__ import division
import numpy.random as npr 
import sys, shlex, subprocess
import commands
import pysam
from Bio.Seq import reverse_complement
import numpy as np


def JuncDic(a,inside=25,outside=75):
    #all the coordinates are 1-indexed and in [,] manner in this function
    dic={}
    with open(a,"r") as f:
        for x in f:
            t=x.strip()
            gn=t.split("@")[1]
            ch=t.split("@")[2].split(":")[0]
            st=t.split("@")[2].split(":")[1]
            U_end   =int(t.split("@")[2].split(":")[2].split("^")[0])-1
            X_start =int(t.split("@")[2].split("^")[1].split("&")[0])+1
            X_end   =int(t.split("@")[2].split("&")[1].split("^")[0])-1
            D_start =int(t.split("@")[2].split("|")[0].split("^")[-1])+1
            if gn in dic:
                if t in dic[gn]:
                    print "Error: Duplicated spliceID"
                    break
            else:
                dic[gn]={}
            dic[gn][t]={"junc":[U_end,X_start,X_end,D_start],
                        "region":[(U_end-inside+1,U_end+outside),(X_start-outside,X_start+inside-1),(X_end-inside+1,X_end+outside),(D_start-outside,D_start+inside-1)]}

    return dic


def is_influenced(x1,x2,y1,y2):
    return max(x1,y1) <= min(x2,y2)
    

def getSeqIfInfluenced(this,tid,allseqs,vcf_info,ai_cons,inside=25,outside=75):
    #all the coordinates are 0-indexed and in [,) manner in this function; e.g. BED format
    returnValue=(0,"","")
    chrom,site,rs,wt,mu,dot1,dot2,attr=this.strip().split("\t")
    mp=int(site)-1
    infl_len=max(len(mu),len(wt))
    infl_offset=len(mu)-len(wt)
    mp_start=mp
    mp_end=mp+infl_len

    strand = tid.split("@")[2].split(":")[1]
    intron1_start_bed  =int(tid.split("@")[2].split(":")[2].split("^")[0])-1
    intron1_end_bed    =int(tid.split("@")[2].split("^")[1].split("&")[0])
    intron2_start_bed  =int(tid.split("@")[2].split("&")[1].split("^")[0])-1
    intron2_end_bed    =int(tid.split("@")[2].split("|")[0].split("^")[-1])
    
    coord_new= []
    gone_index= []
    junct_pos=[intron1_start_bed,intron1_end_bed-2,intron2_start_bed,intron2_end_bed-2]
    #for si,ei,seqi in zip([s0,s1,s2,s3],[e0-1,e1-1,e2-1,e3-1],[seq0,seq1,seq2,seq3]):
        #if is_influenced(mp_start,mp_end-1,si,ei):
    for pi,ji in zip(junct_pos,range(len(junct_pos))):
        if is_influenced(pi,pi+1,mp_start,mp_end-1):
            gone_index.append(ji)
        else:
            #upstream to the mutation
            if (pi+1)<mp_start:
                coord_new.append(pi)
            #downstream to the mutation
            else:
                coord_new.append(pi+infl_offset)

        
    if len(gone_index)==0:
        wt_s0   = intron1_start_bed - inside
        wt_e0   = intron1_start_bed + outside
        wt_s1   = intron1_end_bed   - outside
        wt_e1   = intron1_end_bed   + inside
        wt_s2   = intron2_start_bed - inside
        wt_e2   = intron2_start_bed + outside
        wt_s3   = intron2_end_bed   - outside
        wt_e3   = intron2_end_bed   + inside
    
        mu_s0   = coord_new[0]      - inside
        mu_e0   = coord_new[0]      + outside
        mu_s1   = coord_new[1]+2    - outside
        mu_e1   = coord_new[1]+2    + inside
        mu_s2   = coord_new[2]      - inside
        mu_e2   = coord_new[2]      + outside
        mu_s3   = coord_new[3]+2    - outside
        mu_e3   = coord_new[3]+2    + inside
        
        wt_fa=allseqs[chrom]
        wt_seq0=wt_fa[wt_s0:wt_e0]
        wt_seq1=wt_fa[wt_s1:wt_e1]
        wt_seq2=wt_fa[wt_s2:wt_e2]
        wt_seq3=wt_fa[wt_s3:wt_e3]
        
        mu_fa=wt_fa[0:mp_start]+mu+wt_fa[mp_start+len(wt):]
        mu_seq0=mu_fa[mu_s0:mu_e0]
        mu_seq1=mu_fa[mu_s1:mu_e1]
        mu_seq2=mu_fa[mu_s2:mu_e2]
        mu_seq3=mu_fa[mu_s3:mu_e3]
    
        if strand == "+":
            wt_seq = wt_seq0 + wt_seq1 + wt_seq2 + wt_seq3
            mu_seq = mu_seq0 + mu_seq1 + mu_seq2 + mu_seq3
            anchors = [intron1_start_bed-1,intron1_end_bed,intron2_start_bed-1,intron2_end_bed]    
        else:
            wt_seq0=reverse_complement(wt_seq0)
            wt_seq1=reverse_complement(wt_seq1)
            wt_seq2=reverse_complement(wt_seq2)
            wt_seq3=reverse_complement(wt_seq3)
            mu_seq0=reverse_complement(mu_seq0)
            mu_seq1=reverse_complement(mu_seq1)
            mu_seq2=reverse_complement(mu_seq2)
            mu_seq3=reverse_complement(mu_seq3)
            wt_seq = wt_seq3 + wt_seq2 + wt_seq1 + wt_seq0
            mu_seq = mu_seq3 + mu_seq2 + mu_seq1 + mu_seq0
            anchors = [intron1_start_bed-1,intron1_end_bed,intron2_start_bed-1,intron2_end_bed]
           
        if ("splice_donor_variant" in vcf_info)==False and ("splice_acceptor_variant" in vcf_info)==False and ("frameshift_variant" in vcf_info)==False and ("nonsense" in vcf_info)==False:
            mc_res = ",".join([mc for mc in ["missense_variant","synonymous_variant","intron_variant"] if mc in vcf_info])
        elif ("splice_donor_variant" in vcf_info) or ("splice_acceptor_variant" in vcf_info):
            mc_res = "junctGone_clinvar"
        elif ("frameshift_variant" in vcf_info):
            mc_res = "frameshift"
        elif ("nonsense" in vcf_info):
            if (mp_start >= intron1_end_bed) and (mp_start < intron2_start_bed):
                mc_res = "nonsense_mEX"
            else:
                mc_res = "nonsense_elsewhere"
                
        #anchor points are the coordinates of the most outside nucleotide of the middle exon
        dis2anchors  = [mp-anchor for anchor in anchors]
        dis2anchors1 = [abs(dis2anchor) for dis2anchor in dis2anchors]
        dis_index = dis2anchors1.index(min(dis2anchors1))
        dis_final = dis2anchors[dis_index]
        
        if strand=="+":
            if wt_seq!=mu_seq:
                returnValue= ("%s|Junction_%d|%d|%s|%s|seqDiff"%(strand,dis_index+1,dis_final,ai_cons,mc_res),wt_seq,mu_seq)
            else:
                returnValue= ("%s|Junction_%d|%d|%s|%s|seqSame"%(strand,dis_index+1,dis_final,ai_cons,mc_res),wt_seq,mu_seq)
        else:
            dis_final = dis_final * (-1)
            if wt_seq!=mu_seq:
                returnValue= ("%s|Junction_%d|%d|%s|%s|seqDiff"%(strand,4-dis_index,dis_final,ai_cons,mc_res),wt_seq,mu_seq)
            else:
                returnValue= ("%s|Junction_%d|%d|%s|%s|seqSame"%(strand,4-dis_index,dis_final,ai_cons,mc_res),wt_seq,mu_seq)
    
    else:
        mc_res="junctGone_novel"
        wt_seq=mu_seq="A"*400
        if strand=="+":
            gone_index=[gi+1 for gi in gone_index]
        else:
            gone_index=[4-gi for gi in gone_index]
        returnValue= ("%s|Junction_%d|%f|%s|%s|seqNA"%(strand,min(gone_index),np.nan,ai_cons,mc_res),wt_seq,mu_seq)
    
    return returnValue


def overlap(a,DicJunc,DicFA):
    #all the coordinates are 1-indexed and in [,] manner in this function
    res={"overlap4vcf":[],"overlap4Basset":[]}
    hit=0
    BassetLines=[]
    ch,pos,clv,ref,alt,null1,null2,info=a.strip().split("\t")
    mp=int(pos)-1
    infl_len=max(len(alt),len(ref))
    mp_start=mp
    mp_end=mp+infl_len-1
    info_main=info.strip().split("SpliceAI=")[0]
    info_ai=""
    preds=info.strip().split("SpliceAI=")[-1].split(",")
    for pred in preds:
        allele_alt,gene_symbol,acceptor_gain,acceptor_loss,donor_gain,donor_loss,offset_acceptor_gain,offset_acceptor_loss,offset_donor_gain,offset_donor_loss = pred.split("|")
        acceptor_loss = float(acceptor_loss) * (-1)
        donor_loss = float(donor_loss) * (-1)
        acceptor_gain = float(acceptor_gain)
        donor_gain = float(donor_gain)
        pos_acceptor_loss = int(pos)+int(offset_acceptor_loss)
        pos_donor_loss = int(pos)+int(offset_donor_loss)
        pos_acceptor_gain = int(pos)+int(offset_acceptor_gain)
        pos_donor_gain = int(pos)+int(offset_donor_gain)
        ai_scores={pos_acceptor_gain:acceptor_gain,pos_acceptor_loss:acceptor_loss,pos_donor_gain:donor_gain,pos_donor_loss:donor_loss}
        if gene_symbol in DicJunc:
            for x in DicJunc[gene_symbol]:
                #if a splicing-changing site overlapped with a triplet junction with a significant value
                #overlap_res = [(k in [pos_acceptor_gain, pos_acceptor_loss, pos_donor_gain, pos_donor_loss]) for k in DicJunc[gene_symbol][x]["junc"]]
                #ai_res      = [[acceptor_gain, acceptor_loss, donor_gain, donor_loss][ri] if r else 0 for r,ri in zip(overlap_res,range(len(overlap_res)))]
                ai_res      = [ai_scores[k] if k in list(ai_scores.keys()) else 0 for k in DicJunc[gene_symbol][x]["junc"]]
                st = x.split(":")[1]
                #if ("splice_donor_variant" in info)==False and ("splice_acceptor_variant" in info)==False:   
                #exon skipping due to influence on 3 prime end of middle exon
                if (st=="+" and ai_res[2]<=-0.2) or (st=="-" and ai_res[1]<=-0.2):
                    hit+=1  
                    cons="|exon-skipping"
                    pred_new=pred+cons
                    info_new=info_main+pred_new
                    rel_dis,ref_seq,alt_seq=getSeqIfInfluenced(a,x,DicFA,info,"mEXdrop3")
                    res["overlap4Basset"].append("\t".join([ch,pos,clv,ref,alt,x,rel_dis,info_new,ref_seq,alt_seq])+"\n")
                    res["overlap4vcf"].append("\t".join([ch,pos,clv,ref,alt,x,rel_dis,info_new,ref_seq,alt_seq])+"\n")
                    #with_infl,rel_dis,ref_seq,alt_seq=getSeqIfInfluenced(a,x,DicFA,info,"mEXdrop3")
                    #if with_infl==True and ("frameshift_variant" in info)==False:
                    #    res["overlap4Basset"].append("\t".join([ch,pos,clv,ref,alt,x,rel_dis,info_new,ref_seq,alt_seq])+"\n")
                
                #exon skipping due to influence on 5 prime end of middle exon
                elif (st=="+" and ai_res[1]<=-0.2) or (st=="-" and ai_res[2]<=-0.2):
                    hit+=1  
                    cons="|exon-skipping"
                    pred_new=pred+cons
                    info_new=info_main+pred_new
                    rel_dis,ref_seq,alt_seq=getSeqIfInfluenced(a,x,DicFA,info,"mEXdrop5")
                    res["overlap4Basset"].append("\t".join([ch,pos,clv,ref,alt,x,rel_dis,info_new,ref_seq,alt_seq])+"\n")
                    res["overlap4vcf"].append("\t".join([ch,pos,clv,ref,alt,x,rel_dis,info_new,ref_seq,alt_seq])+"\n")
                    #with_infl,rel_dis,ref_seq,alt_seq=getSeqIfInfluenced(a,x,DicFA,info,"mEXdrop5")
                    #if with_infl==True and ("frameshift_variant" in info)==False:
                    #    res["overlap4Basset"].append("\t".join([ch,pos,clv,ref,alt,x,rel_dis,info_new,ref_seq,alt_seq])+"\n")                     
                
                #more exon inclusion
                elif min(ai_res[1],ai_res[2])>=0 and max(ai_res[1],ai_res[2])>=0.2:
                    hit+=1  
                    cons="|exon-inclusion"
                    pred_new=pred+cons
                    info_new=info_main+pred_new
                    rel_dis,ref_seq,alt_seq=getSeqIfInfluenced(a,x,DicFA,info,"mEXmore")
                    res["overlap4Basset"].append("\t".join([ch,pos,clv,ref,alt,x,rel_dis,info_new,ref_seq,alt_seq])+"\n")
                    res["overlap4vcf"].append("\t".join([ch,pos,clv,ref,alt,x,rel_dis,info_new,ref_seq,alt_seq])+"\n")
                    #with_infl,rel_dis,ref_seq,alt_seq=getSeqIfInfluenced(a,x,DicFA,info,"mEXmore")
                    #if with_infl==True and ("frameshift_variant" in info)==False:
                    #    res["overlap4Basset"].append("\t".join([ch,pos,clv,ref,alt,x,rel_dis,info_new,ref_seq,alt_seq])+"\n")
                      
                elif abs(ai_res[0])>=0.2 or abs(ai_res[3])>=0.2:
                    hit+=1
                    cons="|other-exon-hit"
                    pred_new=pred+cons
                    info_new=info_main+pred_new
                    rel_dis,ref_seq,alt_seq=getSeqIfInfluenced(a,x,DicFA,info,"otherEXchange")
                    res["overlap4vcf"].append("\t".join([ch,pos,clv,ref,alt,x,rel_dis,info_new,ref_seq,alt_seq])+"\n")
   
                #else:
                #    hit+=1
                #    cons="totally-ruined"
                        
    #if hit > 0:
    #    res["overlap4vcf"].append(a)

    return res
        

def curate(a,b):
    #Find donor lost no less than 0.2 in spliceAI prediction
    with open(a,"r") as f:
        ai=f.readlines()
    ai02=[]
    for x in ai:
        hit_ai02=0
        if (x[0] != "#") and ("SpliceAI=" in x):
            preds = x.strip().split("SpliceAI=")[-1].split(",")
            for pred in preds:
                allele_alt,gene_symbol,acceptor_gain,acceptor_loss,donor_gain,donor_loss,offset_acceptor_gain,offset_acceptor_loss,offset_donor_gain,offset_donor_loss = pred.split("|")
                if float(donor_loss) >= 0.2 or float(acceptor_loss) >= 0.2 or float(donor_gain) >= 0.2 or float(acceptor_gain) >= 0.2:
                    hit_ai02=1
        if hit_ai02>0:
            ai02.append(x)
    f1=open(b+a.split("/")[-1].split(".vcf")[0]+"_ds02.vcf","w")
    f1.writelines(ai02)
    f1.close()

    #Further find pathogenic mutations annotated by ClinVAR
    ai02_pathogenic=[]
    for x in ai02:
        if ("Pathogenic" in x) or ("Likely_pathogenic" in x):
            ai02_pathogenic.append(x)
    f2=open(b+a.split("/")[-1].split(".vcf")[0]+"_ds02_pathogenic.vcf","w")
    f2.writelines(ai02_pathogenic)
    f2.close()
    
    #Further find mutations within our 400bp training region along exon triplets
    #Also convert the final VCF to a format that can be used as the input for basset_sad_dadi.py
    sys.stdout.write("Building junction dictionary...")
    juncts=JuncDic("../Input_data/allSpliceIDs.txt")
    sys.stdout.write("done.\n")
    genome = pysam.Fastafile("../Input_data/GRCh37.fa")
    chrs=map(str,list(range(1,23)))+["X","Y"]
    sys.stdout.write("Reading chromosome sequence...")
    chrseqs={chr:genome.fetch(chr) for chr in chrs}
    sys.stdout.write("done.\n")
    
    #f3=open(b+a.split("/")[-1].split(".vcf")[0]+"_ds02_pathogenic_overlapTriplets.vcf","w")
    #f4=open(b+a.split("/")[-1].split(".vcf")[0]+"_ds02_pathogenic_overlapTriplets_4BassetSad.vcf","w")
    f3=open(b+a.split("/")[-1].split(".vcf")[0]+"_ds02_pathogenic_junctions.vcf","w")
    f4=open(b+a.split("/")[-1].split(".vcf")[0]+"_ds02_pathogenic_junctions_mEX_4BassetSad.vcf","w")
    for x in ai02_pathogenic:
        resDic=overlap(x,juncts,chrseqs)
        f3.writelines(resDic["overlap4vcf"])
        f4.writelines(resDic["overlap4Basset"])
    f3.close()
    f4.close()

    return 1


curate("../Output_data/clinvar_20190325_SNV_INDEL1000_mainChr_headerLength_spliceAIpredicted.vcf", "../Output_data/")