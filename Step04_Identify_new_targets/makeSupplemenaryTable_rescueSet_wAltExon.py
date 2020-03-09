#!/usr/bin/env python

def make(a, b, b1):
    #with open("/data/talkowski/dg520/ref/importantGenes/constraint_2018_oe_ub_0.15.txt","r") as f:
    #    ref1=f.readlines()
    #ref0=[x.strip() for x in ref1]
    #res=[]
    with open(a,"r") as f:
        d=f.readlines()
    f1=open(b1,"w")
    f1.write("ClinVAR_ID\tAllele_ID\tReference\tAlteration\tgnomAD_Freq\tDisease\tClinVAR_Molecular_Consquence\tSpliceAI_Prediction\tEnSemblID\tGeneSymbol\tChromosome\tStrand\tIntron1\tExon\tIntron2\tRescue\tPubMed_ID\tInternal_ID\n")
    for x in d[1:]:
        cid,sid,score,ref,alt,target,ref_pred,alt_pred,sad,patho,mc,aid,dis,dis_db,pid,num,freq=x.strip().split("\t")
        eid,gs,info=sid.split("@")
        ch,st,coord=info.split(":")
        in1_start=coord.split("^")[0]
        in1_end=coord.split("^")[1].split("&")[0]
        in2_start=coord.split("&")[1].split("^")[0]
        in2_end=coord.split("&")[1].split("^")[1].split("|")[0]
        ex_start=str(int(in1_end)+1)
        ex_end=str(int(in2_start)-1)
        if target=="t0":
            cls="Rescue_by_inclusion_response"
            ai="splicing_loss"
        elif target=="t1":
            cls="Rescue_by_exclusion_response"
            if "nonsense" in score:
                ai="splicing_loss"
            else:
                ai="splicing_gain"
        f1.write("\t".join([cid,aid,ref,alt,freq,dis,mc,ai,eid,gs,ch,st,in1_start+"-"+in1_end,ex_start+"-"+ex_end,in2_start+"-"+in2_end,cls,pid,sid])+"\n")
        #if gs in ref0:
            #mc=mc.split("|")[1].split("_variant")[0]
        #    content="\t".join([cid,gs,ref,alt,dis,mc,ai,cls])+"\n"
        #    if (content in res) ==False:
        #        res.append(content)        
    f1.close()
    
    #f2=open(b2,"w")
    #f2.write("ClinVAR_ID\tGeneSymbol\tReference\tAlteration\tDisease\tClinVAR_Molecular_Consquence\tSpliceAI_Prediction\tPubMed_ID\n")
    #f2.writelines(res)
    #f2.close()

    with open(b, "r") as f:
        d2 = f.readlines()
    
    with open(b1, "r") as f:
        d1 = f.readlines()

    f2=open(b1[:-4]+"_wAltExon.txt","w")
    f2.write(d1[0])
    for x in d1[1:]:
        t = x.split("\t")
        if ("Rescue_by_exclusion_response" in t):
            if (t[-1] in d2):
                f2.write(x)
        else:
            f2.write(x)
    f2.close()

    return 1


make("../Output_data/sad_table_rescue_alleleFreq.txt",
     "../Input_data/alternative_exon_tripletIDs.txt",
     "../Output_data/sad_table_rescue_alleleFreq_4supplementaryTable_wAltExon.txt")