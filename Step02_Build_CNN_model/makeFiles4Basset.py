#!/usr/bin/env python

from __future__ import division
import numpy.random as npr 
import pysam
from Bio.Seq import reverse_complement


BASSET_FOLDER = ""

def makeBed(a,cl,inside=25,outside=75,reg=2):
    d=[]
    genome = pysam.Fastafile("../Input_data/GRCh37.fa")    
    with open(a,"r") as f:
        d0=f.readlines()
    for x in d0:
        t=x.strip()
        ch=t.split("@")[2].split(":")[0]
        st=t.split("@")[2].split(":")[1]
        end0  =int(t.split("@")[2].split(":")[2].split("^")[0])-1
        start =int(t.split("@")[2].split("^")[1].split("&")[0])
        end   =int(t.split("@")[2].split("&")[1].split("^")[0])-1
        start1=int(t.split("@")[2].split("|")[0].split("^")[-1])
        s0   = end0-inside
        e0   = end0+outside
        s1   = start-outside
        e1   = start+inside
        s2   = end-inside
        e2   = end+outside
        s3   = start1-outside
        e3   = start1+inside

        bedname=t+"::"+cl
        seq0=genome.fetch(ch,s0,e0) if st =="+" else reverse_complement(genome.fetch(ch,s0,e0))
        seq1=genome.fetch(ch,s1,e1) if st =="+" else reverse_complement(genome.fetch(ch,s1,e1))
        seq2=genome.fetch(ch,s2,e2) if st =="+" else reverse_complement(genome.fetch(ch,s2,e2))
        seq3=genome.fetch(ch,s3,e3) if st =="+" else reverse_complement(genome.fetch(ch,s3,e3))
        if reg==4:
            if st=="+":
                d.append(bedname+"\t"+seq0+seq1+seq2+seq3)
            else:
                d.append(bedname+"\t"+seq3+seq2+seq1+seq0)
        elif reg==2:
            if st=="+":
                d.append(bedname+"\t"+seq1+seq2)
            else:
                d.append(bedname+"\t"+seq2+seq1)

    print "Sequences extracted."
    
    return d
    
    
def makeBassetPre(a,limit,pre):
    dic={}
    sz=[]
    cl=[]
    n=0
    k=0
    for d in a:
        sz.append(len(d))
        for i in d:
            t=i.strip().split("\t")
            spliceID,target=t[0].split("::")
            if (target in cl)==False:
                cl.append(target)
            header=spliceID+"::"+target
            seq=t[1]
            dic[n]={"name":header,"seq":seq,"class":["0"]*len(a)}
            dic[n]["class"][k]="1"
            n+=1
        k+=1
        print "Sub-dictionary built."        	
    print "All dictionary built:",len(dic),"entries."	

    train=[]
    valid=[]
    test =[]
    remain=[]
    ini=0
    for x in sz:
        if x <= limit:
            npr.seed(122)
            rnd=npr.choice(xrange(ini,ini+x),x,replace=False)
            valid+=list(rnd[-int(round(x/5)):])
            test+=list(rnd[:-int(round(x/5))][-int(round(x*0.1)):])
            train+=list(rnd[:-int(round(x/5))][:-int(round(x*0.1))])
            ini+=x
        else:
            x0=limit
            npr.seed(122)
            rnd0=npr.choice(xrange(ini,ini+x),x,replace=False)
            rnd=rnd0[:x0]
            valid+=list(rnd[-int(round(x0/5)):])
            test+=list(rnd[:-int(round(x0/5))][-int(round(x0*0.1)):])
            train+=list(rnd[:-int(round(x0/5))][:-int(round(x0*0.1))])
            remain+=list(rnd0[x0:])
            ini+=x
            
    npr.seed(122)
    train=npr.permutation(train)
    npr.seed(122)
    valid=npr.permutation(valid)
    npr.seed(122)
    test=npr.permutation(test)
    npr.seed(122)
    remain=npr.permutation(remain)
    all=list(train)+list(valid)+list(test)
    lft=list(remain)
    print "Data separated:"
    print "Training:",len(train)
    print "Validation:",len(valid)
    print "Test:",len(test)
    print "Left:",len(lft)

    #write act and fa file for Basset H5 generation
    print "Writing to output..."
    f1=open(pre+".fa","w")
    f2=open(pre+"_act.txt","w")
    f2.write("\t"+"\t".join(cl)+"\n")
    for x in all:
        f1.write(">"+dic[x]["name"]+"\n")
        f1.write(dic[x]["seq"]+"\n")
        f2.write(dic[x]["name"]+"\t"+"\t".join(dic[x]["class"])+"\n")
    f1.close()
    f2.close()

    f3=open(pre+"_dataSplit.txt","w")
    f3.write("Train:"+"\t"+str(len(train))+"\n")
    f3.write("Valid:"+"\t"+str(len(valid))+"\n")
    f3.write("Test:"+"\t"+str(len(test))+"\n")
    f3.close()
        
    f4=open(pre+"_left.4Test.fa","w")
    f5=open(pre+"_left.names.txt","w")
    for x in lft:
        f4.write(">"+dic[x]["name"]+"\n")
        f4.write(dic[x]["seq"]+"\n")
        f5.write(dic[x]["name"]+"\n")
    f4.close()
    f5.close()
    
    f6=open(pre+"_test0.4Test.fa","w")
    f7=open(pre+"_test0.names.txt","w")
    for x in test:
        f6.write(">"+dic[x]["name"]+"\n")
        f7.write(dic[x]["name"]+"\n")
        f6.write(dic[x]["seq"]+"\n")
    f6.close()
    f7.close()

    print "Making H5 file for Basset..."
    cmd= BASSET_FOLDER + "/src/seq_hdf5.py -c -v " + str(len(valid)) + " -t " + str(len(test)) + " " + pre + ".fa " + pre + "_act.txt " + pre + "_learn.h5"
    sta,out=commands.getstatusoutput(cmd)
    print out
    #cmd="/data/talkowski/dg520/projects/Basset/src/seq_hdf5.py -c -t " + str(len(remain)) + " " + pre + ".fa " + pre + "_act.txt " + pre + "_testRemain.h5"
    #sta,out=commands.getstatusoutput(cmd)
    #print out
    
    return 1
        
def doit(a,maxi,piece,outp):
    p=[]
    for x in a:
        p.append(makeBed(x[0],x[1],inside=25,outside=75,reg=piece))
    makeBassetPre(p,maxi,outp)
    
    return 1
    
#Main

doit([("../Output_data/up01tripletNames_pos_fdr01.txt","Inclusion"),
      ("../Output_data/dn01tripletNames_pos_fdr01.txt","Exclusion"),
      ("../Output_data/nc01tripletNames_all_382.txt","Stable")],
      700,
      4,
      "../Output_data/3C")
