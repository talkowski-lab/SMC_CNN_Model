#!/usr/bin/env python

from optparse import OptionParser
import sys, copy, os, pdb, random, shutil, subprocess, time
import h5py
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import seaborn as sns
from sklearn import preprocessing
import dna_io
from itertools import product, combinations
import scipy.stats as stats

os.environ['BASSETDIR']=""
################################################################################
# basset_motifs.py
#
# Collect statistics and make plots to explore the first convolution layer
# of the given model using the given sequences.
################################################################################

#sns.set()

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <model_file> <input_file>'
    parser = OptionParser(usage)
    parser.add_option('-a', dest='act_t', default=0.0, type='float', help='Activation threshold (as proportion of max) to consider for PWM [Default: %default]')
    parser.add_option('-c', dest='conv_hdf5_file', default=None, help='Pre-computed model output as HDF5.')
    parser.add_option('-o', dest='out_dir', default='.')
    parser.add_option('-p', dest='plot_logo', default=False, action='store_true', help='Plot LOGO (in probability) for weights of each filter [Default: %default]')
    parser.add_option('-t', dest='test_set', default=False, action='store_true', help='Generate positional convolution heatmap for each sequence in the test set [Default: %default]')
    (options,args) = parser.parse_args()

    if (options.conv_hdf5_file is None) and (len(args) != 2):
        parser.error('Must provide Basset model file and test data in HDF5 format.')
    elif options.conv_hdf5_file is None:
        model_file = args[0]
        input_file = args[1]
        
        if not os.path.isdir(options.out_dir):
            os.mkdir(options.out_dir)

        #################################################################
        # load HDF5 file
        #################################################################        
        try:
            # input_file is HDF5
            input_hdf5 = input_file
            test_hdf5_in = h5py.File(input_hdf5, 'r')
            seqs_1hot = np.array(test_hdf5_in['test_in'])
            seq_class = np.array(test_hdf5_in['test_out'])
            test_hdf5_in.close()
    
        except:
            # input_file is not HDF5
            print "Input file has to be a HDF5 file."
            sys.exit()

        #################################################################
        # Torch Get Model Layers
        #################################################################
        options.conv_hdf5_file = '%s/convolution_out.h5' % options.out_dir
        #torch_cmd = '%s/src/basset_motifs_predict.lua %s %s %s' % (os.environ['BASSETDIR'],model_file, input_hdf5, options.conv_hdf5_file)
        torch_cmd = '%s/src/basset_motifs_positional_infl.lua -batch_size %d %s %s %s' % (os.environ['BASSETDIR'], seq_class.shape[0], model_file, input_hdf5, options.conv_hdf5_file)
        print torch_cmd
        subprocess.call(torch_cmd, shell=True)

    # load convolution output
    conv_hdf5_in = h5py.File(options.conv_hdf5_file, 'r')
    filter_weights = np.array(conv_hdf5_in['filter_weights'])    #shape (50, 4, 5)
    relu_outs = np.array(conv_hdf5_in['relu_outs'])              #shape (131, 50, 400) this is ReLU outputs 
    norm_beta = np.array(conv_hdf5_in['norm_beta'])       #shape (50, )
    norm_gamma = np.array(conv_hdf5_in['norm_gamma'])     #shape (50, )
    norm_mean = np.array(conv_hdf5_in['norm_mean'])       #shape (50, )
    norm_var = np.array(conv_hdf5_in['norm_var'])         #shape (50, )
    pos_infl_entropy = np.array(conv_hdf5_in['filter_infl_entropy'])   #shape (50, 400)
    pos_infl_targets = np.array(conv_hdf5_in['filter_infl_targets'])   #shape (50, 3, 400)
    original_entropy = np.array(conv_hdf5_in['original_entropy'])[0]   #single value
    conv_hdf5_in.close()

    # store useful variables
    num_filters = filter_weights.shape[0]
    filter_size = filter_weights.shape[2]

    #################################################################
    # individual filter plots
    #################################################################
    # also save information contents
    filters_ic = []
    #meme_out = meme_intro('%s/filters_meme.txt'%options.out_dir, seqs)
    
    # generate all possible kmers
    kmers=np.array(list(product([[0,0,0,1], [0,0,1,0], [0,1,0,0],[1,0,0,0]], repeat=filter_size)))
    kmers=np.swapaxes(kmers,1,2) #shape (50,4,5)
    
    #For the paper
    #plot_pos_infl_entropy_heatmap(pos_infl_entropy, filter_size, [1,9,10,18,21,22,25,27,29,37,47,49], "RdBu_r", '%s/positional_infl_entropy_heatmap.pdf' % options.out_dir)
    
    #For the patent
    motif_index1 = [18, 25, 49, 9, 10, 21, 29, 47, 1, 22, 27, 37]
    plot_pos_infl_entropy_heatmap(pos_infl_entropy, filter_size, motif_index1, "RdBu_r", '%s/positional_infl_entropy_heatmap_inclusion_all.pdf' % options.out_dir)
    #motif_index2 = [41,5,20,2,3,40,43,42,44,38,47,29,12,31,21,9,10]
    #plot_pos_infl_entropy_heatmap(pos_infl_entropy, filter_size, motif_index2, "RdBu_r", '%s/positional_infl_entropy_heatmap_exclusion_all.pdf' % options.out_dir)
    #motif_index3 = [37,16,1,22,48,17,27,33,6,46]
    #plot_pos_infl_entropy_heatmap(pos_infl_entropy, filter_size, motif_index3, "RdBu_r", '%s/positional_infl_entropy_heatmap_unchanged_all.pdf' % options.out_dir)

    #plot_pos_pval_heatmap(relu_outs, filter_size, seq_class, [18,25,49], ["Inclusion","Stable"], "Reds", '%s/positional_IvU_heatmap.pdf' % options.out_dir)
    #plot_pos_pval_heatmap(relu_outs, filter_size, seq_class, [18,25,49], ["Exclusion","Stable"], "Blues", '%s/positional_EvU_heatmap.pdf' % options.out_dir)

    #for f in range(num_filters) :
    #    print 'Filter %d' % f
        
        # plot LOGO (in probability) for weights of this filter
        #plot_weight_logo(filter_weights[f,:,:], kmers, norm_beta[f], norm_gamma[f], norm_mean[f], norm_var[f], '%s/filter%d'%(options.out_dir,f), maxpct_t=options.act_t)
        
        # plot positional heatmap
        #set seaborn global params
        
    #    plot_pos_heatmap_single(relu_outs[:,f,:], filter_size, seq_class, '%s/filter%d_positional_heatmap.pdf' % (options.out_dir,f))

        # write possum motif file
        #filter_possum(filter_weights[f,:,:], 'filter%d'%f, '%s/filter%d_possum.txt'%(options.out_dir,f), options.trim_filters)

        

        # make a PWM for the filter
        #filter_pwm, nsites = make_filter_pwm('%s/filter%d_logo.fa'%(options.out_dir,f))

        #if nsites < 10:
            # no information
        #    filters_ic.append(0)
        #else:
            # compute and save information content
        #    filters_ic.append(info_content(filter_pwm))

            # add to the meme motif file
        #    meme_add(meme_out, f, filter_pwm, nsites, options.trim_filters)

    #meme_out.close()


################################################################################
# plot_pos_pval_heatmap
#
# Plot positional heatmap for a single filter based on ReLU output, grouped by treatment response
#
# Input
#  ReLU_out_matrix: np.array of the ReLU output
#  out_pdf:
################################################################################
def plot_pos_infl_entropy_heatmap(mat, filter_size, filter_index, color, out_prefix, seg=4):
    #figure out gap size
    gap = filter_size - 1
    
    #figure out padding size
    pad_side = int(gap // 2)
    
    #figure out interval between two sub-regions
    subinterval = int(mat.shape[-1] / seg)
    
    #figure out sub-region length
    sublen = subinterval - gap
    
    #figure out range of each sub-region
    subs=[]
    for iseg in range(seg):
        subs.append(range(pad_side+subinterval*iseg,pad_side+subinterval*iseg+sublen))
    
    
    p_tensor = mat[filter_index,:]
        
    #plot heatmap
    yticks=["Motif%s" % x for x in filter_index]
    fig = plt.figure(figsize=(8.27, 11.69/3))
    gs =  gridspec.GridSpec(3, seg)
    gs.update(left=0.02, right=0.98, top=0.98, bottom=0.02, wspace=0.02, hspace=0.02)
    for sub, g in zip(subs,range(seg)):
        if g==0 or g==2:
            xticks=list(range(-24,-24+sublen))
        else:
            xticks=list(range(-75,-75+sublen))
        if g>0:
            yticks=False
        ax1=plt.subplot(gs[0,g])
        sns.heatmap(p_tensor[:,sub], linecolor="black", cmap=color, linewidths=0.1, xticklabels=xticks, yticklabels=yticks, cbar=True, cbar_kws={"orientation": "horizontal"}, ax=ax1,vmin=-0.01, vmax=0.01)
        
    plt.savefig(out_prefix,dpi=300)    
    plt.close(fig)





################################################################################
# plot_pos_pval_heatmap
#
# Plot positional heatmap for a single filter based on ReLU output, grouped by treatment response
#
# Input
#  ReLU_out_matrix: np.array of the ReLU output
#  out_pdf:
################################################################################
def plot_pos_pval_heatmap(relu_outs_3D, filter_size, class_vec, filter_index, comparison, color, out_prefix, seg=4):
    #check if the comparison is and only is between two levels
    if len(comparison)!=2:
        print "'comparison' argument has to be a list of length two, where the second one defines the baseline level."
        sys.exit()
    
    #figure out gap size
    gap = filter_size - 1
    
    #figure out padding size
    pad_side = int(gap // 2)
    
    #figure out interval between two sub-regions
    subinterval = int(relu_outs_3D.shape[-1] / seg)
    
    #figure out sub-region length
    sublen = subinterval - gap
    
    #figure out range of each sub-region
    subs=[]
    for iseg in range(seg):
        subs.append(range(pad_side+subinterval*iseg,pad_side+subinterval*iseg+sublen))
    
    #define the total number of test, which is required to correct p values
    cate=list(set(class_vec))
    total_cmps = list(combinations(range(len(cate)),2))
    total_test = len(total_cmps)
    
    #figure out which sequences are in the corresponding groups
    seq_index1=[i for i in range(len(class_vec)) if class_vec[i] == comparison[0]]
    seq_index2=[i for i in range(len(class_vec)) if class_vec[i] == comparison[1]]

    #ReLU 3D tensor for sequences in two classes
    relu_outs1=relu_outs_3D[seq_index1,:,:][:,filter_index,:]
    relu_outs2=relu_outs_3D[seq_index2,:,:][:,filter_index,:]
    
    #Initiate p value tensor
    p_tensor=np.zeros((len(filter_index),relu_outs_3D.shape[-1]),dtype=np.float32)
    
    #1-way ANOVA for each position
    print comparison[0],"vs",comparison[1]
    for mi in range(len(filter_index)):        
        for ri in range(p_tensor.shape[-1]):
            t,p=stats.ttest_ind(relu_outs1[:,mi,ri],relu_outs2[:,mi,ri])
            q = p*total_test
            if np.isnan(p):
                q=1
            elif q >= 0.05:
                q=1
            p_tensor[mi,ri]= np.log10(q)*(-1)
        
    #plot heatmap
    yticks=["Motif%s" % x for x in filter_index]
    
    fig = plt.figure(figsize=(8.27, 11.69/12))
    gs =  gridspec.GridSpec(1, seg)
    gs.update(left=0.02, right=0.98, top=0.98, bottom=0.02, wspace=0.02, hspace=0.02)
    for sub, g in zip(subs,range(seg)):
        if g==0 or g==2:
            xticks=list(range(-24,-24+sublen))
        else:
            xticks=list(range(-75,-75+sublen))
        if g>1:
            yticks=False
        ax1=plt.subplot(gs[0:2,g])
        sns.heatmap(p_tensor[:,sub], linecolor="black", cmap=color, linewidths=0.1, xticklabels=xticks, yticklabels=yticks, cbar=True, cbar_kws={"orientation": "horizontal"}, ax=ax1)
        
    plt.savefig(out_prefix,dpi=300)    
    plt.close(fig)
    

################################################################################
# plot_weight_logo
#
# Plot a weblogo of the filter's weight of the model
#
# Input
#  param_matrix: np.array of the filter's parameter matrix
#  out_pdf:
################################################################################
def plot_weight_logo(filter_weights, seq_tensor, norm_beta, norm_gamma, norm_mean, norm_var, out_prefix, maxpct_t=0.0):
    #Get convolution
    kmer_conv=np.tensordot(filter_weights,seq_tensor,axes=((0,1),(1,2))) #result in a 1D array
    
    #normalization using trained parameters
    kmer_conv_norm = (kmer_conv - norm_mean) * norm_beta / np.sqrt(norm_var + 0.00001) + norm_gamma
    
    #ReLU
    kmer_conv_relu = [0 if val <0 else val for val in kmer_conv_norm]
    
    #calculate a ReLU cutoff to plot weblogo
    if maxpct_t == 0:
        relu_act = 0
    else:
        all_outs = np.ravel(kmer_conv_relu)
        all_outs_mean = all_outs.mean()
        all_outs_norm = all_outs - all_outs_mean
        relu_act = maxpct_t * all_outs_norm.max() + all_outs_mean
    
    #Find which kmer pass the cutoff
    kmer_ok_index=[i for i in range(len(kmer_conv_relu)) if kmer_conv_relu[i]>relu_act]
    if len(kmer_ok_index)>0:
        kmer_ok=np.array([seq_tensor[i,:,:] for i in kmer_ok_index]) #shape (50,4,5)
        kmer_ok_4D=np.expand_dims(kmer_ok, axis=2) #shape (50,4,1,5)
        kmer_seq = dna_io.vecs2dna(kmer_ok_4D)
        kmer_fasta=[">"+str(i)+"\n"+seq+"\n" for i,seq in zip(range(len(kmer_seq)),kmer_seq)]
        kmer_fasta_out=open("%s_activated_kmers.fa" % out_prefix,"w")
        kmer_fasta_out.writelines(kmer_fasta)
        kmer_fasta_out.close()
        if relu_act > 0:
            subprocess.call("weblogo -X NO -Y NO -F pdf --resolution 300 --errorbars NO --fineprint '' -C '#CB2026' A A -C '#34459C' C C -C '#FBB116' G G -C '#0C8040' T T < %s_activated_kmers.fa > %s_weight_LOGO.pdf"%(out_prefix, out_prefix), shell=True)
        else:
            subprocess.call("weblogo -U probability -F pdf --resolution 300 --errorbars NO --fineprint '' -C '#CB2026' A A -C '#34459C' C C -C '#FBB116' G G -C '#0C8040' T T < %s_activated_kmers.fa > %s_weight_LOGO.pdf"%(out_prefix, out_prefix), shell=True)

     
################################################################################
# plot_pos_heatmap_single
#
# Plot positional heatmap for a single filter based on ReLU output, grouped by treatment response
#
# Input
#  ReLU_out_matrix: np.array of the ReLU output
#  out_pdf:
################################################################################
def plot_pos_heatmap_single(relu_outs, filter_size, class_vec, out_prefix, byclass=True, seg=4):
    #figure out gap size
    gap = filter_size - 1
    
    #figure out padding size
    pad_side = gap // 2
    
    #figure out interval between two sub-regions
    subinterval = relu_outs.shape[-1] / seg
    
    #figure out sub-region length
    sublen = subinterval - gap
    
    #figure out range of each sub-region
    subs=[]
    for iseg in range(seg):
        subs.append(range(pad_side+subinterval*iseg,pad_side+subinterval*iseg+sublen))
    
    #define class
    cate=list(set(class_vec))
    
    #Initiate plot
    fig = plt.figure(figsize=(8.27, 11.69))
    gs =  gridspec.GridSpec(6, seg)
    gs.update(left=0.02, right=0.98, top=0.58, bottom=0.02, wspace=0.1, hspace=0.2)
    
    
    conv_vec_byclass=[]
    #group relu vectors into classes
    for aclass, i in zip(cate,range(len(cate))):
        conv_vec=[]
        for cand, j in zip(class_vec,range(relu_outs.shape[0])):
            if cand==aclass:
                conv_vec.append(relu_outs[j,:])
        conv_vec=np.array(conv_vec)
            
        #copy into a final list of each class
        conv_vec_byclass.append(conv_vec)
        #conv_vec_colmean=np.mean(conv_vec,axis=0)
        #conv_vec_colmean=np.expand_dims(conv_vec_colmean,0)
        
    #1-way ANOVA
    cmps = list(combinations(range(len(cate)),2))
    #print cmps
    p_tensor=np.zeros((len(cmps),relu_outs.shape[1]),dtype=np.float32)
    for ci in range(len(cmps)):
        print cate[cmps[ci][0]],"vs",cate[cmps[ci][1]]
        p_vec=np.zeros((1,relu_outs.shape[1]),dtype=np.float32)
        for ri in range(p_tensor.shape[1]):
            t,p=stats.ttest_ind(conv_vec_byclass[cmps[ci][0]][:,ri],conv_vec_byclass[cmps[ci][1]][:,ri])
            q = p*len(cmps)
            if np.isnan(p):
                q=1
            elif q >= 0.05:
                q=1
            p_vec[0,ri]= np.log10(q)*(-1)
        p_tensor[ci,:]=p_vec
        
    #plot
    plot_cbar=False
    for sub, g in zip(subs,range(seg)):
        #lineplot
        ax1=plt.subplot(gs[0:2,g])
        for data in conv_vec_byclass:
            if g==0 or g==2:
                junc=list(range(-24,-24+sublen))*data.shape[0]
            else:
                junc=list(range(-75,-75+sublen))*data.shape[0]
            sns.lineplot(x=junc, y=np.ravel(data[:,sub]), ax=ax1)
            if g>0:
                ax1.tick_params(axis="y",labelleft=False)
        sns.despine()
            
        #heatmap
        ax2=plt.subplot(gs[3,g])
        if g==3:
            plot_cbar=False
        sns.heatmap(p_tensor[:,sub], linecolor="white", linewidths=.2, cbar=plot_cbar, xticklabels=False, yticklabels=False, ax=ax2)
        
        
    plt.savefig(out_prefix,dpi=300)    
    plt.close(fig)
        



################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
