#!/usr/bin/env python3

'''
Execute TIP on a Wiggle file format.

TIP is a probabilistic method that quantitatively measures the regulatory potential of a TF to genes, 
based on the proof that the posterior probability of a TF t targeting gene g is proportional to 
the weighted sum of the binding signal of t over all the positions on gene g:  
    
    * wig2weight, by analyzing mapped reads, retrieves a TF’s binding profile (i.e. its binding signals);
    * calscoreonwig, given the TF's binding profile, computes the regulatory score of the TF to every gene g.

Usage:
    tip_wig.py <path_to_experiments> <bamfilename> 
    
    User defined variable:
    - bamfilename: the name of the bam file as it is in ENCODE database
    - genome: the reference genome (GRCh38 by default)
    - annofilename: the name of annotation file GRCh38
    - width: the window size upstream (10000 by default)
    - smooth: logic value, if smooth=True perform smoothing for the TF binding profile (True by default)
    - typeStep: type of wiggle format (variableStep by default) 

Reference: https://doi.org/10.1093/bioinformatics/btr552
'''

from tqdm import tqdm
from statsmodels.stats.multitest import multipletests
import os, argparse, numpy as np, pandas as pd, scipy.stats

def wig2weight(smooth=True, typeStep='variableStep'):

    '''
    The function takes in input mapped reads and retrieves a TF’s binding profile. 
    The mapped reads are stored in a wiggle format with 'variableStep' and must be splitted chr by chr in order to save 
    substantial space in memory. By consulting the reference genome's annotation file, previously filtered to report 
    only wild type protein coding human genes, the function tracks the reads' frequencies that are centred in gene's TSS, 
    for every gene. Both the gene's TF binding signals and the TF’s binding profile which derives from them are finally 
    archived and saved in numpy binary arrays to speed the following I/O operations up.

    Every position i in the TF’s binding profile (vector weights) is built as the ratio between:
        - the sum of all the gene's TF binding signals (gene vectors) at position i   
        - and the sum of all the gene's TF binding signals at every position [-n,+n].
    '''

    print('wig2weight function:\nwidth:\t%s\nsmooth:\t%s\ngenome:\t%s' %(args.width, args.smooth, args.genome))
    
    ## For this script
    chr_len = chr_info['length'][0:-1]
    chr_nam = chr_info['name'][0:-1]

    # read in gene info: the annotationfilename has been prepared in R
    mygene = pd.read_csv(args.annofilename, sep='\t', names=['name', 'chr', 'str', 'sta', 'end'])
   
    if typeStep == 'variableStep':
        for k in tqdm(range(len(chr_nam)), desc='Build_TFprof'):

            filein = '%szig/%s.zig' %(wd, chr_nam[k])
            myw = np.zeros(args.width * 2, dtype=np.int8) # initialize all-zero object myw for TF t, 20001 pos long

            if os.path.exists(filein):

                # read in wiggle file
                data = open(filein, 'r')
                header = data.readline()
                index = chr_nam.index(chr_nam[k])
                read_cov = np.zeros(chr_len[index], dtype=np.int8) # read.cov is an all-zero vector, long as the chr under study

                for line in data:
                    line = line.rstrip().split()
                    start = int(line[0].split('-')[0])
                    stop = int(line[0].split('-')[1])
                    val = int(line[1])
                    for z in range(start,stop + 1):
                        read_cov[z] = val
                
                np.save('%s%s' %(wd, chr_nam[k]), read_cov)
            else: next

            curgene = mygene[mygene['chr'] == chr_nam[k]].reset_index(drop=True)
            for i in range(len(curgene)):
                try:
                    myw = myw + read_cov[curgene.loc[i,'sta'] : curgene.loc[i,'end']]
                except: 
                    print('Error with:\t%s' %(curgene.loc[i,'name']))
                    next

        ## working with smaller numbers?
        myw = myw / len(mygene)
        
        if smooth == True:
            tmp = np.array([i for i in myw])
            for i in range(len(myw)):
                myw[i] = np.mean(tmp[max(i-250,1):min(i+250,len(myw))])
        
        myw = myw / np.sum(myw)
        np.save('%s%s.npy' %(wd, args.myoutFile1), myw)


def calscoreonwig(width=10000, typeStep='variableStep', genome='GRCh38'):

    '''
    Given the TF's binding profile, calculates the regulatory score of the TF to every gene g.
    The function takes in input numpy binary arrays of gene's TF binding signals and the TF’s binding profile and retrieves 
    the regulatory score of the TF for every gene. The distribution of raw scores' frequencies is assumed to be normal, 
    the regulatory scores are transformed into z-scores and the significance for each gene is estimated based on its z-score. 
    The results are given in a tsv table.

    Every gene's regulatory score is computed as a dot product between:
        - the TF's binding profile  
        - and gene's TF binding signals.
    '''

    print('\ncalscoreonwig function:\nwidth:\t%s\ngenome:\t%s' %(args.width, args.genome))

    ## For this script
    chr_len = chr_info['length'][0:-1]
    chr_nam = chr_info['name'][0:-1]

    mygene = pd.read_csv(args.annofilename, sep='\t', names=['name', 'chr', 'str', 'sta', 'end'])

    if typeStep == 'variableStep':
        mysco = np.zeros(len(mygene))
        gname = ['' for i in range(len(mygene))]
        count = 0

        for k in tqdm(range(len(chr_nam)), desc='Comp_Gene-scores'):

            if os.path.exists('%s%s.npy' %(wd, chr_nam[k])):
                chr_signals = np.load('%s%s.npy' %(wd, chr_nam[k]))
                index = chr_nam.index(chr_nam[k])
                curgene = mygene[mygene['chr'] == chr_nam[k]].reset_index(drop=True)

                for i in range(len(curgene)):
                    tmp = chr_signals[max(1, curgene.loc[i,'sta']) : min(curgene.loc[i,'end'], chr_len[index])]
                    b1 = 1 - (curgene.loc[i,'sta'])
                    b2 = (curgene.loc[i,'sta']) - chr_len[index]
                    if b1>0:
                        add = [0 for z in range(b1)]
                        tmp = np.concatenate((add, tmp))
                    elif b2>0:
                        add = [0 for z in range(b2)]
                        tmp = np.concatenate((tmp, add))

                    ## the regulatory score is a dot product between tf_weights and signals
                    mysco[count] = np.dot(tmp,tf_weights)
                    gname[count] = curgene.loc[i,'name']
                    count += 1

            else: next

    ### The distribution of raw.scores is assumed to be Gaussian-like, 
    ### therefore one could extract zscore, pvalues and apply Benjamini-Hochberg 
    ### multiple hypotheses test correction with FDR at 0.01
    zscore = (mysco - np.mean(mysco)) / np.std(mysco)
    pvals = scipy.stats.norm.pdf(zscore)
    reject, benjhoc_pvals = multipletests(pvals, alpha=args.alpha, method='fdr_bh', is_sorted=False, returnsorted=False)[0:2]
    
    print('%i hypotheses rejected for alpha=%f with fdr_bh\n' %(sum(reject),args.alpha))
    header = 'g.name\traw.score\tz.score\tp.value\tq.value'
    with open('%s%s.tsv' %(wd, args.myoutFile2), 'w') as scores:
        scores.write('%s\n' %header)
        for i in range(len(gname)):
            scores.write('%s\t%f\t%f\t%f\t%f\n' %(gname[i], mysco[i], zscore[i], pvals[i], benjhoc_pvals[i]))
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='TIP probabilistic method on file in wiggle format')
    parser.add_argument('path', help='provide the path to the directory where the experiment is stored')
    parser.add_argument('bamfilename', help='provide the name of the bam experiment (as it is called in ENCODE)')
    parser.add_argument('annofilename', help='provide the genome annotation file')

    parser.add_argument('-g', '--genome', dest='genome', default='GRCh38',
                        help='the reference genome (by default GRCh38, alternatly hg19)')
    parser.add_argument('-og', '--organism', dest='organism', default='human',
                        help='the organism (by default human, alternatly mouse)')
    parser.add_argument('-ow', '--outfile_weights', dest='myoutFile1', default='weights',
                        help='the name to the weights outfile (by default weights)')
    parser.add_argument('-os', '--outfile_scores', dest='myoutFile2', default='score_refseq',
                        help='the name to the scores outfile (by default score_refseq)')
    parser.add_argument('-w', '--width', dest='width', type=int, default=10000,
                        help='provide the width to look back and forth with respect to gene TSS (by default 10000)')
    parser.add_argument('-a', '--alpha', dest='alpha', type=float, default=0.01,
                        help='provide the FWER (alpha value) for Benjamini/Hochberg p-value correction (by default 0.01)')
    parser.add_argument('-s', '--smooth', dest='smooth', type=bool, default=True,
                        help='activate smooth mode on TF binding profile')   
    parser.add_argument('-r', '--remove', dest='remove', default=False, action='store_true',
                        help='remove useless .npy files once used')                                      
    args = parser.parse_args()

    if args.organism == 'human':
        chr_nam = ['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY', 'chrM']
        if args.genome == 'GRCh38': 
            chr_len = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 
                        145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 
                        101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 
                        156040895, 57227415, 16569]
        elif args.genome == 'hg19':
            chr_len = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 
                        146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 
                        102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 
                        155270560, 59373566, 16571]
    elif args.organism == 'mouse':
        chr_nam = ['chr' + str(i) for i in range(1,20)] + ['chrX', 'chrY', 'chrM']
        if args.genome == 'mm10':
            chr_len = [195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 
                        129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 
                        104043685, 98207768, 94987271, 90702639, 61431566, 171031299, 91744698, 16299]

    chr_info = {'length': chr_len, 'name': chr_nam}

    wd = '%s%s/' %(args.path, args.bamfilename)
    weights = '%s%s.npy' %(wd, args.myoutFile1)
    scores = '%s%s.tsv' %(wd, args.myoutFile2)

    if not os.path.exists(weights):
        wig2weight()
    
    if os.path.exists(weights) and not os.path.exists(scores):
        tf_weights = np.load(weights)
        calscoreonwig()
        if args.remove: 
            for chr in chr_nam:
                if os.path.exists('%s%s.npy' %(wd, chr)): os.remove('%s%s.npy' %(wd, chr))