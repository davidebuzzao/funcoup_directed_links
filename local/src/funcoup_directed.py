#!/usr/bin/env python3
'''
Compute and assign LLR to TF-gene links

FunCoup is a redundancy weighted Bayesian integration framework that uses different types of evidence (HT-data) 
and HQ gold standard data to extract log likelihood ratio (LLR) scores and assign FBS to two potential functionally coupled 
proteins by the integration of all evidences. 

Currently, all the evidences are undirected. The FunCoup_directed project aims at building new networks with directional links.
Potential regulatory links (i.e. links from transcription factors (TFs) to their target genes) are extracted from ChIP-Seq data 
and here are assigned of LLRs. To account for redundancy within same evidence type, a weighting schema is implemented in the framework.
The LLR(a,b)t for specific evidence type t is then calculated as a weighted sum of the individual LLR(a,b)e for each evidence e ∈ E 
of type t evidences:
    * LLRs are ranked by their absolute value in decreasing order; 
    * each LLRe is weighted by the product of the distances (the distance includes the Spearman correlation coefficient between TIP 
      overlapping links' raw.scores for each couple of the E experiments for the same Experiment Target) to each evidence with lower rank.

Usage: funcoup_directed [-h] [-b] [-db DB_NAME] [-og ORGANISM] [-fdr FDR] 
                        [-tr] [-fl] [-tf TF] [-lk LINK LINK] [-llr LLR] 
                        [-w] [-nx] [-grn] [-ens] [-Nens]

optional arguments:
  -h, --help            show this help message and exit
  -b, --build           activate LLRs and ENCODE building mode (by default False)
  -db DB_NAME, --database DB_NAME
                        type of gold standard in use (by default trrust, alternately regnet)
  -og ORGANISM, --organism ORGANISM
                        provide the common name of the organism under study (by default human, alternately mouse)
  -fdr FDR              provide the critical value of FDR (by default 01, alternately 05)
  -tr, --threshold      activate threshold mode and save just links with LLR > X (by default False)
  -fl, --filter         activate filter mode (by default False)
  -tf TF, --TransFact TF
                        provide the name of a TF to get all its links (by default False)
  -lk LINK LINK, --link LINK LINK
                        provide the name of a TF and a gene to look for their link (by default False)
  -llr LLR, --loglikelihood LLR
                        provide a LLR threshold (by default 0)
  -w, --weight          activate FunCoup weighting schema mode (by default False)
  -nx, --netx           activate networkX mode to build a proper specific input (by default False)
  -grn, --CancerGRN     activate CancerGRN mode to build a proper specific input (by default False)
  -ens, --ensembl       activate ensembl mode to map gene symbol to Ensembl identifiers (by default False)
  -Nens, --NOTensembl   activate NOTensembl mode to retrieve unmapped gene symbol to Ensembl identifiers as #tf-link (by default False)

Examples of usage:
    * To retrieve LLRs when PGS is TRRUSTv2, find SP corr. matrix and assign LLRs to links:   
        funcoup_directed -b -fdr 01 -w 
    
    * To retrieve LLRs when PGS is RegNetwork, find SP corr. matrix and assign LLRs to links:   
        funcoup_directed -b -db regnet -fdr 05 -w
    
    * To retrieve all FC_directed links:
        funcoup_directed -fdr 05 -w -ens

    * To retrieve only those links with llr bigger than a specified threshold:
        funcoup_directed -fdr 05 -w -ens -tr -llr 2

Reference: http://www.genome.org/cgi/doi/10.1101/gr.087528.108. 
'''

import argparse, numpy as np, pandas as pd, os, pickle, scipy.stats
from tqdm import tqdm

mor_type = ['Activation',
            'Repression',
            'Unknown']
 
def normalize_scores(exp):
    scores_file = 'data/evidences/encode/%s/chip_seq/scores/%s/%s/score.%s.sort.tsv' %(args.organism, args.fdr, exp, args.fdr)
    raw_scores = []
    with open(scores_file) as scores: 
        header = scores.readline()
        for raw_score in scores:
            raw_scores.append(float(raw_score.rstrip().split()[1]))
    return(np.array(raw_scores)/np.max(raw_scores))

def overlap_scores(exp1,exp2):
    '''
    The function takes in input two ChIP-Seq experiments' File Accessions for the same Experiment Target,
    extract gene targets' names and raw.scores assigned by TIP model and returns two lists of ordered floats 
    corresponding to overlapping links' respective raw.scores.

    The overlap of TIP raw.scores is chosen in order to compute afterwards a correlation coefficient between two 
    same sized vectors of numbers. 
    '''

    scores_file1 = 'data/evidences/encode/%s/chip_seq/scores/%s/%s/score.%s.sort.tsv' %(args.organism, args.fdr, exp1, args.fdr)
    scores_file2 = 'data/evidences/encode/%s/chip_seq/scores/%s/%s/score.%s.sort.tsv' %(args.organism, args.fdr, exp2, args.fdr)

    with open(scores_file1) as scores1, open(scores_file2) as scores2: 
        exp1_scores,exp2_scores  = {},{}
        header1, header2 = scores1.readline(), scores2.readline()
        for raw_score1 in scores1: 
            raw_score1 = raw_score1.rstrip().split()
            exp1_scores[raw_score1[0]] = exp1_scores.get(raw_score1[0], float(raw_score1[1]))
        for raw_score2 in scores2: 
            raw_score2 = raw_score2.rstrip().split()
            exp2_scores[raw_score2[0]] = exp2_scores.get(raw_score2[0], float(raw_score2[1]))
    exp_overlap = set(list(exp1_scores.keys())) & set(list(exp2_scores.keys()))

    clean_exp1, clean_exp2 = [], []
    for exp in exp_overlap:
        clean_exp1.append(exp1_scores[exp]); clean_exp2.append(exp2_scores[exp])
    return(clean_exp1,clean_exp2)


def compute_correlation(list_input, normalize=False, padding=False):
    '''
    Calculate a Spearman correlation coefficient with associated p-value.

    The Spearman rank-order Correlation Coefficient (SCC) is a nonparametric measure of the monotonicity of the relationship 
    between two datasets. Unlike the Pearson correlation, the Spearman correlation does not assume that both datasets 
    are normally distributed. Like other correlation coefficients, this one varies between -1 and +1 with 0 implying no correlation. 
    Correlations of -1 or +1 imply an exact monotonic relationship. Positive correlations imply that as x increases, so does y. 
    Negative correlations imply that as x increases, y decreases.
    
    The SCC between TIP raw.scores is computed for each couple of N experiments for the same Experiment Target. 
    The function in use is scipy.stats.spearmanr(X,Y) where X,Y are two equal sized vectors.
        *   here the SCC is computed on the raw.scores TF-gene links overlap, therefore there's no need of padding vectors;
        *   if padding = True: the size of the two vectors is supposed not to coincide, then a padding procedure is implemented to make 
            the size of the shortest one be as long as the longest (by the use of MAX or MEAN value of the shortest vector).

    The matrix storing SCCs per TF-gene link looks like:

          exp1 | exp2 | exp3 | [..] | expN
    -----|-----|------|------|------|------
    exp1 |  1. | 0.5  | 0.6  |   /  |  -1
    -----|-----|------|------|------|------
    exp2 | 0.  |  1.  | 0.4  |   /  |  0.7
    -----|-----|------|------|------|------
    exp3 | 0.  |  0.  |  1.  |   /  |  0.7
    -----|-----|------|------|------|------
    [..] |  /  |   /  |   /  |   /  |  /
    -----|-----|------|------|------|------
    expN | 0.  |  0.  |  0.  |   /  |  1.

    '''
    
    exp = [i[0] for i in list_input]
    corr_matrix = np.zeros((len(exp),len(exp))) ## symmetrix matrix
    if len(exp) > 1:
        for i in range(len(corr_matrix)):
            exp1 = exp[i] ## file1 identifier
            if normalize: exp1 = normalize_scores(exp[i])
            for j in range(len(corr_matrix[i])):
                if corr_matrix[j][i] == 0:
                    exp2 = exp[j] ## file2 identifier
                    if normalize: exp2 = normalize_scores(exp[j])
                    if normalize and padding:
                        diff = len(exp1) - len(exp2)
                        if diff < 0: exp1 = np.concatenate((exp1, np.repeat(np.max(exp1),abs(diff))), axis=None)
                        elif diff > 0: exp2 = np.concatenate((exp2, np.repeat(np.max(exp2),abs(diff))), axis=None)
                        corr_matrix[i][j] = scipy.stats.spearmanr(exp1,exp2)[0]
                    else:
                        clean_exp1,clean_exp2 = overlap_scores(exp1,exp2)
                        corr_matrix[i][j] = scipy.stats.spearmanr(clean_exp1,clean_exp2)[0]

    else: corr_matrix[0][0] = 1.

    return(corr_matrix)

def compute_LLR(save=False):
    pos_gd = pd.read_pickle('data/stat/%s/%s/%s/pos_gd.fl-%s.pkl' %(args.organism, args.db_name, args.fdr, args.filter))
    neg_gd = pd.read_pickle('data/stat/%s/%s/neg_gd.pkl' %(args.organism, args.db_name))
    LLR = pd.DataFrame(np.log(pos_gd['rel_freq']/neg_gd['rel_freq'])).rename(columns={'rel_freq': 'LLR'})
    
    if save:
        LLR.to_pickle('data/stat/%s/%s/%s/LLR.fl-%s.pkl' %(args.organism, args.db_name, args.fdr, args.filter))
        return(LLR)

def build_encode(save=False):
    '''
    The function builds and returns a dictionary like: 
    {
        'experiments': {
            TF1: [exp1,exp2,expM],
            TF2: [exp1,exp2,expM],
            TFX: [exp1,exp2,expM]
                    },
        'links': {
            TF1: [gene1,gene2,geneN],
            TF2: [gene1,gene2,geneN],
            TFX: [gene1,gene2,geneN]
        }
    }
    '''
    encode = {'experiments':{}, 'links':{}}
    for exp in exp_list:
        metadata_file = 'data/evidences/encode/%s/chip_seq/scores/%s/%s/metadata.tsv' %(args.organism, args.fdr, exp)
        with open(metadata_file) as metadata: 
            tf = metadata.readlines()[-1].split()[7].split('-')[0]
            if args.organism == 'mouse': tf = tf[0] + tf[1:].lower() 
        encode['experiments'][tf] = encode['experiments'].get(tf,[])
        encode['experiments'][tf] += [exp]

    for tf in encode['experiments']:
        links = []
        for exp in encode['experiments'][tf]:
            scores_file = 'data/evidences/encode/%s/chip_seq/scores/%s/%s/score.%s.sort.tsv' %(args.organism, args.fdr, exp, args.fdr)
            with open(scores_file) as scores: 
                header = scores.readline()
                for link in scores:
                    links.append(link.rstrip().split()[0])
        encode['links'][tf] = encode.get(tf,list(set(links)))

    if save:
        with open('data/evidences/encode/%s/encode.pkl' %(args.organism), 'wb') as fout:
            pickle.dump(encode, fout)
    return(encode)

def sort_LLR(save=False):
    '''
    The function builds and returns a dictionary like: 
    {
        'TF1': {
            gene1: [(exp1,LLR1),(exp2,LLR2),(expN,LLRN)],
            gene2: [(exp1,LLR1),(exp2,LLR2),(expN,LLRN)],
            geneX: [(exp1,LLR1),(exp2,LLR2),(expN,LLRN)]
        }
        'TFX': {
            gene1: [(exp1,LLR1),(exp2,LLR2),(expN,LLRN)],
            gene2: [(exp1,LLR1),(exp2,LLR2),(expN,LLRN)],
            geneX: [(exp1,LLR1),(exp2,LLR2),(expN,LLRN)]
        }
    }
    where every TF-gene's list stores tuples with ChIP-Seq Experiment's File Accessions and 
    LLRs for the respective link in deacreasing LLR order. 
    '''
    sorted_LLR = {}
    for tf in ENCODE['links']:
        for gene in ENCODE['links'][tf]:
            for exp in ENCODE['experiments'][tf]:
                filein = open('data/evidences/encode/%s/chip_seq/scores/%s/%s/score.%s.sort.pkl' %(args.organism, args.fdr, exp, args.fdr), 'rb')
                score_bins = pickle.load(filein)
                sorted_LLR[tf] = sorted_LLR.get(tf,{})
                if gene in [x for v in score_bins.values() for x in v]: ## list all the values and check if there's gene 
                    chunk = list({k:v for k,v in score_bins.items() if gene in v}.keys())[0] ## take correspondend key for gene 
                    sorted_LLR[tf][gene] = sorted_LLR[tf].get(gene,[])
                    sorted_LLR[tf][gene].append((exp, LLR['LLR'][chunk]))
             
                try: sorted_LLR[tf][gene].sort(key=lambda tup: tup[1], reverse=True)
                except: pass
    
    if save:
        with open('data/stat/%s/%s/%s/sorted_LLR.fl-%s.pkl' %(args.organism, args.db_name, args.fdr, args.filter), 'wb') as fout:
            pickle.dump(sorted_LLR, fout)
    return(sorted_LLR)


def naïve_funcoup(save=False):
    '''
    The naïve FunCoup schema doesn't include a weighting procedure when integrates 
    LLRs for the same Experiment Target, (same type of evidence, different experiment). 
    '''
    fc = {}
    for tf in ENCODE['experiments']:
        for exp in ENCODE['experiments'][tf]:
            filein = open('data/evidences/encode/%s/chip_seq/scores/%s/%s/score.%s.sort.pkl' %(args.organism, args.fdr, exp, args.fdr), 'rb')
            score_bins = pickle.load(filein)
            fc[tf] = fc.get(tf,{})
            for ck in score_bins:
                for gene in score_bins[ck]:
                    fc[tf][gene] = fc[tf].get(gene,0.)
                    fc[tf][gene] += LLR['LLR'][ck]
                    try:
                        for mor in database[tf]:
                            if gene in database[tf][mor]: fc[tf][gene].append(mor)
                    except: pass
    if save:
        with open('data/stat/%s/%s/%s/naive_funcoup.fl-%s.pkl' %(args.organism, args.db_name, args.fdr, args.filter), 'wb') as fout:
            pickle.dump(fc, fout)
        return(fc)

def weighted_funcoup(save_LLR=False, save_corr=False, save_fc=False, alpha=0.7):
    '''
    FunCoup4.1 implements in its Bayesian framework a weighting procedure when integrates 
    LLRs for the same Experiment Target (same type of evidence, different experiment)
    to account for redundancy within the evidences dataset.
    
    Here, for every interaction TF-gene:
        * tup(Experiment, LLR) are sorted in decreasing order by calling the function sort_LLR()
        * the Spearman correlation coefficient is computed between every tup(Experiment, LLR) 
          by calling the function compute_correlation() on:
            - normalized TIP raw.scores by calling the function normalize_scores()
            - padded TIP raw.scores with the max score (when normalized, by adding 1s)
        * LLR(TF-gene) for the specific evidence type t (ChIP-Seq data in this case) is calculated 
          as the weighted sum of the individual LLR(TF-gene)e for each experiment e of type t as follows:
            LLR(a,b)_{t} = \sum_{e} LLR(a,b)_{e} \prod_{k<e} d_{ek}  with  d_{ek} = α(1-max(0,r_{ek}))
    '''
    ## sort LLR session
    if save_LLR: sorted_LLR = sort_LLR(save=True)
    else:
        fin = open('data/stat/%s/%s/%s/sorted_LLR.fl-%s.pkl' %(args.organism, args.db_name, args.fdr, args.filter), 'rb')
        sorted_LLR = pickle.load(fin)
    tfs = list(sorted_LLR.keys())

    ## Spearman correlation session
    if save_corr:
        correlation = {}
        for i in tqdm(range(len(tfs)), desc='compute_correlation'):
            correlation[tfs[i]] = correlation.get(tfs[i],{})
            for gene in sorted_LLR[tfs[i]]:
                correlation[tfs[i]][gene] = correlation[tfs[i]].get(gene)
                correlation[tfs[i]][gene] = compute_correlation(sorted_LLR[tfs[i]][gene])

        with open('data/stat/%s/%s/%s/SPcorrLLR.fl-%s.pkl' %(args.organism, args.db_name, args.fdr, args.filter), 'wb') as fout:
            pickle.dump(correlation, fout)
    else:
        fin = open('data/stat/%s/%s/%s/SPcorrLLR.fl-%s.pkl' %(args.organism, args.db_name, args.fdr, args.filter), 'rb')
        correlation = pickle.load(fin)

    # Now we compute LLR of every TF-gene by taking into account the Spearman 
    # correlation coefficient r_ek as in the formula: d_ek = α(1-max(0,r_ek))
    fc = {}
    for i in tqdm(range(len(tfs)), desc='weighted_FunCoup'):
        fc[tfs[i]] = fc.get(tfs[i],{})

        for gene in sorted_LLR[tfs[i]]:
            fc[tfs[i]][gene] = fc[tfs[i]].get(gene,0)
            N, start = len(sorted_LLR[tfs[i]][gene]) - 1, len(sorted_LLR[tfs[i]][gene]) - 2
            for col in range(N,0,-1):
                distance = 1.
                for row in range(start,0,-1):
                    distance *= alpha * (1 - max(0,float(correlation[tfs[i]][gene][row][col])))
                start -= 1 
                fc[tfs[i]][gene] += sorted_LLR[tfs[i]][gene][col][1] * distance
            fc[tfs[i]][gene] += sorted_LLR[tfs[i]][gene][0][1]
    
    if save_fc:
        with open('data/stat/%s/%s/%s/weighted_funcoup.fl-%s.pkl' %(args.organism, args.db_name, args.fdr, args.filter), 'wb') as fout:
            if args.netx:
                weighted_fc = {}
                for key in fc.keys():
                    for val in fc[key]:
                        if  fc[key][val] > 2:
                            weighted_fc[key] = weighted_fc.get(key,{})
                            weighted_fc[key][val] = weighted_fc[key].get(val,{'weight':0})
                            weighted_fc[key][val]['weight'] = fc[key][val]
                pickle.dump(weighted_fc, fout)
            else:
                pickle.dump(fc, fout)

        return(fc)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute and assign LLR to TF-gene links')
    parser.add_argument('-b','--build', dest='build', default=False, action='store_true',
                        help='activate LLRs and ENCODE building mode (by default False)')
    parser.add_argument('-db', '--database', dest='db_name', default='trrust',
                        help='type of gold standard in use (by default trrust, alternately regnet)')
    parser.add_argument('-og','--organism', dest='organism', default='human',
                        help='provide the common name of the organism under study (by default human, alternately mouse)')
    parser.add_argument('-fdr', dest='fdr', default='01',
                        help='provide the value alpha of FDR (by default 01, alternately 05)')
    parser.add_argument('-tr','--threshold', dest='threshold', default=False, action='store_true',
                        help='activate threshold mode and save just links with LLR > X (by default False)')
    parser.add_argument('-fl','--filter', dest='filter', default=False, action='store_true',
                        help='activate filter mode (by default False)')
    parser.add_argument('-tf','--TransFact', dest='tf', default=False,
                        help='provide the name of a TF to get all its links (by default False)')
    parser.add_argument('-lk','--link', dest='link', default=False, nargs=2,
                        help='provide the name of a TF and a gene to look for their link (by default False)')
    parser.add_argument('-llr','--loglikelihood', dest='llr', type=float, default=0,
                        help='provide a LLR threshold (by default 0)')
    parser.add_argument('-w','--weight', dest='weight', default=False, action='store_true',
                        help='activate FunCoup weighting schema mode (by default False)')
    parser.add_argument('-nx','--netx', dest='netx', default=False, action='store_true',
                        help='activate networkX mode to build a proper input (by default False)')
    parser.add_argument('-grn','--CancerGRN', dest='CancerGRN', default=False, action='store_true',
                        help='activate CancerGRN mode to build a proper input (by default False)')
    parser.add_argument('-ens','--ensembl', dest='ensembl', default=False, action='store_true',
                        help='activate ensembl mode to map gene symbol to Ensembl identifiers (by default False)')
    parser.add_argument('-Nens','--NOTensembl', dest='NOTensembl', default=False, action='store_true',
                        help='activate NOTensembl mode to retrieve unmapped gene symbol to Ensembl identifiers as #tf-link (by default False)')

    args = parser.parse_args()
    path_to_scores = 'data/evidences/encode/%s/chip_seq/scores/%s' %(args.organism, args.fdr)

    ## Upload the PGS 
    db_file = 'data/db/trrust/%s/trrust_db.pkl' %args.organism
    if args.db_name == 'regnet':
        db_file = 'data/db/regnet/%s/regnet_db.pkl' %args.organism
    with open(db_file, 'rb') as db:
        database = pickle.load(db) 

    ## This session is dedicated to build from scratch LLR, ENCODE and FC dataset
    if args.build:
        exp_list = [f for f in os.listdir(path_to_scores) if not f.startswith('.')]
        LLR = compute_LLR(save=True)
        ENCODE = build_encode(save=True)

        if args.weight: 
            FC = weighted_funcoup(save_LLR=True, save_corr=True, save_fc=True)
        else: FC = naïve_funcoup(save=True)
   
    ## This session is dedicated to extract from FC dataset data formatted as TF-gene with respective LLRs
    else:
        if args.weight: filein = open('data/stat/%s/%s/%s/weighted_funcoup.fl-%s.pkl' %(args.organism, args.db_name, args.fdr, args.filter), 'rb')
        else: filein = open('data/stat/%s/%s/%s/naive_funcoup.fl-%s.pkl' %(args.organism, args.db_name, args.fdr, args.filter), 'rb')
        FC = pickle.load(filein)

        if args.tf:
            for g in FC[args.tf]:
                if FC[args.tf][g][0] > args.fbs:
                    print('%s\t%s\t%0.2f' %(args.tf, g, FC[args.tf][g]))
        elif args.link:
            tf, gene = args.link[0], args.link[1]
            if gene in FC[tf]:
                print('%s-%s\t%0.2f' %(tf, gene, FC[tf][gene]))
            else:
                print('There is no %s-%s link according to FunCoup framework' %(tf, gene))
        else:
            if args.ensembl: 
                FC41map_file = open('data/db/fc4.1/%s/FC4.1_H.map.txt' %args.organism) #open('data/evidences/encode/%s/biomart_ensembl100-tfgene.txt' %args.organism) 
                header = FC41map_file.readline()
                FC41map = {}
                for line in FC41map_file:
                    line = line.rstrip().split()
                    FC41map[line[1]] = FC41map.get(line[1],line[0]) #if mart_export.txt --> tf:line[2], g:line[0]

            if args.ensembl: 
                fout1 = open('data/stat/%s/%s/%s/weighted_funcoup.ensemblID.fl-%s.tsv' %(args.organism, args.db_name, args.fdr, args.filter), 'w')
                if args.NOTensembl: fout1 = open('data/stat/%s/%s/%s/weighted_funcoup.fl-%s.notEns.tsv' %(args.organism, args.db_name, args.fdr, args.filter), 'w')
                if args.threshold: 
                    fout1 = open('data/stat/%s/%s/%s/weighted_funcoup.llr%i.ensemblID.fl-%s.tsv' %(args.organism, args.db_name, args.fdr, args.llr, args.filter), 'w')
                    if args.NOTensembl: fout1 = open('data/stat/%s/%s/%s/weighted_funcoup.llr%i.fl-%s.notEns.tsv' %(args.organism, args.db_name, args.fdr, args.llr, args.filter), 'w')
                    fout2 = open('data/stat/%s/%s/%s/weighted_funcoup.llr%i.ensemblID.fl-%s.pkl' %(args.organism, args.db_name, args.fdr, args.llr, args.filter), 'wb')
            else:
                fout1 = open('data/stat/%s/%s/%s/weighted_funcoup.fl-%s.tsv' %(args.organism, args.db_name, args.fdr, args.filter), 'w')
                if args.threshold: 
                    fout1 = open('data/stat/%s/%s/%s/weighted_funcoup.llr%i.fl-%s.tsv' %(args.organism, args.db_name, args.fdr, args.llr, args.filter), 'w')
                    fout2 = open('data/stat/%s/%s/%s/weighted_funcoup.llr%i.fl-%s.pkl' %(args.organism, args.db_name, args.fdr, args.llr, args.filter), 'wb')

            weighted_fc = {}
            if args.CancerGRN: maximum=0
            for tf in FC:
                for g in FC[tf]:
                    if args.threshold:
                        if FC[tf][g] > args.llr:
                            weighted_fc[tf] = weighted_fc.get(tf,{})
                            weighted_fc[tf][g] = weighted_fc[tf].get(g,FC[tf][g])
                            if args.CancerGRN and maximum < FC[tf][g]: maximum = FC[tf][g]
                            if args.ensembl: 
                                try: tf2, g2 = FC41map[tf], FC41map[g]
                                except: 
                                    if args.NOTensembl: fout1.write('#%s\t%s\n' %(tf2, g2))
                                    else: pass
                                else: fout1.write('%s\t%s\t%0.3f\n' %(tf2, g2, weighted_fc[tf][g]))
                            else: fout1.write('%s\t%s\t%0.3f\n' %(tf, g, weighted_fc[tf][g]))
                    else:
                        weighted_fc[tf] = weighted_fc.get(tf,{})
                        weighted_fc[tf][g] = weighted_fc[tf].get(g,FC[tf][g])
                        if args.ensembl: 
                            try: tf2, g2 = FC41map[tf], FC41map[g]
                            except: 
                                if args.NOTensembl: fout1.write('#%s\t%s\n' %(tf2, g2))
                                else: pass
                            else: fout1.write('%s\t%s\t%0.3f\n' %(tf2, g2, weighted_fc[tf][g]))
                        else: fout1.write('%s\t%s\t%0.3f\n' %(tf, g, weighted_fc[tf][g]))
            
            fout1.close()
            
            if args.threshold: 
                pickle.dump(weighted_fc, fout2)
                fout2.close()

            if args.CancerGRN:
                for tf in weighted_fc:
                    for g in weighted_fc[tf]:
                        print('%s\t%s\t1\t1\t1\t%f' %(tf, g, weighted_fc[tf][g]/maximum))
