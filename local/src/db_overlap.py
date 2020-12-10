#!/usr/bin/env python3

'''
Compute overlap of nodes and links between databases'network.

Program usage: db_overlap [-h] [-db DB_NAME] [-og ORGANISM] [-ov TARGET_OVERLAP] 
                            [-fdr FDR] [-r] [-t] [-e] [-ng] [-p] [-sp]

Optional arguments:
    -h, --help            show this help message and exit
    -db DB_NAME, --database DB_NAME
                            type of gold standard in use (by default trrust, alternately regnet)
    -og ORGANISM, --organism ORGANISM
                            the organism (by default human, alternatly mouse)
    -ov TARGET_OVERLAP, --orverlap TARGET_OVERLAP
                            provide the overlaps target (by default tf, alternatly link)
    -fdr FDR              provide the value alpha of FDR (by default 01, alternately 05)
    -r, --regnet          activate regnet mode (by default False)
    -t, --trrust          activate trrust mode (by default False)
    -e, --encode          activate encode mode (by default False)
    -ng, --negold         activate statistics on negative golden standard (by default False)
    -p, --plot            activate plot mode (by default False)
    -sp, --saveplot       activate save plot mode (by default False)

Examples of usage: 
    * db_overlap  -db trrust -og human -ov link -sp
    * db_overlap -db regnet -og mouse -ov tf -p -fdr 05

'''

import argparse, numpy as np, os, pickle
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3, venn2_unweighted, venn3_circles

def tf_overlap():

    f1 = open('data/evidences/encode/%s/tf.id' %args.organism)
    encode = set(f1.read().splitlines())

    if args.regnet and args.encode and not args.trrust:
        f2 = open('data/db/regnet/%s/tf.id' %args.organism)
        regnet = set(f2.read().splitlines())
        diagram = venn2([regnet,encode], set_labels=('RegNetwork', 'ENCODE'), set_colors=('purple', 'skyblue'), alpha = 0.7)
        plt.title('TF overlap between RegNetwork and ENCODE')
        plt.show()

    elif args.trrust and args.encode and not args.regnet:
        f3 = open('data/db/trrust/%s/tf.id' %args.organism)
        trrust = set(f3.read().splitlines())
        plt.close()
        diagram2 = venn2([trrust,encode], set_labels=('TRRUSTv2', 'ENCODE'), set_colors=('purple', 'skyblue'), alpha = 0.7)
        plt.title('TF overlap between TRRUST and ENCODE')
        plt.show()

    elif args.regnet and args.trrust and args.encode:
        f2 = open('data/db/regnet/%s/tf.id' %args.organism)
        regnet = set(f2.read().splitlines())
        f3 = open('data/db/trrust/%s/tf.id' %args.organism)
        trrust = set(f3.read().splitlines())
        diagram3 = venn3([regnet,trrust,encode], set_labels=('RegNetwork', 'TRRUSTv2', 'ENCODE'), set_colors=('purple', 'skyblue', 'blue'), alpha = 0.7)
        plt.title('TF overlap between RegNetwork, TRRUST and ENCODE')
        plt.show(diagram3)

def links_overlap():
    
    ## Map experiment accession file(s) to target (TF)
    exp_list = [f for f in os.listdir('data/evidences/encode/%s/chip_seq/scores/%s/' %(args.organism, args.fdr))]
    tf_exp = {}
    for exp1 in exp_list:
        metadata_file = 'data/evidences/encode/%s/chip_seq/scores/%s/%s/metadata.tsv' %(args.organism, args.fdr, exp1)
        with open(metadata_file) as metadata: 
            tf = metadata.readlines()[-1].split()[7].split('-')[0]
            if args.organism == 'mouse': tf = tf[0] + tf[1:].lower() 
        tf_exp[tf] = tf_exp.get(tf,[])
        tf_exp[tf] += [exp1]

    ## Build ENCODE TF-gene dictionary 
    mor_type = ['Activation', 'Repression', 'Unknown']
    ENCODE = dict((tf,{'Unknown':[]}) for tf in tf_exp)

    for tf1 in ENCODE:
        links = []
        for exp2 in tf_exp[tf1]:
            scores_file = 'data/evidences/encode/%s/chip_seq/scores/%s/%s/score.%s.sort.tsv' %(args.organism, args.fdr, exp2, args.fdr)
            with open(scores_file) as scores: 
                for link in scores:
                    links.append(link.rstrip().split()[0])
        ENCODE[tf1]['Unknown'] += list(set(links))

    ## Import PSG dataset
    db_file = 'data/db/trrust/%s/trrust_db.pkl' %args.organism
    if args.neg_gold: db_file = 'data/db/trrust/%s/neg_gold.pkl' %args.organism
    if args.db_name == 'regnet':
        db_file = 'data/db/regnet/%s/regnet_db.pkl' %args.organism
        if args.neg_gold: db_file = 'data/db/regnet/%s/neg_gold.pkl' %args.organism
    
    with open(db_file, 'rb') as db: PSG = pickle.load(db) 

    ## Set operations between ENCODE and PSG 
    labels = ['inter', 'diff1', 'diff2'] #inter:PSG & ENCODE, diff1:PSG - ENCODE, diff2:ENCODE - PSG
    overlap = dict((tf2, dict((lab1, []) for lab1 in labels)) for tf2 in PSG)
    TF_inter = list(set(ENCODE.keys()) & set(PSG.keys())) #shared TFs between ENCODE and PSG
    TF_diff1 = list(set(set(PSG.keys()) - set(ENCODE.keys()))) #TF uniquely present in PSG
    TF_diff2 = list(set(ENCODE.keys()) - set(PSG.keys())) #TF uniquely present in ENCODE
    inter, diff1, diff2 = 0,0,0
    for tf3 in TF_inter:
        overlap[tf3]['inter'] = list(set(PSG[tf3]['Unknown']) & set(ENCODE[tf3]['Unknown']))
        overlap[tf3]['diff1'] = list(set(PSG[tf3]['Unknown']) - set(ENCODE[tf3]['Unknown']))
        overlap[tf3]['diff2'] = list(set(ENCODE[tf3]['Unknown']) - set(PSG[tf3]['Unknown']))
        inter, diff1, diff2 = inter+len(overlap[tf3]['inter']), diff1+len(overlap[tf3]['diff1']), diff2+len(overlap[tf3]['diff2'])
    for tf4 in TF_diff1:
        diff1 += len(set(PSG[tf4]['Unknown']))
    for tf5 in TF_diff2:
        diff2 += len(set(ENCODE[tf5]['Unknown']))

    print(inter, diff1, diff2)
    if args.plot or args.save_plot:
        label = 'TRRUSTv2'
        if args.db_name == 'regnet': label = 'RegNetwork'
        if args.neg_gold: label = 'NGS'
        diagram = venn2_unweighted(subsets=(diff1,diff2,inter), set_labels=(label, 'ENCODE'), set_colors=('purple', 'skyblue'), alpha = 0.7)
        plt.title('TF-gene overlap between %s and ENCODE' %(label))
        if args.plot: plt.show()
        elif args.save_plot: plt.savefig('data/stat/%s/%s/%s/links.ENCODE-%s.png' %( args.organism, args.db_name, args.fdr, label))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute overlap of nodes and links between databases\'network')

    parser.add_argument('-db', '--database', dest='db_name', default='trrust',
                        help='type of gold standard in use (by default trrust, alternately regnet)')
    parser.add_argument('-og', '--organism', dest='organism', default='human',
                        help='the organism (by default human, alternatly mouse)')
    parser.add_argument('-ov', '--orverlap', dest='target_overlap', default='tf',
                        help='provide the overlaps target (by default tf, alternatly link)')
    parser.add_argument('-fdr', dest='fdr', default='01',
                        help='provide the value alpha of FDR (by default 01, alternately 05)')
    parser.add_argument('-r', '--regnet', dest='regnet', default=False, action='store_true',
                        help='activate regnet mode (by default False)') 
    parser.add_argument('-t', '--trrust', dest='trrust', default=False, action='store_true',
                        help='activate trrust mode (by default False)') 
    parser.add_argument('-e', '--encode', dest='encode', default=False, action='store_true',
                        help='activate encode mode (by default False)') 
    parser.add_argument('-ng','--negold', dest='neg_gold', default=False, action='store_true',
                        help='activate statistics on negative golden standard (by default False)')
    parser.add_argument('-p','--plot', dest='plot', default=False, action='store_true',
                        help='activate plot mode (by default False)')
    parser.add_argument('-sp','--saveplot', dest='save_plot', default=False, action='store_true',
                        help='activate save plot mode (by default False)')

    args = parser.parse_args()

    if args.target_overlap == 'tf':
        tf_overlap()
    elif args.target_overlap == 'link':
        links_overlap()