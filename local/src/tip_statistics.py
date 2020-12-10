#!/usr/bin/env python3

'''
Compute some naïve statistics from TIP execution's results.

The FunCoup program splits every evidence's distribution of continous scores in bins by the execution
of a discretization algorithm. Through a training set, each bin is assigned probabilistic scores (likelihood ratios), 
which are calculated by dividing the evidence occurrence in the positive training set with its background frequency.

By now, if we want to verify the accuracy with which the interactions manually curated by authors of TF-gene databases
(like TRRUSTv2 or RegNetwork) are included when the TIP probabilistic method is executed to quantitatively measure the regulatory 
potential of a TF to genes, we can:
 - divide each list of scores of resulting interactions into a bunch of chunks of equal size (10 would be fine for FC),
 - sort them by decreasing significance (qvalue) and 
 - control the frequency with which they are found in the 1st chunk, in the 2nd chunk, in the 3rd chunk and so on and so forth. 

Program usage: tip_statistics [-h] [-db DB_NAME] [-ck NUM_CHUNKS] [-fdr FDR] [-fr FREQUENCY] 
                                [-ss STAT] [-og ORGANISM] [-p] [-sp] [-fl] [-ng]

Optional arguments:
  -h, --help            show this help message and exit
  -db DB_NAME, --database DB_NAME
                        type of PGS in use (by default trrust, alternately regnet)
  -ck NUM_CHUNKS, --chunks NUM_CHUNKS
                        an integer to split the dataset (by default 10)
  -fdr FDR              provide the critical value of FDR (by default 01, alternately 05)
  -fr FREQUENCY, --freq FREQUENCY
                        to compute absolute or relative frequency statistics (by default relative, alternately absolute)
  -ss STAT, --stat STAT
                        single, average or both type of statistics to do (by default average, alternately single/both)
  -og ORGANISM, --organism ORGANISM
                        provide the common name of the organism under study (by default human, alternately mouse)
  -p, --plot            activate plot mode (by default False)
  -sp, --saveplot       activate save plot mode (by default False)
  -fl, --filter         activate filter mode (by default False)
  -ng, --negold         activate statistics on NGS (by default False)

Examples of usage: 
    * To get stats on evidences with qvalue < 0.01 when PGS is TRRUSTv2 and save barplot:
        tip_statistics  -sp

    * To get stats on evidences with qvalue < 0.05 when PGS is RegNetwork, activate filter mode
     and save barplot:
        tip_statistics -db regnet -fdr 05 -fl -sp
'''
import os, pickle, argparse, pandas as pd, numpy as np, matplotlib.pyplot as plt, matplotlib.ticker as mtick

mor_type = ['Activation',
            'Repression',
            'Unknown']

def check_bin(exp_path, exp_name, tf, chunk_df, db, lack_info=False):

    '''
    The function takes in input a list of potential interactions' scores for TF-gene couples that are sorted 
    in p-value ascending order, it spearates them in 10 equal sized bins, check in which chunk the manually curated 
    interactions for the TF under study (golden standard) have been found and store the occurencies in a pandas 
    dataframe.   
    
    If lack_info is set to True, the function returns also information about:
        - TF for which no manually curated interaction is present (lack_info['no_tf'])
        - experiment for which no manually curated interaction is present (lack_info['no_exp'])
        - gene for which no curated TF-gene within the potential interactions has been found (lack_info['no_gene'])
    '''

    chunks = list(chunk_df.index)
    num_chunks = len(chunks)
    num_lines = sum(1 for line in open(exp_path)) - 1

    with open(exp_path) as exp:
        if os.path.exists('%s/%s/score.%s.sort.pkl' %(path_to_scores, exp_name, args.fdr)): 
            fin = open('%s/%s/score.%s.sort.pkl' %(path_to_scores, exp_name, args.fdr), 'rb')
            exp_bins = pickle.load(fin)
        else: 
            exp_bins = dict([(ck, []) for ck in chunks])
            count = 1
            for ck in chunks:
                while count <= (num_lines/num_chunks) * (chunks.index(ck) +1) and count <= num_lines:
                    exp_bins[ck].append(exp.readline().split()[0])
                    count += 1

        for mor in mor_type:
            try: ## if tf is not in db it raises key error and both tf,exp name are saved
                for target in db[tf][mor]:
                    count = 1
                    for ck in chunks:
                        if target in exp_bins[ck]:
                            chunk_df['abs_freq'][ck] += 1
                            break
                        else: count += 1
                    if count > num_chunks and lack_info and target not in lack_info['no_gene']:
                        lack_info['no_gene'].append(target)
            except:
                if lack_info and tf not in lack_info['no_tf']: 
                    lack_info['no_tf'].append(tf)
                    lack_info['no_exp'].append(exp_name)
    
    if not os.path.exists('%s/%s/score.%s.sort.pkl' %(path_to_scores, exp_name, args.fdr)): 
            fout = open('%s/%s/score.%s.sort.pkl' %(path_to_scores, exp_name, args.fdr),'wb')
            pickle.dump(exp_bins, fout) 
            
    if lack_info:
        return(chunk_df, lack_info)
    else: 
        return(chunk_df)

def barplot(dataframe, db_name='trrust', freq_type='relative', color='Greys_r', experiment=None, title=None):

    '''
    Plot dataframes in input. 
    If args.save_plot==True, save plots. 
    '''

    df = dataframe.plot(kind='bar', rot=0, colormap=color, edgecolor='black')   

    if freq_type == 'relative':
        df.set(xlabel="Equal size chunks", ylabel="Relative Frequency", title=title)
        df.yaxis.set_major_formatter(mtick.PercentFormatter())

    elif freq_type == 'absolute':
        df.set(xlabel="Equal size chunks", ylabel="Absolute Frequency", title=title)

    elif freq_type == 'both':
        df.set(xlabel="Equal size chunks", ylabel="Frequency", title=title)

    if args.save_plot and freq_type== 'relative':
        axes = plt.gca()
        axes.set_ylim([None,25])
        if args.neg_gold: plt.savefig('data/stat/%s/%s/neg_gold.png' %(args.organism, args.db_name), dpi=300, box_inches='tight')       
        else: plt.savefig('data/stat/%s/%s/%s/pos_gd.fl-%s.png' %(args.organism, args.db_name, args.fdr, args.filter))
        plt.close() ## otherwise plots are superimposed
    elif args.save_plot and freq_type== 'absolute':
        plt.savefig('data/stat/%s/%s/%s/single/%s.png' %(args.organism, args.db_name, args.fdr, experiment))
        plt.close() ## otherwise plots are superimposed
    else:
        axes = plt.gca()
        axes.set_ylim([None,25])
        plt.show()
        plt.close()

def counts(exp_input):

    '''
    Provided a list of experiment's names, the function has access to experiments' metadata.tsv
    and score.01.sort.tsv and gives back:
        * a list with different TF the experiments are pointing to
        * the number of TF-gene evidences graph's links 
    '''

    exp_list = [f for f in os.listdir(path_to_scores) if f in exp_input]
    tf_list = []
    links = {}
    for exp in exp_list:
        metadata_file = '%s/%s/metadata.tsv' %(path_to_scores, exp)
        scores_file = '%s/%s/score.%s.sort.tsv' %(path_to_scores, exp, args.fdr)
        with open(metadata_file) as metadata, open(scores_file) as scores:
            tf = metadata.readlines()[-1].split()[7].split('-')[0]
            if args.organism == 'mouse': tf = tf[0] + tf[1:].lower() 
            if tf not in tf_list: tf_list.append(tf)
            links[tf] = links.get(tf, [])
            for link in scores:
                link = link.rstrip().split()[0]
                if link not in links[tf]:
                    links[tf].append(link)
    num_links = 0
    for tf in links:
        num_links += len(links[tf])
    
    return(tf_list, num_links)

def main(database, color_plot='Greys_r'):

    chunks = ['ck' + str(i) for i in range(1, args.num_chunks + 1)]
    exp_list = [f for f in os.listdir(path_to_scores) if not f.startswith('.')]
    lack_info = {'no_tf': [], 'no_exp': [], 'no_gene': []}
    stat_info = {'exp':{'in_use':[], 'overall':[]}}

    if args.stat == 'single':
        for exp in exp_list:
            chunk_arr = np.zeros((args.num_chunks,2), dtype='float')

            if not os.path.exists('data/stat/%s/%s/%s/single/%s.png' %(args.organism, args.db_name, args.fdr, exp)): # --> True/False
                chunk_df = pd.DataFrame(data=chunk_arr, index=[ck for ck in chunks], columns=['abs_freq', 'rel_freq'])
                
                metadata_file = '%s/%s/metadata.tsv' %(path_to_scores, exp)
                exp_path = '%s/%s/score.%s.sort.tsv' %(path_to_scores, exp, args.fdr)
                with open(metadata_file) as metadata:
                    tf = metadata.readlines()[-1].split()[7].split('-')[0]
                    if args.organism == 'mouse': tf = tf[0] + tf[1:].lower() 
                chunk_df = check_bin(exp_path, exp, tf, chunk_df, database)

                if sum(chunk_df['abs_freq']) > 0:
                    for i in range(args.num_chunks):
                        chunk_df['rel_freq'][i] = (chunk_df['abs_freq'][i]/sum(chunk_df['abs_freq']))*100
                        stat_info['exp']['in_use'].append(exp) 

                    stat_info['exp']['overall'].append(exp) 
                    
                    if args.plot or args.save_plot: #nd sum(chunk_df['abs_freq'][:args.num_chunks//2]) > sum(chunk_df['abs_freq'][args.num_chunks//2:]):
                        barplot(chunk_df['abs_freq'], freq_type='absolute', color=color_plot, experiment=exp, title=tf)

    elif args.stat == 'average':
        exp_count = 0
        chunk_arr = np.zeros((args.num_chunks,2), dtype='float')
        chunk_df = pd.DataFrame(data=chunk_arr, index=[ck for ck in chunks], columns=['abs_freq', 'rel_freq'])

        if args.filter:
            listatf = []
            for exp in exp_list:
                tmp_arr = np.zeros((args.num_chunks,2), dtype='float')
                tmp_df = pd.DataFrame(data=tmp_arr, index=[ck for ck in chunks], columns=['abs_freq', 'rel_freq'])
                metadata_file = '%s/%s/metadata.tsv' %(path_to_scores, exp)
                exp_path = '%s/%s/score.%s.sort.tsv' %(path_to_scores, exp, args.fdr)
                with open(metadata_file) as metadata:
                    tf = metadata.readlines()[-1].split()[7].split('-')[0]
                    if args.organism == 'mouse': tf = tf[0] + tf[1:].lower() 
                    tmp_df, lack_info = check_bin(exp_path, exp, tf, tmp_df, database, lack_info=lack_info)

                    if sum(tmp_df['abs_freq']) > 0 and tf not in listatf: #and sum(tmp_df['abs_freq'][:(args.num_chunks//2)]) > sum(tmp_df['abs_freq'][(args.num_chunks//2):]):
                        chunk_df += tmp_df
                        exp_count += 1
                        stat_info['exp']['in_use'].append(exp) 
                        listatf.append(tf)
                    stat_info['exp']['overall'].append(exp) 
        else:
            for exp in exp_list:
                tmp_arr = np.zeros((args.num_chunks,2), dtype='float')
                tmp_df = pd.DataFrame(data=tmp_arr, index=[ck for ck in chunks], columns=['abs_freq', 'rel_freq'])
                metadata_file = '%s/%s/metadata.tsv' %(path_to_scores, exp)
                exp_path = '%s/%s/score.%s.sort.tsv' %(path_to_scores, exp, args.fdr)
                with open(metadata_file) as metadata:
                    tf = metadata.readlines()[-1].split()[7].split('-')[0]
                    if args.organism == 'mouse': tf = tf[0] + tf[1:].lower()

                    tmp_df, lack_info = check_bin(exp_path, exp, tf, tmp_df, database, lack_info=lack_info)                    
                    if sum(tmp_df['abs_freq']) > 0: 
                        chunk_df += tmp_df
                        exp_count += 1
                        stat_info['exp']['in_use'].append(exp) 

                    stat_info['exp']['overall'].append(exp) 

        for i in range(args.num_chunks):
            chunk_df['rel_freq'][i] = (chunk_df['abs_freq'][i]/sum(chunk_df['abs_freq']))*100
            
        if lack_info:
            with open('data/stat/%s/%s/%s/no_tf.txt' %(args.organism, args.db_name, args.fdr), 'a') as no_tf:
                for info in lack_info['no_exp']:
                    no_tf.write('%s\n' %info)
             
        if args.plot or args.save_plot:
            print(chunk_df)
            if args.neg_gold: chunk_df.to_pickle('data/stat/%s/%s/neg_gd.pkl' %(args.organism, args.db_name))
            else: chunk_df.to_pickle('data/stat/%s/%s/%s/pos_gd.fl-%s.pkl' %(args.organism, args.db_name, args.fdr, args.filter))

            title = 'TIP naïve statistics (filter=%s)' %args.filter
            if args.frequency == 'absolute': barplot(chunk_df['abs_freq'], freq_type='absolute', db_name=args.db_name, color=color_plot, title=title)
            elif args.frequency == 'relative': barplot(chunk_df['rel_freq'], freq_type='relative', db_name=args.db_name, color=color_plot, title=title)
            elif args.frequency == 'both': barplot(chunk_df, freq_type='both', color=color_plot, title=title)
            
    elif args.stat == 'both':
        main(path_to_scores, database, freq_type='absolute', color_plot='Greys_r', type_study='single')
        main(path_to_scores, database, freq_type='relative', color_plot='Greys_r', type_study='average')

    return(stat_info) #lack_info


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='TIP naïve statistics')

    parser.add_argument('-db', '--database', dest='db_name', default='trrust',
                        help='type of gold standard in use (by default trrust, alternately regnet)')
    parser.add_argument('-ck', '--chunks', dest='num_chunks', type=int, default=10,
                        help='an integer to split the dataset (by default 10)')
    parser.add_argument('-fdr', dest='fdr', default='01',
                        help='provide the value alpha of FDR (by default 01, alternately 05)')
    parser.add_argument('-fr', '--freq', dest='frequency', default='relative',
                        help='absolute or relative frequency to finally display (by default relative, alternately absolute)')
    parser.add_argument('-ss','--stat', dest='stat', default='average',
                        help='single, average or both type of statistics to do (by default average, alternately single/both)')
    parser.add_argument('-og','--organism', dest='organism', default='human',
                        help='provide the common name of the organism under study (by default human, alternately mouse)')
    parser.add_argument('-p','--plot', dest='plot', default=False, action='store_true',
                        help='activate plot mode (by default False)')
    parser.add_argument('-sp','--saveplot', dest='save_plot', default=False, action='store_true',
                        help='activate save plot mode (by default False)')
    parser.add_argument('-fl','--filter', dest='filter', default=False, action='store_true',
                        help='activate filter mode (by default False)')
    parser.add_argument('-ng','--negold', dest='neg_gold', default=False, action='store_true',
                        help='activate statistics on negative golden standard (by default False)')

    args = parser.parse_args()

    ## Define path to PGS or NGS
    path_to_scores = 'data/evidences/encode/%s/chip_seq/scores/%s' %(args.organism, args.fdr)
    db_file = 'data/db/trrust/%s/pos_gold.pkl' %args.organism
    if args.neg_gold: db_file = 'data/db/trrust/%s/neg_gold.pkl' %args.organism
    if args.db_name == 'regnet':
        db_file = 'data/db/regnet/%s/pos_gold.pkl' %args.organism
        if args.neg_gold: db_file = 'data/db/regnet/%s/neg_gold.pkl' %args.organism
    ## Check if PGS/NGS exists
    try:
        with open(db_file, 'rb') as db:
            database = pickle.load(db) 
    except:
        print('provide the db of manually curated TF interactions as a python pickled dictionary;\n\
        if not available, make use of Dataset class to build it')
        raise SystemExit
    else:
        ## Execute main() function and plot stats' results
        print(db_file)
        stat_info = main(database, color_plot='Greys_r')

        ## Print to stdout numerical stats' results 
        if not args.plot and not args.save_plot and args.stat == 'average':
            tf_in_use, link_in_use = counts(stat_info['exp']['in_use'])
            tf_overall, link_overall = counts(stat_info['exp']['overall'])
            a = np.zeros((3,3))
            statistics = pd.DataFrame(a, index=['Exp','TF', 'links'], columns=['In use', 'Overall', '% In use/Overall'])

            statistics['In use']['Exp'], statistics['In use']['TF'], statistics['In use']['links'] = len(stat_info['exp']['in_use']), len(tf_in_use), link_in_use
            statistics['Overall']['Exp'], statistics['Overall']['TF'], statistics['Overall']['links'] = len(stat_info['exp']['overall']), len(tf_overall),  link_overall
            statistics['% In use/Overall']['Exp'], statistics['% In use/Overall']['TF'], statistics['% In use/Overall']['links'] = (len(stat_info['exp']['in_use'])/len(stat_info['exp']['overall'])*100), (len(tf_in_use)/len(tf_overall))*100, (link_in_use/link_overall)*100
            print(statistics)
            if args.neg_gold: statistics.to_pickle('data/stat/%s/%s/neg_gd.stat.pkl' %(args.organism, args.db_name))
            else: statistics.to_pickle('data/stat/%s/%s/%s/pos_gd.fl-%s.stat.pkl' %(args.organism, args.db_name, args.fdr, args.filter))

            ## save tf and exp in use
            f1 = open('data/stat/%s/%s/%s/tf.fl-%s.txt' %(args.organism, args.db_name, args.fdr, args.filter), 'w')
            f2 = open('data/stat/%s/%s/%s/exp.fl-%s.txt' %(args.organism, args.db_name, args.fdr, args.filter), 'w')
            for tf in tf_in_use:
                f1.write('%s\n' %tf)
            for exp in stat_info['exp']['in_use']:
                f2.write('%s\n' %exp)