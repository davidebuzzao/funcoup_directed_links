#!/usr/bin/env python3

'''
Use the following class Golden Standards to build positive and negative gold standards (PGS/NGS).

Gold standard datasets are in Biology collections of data for which experimental evidence exists 
and has been manually curated. For the purpose of the funcoup_directed, the gold standard to look for becomes 
a dataset of experimentally-verified regulatory interactions between TFs and target genes. 
In order to train a naïve Bayesian network, one requires both datasets with positive and negative examples, 
whose data quality will determine the resulting performance. 

PGS
The positive gold standard is s set of known FC. For the project two main palusible GS have been investigated:
TRRUSTv2[1] and RegNetwork[2].

Program Usage:
    from local.src.Dataset import GoldenStandard
    GoldenStandard()
        .build('data/db/trrust/human/trrust_rawdata.human.tsv)
        .save()
    or 
    GoldenStandard(db_name='regnet')
        .build('data/db/regnet/human/tf-gene.human.tsv')
        .save()

NGS
In FunCoup it is assumed that interactions that are not in the PGS can be considered as NGS. 

In this simplified case:
    1.  picking a number of random interactions between genes would not bring to a realistic dataset, 
        being transcription factors an exclusive portion of them.
    2.  for a regulatory network, where links connect TFs to genes with direction, we can distinguish 
        ingoing and outgoing connectivity for every node and it’s proved scale freeness is still valid 
        only when we’re looking at outgoing connections [1].

That’s why one could instead:
    *   choose a list of known TFs in the respective genome;
    *   pick for each TF a bunch of random connections* from all the genes; eliminate TF-gene interactions 
        found in the PGS.
    
Program Usage:
    from local.src.Dataset import GoldenStandard
    GoldenStandard(setype='negative')
        .build()
        .save()
    or 
    GoldenStandard(setype='negative', db_name='regnet')
        .build()
        .save()

References:
[1] Han, H., Cho, J. W., Lee, S., Yun, A., Kim, H., Bae, D., ... Lee, I. (2018). 
    TRRUST v2: An expanded reference database of human and mouse transcriptional 
    regulatory interactions. Nucleic Acids Research, 46(D1), D380–D386. https://doi.org/10.1093/nar/gkx1013
[2] Liu, Z. P., Wu, C., Miao, H., & Wu, H. (2015). 
    RegNetwork: An integrated database of transcriptional and post-transcriptional regulatory 
    networks in human and mouse. Database, 2015, 1–12. https:// doi.org/10.1093/database/bav095

'''

from sys import argv
import numpy as np, pickle, random

class GoldenStandard():
    '''
    The class Dataset is intended to be used to build a dicitonary 
    contatining information about TF-gene logic interactions and, 
    according to availability, mode of regulation (activation, repression, unknown). 
    '''

    def __init__(self, setype='positive', db_name='trrust', organism='human'):
        self.mor_type = ['Activation', 'Repression', 'Unknown']
        self.setype = setype
        self.db_name = db_name
        self.organism = organism

        if setype == 'positive': 
            f = open('data/db/%s/%s/tf.id' %(db_name, organism))
            id_list = f.read().splitlines()
            self.dataset = dict((id, dict((mor, []) for mor in self.mor_type)) for id in id_list)

        elif setype == 'negative':
            evidences_id = open('data/evidences/encode/%s/tf.id' %organism) 
            self.dataset = dict((id.rstrip(), {'Unknown',[]}) for id in evidences_id)
        else:
            self.dataset = None

    def build(self, links=None, mode_reg=False):
        self.mode_reg = mode_reg
        if self.setype == 'positive':
            try: links is not None
            except:
                print('Provide TF-gene links dataset')
                raise SystemExit
            else: self._build_positive(links)   

        elif self.setype == 'negative':
            self._build_negative()

        return(self)

    def _build_positive(self, links):
        ############################################################################
        ## Simply open the dataset with curated links, parse the file and store
        ## in a key-value fashion TF-gene interactions. Here the program runs only on 
        ## TRRUSTv2 and RegNetwork datasets.
        with open(links) as data:
            if self.db_name == 'trrust':
                for line in data:
                    line = line.rstrip().split()
                    tf, gene, mor = line[0], line[1], line[2]
                    self.dataset[tf] = self.dataset.get(tf,True)
                    if self.mode_reg:
                        self.dataset[tf][mor].append(gene)
                    else:
                        if gene not in self.dataset[tf]['Unknown']:
                            self.dataset[tf]['Unknown'].append(gene)

            elif self.db_name == 'regnet':
                for line in data:
                    line = line.rstrip().split()
                    tf, gene = line[0], line[1]
                    self.dataset[tf] = self.dataset.get(tf,True)
                    self.dataset[tf]['Unknown'].append(gene)
        
        return(self)

    def _build_negative(self):
        ############################################################################
        ## Load the positive gold standard in use 
        pos_gold = self._load('data/db/%s/%s/pos_gold.pkl' %(self.db_name, self.organism))
        
        ############################################################################
        ## Count the average and std of links per TF frequency distribution in the pos gs
        lenghts = []
        for i in pos_gold: 
            links = 0
            for m in pos_gold[i]: 
                for j in pos_gold[i][m]: 
                    links += 1 
            lenghts.append(links)
        lenghts_corrected = [lenghts[i] for i in range(len(lenghts)) if lenghts[i] < 200] # set a number to exlcude TF with > 200 links
        average_length, stdev = float(np.mean(lenghts_corrected)), float(np.std(lenghts_corrected))
        print('PGS num of links\nAverage: %i\nStd: %i' %(round(average_length),round(stdev)))
        
        ##########################################################################################
        ## Build the negative gold standard with available TF in ENCODE randomly connecting  
        ## a number of never seen interactions equivalent to the positive gold standard's 
        ## (average + stdev/2) number of links per TF
        evidences_id = open('data/evidences/encode/%s/tf.id' %self.organism) 
        neg_gold_tf = list(self.dataset.keys())
        pos_gold_tf = list(pos_gold.keys())

        data = open('data/annotations/%s/gene.id' %self.organism)
        gene_list = data.read().splitlines()

        random.shuffle(neg_gold_tf) 
        for i in range(len(neg_gold_tf)):
            random.shuffle(gene_list) 
            self.dataset[neg_gold_tf[i]['Unknown']] += gene_list[:round(average_length+(stdev//2))]

        links = []
        for tf in pos_gold: 
            for m in self.mor_type: 
                for j in pos_gold[tf][m]: 
                    try:  
                        if j in self.dataset[tf]['Unknown']:
                            index = self.dataset[tf]['Unknown'].index(j) 
                            del self.dataset[tf]['Unknown'][index] 
                    except: next 
            if tf in self.dataset: links.append(self.dataset[tf]['Unknown'])
        return(self)

    def save(self) -> object:
        if self.setype == 'positive': fout_name = 'pos_gold'
        elif self.setype == 'negative': fout_name = 'neg_gold' 
        with open('data/db/%s/%s/%s.pkl' %(self.db_name, self.organism, fout_name), 'wb') as fileout:
            pickle.dump(self.dataset, fileout)
        return(self)

    def _load(self, path=False):
        try: path != False
        except: 
            print('Method usage: obj.load(path=X)')
            raise SystemExit
        else:
            with open(path, 'rb') as filein:
                return(pickle.load(filein))
        
    def fetch_dict(self) -> dict:
        return(self.dataset)

    def __len__(self):
        return len(self.dataset)