#!/usr/bin/env python3

'''
Reworks the content of a Wiggle file to save space.

A file in Wiggle format can stores a lot of bits and I/O operations can be really RAM consuming.
The scripts is intended:
    *   first, to split a file in wiggle format that is sorted in order of chromosomes
        and is composed of header and depth scores at variable steps;
    *   second, to rebuild every chromosome's file in wiggle format in a way it contains 
        the same information but in a condensed fashion.

Program usage: compress_wiggle [-h] [-og ORGANISM] [-r] experiments wigfilename

positional arguments:
  experiments           provide the path to the directory where the experiment is stored
  wigfilename           provide the file accession

optional arguments:
  -h, --help            show this help message and exit
  -og ORGANISM, --organism ORGANISM
                        the organism (by default human, alternatly mouse)
  -r, --remove          remove useless wiggle files once compressed
'''

import os, argparse
from tqdm import tqdm

def split_chr():
    
    '''
    The function takes in input a .wig file that is sorted in order of chromosomes
    and is composed of header and depth scores at variable steps.
    '''
    wigfile = '%s%s/%s.wig' %(args.experiments, args.wigfilename, args.wigfilename)
    num_lines = sum(1 for line in open(wigfile))

    with open(wigfile) as filein:
        header = filein.readline() ## header's not interesting
        line = filein.readline()
        num_line = 1
        while num_line < num_lines:
            chr = line.rstrip().split()[1].split('=')[1]
            if chr not in chr_nam:
                line, num_line = filein.readline(), num_line + 1
                while not line.startswith('variableStep') and line:
                    line = filein.readline()
                    num_line += 1
            else:
                with open('%s%s/%s.wig' %(args.experiments, args.wigfilename, chr), 'w') as fileout:
                    fileout.write('variableStep chrom=%s\n' %chr)
                    line, num_line = filein.readline(), num_line + 1
                    while not line.startswith('variableStep') and line:
                        fileout.write(line)
                        line = filein.readline()
                        num_line += 1

def compress_wig(wigfilename, zigfilename):

    '''
    The function takes in input a .chr*.wig file and reorganize the information
    so that it contains exactly the same information but in a ~1/20 of bits fashion.
    '''

    num_lines = sum(1 for line in open(wigfilename))

    if num_lines > 1:
        with open(wigfilename) as filein, open(zigfilename,'w') as fileout:
            line = filein.readline()
            fileout.write(line)
            
            num_l = 1

            line = filein.readline().split()
            start, pos, val, new_val = int(line[0]), int(line[0]), line[1], line[1]
            all_pos = [start]

            while num_l < num_lines - 2:
                count = 0
                while (pos - (start + count)) <= 1 and new_val == val and num_l <= num_lines - 2:
                    line = filein.readline().rstrip().split()
                    pos = int(line[0])
                    new_val = line[1]
                    all_pos.append(pos)
                    count += 1
                    num_l += 1

                end = str(all_pos[-2])
                fileout.write('%s-%s %s\n' %(start,end,val))

                start, val = all_pos[-1], new_val
                all_pos = [start]
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reworks the content of a Wiggle file to save space')
    parser.add_argument('experiments', help='provide the path to the directory where the experiment is stored')
    parser.add_argument('wigfilename', help='provide the file accession')

    parser.add_argument('-og', '--organism', dest='organism', default='human',
                        help='the organism (by default human, alternatly mouse)')
    parser.add_argument('-r', '--remove', dest='removing', default=False, action='store_true',
                        help='remove useless wiggle files once compressed') 
    args = parser.parse_args()

    if args.organism == 'human':
        chr_nam = ['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY', 'chrM']
    elif args.organism == 'mouse':
        chr_nam = ['chr' + str(i) for i in range(1,20)] + ['chrX', 'chrY', 'chrM']

    
    split_chr()

    for chr in tqdm(range(len(chr_nam)), desc='Compress_chr'):
        filein = '%s%s/%s.wig' %(args.experiments, args.wigfilename, chr_nam[chr])
        if os.path.exists(filein): # --> True/False
            fileout = '%s%s/zig/%s.zig' %(args.experiments, args.wigfilename, chr_nam[chr])
            if not os.path.exists(fileout):
                compress_wig(filein, fileout) 
            if args.removing: os.remove(filein)
    if args.removing: os.remove('%s%s/%s.wig' %(args.experiments, args.wigfilename, args.wigfilename))