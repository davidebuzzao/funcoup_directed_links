#!/usr/bin/env python3

'''
Get in touch with ENCODE-API.

ENCODE stores metadata in standard JavaScript Object Notation (JSON) format.
This script interact with the database through an industry-standard, Hypertext-Transfer-Protocol-(HTTP)-based, 
Representational-state-transfer (RESTful) application programming interface (API). It uses the libraries  
to handle the network connection and parse the objects returned (requests and json for Python).

Usage: encode_api [-h] [-f FORMAT] [-a ASSEMBLY] files

Positional arguments:
  files                 provide the text file with urls to ENCODE files in different formats

Optional arguments:
  -h, --help            show this help message and exit
  -f FORMAT, --format FORMAT
                        provide the format of the files under study (by default bam, alternately fastq.gz, bed.gz, bigBed, bigWig)
  -a ASSEMBLY, --assembly ASSEMBLY
                        provide the assembly in use (by default GRCh38, alternately mm10)

Example of usage:
    * encode_api data/evidences/encode/human/files.txt

Reference: https://www.encodeproject.org/help/rest-api/
'''

import argparse
import requests, json
    
def get_ENCODE(url, experiment_urls, extension='bam'):

    # Force return from the server in JSON format
    headers = {'accept': 'application/json'}
    # GET the object
    response = requests.get(url, headers=headers)
    # Extract the JSON response as a Python dictionary
    data = response.json()

    # Print the data columns available
    #print(data['columns'])

    print('Number of experiments: %s' %data['total'])
    print('Experiment accession\tExperiment target\tAssay\tFile accession')

    # Extract experiment_accession, target, assay_title and specific formatted file
    for i in range(len(data['@graph'])):
        experiment_accession = data['@graph'][i]['accession']
        target = data['@graph'][i]['target']['label'] ## check if label != genes[symbol]
        assay_title = data['@graph'][i]['assay_title']
        files = [f['@id'].split('/')[2] for f in data['@graph'][i]['files']] ## is there a way to distinguish BAM formatted files?
        desired_files = [d for d in files if 'https://www.encodeproject.org/files/%s/@@download/%s.%s' %(d, d, extension) in experiment_urls]
        
        print('%s\t%s\t%s\t%s' %(experiment_accession, target, assay_title, ','.join(desired_files)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get in touch with ENCODE')
    parser.add_argument('files', help='provide the text file with urls to ENCODE files in different formats')
    parser.add_argument('-f', '--format', dest='format', default='bam',
                        help='provide the format of the files under study (by default bam, alternately fastq.gz, bed.gz, bigBed, bigWig)')
    parser.add_argument('-a', '--assembly', dest='assembly', default='GRCh38',
                        help='provide the assembly in use (by default GRCh38, alternately mm10)')
    
    args = parser.parse_args()

    # This URL locates the ENCODE human TF ChIP-seq experiments
    url = ('https://www.encodeproject.org'
    '/search/?type=Experiment'
    '&status=released'
    '&assay_title=TF+ChIP-seq'
    '&assembly=%s'
    '&assay_title=total+RNA-seq'
    '&target.investigated_as=transcription+factor'
    '&files.file_type=%s'
    '&limit=all'
    '&frame=obejct' %(args.assembly, args.format))

    # Open the text file with links to metadata and to all ENCODE files 
    filein = open(args.files)
    metadata = filein.readline()

    exp_urls = []
    for line in filein:
        if line.rstrip().split('.')[-1] == args.format: 
            exp_urls.append(line.rstrip())

    get_ENCODE(url, exp_urls, extension=args.format)
