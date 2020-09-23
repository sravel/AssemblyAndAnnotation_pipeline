#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
# @author Florian CHARRIAT



##################################################
## Modules
##################################################
#Import MODULES_SEB
import sys, os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"
sys.path.insert(1,current_dir+'../modules/')

## Python modules
import argparse, collections
from time import localtime, strftime
from module_Flo import sort_human
import pandas as pd
## BIO Python modules
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='26/06/2017'


if __name__ == "__main__":

    ######################################################### Argument ######################################################
    # Initializations
    start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())

    # Parameters recovery
    parser = argparse.ArgumentParser(prog=__file__, description='''This Programme rename genes name''')
    parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
                        'display '+__file__+' version number and exit')
    filesreq = parser.add_argument_group('Input  infos for running')

    filesreq.add_argument('-a', '--path_aln',  type=str, required=True, dest = 'path_aln',
                          help = 'path to directory which contains all alignement')

    filesreq.add_argument('-b', '--path_busco', type=str, required=True, dest='path_busco',
                          help='path to directory which contains all busco file')

    filesreq.add_argument('-db', '--db_busco', type=str, required=True, dest='db_busco',
                          help='name of the busco data base used')

    filesreq.add_argument('-p', '--prefix', type=str, required=True, dest='prefix',
                          help='prefix of the output file')

    args = parser.parse_args()
    path_aln = args.path_aln
    path_busco = args.path_busco
    db_busco = args.db_busco
    prefix = args.prefix

################ MAIN ##########################################################################
    dico_mapping = collections.defaultdict(dict)
    for sample in os.listdir(path_aln):
        dico_mapping[sample] = collections.defaultdict(dict)
        for rnaseq in os.listdir(f'{path_aln}{sample}'):
            if rnaseq != 'bamList':
                with open(f'{path_aln}{sample}/{rnaseq}/Summary_alignement.txt', 'r') as f:
                    for line in f:
                        if 'overall alignment rate' in line :
                            percent = line.split('%')[0]
                            dico_mapping[sample][rnaseq] = float(percent)

    # Transform dico to data frame with pandas
    data = pd.DataFrame(dico_mapping)
    data = data.transpose()
    data.to_csv(f'{prefix}_mapping.txt',sep='\t')

    dico_busco = collections.defaultdict(dict)

    for sample in os.listdir(path_busco):
        if sample != 'busco_downloads' :
            file = f'{path_busco}{sample}/short_summary.specific.{db_busco}.{sample}.txt'
            with open(file,'r') as f :
                for line in f :
                    if 'number of BUSCOs:' in line :
                        number = line.split('number of BUSCOs:')[-1].split(')')[0]
                        dico_busco[sample]['number'] = int(number)
                    if 'Complete BUSCOs' in line :
                        percent = round((int(line.strip().split()[0])/dico_busco[sample]['number'])*100,2)
                        dico_busco[sample]['complete'] = percent
                    if 'Complete and single-copy BUSCOs' in line :
                        percent = round((int(line.strip().split()[0])/dico_busco[sample]['number'])*100,2)
                        dico_busco[sample]['SCO'] = percent
                    if 'Complete and duplicated BUSCOs' in line :
                        percent = round((int(line.strip().split()[0])/dico_busco[sample]['number'])*100,2)
                        dico_busco[sample]['MCO'] = percent
                    if 'Fragmented BUSCOs' in line :
                        percent = round((int(line.strip().split()[0])/dico_busco[sample]['number'])*100,2)
                        dico_busco[sample]['Frag'] = percent
                    if 'Missing BUSCOs' in line :
                        percent = round((int(line.strip().split()[0])/dico_busco[sample]['number'])*100,2)
                        dico_busco[sample]['missing'] = percent

    data = pd.DataFrame(dico_busco)
    data = data.transpose()
    data.to_csv(f'{prefix}_busco.txt',sep='\t')