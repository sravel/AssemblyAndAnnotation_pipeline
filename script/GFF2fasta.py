#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
# @package mergeBraker_augustus.py
# @author Florian Charriat
########## Module ###############
## Python modules
import argparse, os, sys, gzip , collections

# Import BioPython

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#Import module_Flo
from module_Flo import verifDir, createDir , form, isFasta, recupId ,verifFichier,fasta2dict,sort_human


if __name__ == "__main__":

	version="0.1" 
	
############ Argparse #####################
	parser = argparse.ArgumentParser(prog=__file__, description='''
		This program is used to merged Augustus ouput with braker output. This program is used in Braker_pipeline. Augustus output must be at gff3 format.''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
'display mergeBraker_augustus version number and exit')


	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-g', '--gff',type = str, default = 'None', dest = 'gff', help = 'Path of the gff file')
	filesreq.add_argument('-f', '--fasta',type = str, required=True, dest = 'fasta', help = 'Path of the fasta file')
	filesreq.add_argument('-p', '--prefix',type = str, required=True, dest = 'prefix', help = 'Prefix of the output files')


	
######### Recuperation arguments ###########

	args = parser.parse_args()
	gff = os.path.abspath(args.gff)
	fasta = os.path.abspath(args.fasta)
	prefix = os.path.abspath(args.prefix)

########### Gestion directory ##############
	verifFichier(gff)
	verifFichier(fasta)
	nameGFF = gff.split('/')[-1].replace('.gz','')


####################### main #################
	# Create Variable  for the start and stop codon, this variable help for evaluate the annotation
	start_codon = 'ATG'
	stop_codon = ['TGA','TAA','TAG']
	# Create dictionary for cds and gene information (start, stop and strand)
	dico_cds = collections.defaultdict(dict)
	dico_gene = collections.defaultdict(list)
	augustus = False
	fasta_dico = fasta2dict(fasta)
	liste_scaff = fasta_dico.keys()
	# Initiate the dictionary for cds
	for elt in liste_scaff :
		dico_cds[elt] = collections.defaultdict(list)
	# Open gff for retrieve all information for cds and gene
	with open(gff,'rt') as gff_file :
		for line in gff_file :
			if '# This output was generated with AUGUSTUS' in line :
				augustus = True
			if line[0] != '#' :
				tabLine = line.split('\t')
				type = tabLine[2]
				# Parse gene information
				if type == 'gene' :
				 	id = tabLine[-1].strip().split(';')[0].replace('ID=','')
				 	scaff = tabLine[0]
				 	start = tabLine[3]
				 	end = tabLine[4]
				 	sens = tabLine[6]
				 	dico_gene[id] = [scaff,start,end]
				#Parse transcript information for retrieve ID
				if type == 'mRNA' or type == 'transcript' or type == 'rRNA' or type == 'tRNA':
					# id = tabLine[-1].replace('\n','').split('=')[1].split(';')[0]
					id = tabLine[-1].strip().split(';')[0].replace('ID=','')
					scaff = f'{tabLine[0]}'
				#Parse CDS information
				if type == 'CDS' :
					start = tabLine[3]
					end = tabLine[4]
					sens = tabLine[6]
					tools = tabLine[1]
					if scaff in liste_scaff :
						dico_cds[f'{scaff}'][id].append([start,end,sens])
					else :
						print(f'{scaff} doen"t exist')

	# Initiate gene_false variable for count the number of gene without
	gene_false = 0
	# Initiate nb variable for count the number of gene
	nb = 0
	#Write 3 files for gene, cds and proteine fasta file
	with open(f'{prefix}_cds.fasta','w') as f ,\
			open(f'{prefix}_protein.fasta','w') as f_prot,\
			open(f'{prefix}_gene.fasta','w') as f_gene :
		# loop on scaffolds of the gff
		for scaffold in dico_cds.keys() :
			# Create variable seq which contain the scaffold sequence
			seq = str(fasta_dico[scaffold].seq)
			dico_scaff = dico_cds[scaffold]
			# Processes IDs from the cds dictionary
			for ID in dico_scaff.keys() :
				seq_cds = ""
				nb += 1
				for cds_pos in dico_scaff[ID]:
					start_mRNA =int(cds_pos[0]) - 1
					break
				for cds_pos in dico_scaff[ID]:
					start = int(cds_pos[0]) - 1
					stop = int(cds_pos[1])
					seq_cds = seq_cds + seq[start:stop]
					sens = cds_pos[2]
				if augustus and sens == '+' :
					seq_cds = seq_cds + seq[stop:stop+3]
				elif augustus and sens == '-' :
					seq_cds =  seq[start_mRNA-3:start_mRNA] + seq_cds
				seq_cds = Seq(seq_cds)
				gene_ID = ID.split('.')[0].replace('T0','')
				seq_gene = Seq(seq[int(dico_gene[gene_ID][1]):int(dico_gene[gene_ID][2])])
				if sens == '-' :
					seq_cds = seq_cds.reverse_complement()
				seq_cds = str(seq_cds)
				record = SeqRecord(Seq(seq_cds), id=f'{ID}', name=f'{ID}', description='')
				SeqIO.write(record, f, "fasta")
				record = SeqRecord(Seq(seq_cds).translate(), id=f'{ID}', name=f'{ID}', description='')
				SeqIO.write(record, f_prot, "fasta")
				record = SeqRecord(seq_gene, id=f'{ID}', name=f'{ID}', description='')
				SeqIO.write(record, f_gene, "fasta")
				if seq_cds[0:3].upper() != start_codon or seq_cds[-3:].upper() not in stop_codon :
					gene_false += 1
	print(f'{nameGFF}  :  {gene_false} / {nb} have not good start or stop codon' )





			
							
				
				
				
				
			
				
				
				
				
				
							
					
									
									
