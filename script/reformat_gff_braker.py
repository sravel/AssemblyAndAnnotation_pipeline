#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
# @package renameAndAddParent.py
# @author Sebastien Ravel



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
## BIO Python modules
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

##################################################
## Variables Globales
version="0.1"
VERSION_DATE='26/06/2017'


##################################################
## Function
def reformat_gff_braker(gff_input,gff_output):
	'''
	This function reformat braker.gtf to a standard gff file whith new gene ID
	Parameters :
		gff_input (str) : Path of the gff file to reformat
		gff_output (str) : Path of the output gff

	'''
	#Creat dict which contains all gff's information
	dico = collections.defaultdict(dict)
	dico_order = collections.defaultdict(dict)
	with open(gff_input,'r') as gff_file :
		for line in gff_file :
			scaffold,tools,type,start,end,frame,sens,unknow,ID = line.strip().split('\t')
			if type != 'gene' and tools == 'AUGUSTUS':
				ID = ID.split('"')[1].split('_')[-1].split('.')[0]
			if type != 'gene' and tools == 'GeneMark.hmm':
				ID = ID.split('"')[1].split('_')[-2]
			if ID not in dico.keys() :
				# Initiate a dico order for order the new gff file with scaffold and start position
				dico_order[scaffold][start] = ID
				index_cds = index_exon = index_intron = 1
			if type == 'intron' or type == 'INTRON' :
				type = f'{type}_{index_intron}'
				index_intron += 1
			if type == 'CDS' or type == 'cds' :
				type = f'{type}_{index_cds}'
				index_cds += 1
			if type == "exon" :
				type = f'{type}_{index_exon}'
				index_exon += 1
			dico[ID][type] = {"scaffold":scaffold,"start" : start,'end':end,"tool" : tools,'frame' : frame,
										"sens" : sens,"unknow" : unknow}
	################################################################################################
	# Write the new fasta, and create line missing of the first gff
	## The two first loops use the dico order variable for obtain the good order in the new gff file
	#for rename ID in the good order
	i = 0
	with open(gff_output,'w') as output_file :
		output_file.write('##gff-version 3\n')
		for scaffold in sorted(dico_order, key = sort_human) :
			for position in sorted(dico_order[scaffold],key = sort_human) :
				i += 1
				# Retieve ID for the dico order for write the good information
				ID = dico_order[scaffold][position]
				########################################################################################
				### In this part, I add the information missing in the old gff. For every genomics type
				### (gene, trancript, start_codon etc...) I use the CDS_1 information (except for start and end position)
				if 'start_codon' not in dico[ID] :
					if dico[ID]["CDS_1"]['sens'] == '+' :
						# For positive strand I search start codon position in the min position in all genomics type of the current ID
						start = min([int(dico[ID][type]["start"]) for type in dico[ID] ])
						end = start + 2
					elif dico[ID]["CDS_1"]['sens'] == '-' :
						# For negative strand I search start codon position in the max position in all genomics type of the current ID
						end = max([int(dico[ID][type]["end"]) for type in dico[ID] ])
						start = end - 2
					# Add start codon information at dico
					dico[ID]['start_codon'] = {"scaffold": dico[ID]["CDS_1"]['scaffold'],
										"start": start,
										"end": end,
										"tool": dico[ID]["CDS_1"]['tool'],
										"frame": dico[ID]["CDS_1"]['frame'],
										"sens": dico[ID]["CDS_1"]['sens'],
										"unknow": dico[ID]["CDS_1"]['unknow']}

				if 'stop_codon' not in dico[ID] :
					if dico[ID]["CDS_1"]['sens'] == '-' :
						# For  negative strand  I search stop codon position in the min position in all genomics type of the current ID
						start = min([int(dico[ID][type]["start"]) for type in dico[ID] ])
						end = start + 2
					elif dico[ID]["CDS_1"]['sens'] == '+' :
						# For positive strand,  I search stop codon position in the max position in all genomics type of the current ID
						end = max([int(dico[ID][type]["end"]) for type in dico[ID] ])
						start = end - 2
					# Add stop codon information at dico
					dico[ID]['stop_codon'] = {"scaffold": dico[ID]["CDS_1"]['scaffold'],
										"start":start,
										"end": end,
										"tool": dico[ID]["CDS_1"]['tool'],
										"frame": dico[ID]["CDS_1"]['frame'],
										"sens": dico[ID]["CDS_1"]['sens'],
										"unknow": dico[ID]["CDS_1"]['unknow']}

				if 'gene' not in dico[ID] :
					# For the start position of the gene, I search the min position in all genomics type of the current ID
					start = min([int(dico[ID][type]["start"]) for type in dico[ID] ])
					# For the end position of the gene, I search the max position in all genomics type of the current ID
					end = max([int(dico[ID][type]["end"]) for type in dico[ID] ])
					# Add gene information at dico
					dico[ID]['gene'] = {"scaffold": dico[ID]["CDS_1"]['scaffold'],
										"start" : start,
										"end" : end,
										"tool": dico[ID]["CDS_1"]['tool'],
										"frame": dico[ID]["CDS_1"]['frame'],
										"sens": dico[ID]["CDS_1"]['sens'],
										"unknow": dico[ID]["CDS_1"]['unknow']}

				if 'transcript' not in dico[ID] or 'mRNA' not in dico[ID] :
					# For the start position of the transcript, I search the min position in all genomics type of the current ID
					start = min([int(dico[ID][type]["start"]) for type in dico[ID] ])
					# For the start position of the transcript, I search the min position in all genomics type of the current ID
					end = max([int(dico[ID][type]["end"]) for type in dico[ID] ])
					dico[ID]['transcript'] = {"scaffold": dico[ID]["CDS_1"]['scaffold'],
										"start": start,
										"end": end,
										"tool": dico[ID]["CDS_1"]['tool'],
										"frame": dico[ID]["CDS_1"]['frame'],
										"sens": dico[ID]["CDS_1"]['sens'],
										"unknow": dico[ID]["CDS_1"]['unknow']}

				###########################################################################################
				#In this part, the script write the new content of the gff
				## In first the gene information
				gene = dico[ID]['gene']
				output_file.write(f'{gene["scaffold"]}\t'
					  f'BRAKER\t'
					  f'gene\t'
					  f'{gene["start"]}\t'
					  f'{gene["end"]}\t'
					  f'{gene["frame"]}\t'
					  f'{gene["sens"]}\t'
					  f'{gene["unknow"]}\t'
					  f'ID=g{i}\n')

				## Then the transcript information
				transcript = dico[ID]['transcript']
				output_file.write(f'{transcript["scaffold"]}\t'
					  f'BRAKER\t'
					  f'transcript\t'
					  f'{transcript["start"]}\t'
					  f'{transcript["end"]}\t'
					  f'{transcript["frame"]}\t'
					  f'{transcript["sens"]}\t'
					  f'{transcript["unknow"]}\t'
					  f'ID=g{i}.t1; Parent=g{i}\n')

				## Then the first coddon information (start codon for strand + and stop condon for strand -)
				if dico[ID]['transcript']['sens'] == '+' :
					codon_first = 'start_codon'
					codon_last = 'stop_codon'
				elif dico[ID]['transcript']['sens'] == '-' :
					codon_first = 'stop_codon'
					codon_last = 'start_codon'
				start = dico[ID][codon_first]
				output_file.write(f'{start["scaffold"]}\t'
					  f'BRAKER\t'
					  f'{codon_first}\t'
					  f'{start["start"]}\t'
					  f'{start["end"]}\t'
					  f'{start["frame"]}\t'
					  f'{start["sens"]}\t'
					  f'{start["unknow"]}\t'
					  f'ID=g{i}.t1; Parent=g{i}\n')
				# Then we write the CDS, exon and cds in the good order
				for index in sorted([int(i.split('_')[1]) for i in dico[ID] if 'CDS' in i], key = sort_human) :
					CDS = dico[ID][f'CDS_{index}']
					output_file.write(f'{CDS["scaffold"]}\t'
						  f'BRAKER\t'
						  f'CDS\t'
						  f'{CDS["start"]}\t'
						  f'{CDS["end"]}\t'
						  f'{CDS["frame"]}\t'
						  f'{CDS["sens"]}\t'
						  f'{CDS["unknow"]}\t'
						  f'ID=g{i}.t1.cds; Parent=g{i}t1\n')
					if f'exon_{index}' in dico[ID] :
						exon = dico[ID][f'exon_{index}']
					else :
						exon = dico[ID][f'CDS_{index}']
					output_file.write(f'{exon["scaffold"]}\t'
						  f'BRAKER\t'
						  f'exon\t'
						  f'{exon["start"]}\t'
						  f'{exon["end"]}\t'
						  f'{exon["frame"]}\t'
						  f'{exon["sens"]}\t'
						  f'{exon["unknow"]}\t'
						  f'ID=g{i}.t1.exon; Parent=g{i}t1\n')
					if f'intron_{index}' in dico[ID] :
						intron = dico[ID][f'intron_{index}']
						output_file.write(f'{intron["scaffold"]}\t'
							  f'BRAKER\t'
							  f'intron\t'
							  f'{intron["start"]}\t'
							  f'{intron["end"]}\t'
							  f'{intron["frame"]}\t'
							  f'{intron["sens"]}\t'
							  f'{intron["unknow"]}\t'
							  f'Parent=g{i}\n')
				# And for the end, the last codon is write
				end = dico[ID][codon_last]
				output_file.write(f'{start["scaffold"]}\t'
					  f'BRAKER\t'
					  f'{codon_last}\t'
					  f'{start["start"]}\t'
					  f'{start["end"]}\t'
					  f'{start["frame"]}\t'
					  f'{start["sens"]}\t'
					  f'{start["unknow"]}\t'
					  f'ID=g{i}.t1; Parent=g{i}\n')
##################################################
## Main code
##################################################
if __name__ == "__main__":

	# Initializations
	start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())

	# Parameters recovery
	parser = argparse.ArgumentParser(prog=__file__, description='''This Programme rename genes name''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display '+__file__+' version number and exit')
	#parser.add_argument('-dd', '--debug',choices=("False","True"), dest='debug', help='enter verbose/debug mode', default = "False")

	filesreq = parser.add_argument_group('Input  infos for running')
	filesreq.add_argument('-g', '--gff', metavar="<path/to/gff>", type=str, required=True, dest = 'gffFileIn', help = 'path to gff file')
	filesreq.add_argument('-o', '--out', metavar="<path/to/output/file>", required=True, dest = 'outDir', help = 'Path of output file')


	# Check parameters
	args = parser.parse_args()


	#Welcome message
	print("#################################################################")
	print("#              Welcome in renameGFF (Version " + version + ")               #")
	print("#################################################################\n\n")


	gffFileIn = args.gffFileIn
	ouputDir = args.outDir
	reformat_gff_braker(gffFileIn,ouputDir)
	################################################################################################
