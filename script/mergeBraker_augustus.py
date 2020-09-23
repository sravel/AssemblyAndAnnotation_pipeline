#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
# @package mergeBraker_augustus.py
# @author Florian Charriat

"""
	The mergeBraker_augustus script
	===============================
	:author: Charriat Florian
	:contact: florian.charriat@inra.fr
	:date: 11/06/2018
	:version: 0.1

	Script description
	------------------

	This program is used to merged Augustus ouput with braker output. This program is used in Braker_pipeline. Augustus output must be at gff3 format.
	
	Example
	-------

	>>> mergeBraker_augustus.py  -a /work/augustus.gff3 -b /work/braker.gff3 -outFile /work/result


	Help Programm
	-------------

	optional arguments:
		- \-h, --help
						show this help 	message and exit
		- \-v, --version
						display mergeBraker_augustus.py version number and exit
						
	Input mandatory infos for running:
								
		- \-a <path/to/augustus/output/file>, --file <path/to/augustus/output/file>
						path of the augustus output file
		- \-b <path/to/braker/output/file>, --file <path/to/braker/output/file>
						path of the braker output file
		- \-o <path/to/output/directory>, --outFile <path/to/output/directory>
						path and name of the output directory

"""


########## Module ###############
## Python modules
import argparse, os, sys

# Import BioPython

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from reformat_gff_braker import reformat_gff_braker
from renameGFF import renameGFF

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
	filesreq.add_argument('-a', '--augustus',type = str, default = 'None', dest = 'augustus', help = 'Path of the augustus output file')
	filesreq.add_argument('-b', '--braker',type = str,  required=True, dest = 'BRAKER', help = 'Path of the braker output file')
	filesreq.add_argument('-o', '--outFile',type = str, required=True, dest = 'outFile', help = 'Path of the output directory')

	
######### Recuperation arguments ###########

	args = parser.parse_args()
	brakerFile = os.path.abspath(args.BRAKER)
	augustusFile = os.path.abspath(args.augustus)
	outFile= os.path.abspath(args.outFile)	
########### Gestion directory ##############

	name = augustusFile.split('/')[-1].replace('.gff3','')
	brakerPath = brakerFile.replace(brakerFile.split('/')[-1],'')
	augustusPath = augustusFile.replace(augustusFile.split('/')[-1],'')


############### start message ########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("       Welcome in mergeBraker_augustus   (Version " + version + ")      ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold')+'\n')

####################### main #################"
	nouveauGene = []
	listeGene = []
	liste = []
	liste1 = []
	nbAjout = 0
	if nbAjout != 0 :
		print('Nombre de gene identique : '+str(nbAjout))
	nbAjout = 0
	## Launch reformat function for obtain a standard gff3 file
	reformat_gff_braker(gff_input = brakerFile,gff_output = brakerFile.replace("braker","new_braker"))
	print(name +' in process')
	## Open augustus file for obtain all gene information in list. This list help for comparate gene content
	## between braker and augutus
	with open(augustusFile,'r') as augustus_file :
		for lineA in augustus_file :
			if lineA[0] != '#' :
				listeGene.append(lineA)
			elif 'end gene' in lineA :
				liste.append(listeGene)
				liste1.append(listeGene)
				listeGene = []
	print('Liste Augustuse done')
	print('Nombre de gene dans le fichier augustus : '+ str(len(liste)))
	# Open Braker and compare the gene content for remove augustus gene present in braker gff
	## Loop on the gene content of augustus gff
	for elt in liste :
		gene = elt[0]
		scaffoldA = gene.split('\t')[0]
		## Open braker for compare the gene content with braker information
		with  open(brakerFile.replace("braker", "new_braker"), 'r') as braker_file:
			for lineB in braker_file :
				scaffoldB = lineB.split('\t')[0]
				# See if the scaffolds are same in augustus and braker gene information
				if scaffoldA == scaffoldB :
					scaffold = 'done'
					typeB = lineB.split('\t')[2]
					if typeB == 'gene' :
						startA = int(gene.split('\t')[3])
						endA = int(gene.split('\t')[4])
						startB = int(lineB.split('\t')[3])
						endB = int(lineB.split('\t')[4])
						id = lineB.split('\t')[-1].replace('\n','')
						### Search braker and augustus gene overlap
						if (startB <= endA <= endB) or (startB <= startA <= endB) or (startA <= endB <= endA) or (startA <= startB <= endA):
							nbAjout +=1
							# remove elt on liste 1 because the gene augustus is already in braker gff file
							liste1.remove(elt)
							break

	# Retrieve the last gene ID for add augustus gene
	with  open(brakerFile.replace("braker", "new_braker"), 'r') as braker_file:
		numberGeneBraker = braker_file.readlines()[-1].split('=g')[-1].split('.t')[0]
	print('Comparaison done')
	# Rename braker gff for obtain the same ID in the new gff merge
	renameGFF(gff_input = brakerFile.replace("braker","new_braker") , strainName = name, tools = 'BRAKER', num = 1, gff_output = f'{brakerPath}{name}_braker.gff3')
	print('rename Braker file done')
	print('Nombre de genes a ajouter : '+ str(len(liste1)))
	# Write the new augutus file with only the gene not present in braker gff
	f = open(augustusPath+name+'_prov.gff3','w')
	for elt in liste1:
		for sousElt in elt :
			f.write(sousElt)
	f.close()
	print('Write new gff3 for augustus : Done')
	# Rename augustus gff for obtain the same ID in the new gff merge
	renameGFF(gff_input = f'{augustusPath}{name}_prov.gff3' , strainName = name, tools = 'AUGUSTUS', num = (int(numberGeneBraker)+1), gff_output = f'{augustusPath}{name}_augustus.gff3')
	os.system('rm %s'%(augustusPath+name+'_prov.gff3'))
	# Concat the two gff
	os.system('cat %s%s_braker.gff3 %s%s_augustus.gff3 > %s'%(brakerPath,name,augustusPath,name,outFile))
			
				
############## end message ###########################

	print(form("\n\t---------------------------------------------------------",'yellow','bold'))
	print("\t"+form("|",'yellow','bold')+form("                    End of execution                   ",type='bold')+form("|",'yellow','bold'))
	print(form("\t---------------------------------------------------------",'yellow','bold'))




			
							
				
				
				
				
			
				
				
				
				
				
							
					
									
									
