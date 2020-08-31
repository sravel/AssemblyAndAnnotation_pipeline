#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
# @package module_config.py
# @author Florian Charriat


########## Description ##########
# This module is used for control every arguments of the config file for the Annotation pipeline

########## Module ###############
import os, time

# TODO : Maybe add in verif_fasta_file if the file is a fasta (use biopython)
######### Fonction #############

# Data directory for fastq/fasta and RNAseq
def verif_directory(path_directory_fasta,extension) :
    '''
    This function checks if the directory exists as well as
    the presence of file with the correct extension
    Parameters :
		path_directory_fasta (str) : path of the fasta/fastq directory give in input
		extension (str) : extention of files traited (fasta, fastq, fastq.gz etc ....)
    '''
    # Add the '/' caracter if the directory had none
    if path_directory_fasta.endswith('/') == False:
        path_directory_fasta = path_directory_fasta + '/'
    # Checks if the directory exist
    if os.path.isdir(path_directory_fasta) == False :
        raise ValueError(form(f"ERROR : The directory '{path_directory_fasta}' is not valid path,"
                              f" please check if your directory exists",'red', 'bold'))
    # Cheks if files with the good extention are present in the directory
    for elt in os.listdir(path_directory_fasta) :
        if elt.endswith(extension) :
            return path_directory_fasta
    raise ValueError(form(
        f"ERROR : They have not {extension} file in your input directory",'red', 'bold'))

# Protein fasta for annotation with exonerate and ET data base for repeatMasker
def verif_fasta_file(path_fasta) :
    '''
    This function checks if the fasta file exists and if is a fasta file
    Parameters :
        path_fasta (str) : path of the fasta protein file
    '''
    liste_fasta_extension = ['fasta','fa','faa','fst']
    if os.path.isfile(path_directory_fasta) == False :
        raise ValueError(form(f"ERROR : The path '{path_directory_fasta}' is not valid ,"
                              f" please check if your file exists",'red', 'bold'))
    if path_fasta.split('.')[-1]  in liste_fasta = ['fasta','fa','faa','fst'] or \
            path_fasta.split('.')[-2] in liste_fasta = ['fasta', 'fa', 'faa', 'fst'] and path_fasta.split('.')[-1] == 'gz':
        return True
    else :
        raise ValueError(form(
            f"ERROR : You file {path_fasta} is not a fasta file (the extension don't corresponding to a fasta file)", 'red', 'bold'))

# Species for augustus config
def verif_species(species_id,path_augustus_config) :
    '''
    This function checks if the species ID is in data set of augustus
    Parameters :
        species_id (str) : name of species which be used by augustus
        path_augustus_config (str) : path of the augustus config directory used by Augustus and braker
    '''
    liste_id_species = os.listdir(f'{path_augustus_config}species/')
    if species_id not in liste_id_species :
        raise ValueError(form(
            f"ERROR : You species id {path_fasta} isn't a id used by augustus, please change your id or "
            f"add you training set in {path_augustus_config}species directory",
            'red', 'bold'))

# Path of gm_key for Braker tools (geneMark licence)
def verif_gm_key(path_gm_key) :
    '''
    This function cheks if gm key file exist and if the date of creation of the file is less than the
    end of the licences (300 days)
    Parameters :
        path_gm_key (str) : path of the gm key (licence for genemark)
    '''
    if os.path.isfile(path_gm_key) == False:
        raise ValueError(form(f"ERROR : Missing licence for genemark, please get the gm_key at "
                              f"{path_gm_key}", 'red', 'bold'))
    file_creation_date = time.ctime(os.path.getctime(path_gm_key))
    date_last_modification = (time.time() - os.path.getmtime(path_gm_key) )/(60*60*24)
    if date_last_modification > 300 :
        raise ValueError(form(f"ERROR : Too old genemark-ES licence, please download an new licence (gm_key)"
                              f"at http://exon.gatech.edu/GeneMark/license_download.cgi", 'red', 'bold'))
# Output directory
def make_output_directory(path_directory) :
    '''
    This function create the output directory for pipeline
    Parameters :
        path_directory (str) : output directory path to create
    '''
    # Add the '/' caracter if the directory had none
    if path_directory.endswith('/') == False:
        path_directory = path_directory + '/'
    # Create the directory if didn't exist
    if not os.path.exists(path_directory):
        os.makedirs(path_directory)