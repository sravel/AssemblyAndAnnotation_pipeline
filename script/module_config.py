#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
# @package module_config.py
# @author Florian Charriat


########## Description ##########
# This module is used for control every arguments of the config file for the Annotation pipeline

########## Module ###############
import os, time

#TODO : Maybe add in verif_fasta_file if the file is a fasta (use biopython)
######### Fonction #############
def form(text, col='white', type='none'):
    '''
    Used to format the texts displayed on the terminal.
    :Parameters:
         text (str) : Text to format
         col (str) : The desired color between the colors red, green, yellow, orange, blue and purple
         type (str/list)  : bold, underline, blind et highligth
    '''
    W = '\033[0'  # white (normal)
    R = '\033[31'  # red
    G = '\033[32'  # green
    Y = '\033[33'  # yellow
    O = '\033[33'  # orange
    B = '\033[34'  # blue
    P = '\033[35'  # purple
    end = '\033[0m'  # white (normal)
    Bold = ';1'
    underline = ';4'
    blind = ';5'
    highlight = ';7'
    text = 'm' + text
    if 'bold' in type:
        text = Bold + text
    if 'underline' in type:
        text = underline + text
    if 'highlight' in type:
        text = blind + text
    if 'highlight' in type:
        text = highlight + text
    if col == 'red':
        return R + text + end
    elif col == 'white':
        return W + text + end
    elif col == 'green':
        return G + text + end
    elif col == 'yellow':
        return Y + text + end
    elif col == 'orange':
        return O + text + end
    elif col == 'blue':
        return B + text + end
    elif col == 'purple':
        return P + text + end

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
    if os.path.isfile(path_fasta) == False :
        raise ValueError(form(f"ERROR : The path '{path_fasta}' is not valid , please check if your file exists",'red', 'bold'))
    if path_fasta.split('.')[-1]  in liste_fasta_extension or path_fasta.split('.')[-2] in liste_fasta_extension and path_fasta.split('.')[-1] == 'gz':
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
    return path_directory

def verif_busco_id(ids) :
    '''
    This function cheks if the busco id is in data base
    '''
    liste_db_busco = ['bacteria_odb10','acidobacteria_odb10','actinobacteria_phylum_odb10',
                  'actinobacteria_class_odb10','corynebacteriales_odb10','micrococcales_odb10',
                  'propionibacteriales_odb10','streptomycetales_odb10','streptosporangiales_odb10',
                  'coriobacteriia_odb10','coriobacteriales_odb10','aquificae_odb10','bacteroideteschlorobi_group_odb10',
                  'bacteroidetes_odb10','bacteroidia_odb10','bacteroidales_odb10','cytophagia_odb10','cytophagales_odb10',
                  'flavobacteriia_odb10','flavobacteriales_odb10','sphingobacteriia_odb10','chlorobi_odb10',
                  'chlamydiae_odb10','chloroflexi_odb10','cyanobacteria_odb10','chroococcales_odb10','nostocales_odb10',
                  'oscillatoriales_odb10','synechococcales_odb10','firmicutes_odb10','bacilli_odb10','bacillales_odb10',
                  'lactobacillales_odb10','clostridia_odb10','clostridiales_odb10','thermoanaerobacterales_odb10',
                  'selenomonadales_odb10','tissierellia_odb10','tissierellales_odb10','fusobacteria_odb10',
                  'fusobacteriales_odb10','planctomycetes_odb10','proteobacteria_odb10','alphaproteobacteria_odb10',
                  'rhizobiales_odb10','rhizobiumagrobacterium_group_odb10','rhodobacterales_odb10','rhodospirillales_odb10',
                  'rickettsiales_odb10','sphingomonadales_odb10','betaproteobacteria_odb10','burkholderiales_odb10',
                  'neisseriales_odb10','nitrosomonadales_odb10','deltaepsilonsubdivisions_odb10','deltaproteobacteria_odb10',
                  'desulfobacterales_odb10','desulfovibrionales_odb10','desulfuromonadales_odb10','epsilonproteobacteria_odb10',
                  'campylobacterales_odb10','gammaproteobacteria_odb10','alteromonadales_odb10','cellvibrionales_odb10',
                  'chromatiales_odb10','enterobacterales_odb10','legionellales_odb10','oceanospirillales_odb10','pasteurellales_odb10',
                  'pseudomonadales_odb10','thiotrichales_odb10','vibrionales_odb10','xanthomonadales_odb10','spirochaetes_odb10',
                  'spirochaetia_odb10','spirochaetales_odb10','synergistetes_odb10','tenericutes_odb10','mollicutes_odb10',
                  'entomoplasmatales_odb10','mycoplasmatales_odb10','thermotogae_odb10','verrucomicrobia_odb10','archaea_odb10',
                  'thaumarchaeota_odb10','thermoprotei_odb10','thermoproteales_odb10','sulfolobales_odb10','desulfurococcales_odb10',
                  'euryarchaeota_odb10','thermoplasmata_odb10','methanococcales_odb10','methanobacteria_odb10','methanomicrobia_odb10',
                  'methanomicrobiales_odb10','halobacteria_odb10','halobacteriales_odb10','natrialbales_odb10','haloferacales_odb10',
                  'eukaryota_odb10','alveolata_odb10','apicomplexa_odb10','aconoidasida_odb10','plasmodium_odb10','coccidia_odb10',
                  'euglenozoa_odb10','fungi_odb10','ascomycota_odb10','dothideomycetes_odb10','capnodiales_odb10','pleosporales_odb10',
                  'eurotiomycetes_odb10','chaetothyriales_odb10','eurotiales_odb10','onygenales_odb10','leotiomycetes_odb10',
                  'helotiales_odb10','saccharomycetes_odb10','sordariomycetes_odb10','glomerellales_odb10','hypocreales_odb10',
                  'basidiomycota_odb10','agaricomycetes_odb10','agaricales_odb10','boletales_odb10','polyporales_odb10',
                  'tremellomycetes_odb10','microsporidia_odb10','mucoromycota_odb10','mucorales_odb10','metazoa_odb10',
                  'arthropoda_odb10','arachnida_odb10','insecta_odb10','endopterygota_odb10','diptera_odb10','hymenoptera_odb10',
                  'lepidoptera_odb10','hemiptera_odb10','mollusca_odb10','nematoda_odb10','vertebrata_odb10','actinopterygii_odb10',
                  'cyprinodontiformes_odb10','tetrapoda_odb10','mammalia_odb10','eutheria_odb10','euarchontoglires_odb10',
                  'glires_odb10','primates_odb10','laurasiatheria_odb10','carnivora_odb10','cetartiodactyla_odb10',
                  'sauropsida_odb10','aves_odb10','passeriformes_odb10','stramenopiles_odb10','viridiplantae_odb10',
                  'chlorophyta_odb10','embryophyta_odb10','liliopsida_odb10','poales_odb10','eudicots_odb10',
                  'brassicales_odb10','fabales_odb10','solanales_odb10']
    if ids not in liste_db_busco :
        raise ValueError(
            form(f"ERROR : {id} is not inclued in the current busco data base, please cheks on busco data base", 'red', 'bold'))