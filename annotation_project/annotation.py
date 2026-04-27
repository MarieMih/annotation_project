"""
Main part of annotation.
"""

import sys
import subprocess
import logging
import os
sys.path.append(os.path.dirname(__file__))
from divide_tsv import divide_tsv
from divide_fasta import divide_fasta
from catch_ids import catch_ids
from correcting_gff import correcting_gff
from converting_to_gtf import convert_gff_to_gtf
from helpers import check_file_exists
import common_variables


def annotation(start_file):
    """
    Modification of Bakta output. Get V350045701_L04_26_1.tsv
    """

    log = logging.getLogger('ANNOTATION')

    check_file_exists(start_file)

    prefix_of_sample                          = start_file.rpartition('.')[0]
    uniref100_data                            = prefix_of_sample + '_uniref100.tsv'
    uniref100_upimapi_search_input_file       = prefix_of_sample + '_uniref100_uniref100_ids.csv'
    uniref100_upimapi_search_output_directory = prefix_of_sample + '_upimapi_ref2ref'
    kb_upimapi_output_directory               = os.path.join(uniref100_upimapi_search_output_directory, 'uniprotkb')
    file_for_converting                       = os.path.dirname(start_file)

    try:
        divide_tsv(start_file)
    except:
        log.error('Wrong genome file format!', exc_info=True)
        sys.exit('ERROR: wrong genome file format!')

###### begin - block for records with UniRef100 and without UserProtein
    try:
        extract_uniref(uniref100_data)
    except:
        log.error('extract_uniref error!', exc_info=True)
        sys.exit('ERROR: extract_uniref failed!')

    if not os.path.exists(kb_upimapi_output_directory):
        os.makedirs(kb_upimapi_output_directory)
    
    subprocess.run(['upimapi',
                    '-i', uniref100_upimapi_search_input_file,
                    '-o', kb_upimapi_output_directory,
                    '--from-db', 'UniProtKB AC/ID',
                    '--to-db', 'UniProtKB',
                    '--columns', "Entry&Entry Name&Gene Names&Protein names&EC number&Function [CC]&Pathway&Keywords&Protein existence&Gene Ontology (GO)&Protein families&Taxonomic lineage&Taxonomic lineage (Ids)&Taxonomic lineage IDs (SPECIES)&Taxonomic lineage (SPECIES)&Organism&Organism (ID)&BioCyc&BRENDA&CDD&eggNOG&Ensembl&InterPro&KEGG&Pfam&Reactome&RefSeq&UniPathway",
                    '-t', '1'],
                    check=True)
###### end - block for records with UniRef100 and without UserProtein

    try:
        file_annotation_gff = correcting_gff(file_for_converting)
    except:
        log.error('correcting gff error!', exc_info=True)
        sys.exit('ERROR: correcting_gff failed!')

    try:
        divide_fasta_res = divide_fasta(start_file)
    except:
        log.error('divide fasta error!', exc_info=True)
        sys.exit('ERROR: divide fasta failed!')

