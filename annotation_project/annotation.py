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
from extract_uniref import extract_uniref
from catch_ids import catch_ids
from correcting_gff import correcting_gff
from converting_to_gtf import convert_gff_to_gtf
import common_variables


def annotation(start_file):
    """
    Modification of Bakta output. Get V350045701_L04_26_1.tsv
    """

    log = logging.getLogger('ANNOTATION')

    pref = start_file.rpartition('.')[0]
    uniref100_data = pref + '_uniref100.tsv'
    uniref100_upimapi_search_input_file = pref + '_uniref100_uniref100_ids.csv'
    uniref100_upimapi_search_output_directory = pref + '_upimapi_ref2ref'
    kb_upimapi_output_directory = uniref100_upimapi_search_output_directory + '/uniprotkb'
    kb_tsv = pref + '_cds_sorf.tsv'
    kb_faa = pref + ".faa"
    file_for_converting = os.path.split(start_file)[0]

###### begin - block of creating files for process
    if os.path.exists(start_file):
        print(f'The file {start_file} exists')
    else:
        print(f'The file {start_file} does not exist')
        exit()

    try:
        divide_tsv(start_file)
    except:
        log.error('Wrong genome file format!', exc_info=True)
        sys.exit('ERROR: wrong genome file format!')
###### end - block of creating files for process

###### begin - block for records with UniRef100 and without UserProtein
    try:
        extract_uniref(uniref100_data)
    except:
        log.error('extract_uniref error!', exc_info=True)
        sys.exit('ERROR: extract_uniref failed!')

    if not os.path.exists(kb_upimapi_output_directory):
        os.makedirs(kb_upimapi_output_directory)
    
    result_upimapi_ref2kb = subprocess.run(['upimapi',
                                            '-i', uniref100_upimapi_search_input_file,
                                            '-o', kb_upimapi_output_directory,
                                            '--from-db', 'UniProtKB AC/ID',
                                            '--to-db', 'UniProtKB',
                                            '-rd', common_variables.UPIMAPI_RESOURCES,
                                            '-t', '1'],
                                           check=True)
###### end - block for records with UniRef100 and without UserProtein

# ######## begin - block for records without UniRef100 and without UserProtein
#     try:
#         catch_ids(kb_tsv, kb_faa)
#     except:
#         log.error('catch_ids error!', exc_info=True)
#         sys.exit('ERROR: catch_ids failed!')
# ######## end - block for records without UniRef100 and without UserProtein

######## begin - block creating gff
    try:
        file_annotation_gff = correcting_gff(file_for_converting)
    except:
        log.error('correcting gff error!', exc_info=True)
        sys.exit('ERROR: correcting_gff failed!')
######## end - block creating gff

###### begin - block of creating fasta files for process
    try:
        divide_fasta_res = divide_fasta(start_file)
    except:
        log.error('divide fasta error!', exc_info=True)
        sys.exit('ERROR: divide fasta failed!')
###### begin - block of creating fasta files for process
