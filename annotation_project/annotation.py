"""
Main part of annotation.
"""

import sys
import subprocess
import logging
import os
from divide_tsv import divide_tsv
from extract_uniref import extract_uniref
from ref2kb import converting_uniref_to_uniprotkb
from catch_ids import catch_ids
from correcting_gff import correcting_gff
from converting_to_gtf import convert_gff_to_gtf

UPIMAPI_RESOURCES = "/storage/data1/marmi/upimapi_databases"


def annotation(start_file):
    """
    Modification of Bakta output.
    """
    # gets tsv file!!

    log = logging.getLogger('ANNOTATION')

    pref = start_file.rpartition('.')[0]
    uniref100_data = pref + '_uniref100.tsv'
    uniref100_upimapi_search_input_file = pref + '_uniref100_uniref100_ids.csv'
    uniref100_upimapi_search_output_directory = pref + '_upimapi_ref2ref'
    kb_representative_data = uniref100_upimapi_search_output_directory + '/uniprotinfo.tsv'
    kb_upimapi_input_file = uniref100_upimapi_search_output_directory + '/uniprotinfo_uniref_representative_ids.csv'
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

    result_upimapi_ref2ref = subprocess.run(['upimapi',
                                             '-i', uniref100_upimapi_search_input_file,
                                             '-o', uniref100_upimapi_search_output_directory,
                                             '--from-db', 'UniRef100',
                                             '--to-db', 'UniRef100'],
                                            check=True)

    try:
        converting_uniref_to_uniprotkb(kb_representative_data)
    except:
        log.error('error while coverting ids!', exc_info=True)
        sys.exit('ERROR: converting_uniref_to_uniprotkb failed!')

    result_upimapi_ref2kb = subprocess.run(['upimapi',
                                            '-i', kb_upimapi_input_file,
                                            '-o', kb_upimapi_output_directory,
                                            '--from-db', 'UniProtKB AC/ID',
                                            '--to-db', 'UniProtKB'],
                                           check=True)
###### end - block for records with UniRef100 and without UserProtein

# ######## begin - block for records without UniRef100 and without UserProtein
#     try:
#         catch_ids(kb_tsv, kb_faa)
#     except:
#         log.error('catch_ids error!', exc_info=True)
#         sys.exit('ERROR: catch_ids failed!')
# ######## end - block for records without UniRef100 and without UserProtein

    try:
        file_annotation_gff = correcting_gff(file_for_converting)
    except:
        log.error('correcting gff error!', exc_info=True)
        sys.exit('ERROR: correcting_gff failed!')


    try:
        convert_gff_to_gtf(file_annotation_gff)
    except:
        log.error('converting gff error!', exc_info=True)
        sys.exit('ERROR: converting failed!')
