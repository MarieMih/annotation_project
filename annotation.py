import sys
import subprocess
import logging
import os
from divide_tsv import divide_tsv
from extract_uniref import extract_uniref
from ref2kb import converting_uniref_to_uniprotkb
from catch_ids import catch_ids
from correcting_gff import correcting_gff

log = logging.getLogger('ANNOTATION')

start_file = sys.argv[1] # gets tsv file!!
pref = start_file.rpartition('.')[0]
uniref100 = pref + '_uniref100.tsv'
uniref100_input_upimapi = pref + '_uniref100_uniref100_ids.csv'
uniref100_output = pref + '_upimapi_ref2ref'
kb = uniref100_output + '/uniprotinfo.tsv'
kb_input_upimapi = uniref100_output + '/uniprotinfo_uniref_representative_ids.csv'
kb_output  =  uniref100_output + '/uniprotkb'
kb_tsv = pref + '_cds_sorf.tsv'
kb_faa = pref + ".faa"
annotation_input = kb_tsv.rpartition('.')[0]+"_by_bakta_tag.faa"
annotation_output = kb_output + "/annotation"
converting = os.path.split(start_file)[0]

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

try:
  extract_uniref(uniref100)
except:
  log.error('extract_uniref error!', exc_info=True)
  sys.exit('ERROR: extract_uniref failed!')

result_upimapi_ref2ref = subprocess.run(['upimapi', '-i', uniref100_input_upimapi, \
                                         '-o', uniref100_output, \
                                         '--from-db', 'UniRef100', '--to-db', 'UniRef100'], \
                                        check = True)

try:
  converting_uniref_to_uniprotkb(kb)
except:
  log.error('error while coverting ids!', exc_info=True)
  sys.exit('ERROR: converting_uniref_to_uniprotkb failed!')

result_upimapi_ref2kb = subprocess.run(['upimapi', '-i', kb_input_upimapi, \
                                         '-o', kb_output, \
                                         '--from-db', 'UniProtKB AC/ID', '--to-db', 'UniProtKB'], \
                                       check = True)
try:
  catch_ids(kb_tsv, kb_faa)
except:
  log.error('catch_ids error!', exc_info=True)
  sys.exit('ERROR: catch_ids failed!')

result_anno = subprocess.run(['upimapi', '-i', annotation_input, \
                                         '-o', annotation_output, \
                                         '-db', 'uniprot'], \
                                       check = True)

try:
  correcting_gff(converting)
except:
  log.error('correcting gff error!', exc_info=True)
  sys.exit('ERROR: correcting_gff failed!')
