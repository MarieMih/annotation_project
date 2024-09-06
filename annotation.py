import sys
import subprocess
import os

start_file = sys.argv[1]
uniref100 = start_file.rpartition('.')[0] + '_uniref100.tsv'
uniref100_input_upimapi = start_file.rpartition('.')[0] + '_uniref100_uniref100_ids.csv'
uniref100_output = start_file.rpartition('.')[0] + '_upimapi_ref2ref'
kb = uniref100_output + '/uniprotinfo.tsv'
kb_input_upimapi = uniref100_output + '/uniprotinfo_uniref_representative_ids.csv'
kb_output  =  uniref100_output + '/uniprotkb'
kb_tsv = start_file.rpartition('.')[0] + '_cds_sorf.tsv'
kb_faa = start_file.rpartition('.')[0] + ".faa"
annotation_input = kb_tsv.rpartition('.')[0]+"_by_bakta_tag.faa"
annotation_output = kb_output + "/annotation"
converting = os.path.split(start_file)[0]

if os.path.exists(start_file):
    print(f'The file {start_file} exists')
else:
    print(f'The file {start_file} does not exist')
    exit()

result_divide_tsv = subprocess.run(['python', 'divide_tsv.py', start_file], check = True)
result_ref100 = subprocess.run(['python', 'extract_uniref.py', uniref100], check = True)
result_upimapi_ref2ref = subprocess.run(['upimapi', '-i', uniref100_input_upimapi, \
                                         '-o', uniref100_output, \
                                         '--from-db', 'UniRef100', '--to-db', 'UniRef100'], \
                                        check = True)
result_ref2kb = subprocess.run(['python', 'ref2kb.py', kb], check = True)
result_upimapi_ref2kb = subprocess.run(['upimapi', '-i', kb_input_upimapi, \
                                         '-o', kb_output, \
                                         '--from-db', 'UniProtKB AC/ID', '--to-db', 'UniProtKB'], \
                                       check = True)
result_catch_faa = subprocess.run(['python', 'catch_ids.py', kb_tsv, kb_faa], check  = True)
result_anno = subprocess.run(['upimapi', '-i', annotation_input, \
                                         '-o', annotation_output, \
                                         '-db', 'uniprot'], \
                                       check = True)
result_convert = subprocess.run(['python', 'correcting_gff.py', converting], check = True)
