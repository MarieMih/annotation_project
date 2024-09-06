import csv
import sys
import re

# Файлы ввода и вывода
input_file = sys.argv[1]
output_file0 = input_file.rpartition('.')[0]+'_userproteins_only.tsv'
output_file1 = input_file.rpartition('.')[0]+'_userproteins.tsv'
output_file2 = input_file.rpartition('.')[0]+'_uniref100.tsv'
output_file3 = input_file.rpartition('.')[0]+'_semidefined.tsv'
output_file3_1 = input_file.rpartition('.')[0]+'_cds_sorf.tsv'
output_file3_2 = input_file.rpartition('.')[0]+'_rna.tsv'


# Открываем файлы для чтения и записи
with open(input_file, 'r', newline='') as infile, \
     open(output_file0, 'w', newline='') as outfile0, \
     open(output_file1, 'w', newline='') as outfile1, \
     open(output_file2, 'w', newline='') as outfile2, \
     open(output_file3, 'w', newline='') as outfile3, \
     open(output_file3_1, 'w', newline='') as outfile3_1, \
     open(output_file3_2, 'w', newline='') as outfile3_2:
    
    reader = csv.reader(infile, delimiter='\t')
    writer0 = csv.writer(outfile0, delimiter='\t')
    writer1 = csv.writer(outfile1, delimiter='\t')
    writer2 = csv.writer(outfile2, delimiter='\t')
    writer3 = csv.writer(outfile3, delimiter='\t')
    writer3_1 = csv.writer(outfile3_1, delimiter='\t')
    writer3_2 = csv.writer(outfile3_2, delimiter='\t')
    
    all_the_rest = []


    for row in reader:
        row_str = '\t'.join(row)
        if 'UserProtein' in row_str and 'UniRef100' not in row_str:
            writer0.writerow(row)
        elif 'UserProtein' in row_str and 'UniRef100' in row_str:
            uniprotkb = row_str
            entry = re.compile(r"UserProtein:[^|]*\|([^,]*)")
            if isinstance(uniprotkb, str):
                match = re.search(entry, uniprotkb)
                if match:
                    uniprotkb = match.group(1)

            unirefkb = row_str
            entry = re.compile(r"UniRef:UniRef100_([^,]*)")
            if isinstance(unirefkb, str):
                match = re.search(entry, unirefkb)
                if match:
                    unirefkb = match.group(1)

            if unirefkb == uniprotkb:
                writer1.writerow(row)
            else:
                # writer0.writerow(row)
                writer2.writerow(row)
        elif 'UniRef100' in row_str and 'UserProtein' not in row_str:
            writer2.writerow(row)
        else:
            writer3.writerow(row)
            all_the_rest.append(row)

    for row in all_the_rest:
        if len(row) > 1 and (row[1] == 'cds' or row[1] == 'sorf'):
            writer3_1.writerow(row)
        else:
            writer3_2.writerow(row)

print("Файл успешно разделен на шесть файлов.")
