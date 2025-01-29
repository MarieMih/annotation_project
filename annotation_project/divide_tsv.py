import csv
import re


def divide_tsv(input_file):
    """
    This function gets input annotation file from bakta (.tsv)
    and divides into six files.
    """
    pref = input_file.rpartition('.')[0]
    output_file0 = pref + '_userproteins_only.tsv'  # для которых есть UserProtein независимо от uniref100
    output_file1 = pref + '_userproteins.tsv'       # для тех, у которых UserProtein и Uniref100 совпадают
    output_file2 = pref + '_uniref100.tsv'          # с uniref100 и без UserProtein, для которых вытягиваются id Uniprot без проверки
    output_file3 = pref + '_semidefined.tsv'        # все, у кого нет UserProtein и нет uniref100
    output_file3_1 = pref + '_cds_sorf.tsv'         # белки, у которых нет UserProtein и нет uniref100
    output_file3_2 = pref + '_rna.tsv'              #

    with open(input_file, 'r', newline='') as infile,          \
         open(output_file0, 'w', newline='') as outfile0,      \
         open(output_file1, 'w', newline='') as outfile1,      \
         open(output_file2, 'w', newline='') as outfile2,      \
         open(output_file3, 'w', newline='') as outfile3,      \
         open(output_file3_1, 'w', newline='') as outfile3_1,  \
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
                writer0.writerow(row)

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
