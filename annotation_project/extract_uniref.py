import csv
import re


def extract_uniref(input_file):
    """
    Gets uniref100.tsv
    Input: file _uniref100.tsv from bakta with records containing UniRef100 
    Output: _columns.tsv with columns - locus_tag of protein and UniRef100
    Output: _uniref100_ids.csv - all UniRef100 separated by comma for sending to upimapi
    """

    pref = input_file.rpartition('.')[0]
    output_file_tsv = pref + '_columns.tsv'
    output_file_csv = pref + '_uniref100_ids.csv'

    extracted_data = []
    uniref100_list = []

    uniref100_pattern = re.compile(r'UniRef:UniRef100_[A-Za-z0-9]+')

    with open(input_file, 'r', newline='') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:
            if len(row) >= 6:
                sixth_column = row[5]
                last_column = row[-1]

                uniref100_ids = uniref100_pattern.findall(last_column)
                if uniref100_ids:
                    for uniref100_id in uniref100_ids:
                        uniref100_id = uniref100_id.split('UniRef:')[1]
                        extracted_data.append([sixth_column, uniref100_id])
                        uniref100_list.append(uniref100_id)

    with open(output_file_tsv, 'w', newline='') as outfile_tsv:
        writer = csv.writer(outfile_tsv, delimiter='\t')
        writer.writerows(extracted_data)

    with open(output_file_csv, 'w', newline='') as outfile_csv:
        writer = csv.writer(outfile_csv)
        writer.writerow(uniref100_list)

    print("Файлы успешно созданы.")
