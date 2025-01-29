import csv
import ast


def converting_uniref_to_uniprotkb(input_file):
    """
    Gets id of representation members and find uniprotkb ids.
    Input: uniprotinfo.tsv from output of upimapi with information about representation members.
    Output: _uniref_representative_ids.tsv/.csv.
    """

    pref = input_file.rpartition('.')[0]
    output_file_tsv = pref + '_uniref_representative_ids.tsv'
    output_file_csv = pref + '_uniref_representative_ids.csv'

    extracted_data = []
    uniref100_list = []

    with open(input_file, 'r', newline='') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:
            if len(row) >= 2:
                first_column = row[0]
                json_str = row[1]

                try:
                    json_data = ast.literal_eval(json_str)

                    if 'accessions' in json_data['representativeMember']:
                        member_id_value = json_data['representativeMember']['accessions'][0]
                        extracted_data.append([first_column, member_id_value])
                        uniref100_list.append(member_id_value)

                except (SyntaxError, ValueError) as e:
                    print(f"Предупреждение при парсинге JSON в строке: {row}, error: {e}")

    with open(output_file_tsv, 'w', newline='') as outfile_tsv:
        writer = csv.writer(outfile_tsv, delimiter='\t')
        writer.writerows(extracted_data)

    with open(output_file_csv, 'w', newline='') as outfile_csv:
        writer = csv.writer(outfile_csv)
        writer.writerow(uniref100_list)

    print("Файл успешно создан.")
