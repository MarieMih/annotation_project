import csv


def divide_tsv(input_file):
    """
    This function gets input annotation file from bakta (.tsv)
    and divides into six files.
    """
    prefix_of_sample        = input_file.rpartition('.')[0]
    output_file_userprotein = prefix_of_sample + '_userproteins_only.tsv'  # для которых есть UserProtein независимо от uniref100
    # output_file1   = prefix_of_sample + '_userproteins.tsv'       # для тех, у которых UserProtein и Uniref100 совпадают
    output_file_uniref      = prefix_of_sample + '_uniref100.tsv'          # с uniref100 и без UserProtein, для которых вытягиваются id Uniprot без проверки
    output_file_unchar      = prefix_of_sample + '_semidefined.tsv'        # все, у кого нет UserProtein и нет uniref100
    output_file_unchar_cds  = prefix_of_sample + '_cds_sorf.tsv'         # белки, у которых нет UserProtein и нет uniref100
    output_file_unchar_rna  = prefix_of_sample + '_rna.tsv'              #

    with open(input_file,              'r', newline='') as infile,       \
         open(output_file_userprotein, 'w', newline='') as of_user,      \
         open(output_file_uniref,      'w', newline='') as of_unir,      \
         open(output_file_unchar,      'w', newline='') as of_unch,      \
         open(output_file_unchar_cds,  'w', newline='') as of_unch_cds,  \
         open(output_file_unchar_rna,  'w', newline='') as of_unch_rna:

        reader          = csv.reader(infile,      delimiter='\t')
        writer_user     = csv.writer(of_user,     delimiter='\t')
        writer_unir     = csv.writer(of_unir,     delimiter='\t')
        writer_unch     = csv.writer(of_unch,     delimiter='\t')
        writer_unch_cds = csv.writer(of_unch_cds, delimiter='\t')
        writer_unch_rna = csv.writer(of_unch_rna, delimiter='\t')

        all_the_rest = []

        for row in reader:
            row_str = '\t'.join(row)

            if 'UserProtein' in row_str:
                writer_user.writerow(row)
            elif 'UserProtein' not in row_str and 'UniRef' in row_str:
                writer_unir.writerow(row)
            else:
                writer_unch.writerow(row)
                all_the_rest.append(row)

        for row in all_the_rest:
            if len(row) > 1 and (row[1] == 'cds' or row[1] == 'sorf'):
                writer_unch_cds.writerow(row)
            else:
                writer_unch_rna.writerow(row)

    print("Файл успешно разделен по категориям.")
