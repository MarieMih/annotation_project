import csv


def divide_tsv(input_file):
    """
    This function gets input annotation file from bakta (.tsv)
    and divides into different files.
    """
    prefix_of_sample        = input_file.rpartition('.')[0]
    output_file_unchar_cds  = prefix_of_sample + '_cds_sorf.tsv'
    output_file_unchar_rna  = prefix_of_sample + '_others.tsv'

    with open(input_file,              'r', newline='') as infile,       \
         open(output_file_unchar_cds,  'w', newline='') as of_unch_cds,  \
         open(output_file_unchar_rna,  'w', newline='') as of_unch_rna:

        reader          = csv.reader(infile,      delimiter='\t')
        writer_unch_cds = csv.writer(of_unch_cds, delimiter='\t')
        writer_unch_rna = csv.writer(of_unch_rna, delimiter='\t')

        for row in reader:
            if len(row) > 1 and (row[1] == 'cds' or row[1] == 'sorf'):
                writer_unch_cds.writerow(row)
            else:
                writer_unch_rna.writerow(row)

    print("Файл успешно разделен по категориям.")
