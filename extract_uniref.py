import csv
import re

def extract_uniref(input_file):
    """
    Gets uniref100.tsv
    """

    pref = input_file.rpartition('.')[0]
    output_file_tsv = pref + '_columns.tsv'
    output_file_csv = pref + '_uniref100_ids.csv'

    # Списки для хранения данных
    extracted_data = []
    uniref100_list = []

    # Регулярное выражение для поиска UniRef100 идентификаторов
    uniref100_pattern = re.compile(r'UniRef:UniRef100_[A-Za-z0-9]+')

    # Читаем входной файл
    with open(input_file, 'r', newline='') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:
            if len(row) >= 6:  # Проверка, что строка содержит как минимум 6 колонок
                sixth_column = row[5]
                last_column = row[-1]
                
                # Поиск всех вхождений UniRef100 в последней колонке
                uniref100_ids = uniref100_pattern.findall(last_column)
                if uniref100_ids:
                    # Добавляем шестую колонку и все найденные UniRef100 идентификаторы в extracted_data
                    for uniref100_id in uniref100_ids:
                        uniref100_id = uniref100_id.split('UniRef:')[1]
                        extracted_data.append([sixth_column, uniref100_id])
                        uniref100_list.append(uniref100_id)

    # Записываем извлеченные данные в TSV-файл
    with open(output_file_tsv, 'w', newline='') as outfile_tsv:
        writer = csv.writer(outfile_tsv, delimiter='\t')
        writer.writerows(extracted_data)

    # Записываем UniRef100 идентификаторы в CSV-файл
    with open(output_file_csv, 'w', newline='') as outfile_csv:
        writer = csv.writer(outfile_csv)
        writer.writerow(uniref100_list)

    print("Файлы успешно созданы.")
