def find_gaps_with_overlap(gff_file, output_file, pattern_prefix="ZVL_missed"):
    """
    Находит неаннотированные регионы с учетом того, что следующие записи могут начинаться раньше текущего end
    """
    # Читаем все записи
    comment = []
    features = []
    genome = []
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                comment.append(line)
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 5:
                try:
                    start = int(fields[3])
                    end = int(fields[4])
                    features.append({
                        'line': line.strip(),
                        'start': start,
                        'end': end,
                        'seqid': fields[0]
                    })
                except ValueError:
                    continue
            else:
                genome.append(line)
                continue
    
    gaps = []
    gap_counter = 1
    current_max_end = features[0]['end']
    
    for i in range(1, len(features)):
        current_feature = features[i]
        
        if current_feature['start'] > current_max_end:
            gap_start = current_max_end
            gap_end = current_feature['start'] - 1
            
            if gap_start < gap_end:
                gaps.append({
                    'seqid': current_feature['seqid'],
                    'start': gap_start,
                    'end': gap_end,
                    'counter': gap_counter
                })
                gap_counter += 1
        
        if current_feature['end'] > current_max_end:
            current_max_end = current_feature['end']
    
    with open(output_file, 'w') as f_out:

        # Сначала пишем оригинальные features
        for feature in features:
            f_out.write(feature['line'] + '\n')
        
        # Затем добавляем gaps в конец (потом отсортируем)
        for gap in gaps:
            gap_line = create_gap_line(
                gap['seqid'], gap['start'] + 1, gap['end'], pattern_prefix, gap['counter']
            )
            f_out.write(gap_line)
    
    # Пересортируем файл чтобы gaps были в правильных позициях
    sort_gff_file(output_file, comment, genome)
    
    print(f"Найдено {len(gaps)} неаннотированных регионов")

def create_gap_line(seqid, start, end, pattern_prefix, counter):
    """Создает GFF строку для gap'а"""
    attributes = f"ID={pattern_prefix}_{counter};Name=not_annotated_region;locus_tag={pattern_prefix}_{counter}"
    return f"{seqid}\tgap_finder\tnar\t{start}\t{end}\t.\t+\t0\t{attributes}\n"
    

def sort_gff_file(gff_file, comment, genome):
    """Сортирует GFF файл по позиции start"""
    headers = []
    data_lines = []
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                headers.append(line)
            else:
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    try:
                        start = int(fields[3])
                        data_lines.append((start, line))
                    except ValueError:
                        data_lines.append((0, line))
    
    # Сортируем по start
    data_lines.sort(key=lambda x: x[0])
    
    # Записываем обратно
    with open("final.gff3", 'w') as f:
        f.writelines(comment)
        for start, line in data_lines:
            f.write(line)
        f.writelines(genome)

