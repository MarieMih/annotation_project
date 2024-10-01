import sys
import csv
from collections import defaultdict


def process_tsv(input_file):
    counts = defaultdict(lambda: {
        'empty_sixth': 0,
        'non_empty_sixth': 0,
        'gene_uni': set(),
        'user_protein': 0,
        'total_lines': 0,
        'uncharacterized_protein': 0,
        'hypothetical_protein': 0,
        'uniprotkb': 0,
        'uniprotkb_uni': set()
    })

    with open(input_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')

        for row in reader:
            if row[0].startswith('#'):
                continue

            key = row[1]
            counts[key]['total_lines'] += 1

            if row[6] == '':
                counts[key]['empty_sixth'] += 1
            else:
                counts[key]['non_empty_sixth'] += 1
                counts[key]['gene_uni'].add(row[6])

                if "UserProtein" in row[8]:
                    counts[key]['user_protein'] += 1

                if "Uncharacterized protein" in row[7]:
                    counts[key]['uncharacterized_protein'] += 1

            if "hypothetical protein" in row[7]:
                counts[key]['hypothetical_protein'] += 1

            if row[10] != '':
                counts[key]['uniprotkb'] += 1
                counts[key]['uniprotkb_uni'].add(row[10])

    print(f"{sys.argv[1].rpartition('/')[2]}")

    for key, count in sorted(counts.items(), key=lambda item: item[0]):
        print(f"Group {key}:")
        print(f"  Total lines: {count['total_lines']}")
        print(f"  Empty gene names: {count['empty_sixth']}")
        print(f"  Non-empty gene names: {count['non_empty_sixth']}")
        print(f"  Non-empty unique gene names: {len(count['gene_uni'])}")
        print(f"  Non-empty gene names and 'UserProtein': {count['user_protein']}")
        print(f"  Non-empty gene names and 'Uncharacterized protein': {count['uncharacterized_protein']}")
        print(f"  Lines with 'hypothetical protein': {count['hypothetical_protein']}")
        print(f"  UniProtKB: {count['uniprotkb']}")
        print(f"  UniProtKB (unique): {len(count['uniprotkb_uni'])}")
        print()


# input_file = sys.argv[1]
# process_tsv(input_file)
