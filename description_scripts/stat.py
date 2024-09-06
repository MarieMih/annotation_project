import sys
import csv
from collections import defaultdict

def process_tsv(input_file):
    # Initialize dictionaries to store counts grouped by the third column's value
    counts = defaultdict(lambda: {
        'empty_sixth': 0,
        'non_empty_sixth': 0,
        'user_protein': 0,
        'total_lines': 0,
        'uncharacterized_protein': 0,
        'hypothetical_protein': 0,
        'uniprotkb': 0
    })
    
    # Open and read the input TSV file
    with open(input_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        
        for row in reader:
            # Skip header lines starting with #
            if row[0].startswith('#'):
                continue
            
            # Get the value of the second column
            key = row[1]
            
            # Increment the total line count for this key
            counts[key]['total_lines'] += 1
            
            # Check if the seventh column is empty or not
            if row[6] == '':
                counts[key]['empty_sixth'] += 1
            else:
                counts[key]['non_empty_sixth'] += 1
                
                # Check if the ninth column contains "UserProtein"
                if "UserProtein" in row[8]:
                    counts[key]['user_protein'] += 1

                # Check if the eighth column contains "Uncharacterized protein"
                if "Uncharacterized protein" in row[7]:
                    counts[key]['uncharacterized_protein'] += 1
            
            # Check if the eighth column contains "hypothetical protein"
            if "hypothetical protein" in row[7]:
                counts[key]['hypothetical_protein'] += 1

            if row[10] != '':
                counts[key]['uniprotkb'] += 1
    print(f"{sys.argv[1].rpartition('/')[2]}")
    # Print the results sorted(x.items(), key=lambda item: item[1])}
    # for key, count in counts.items():
    for key, count in sorted(counts.items(), key=lambda item: item[0]):
        print(f"Group {key}:")
        print(f"  Total lines: {count['total_lines']}")
        print(f"  Empty gene names: {count['empty_sixth']}")
        print(f"  Non-empty gene names: {count['non_empty_sixth']}")
        print(f"  Non-empty gene names and 'UserProtein': {count['user_protein']}")
        print(f"  Non-empty gene names and 'Uncharacterized protein': {count['uncharacterized_protein']}")
        print(f"  Lines with 'hypothetical protein': {count['hypothetical_protein']}")
        print(f"  UniProtKB: {count['uniprotkb']}")
        print()

# Example usage
input_file = sys.argv[1]  # Replace with your input file path
process_tsv(input_file)
