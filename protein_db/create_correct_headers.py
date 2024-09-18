import re
import sys

if len(sys.argv) != 3:
    print("Usage: python script.py input_file output_file")
    sys.exit(1)

DBNAME = "UniProtKB"

input_file = sys.argv[1]
output_file = sys.argv[2]

# pattern_id = r'>(.*?)\s+'

# pattern_id = r'>[a-zA-Z0-9_-]+\|[a-zA-Z0-9_-]+\|(.*?)\s'
pattern_id = r'>[a-zA-Z0-9_-]+\|(.*?)\|[a-zA-Z0-9_-]+\s'
pattern_product = r'>[a-zA-Z0-9_-]+\|[a-zA-Z0-9_-]+\|[a-zA-Z0-9_]+\s(.*?)\sOS='
pattern_gene = r'GN=(.*?)\s'

with open (output_file, "w") as w:
    with open (input_file, "r") as f:
        for line_file in f:
            match_id = re.search(pattern_id, line_file)
            match_product = re.search(pattern_product, line_file)
            match_gene = re.search(pattern_gene, line_file)
            if match_product:
                substring_product = match_product.group(1)
                if match_gene:
                    substring_gene = match_gene.group(1)
                else:
                    substring_gene = ""
                if match_id:
                    substring_id = f"{DBNAME}" + "|" + match_id.group(1)
                else:
                    substring_id = f"{DBNAME}"
                header_str = ">" + substring_id + " " + substring_gene + "~~~" + substring_product + "~~~" + "\n"
                w.write(header_str)
            else:
                w.write(line_file)
