import re
import requests
import gzip
import os
import shutil
import subprocess
from io import BytesIO
from requests.adapters import HTTPAdapter, Retry
from tqdm import tqdm

re_next_link = re.compile(r'<(.+)>; rel="next"')

retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)
    return None

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers.get("x-total-results", "unknown")
        yield response, total
        batch_url = get_next_link(response.headers)

def generate_url(taxids: list):
    prefix_url = 'https://rest.uniprot.org/uniprotkb/search?compressed=true&format=fasta&query=('
    postfix_url = '+AND+(reviewed%3Atrue))&size=500'
    ref_taxids = ["(taxonomy_id%3A" + str(t) + ")" for t in taxids]
    if len(ref_taxids) > 1:
        tax_str = "+OR+".join(ref_taxids)
        tax_str = "(" + tax_str + ")"
        return prefix_url + tax_str + postfix_url
    elif len(ref_taxids) == 1:
        return prefix_url + ref_taxids[0] + postfix_url
    else:
        return None

def generate_url_trembl(taxids: list):
    prefix_url = 'https://rest.uniprot.org/uniprotkb/search?compressed=true&format=fasta&query=('
    postfix_url = '+AND+(reviewed%3Afalse))&size=500'
    ref_taxids = ["(taxonomy_id%3A" + str(t) + ")" for t in taxids]
    if len(ref_taxids) > 1:
        tax_str = "+OR+".join(ref_taxids)
        tax_str = "(" + tax_str + ")"
        return prefix_url + tax_str + postfix_url
    elif len(ref_taxids) == 1:
        return prefix_url + ref_taxids[0] + postfix_url
    else:
        return None

def check_taxids(taxids: list):
    return True

def correct_short_headers_for_bakta(input_file, output_file, long_format = False, min_identity = str(90), min_query_cov = str(80), min_subject_cov = str(80)):

    DBNAME = "UniProtKB"
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
                    if long_format:
                        header_str = ">" + substring_id + " " + min_identity + "~~~" + min_query_cov + "~~~" + min_subject_cov + "~~~" + substring_gene + "~~~" + substring_product + "~~~" + "\n"
                    else:
                        header_str = ">" + substring_id + " " + substring_gene + "~~~" + substring_product + "~~~" + "\n"
                    w.write(header_str)
                else:
                    w.write(line_file)

def create_trembl_db(PROTEIN_DB_FOLDER_USER = "/storage/data1/marmi/annotation_project/protein_db", tax = ['562']):

    PROTEIN_DB_FOLDER = os.path.join(PROTEIN_DB_FOLDER_USER, "upimapi_trembl_taxid")
    url = generate_url_trembl(tax)

    DB_FILE = os.path.join(PROTEIN_DB_FOLDER, 'uniprot_sequences_' + "_".join(tax) + '_trembl.fasta.gz')
    os.makedirs(PROTEIN_DB_FOLDER)

    with open(DB_FILE, 'ab') as f:
        progress_bar = None
        for response, total in get_batch(url):
            f.write(response.content)
            if progress_bar is None:
                if total.isdigit():
                    progress_bar = tqdm(total=int(total), desc="Downloading", unit=" entry")
                else:
                    progress_bar = tqdm(desc="Downloading", unit=" entry")
            progress_bar.update(500)
        if progress_bar:
            progress_bar.close()

    with open(DB_FILE[:-3], 'wb') as f_out:
        with gzip.open(DB_FILE, 'rb') as f_in:
            shutil.copyfileobj(f_in, f_out)

    if os.path.exists(DB_FILE):
        os.remove(DB_FILE)

def create_protein_db(PROTEIN_DB_FOLDER_USER = "/storage/data1/marmi/annotation_project/protein_db", tax = ['562']):

    PROTEIN_DB_FOLDER = os.path.join(PROTEIN_DB_FOLDER_USER, "usertaxids_colinca")
    url = generate_url(tax)

    DB_FILE = os.path.join(PROTEIN_DB_FOLDER, 'uniprot_sequences_' + "_".join(tax) + '.fasta.gz')
    os.makedirs(PROTEIN_DB_FOLDER)

    with open(DB_FILE, 'ab') as f:
        progress_bar = None
        for response, total in get_batch(url):
            f.write(response.content)
            if progress_bar is None:
                if total.isdigit():
                    progress_bar = tqdm(total=int(total), desc="Downloading", unit=" entry")
                else:
                    progress_bar = tqdm(desc="Downloading", unit=" entry")
            progress_bar.update(500)
        if progress_bar:
            progress_bar.close()

    with open(DB_FILE[:-3], 'wb') as f_out:
        with gzip.open(DB_FILE, 'rb') as f_in:
            shutil.copyfileobj(f_in, f_out)

    if os.path.exists(DB_FILE):
        os.remove(DB_FILE)

    ferr = open(os.path.join(PROTEIN_DB_FOLDER, "mmseqs2.log"), 'w')
    subprocess.run(['mmseqs', 'easy-cluster',
                    DB_FILE[:-3],
                    os.path.join(PROTEIN_DB_FOLDER, 'uniprot_sequences_' + "_".join(tax)),
                    os.path.join(PROTEIN_DB_FOLDER, "tmp"),
                    "--cov-mode", "0",
                    "-c", "0.99",
                    "--min-seq-id", "0.99"],
                check=True, stdout = ferr)
    ferr.close()

    if os.path.exists(DB_FILE[:-3]):
        os.remove(DB_FILE[:-3])

    if os.path.exists(os.path.join(PROTEIN_DB_FOLDER, "tmp")):
        shutil.rmtree(os.path.join(PROTEIN_DB_FOLDER, "tmp"))

    os.remove(os.path.join(PROTEIN_DB_FOLDER, 'uniprot_sequences_' + "_".join(tax) + "_all_seqs.fasta"))
    os.remove(os.path.join(PROTEIN_DB_FOLDER, 'uniprot_sequences_' + "_".join(tax) + "_cluster.tsv"))

    correct_short_headers_for_bakta(os.path.join(PROTEIN_DB_FOLDER, 'uniprot_sequences_' + "_".join(tax) + "_rep_seq.fasta"), os.path.join(PROTEIN_DB_FOLDER, 'uniprot_sequences_' + "_".join(tax) + "_rep.fasta"))
    os.remove(os.path.join(PROTEIN_DB_FOLDER, 'uniprot_sequences_' + "_".join(tax) + "_rep_seq.fasta"))
