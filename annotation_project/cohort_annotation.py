import os
import sys
import shutil
sys.path.append(os.path.dirname(__file__))
import common_variables
from converting_to_gtf import convert_gff_to_gtf
from correct_annotation_files import correct_tsv_file, precorrect_tsv_file
from pangenome_analysis import pangenome_curves, create_presence_absence_matrix, pangenome_tsv, pangenome_fasta
from preparation import bakta_annotation
from make_common_protein_fasta import make_cds_clusters
from make_common_feature_fasta import make_feature_clusters
from metrics.stat import make_stat_file
from helpers import create_directory, create_directory_with_soft_links, check_file_exists
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed


def _annotate_single_fasta(fasta_path):
    name = os.path.split(fasta_path)[1].partition('.')[0]
    print(fasta_path, name)
    print("Bakta prefix: 24 symbols max. Taken latest 20 sym.")
    bakta_annotation(fasta_path, name[-24:])


def annotate_fasta_in_dir(directory, jobs=1):
    files = []
    for filename in os.listdir(directory):
        if filename.endswith('.fasta') or filename.endswith('.fa') or filename.endswith('.fna'):
            files.append(os.path.join(directory, filename))

    if jobs is None or jobs < 1:
        jobs = 1

    if jobs == 1:
        for fasta in files:
            _annotate_single_fasta(fasta)
        return

    with ThreadPoolExecutor(max_workers=jobs) as executor:
        future_to_fasta = {executor.submit(_annotate_single_fasta, fasta): fasta for fasta in files}
        for future in as_completed(future_to_fasta):
            fasta = future_to_fasta[future]
            try:
                future.result()
            except Exception as exc:
                print(f"Error annotating {fasta}: {exc}")


def cohort_annotation(directory, data_line):
    tsvs      = []
    inference = []
    faa       = []
    ffn       = []
    files     = []
    for file_path in os.listdir(directory):
        if file_path.startswith('bakta_annotation') and not file_path.endswith('.log'):
            files.append(os.path.join(directory, file_path))

    for i in files:
        name = os.path.split(i)[-1].replace("bakta_annotation_", "")
        annotation_tsv = os.path.join(i, name[-24:] + ".tsv")
        print(annotation_tsv)

        if check_file_exists(annotation_tsv.replace(".tsv", ".bakta.tsv")) != 0:
            shutil.copy(annotation_tsv, annotation_tsv.replace(".tsv", ".bakta.tsv"))
        precorrect_tsv_file(annotation_tsv)

        tsvs.append(annotation_tsv)
        inference.append(annotation_tsv.replace(".tsv", ".inference.tsv"))
        faa.append(annotation_tsv.replace(".tsv", ".faa"))
        ffn.append(annotation_tsv.replace(".tsv", ".ffn"))
        

    common_pangenome_path = os.path.join(directory, "pangenome_data_" + data_line)
    path_for_tsvs         = os.path.join(common_pangenome_path, "tsvs")
    path_for_infr         = os.path.join(common_pangenome_path, "inferences")
    path_for_faa          = os.path.join(common_pangenome_path, "faa")
    path_for_ffn          = os.path.join(common_pangenome_path, "ffn")
    create_directory(common_pangenome_path)
    create_directory_with_soft_links(tsvs, path_for_tsvs)
    create_directory_with_soft_links(inference, path_for_infr)
    create_directory_with_soft_links(faa, path_for_faa)
    create_directory_with_soft_links(ffn, path_for_ffn)

    """
    Part for CDS and sORF.
    """
    cluster_file = make_cds_clusters(common_variables.TOOL, [tsvtmp.replace(".tsv", ".faa") for tsvtmp in tsvs], common_pangenome_path, common_pangenome_path)
    data = pangenome_existing_feature(path_for_tsvs, common_pangenome_path, tsvs, cluster_file, feature="cds", type="aa", path_for_ffn=path_for_ffn, path_for_faa=path_for_faa)

    print("Start pangenome")
    coresize, pansize = pangenome_curves(data[0])
    print(f"Core: {coresize}, pan: {pansize}")

    """
    Part for ncRNA.
    """

    rna_file = make_feature_clusters(common_variables.TOOL_FOR_RNA, tsvs, common_pangenome_path, common_pangenome_path)
    pangenome_existing_feature(path_for_tsvs, common_pangenome_path, tsvs, rna_file, feature="rna", type="nc", path_for_ffn=path_for_ffn, path_for_faa=path_for_faa)



def pangenome_existing_feature(path_for_tsvs, common_pangenome_path, tsvs, cluster_file, feature="rna", type="aa", path_for_ffn=None, path_for_faa=None):
    matrix_binary  = create_presence_absence_matrix(path_for_tsvs, cluster_file, "binary", common_pangenome_path, f"presence_absence_matrix_{feature}")
    matrix_numeric = create_presence_absence_matrix(path_for_tsvs, cluster_file, "numeric", common_pangenome_path, f"presence_absence_matrix_{feature}")
    matrix         = create_presence_absence_matrix(path_for_tsvs, cluster_file, "locus", common_pangenome_path, f"presence_absence_matrix_{feature}")

    pangenome = pangenome_tsv(path_for_tsvs, cluster_file, matrix, common_pangenome_path, f"pangenome_table_{feature}")

    match type:
        case "aa":
            pangenome_fasta(path_for_faa, pangenome, common_pangenome_path, "faa", f"pangenome_{feature}")
            pangenome_fasta(path_for_ffn, pangenome, common_pangenome_path, "ffn", f"pangenome_{feature}")
        case "nc":
            pangenome_fasta(path_for_ffn, pangenome, common_pangenome_path, "ffn", f"pangenome_{feature}")
        case _:
            pass

    for i in tsvs:
        correct_tsv_file(str(Path(i).resolve()).replace(".tsv", "_pangenome.tsv"), pangenome, cluster_file)

    return [matrix_binary, matrix_numeric, matrix, pangenome]