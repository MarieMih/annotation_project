import os
import sys
import shutil
sys.path.append(os.path.dirname(__file__))
import common_variables
from converting_to_gtf import convert_gff_to_gtf
from correct_annotation_files import correct_tsv_file, precorrect_tsv_file
from pangenome.pangenome_analysis import pangenome_curves, create_presence_absence_matrix, pangenome_tsv, pangenome_fasta
from preparation import bakta_annotation
from make_common_protein_fasta import make_common_protein_fasta
from make_common_nucl_fasta import make_common_rna_fasta
from metrics.stat import make_stat_file
from helpers import create_directory, create_directory_with_soft_links, check_file_exists
from pathlib import Path


def annotate_fasta_in_dir(directory):
    files = []
    for filename in os.listdir(directory):
        if filename.endswith('.fasta') or filename.endswith('.fa') or filename.endswith('.fna'):
            file_path = os.path.join(directory, filename)
            files.append(file_path)

    for i in files:
        name = os.path.split(i)[1].partition('.')[0]
        print(i, name)
        print("Bakta prefix: 24 symbols max. Taken latest 20 sym.")
        bakta_annotation(i, name[-24:])


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
    cluster_file = make_common_protein_fasta(common_variables.TOOL, [tsvtmp.replace(".tsv", ".faa") for tsvtmp in tsvs], common_pangenome_path, common_pangenome_path)
    
    matrix_binary  = create_presence_absence_matrix(path_for_tsvs, cluster_file, "binary", common_pangenome_path)
    matrix_numeric = create_presence_absence_matrix(path_for_tsvs, cluster_file, "numeric", common_pangenome_path)
    matrix         = create_presence_absence_matrix(path_for_tsvs, cluster_file, "locus", common_pangenome_path)

    pangenome = pangenome_tsv(path_for_tsvs, cluster_file, matrix, common_pangenome_path)
    pangenome_fasta(path_for_faa, pangenome, common_pangenome_path, "faa")
    pangenome_fasta(path_for_ffn, pangenome, common_pangenome_path, "ffn")

    for i in tsvs:
        correct_tsv_file(str(Path(i).resolve()), pangenome, cluster_file)

    # for i in tsvs:
    #     fn = i.replace(".tsv", ".gff3")
    #     convert_gff_to_gtf(fn)


    print("Start pangenome")
    coresize, pansize = pangenome_curves(matrix_binary)
    print(f"Core: {coresize}, pan: {pansize}")



    """
    Part for ncRNA.
    """

    rna_file     = make_common_rna_fasta(tsvs, common_pangenome_path)

    matrix_binary_rna  = create_presence_absence_matrix(path_for_tsvs, rna_file, "binary", common_pangenome_path, "presence_absence_matrix_rna")
    matrix_numeric_rna = create_presence_absence_matrix(path_for_tsvs, rna_file, "numeric", common_pangenome_path, "presence_absence_matrix_rna")
    matrix_rna         = create_presence_absence_matrix(path_for_tsvs, rna_file, "locus", common_pangenome_path, "presence_absence_matrix_rna")

    pangenome_rna = pangenome_tsv(path_for_tsvs, rna_file, matrix_rna, common_pangenome_path, "pangenome_table_rna.tsv")
    pangenome_fasta(path_for_ffn, pangenome_rna, common_pangenome_path, "ffn", "pangenome_rna")

    for i in tsvs:
        correct_tsv_file(str(Path(i).resolve()).replace(".tsv", "_pangenome.tsv"), pangenome_rna, rna_file)

    # print("Start stat creation")
    # make_stat_file(common_pangenome_path)
