import os
import sys
import asyncio
import shutil
import subprocess
sys.path.append(os.path.dirname(__file__))
from annotation import annotation
from converting_to_gtf import convert_gff_to_gtf
from correct_annotation_files import correct_annotation_files
from create_union_protein_fasta_from_gffs import create_fasta_file
from pangenome.pangenome_analysis import pangenome_analysis, create_directory_with_soft_links
from preparation import assembly_unicycler_pe, bakta_annotation, send_smth, filtering_fastq_pe
from make_common_protein_fasta import make_common_protein_fasta
from metrics.stat import make_stat_file
from create_protein_trusted_list import create_protein_db, check_taxids


def pipeline_setting():
    """
    Download Bakta, upimapi, userprotein databases.
    Set telegram-bot if needed.
    """
    print("Hello! This is script for set the main databases needed in pipeline.")

    BAKTA_TRUE = ""
    BAKTA_DB_TYPE = ""
    BAKTA_DB = ""

    while (BAKTA_TRUE != "y") and (BAKTA_TRUE != "n"):
        BAKTA_TRUE = input("Do you want to install Bakta DB? (y/n): ")

    if BAKTA_TRUE == "y":
        while (BAKTA_DB_TYPE != "full") and (BAKTA_DB_TYPE != "light"):
            BAKTA_DB_TYPE = input("Do you want to install full (~70 Gb) or light (3 Gb) version of Bakta DB? (full/light): ")
        while not os.path.exists(BAKTA_DB):
            BAKTA_DB = input("Where do you want to install Bakta DB? ")
            if not os.path.exists(BAKTA_DB): print("Folder doesn't exists.")
        subprocess.run(['bakta_db', 'download',
                        "--output", BAKTA_DB,
                        "--type", BAKTA_DB_TYPE],
                        check=True)

    USERPROTEIN_TRUE = ""
    USERPROTEIN_DB_TAXIDS = ""
    USERPROTEIN_DB = ""

    while (USERPROTEIN_TRUE != "y") and (USERPROTEIN_TRUE != "n"):
        USERPROTEIN_TRUE = input("Do you want to create protein-trusted file? (y/n): ")

    if USERPROTEIN_TRUE == "y":
        while (USERPROTEIN_DB_TAXIDS == ""):
            USERPROTEIN_DB_TAXIDS = input("What taxids do you want to include? (separate by one space) ")
            real_taxids = check_taxids(USERPROTEIN_DB_TAXIDS)
            if real_taxids != True:
                print(f"Taxids {real_taxids} aren't correct!")
                USERPROTEIN_DB_TAXIDS = ''
        while not os.path.exists(USERPROTEIN_DB):
            USERPROTEIN_DB = input("Where do you want to download UserProtein DB? ")
            if not os.path.exists(USERPROTEIN_DB): print("Folder doesn't exists.")
        real_taxids = USERPROTEIN_DB_TAXIDS.strip().split(" ")
        create_protein_db(USERPROTEIN_DB, real_taxids)
        print(f"{os.path.join(USERPROTEIN_DB, 'usertaxids_colinca', 'uniprot_sequences_' + '_'.join(real_taxids) + '_rep.fasta')} created!")

    TELEGRAM_BOOL = ""
    TELEGRAM_NEW = ""

    while (TELEGRAM_BOOL != "y") and (TELEGRAM_BOOL != "n"):
        TELEGRAM_BOOL = input("Do you want to receive messages in your own telegram-bot? (y/n): ")
    if TELEGRAM_BOOL == "y":
        while (TELEGRAM_NEW != "y") and (TELEGRAM_NEW != "n"):
            TELEGRAM_NEW = input("Do you want to create new telegram-bot (otherwise will be used the one configured now in system)? (y/n): ")
        if TELEGRAM_NEW == 'y':
            print("Okey! Let's started")
            subprocess.run(['telegram-send', '-c'],
                        check=True)

def pipeline_since_fastq(directory):
    """"
    Full pipeline with filteration, assembling, annotation
    and pangenome analysis.
    """
    files = []
    for filename in os.listdir(directory):
        if filename.endswith('_1.fq.gz'):
            file_path = os.path.join(directory, filename)
            files.append(file_path)

    tsvs = []
    path_for_tsvs = directory + "/matrix_tsv"

    for i in files:
        name = os.path.split(i)[1].partition('.')[0]
        print(i, name)
        read_1, read_2 = filtering_fastq_pe(i)
        assembly = assembly_unicycler_pe(read_1, read_2)
        assembly = os.path.split(i)[0] + "/assembly_" + name + "_sub" + "/assembly.fasta"
        bakta_annotation(assembly, name[-24:])
        annotation_tsv = os.path.split(i)[0] + "/assembly_" + name + "_sub" + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + ".tsv"
        annotation(annotation_tsv)
        annotation_tsv = os.path.split(i)[0] + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + "_extended.tsv"
        tsvs.append(annotation_tsv)

        shutil.copytree(os.path.split(i)[0] + "/assembly_" + name + "_sub" + "/bakta_annotation_" + name[-24:],
                        os.path.split(i)[0] + "/bakta_annotation_" + name[-24:])
        shutil.rmtree(os.path.split(i)[0] + "/assembly_" + name + "_sub" + "/bakta_annotation_" + name[-24:])

    make_common_protein_fasta(tsvs)
    correct_annotation_files(tsvs)

    for i in tsvs:
        fn = i.replace(".tsv", ".gff3")
        convert_gff_to_gtf(fn)

    create_directory_with_soft_links(tsvs, path_for_tsvs)

    print("Start pangenome")
    pangenome_analysis(path_for_tsvs)

    print("Start creating faa")
    create_fasta_file(tsvs, directory)

    print("Start stat creation")
    make_stat_file(path_for_tsvs)

    asyncio.run(send_smth(text=["Ends annotation"]))


def pipeline_assembly_file(file):
    """
    Unfortunately, .fasta files must be in different directories
    which lies in the same one as sample.txt.
    """
    files = []
    directory = os.path.split(file)[0]
    with open(file, "r") as f:
        for i in f:
            files.append(i[:-1])

    tsvs = []
    path_for_tsvs = directory + "/matrix_tsv"

    for i in files:
        name = os.path.split(i)[1].partition('.')[0]
        print(i, name)
        bakta_annotation(i, name[-24:])
        annotation_tsv = os.path.split(i)[0] + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + ".tsv"
        annotation(annotation_tsv)
        annotation_tsv = directory + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + "_extended.tsv"
        tsvs.append(annotation_tsv)

        shutil.copytree(os.path.split(i)[0] + "/bakta_annotation_" + name[-24:],
                        directory + "/bakta_annotation_" + name[-24:])
        shutil.rmtree(os.path.split(i)[0] + "/bakta_annotation_" + name[-24:])

    make_common_protein_fasta(tsvs)
    correct_annotation_files(tsvs)

    for i in tsvs:
        fn = i.replace(".tsv", ".gff3")
        convert_gff_to_gtf(fn)

    create_directory_with_soft_links(tsvs, path_for_tsvs)

    print("Start pangenome")
    pangenome_analysis(path_for_tsvs)

    print("Start creating faa")
    create_fasta_file(tsvs, directory)

    print("Start stat creation")
    make_stat_file(path_for_tsvs)

    asyncio.run(send_smth(text=["Ends annotation"]))


def pipeline_assembly_bakta_only(directory):
    """
    Annotate all .fasta in directory by bakta with custom db.
    """
    files = []
    for filename in os.listdir(directory):
        if filename.endswith('.fasta') or filename.endswith('.fa'):
            file_path = os.path.join(directory, filename)
            files.append(file_path)

    tsvs = []
    target = directory + "/matrix_tsv_bakta"

    for i in files:
        name = os.path.split(i)[1].partition('.')[0]
        print(i, name)
        bakta_annotation(i, name[-24:])
        annotation_tsv = os.path.split(i)[0] + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + ".tsv"
        tsvs.append(annotation_tsv)

    create_directory_with_soft_links(tsvs, target)
    make_stat_file(target)


def pipeline_stat_all_tsv_in_dir(directory):
    """
    Make
    """
    make_stat_file(directory)


def polishing_annotation_for_cohort(directory):
    """
    Polishing for assemblies placed in different directories.
    Works with this pipeline assemblies and annotations only.
    Structure of input directory like:
    assembly_S1_sub/bakta_annotation_S1
    assembly_S2_sub/bakta_annotation_S2
    assembly_S3_sub/bakta_annotation_S3
    ...
    """
    files = []
    for filename in os.listdir(directory):
        if filename.endswith('_sub'):
            file_path = os.path.join(directory, filename)
            files.append(file_path)

    tsvs = []
    target = directory + "/matrix_tsv"

    for i in files:
        name = os.path.split(i)[1]
        name = name.replace("_sub", "").replace("assembly_", "")
        print(i, name)
        if os.path.exists(i + "/bakta_annotation"):
            os.rename(i + "/bakta_annotation", i + "/bakta_annotation_" + name[-24:])
        annotation_tsv = i + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + ".tsv"
        annotation(annotation_tsv)
        annotation_tsv = i + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + "_extended.tsv"
        tsvs.append(annotation_tsv)

    create_directory_with_soft_links(tsvs, target)
    print("Start pangenome")
    pangenome_analysis(target)
    print("Start creating faa")
    create_fasta_file(tsvs, directory)
    print("Start stat creation")
    make_stat_file(target)
    asyncio.run(send_smth(text=["Ends annotation"]))


def polishing_annotation(directory):
    """
    Finding all bakta_annotation_* directories in input directory and reannotate them.
    """
    files = []
    for file_path in os.listdir(directory):
        if file_path.startswith('bakta_annotation') and not file_path.endswith('.log'):
            files.append(os.path.join(directory, file_path))

    tsvs = []
    path_for_tsvs = os.path.join(directory, "matrix_tsv")

    for i in files:
        name = os.path.split(i)[-1].replace("bakta_annotation_", "")
        print(i, name)
        annotation_tsv = os.path.join(i, name[-24:] + ".tsv")
        annotation(annotation_tsv)
        annotation_tsv = os.path.join(i, name[-24:] + "_extended.tsv")
        tsvs.append(annotation_tsv)

    make_common_protein_fasta(tsvs)
    correct_annotation_files(tsvs)

    for i in tsvs:
        fn = i.replace(".tsv", ".gff3")
        convert_gff_to_gtf(fn)

    create_directory_with_soft_links(tsvs, path_for_tsvs)

    print("Start pangenome")
    pangenome_analysis(path_for_tsvs)

    print("Start creating faa")
    create_fasta_file(tsvs, directory)

    print("Start stat creation")
    make_stat_file(path_for_tsvs)

    asyncio.run(send_smth(text=["Ends annotation"]))


def pipeline_assembly(directory):
    """
    Annotate all .fasta or .fa assemblies in input directory.
    """

    files = []
    for filename in os.listdir(directory):
        if filename.endswith('.fasta') or filename.endswith('.fa'):
            file_path = os.path.join(directory, filename)
            files.append(file_path)

    for i in files:
        name = os.path.split(i)[1].partition('.')[0]
        print(i, name)
        bakta_annotation(i, name[-24:])

    files = []
    for file_path in os.listdir(directory):
        if file_path.startswith('bakta_annotation') and not file_path.endswith('.log'):
            files.append(os.path.join(directory, file_path))

    tsvs = []
    path_for_tsvs = os.path.join(directory, "matrix_tsv")

    for i in files:
        name = os.path.split(i)[-1].replace("bakta_annotation_", "")
        print(i, name)
        annotation_tsv = os.path.join(i, name[-24:] + ".tsv")
        annotation(annotation_tsv)
        annotation_tsv = os.path.join(i, name[-24:] + "_extended.tsv")
        tsvs.append(annotation_tsv)

    make_common_protein_fasta(tsvs)
    correct_annotation_files(tsvs)

    for i in tsvs:
        fn = i.replace(".tsv", ".gff3")
        convert_gff_to_gtf(fn)

    create_directory_with_soft_links(tsvs, path_for_tsvs)

    print("Start pangenome")
    pangenome_analysis(path_for_tsvs)

    print("Start creating faa")
    create_fasta_file(tsvs, directory)

    print("Start stat creation")
    make_stat_file(path_for_tsvs)

    asyncio.run(send_smth(text=["Ends annotation"]))
