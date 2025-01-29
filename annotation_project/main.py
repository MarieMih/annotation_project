import os
import asyncio
from preparation import bakta_annotation
from preparation import send_smth
from preparation import filtering_fastq_pe
from preparation import assembly_unicycler_pe
from annotation import annotation
from pangenome.pangenome_analysis import create_directory_with_soft_links
from pangenome.pangenome_analysis import pangenome_analysis
from create_union_protein_fasta_from_gffs import create_fasta_file
from metrics.stat import make_stat_file


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
    target = directory + "/matrix_tsv"

    for i in files:
        name = os.path.split(i)[1].partition('.')[0]
        print(i, name)
        read_1, read_2 = filtering_fastq_pe(i)
        assembly = assembly_unicycler_pe(read_1, read_2)
        assembly = os.path.split(i)[0] + "/assembly_" + name + "_sub" + "/assembly.fasta"
        bakta_annotation(assembly, name[-24:])
        annotation_tsv = os.path.split(i)[0] + "/assembly_" + name + "_sub" + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + ".tsv"
        annotation(annotation_tsv)
        annotation_tsv = os.path.split(i)[0] + "/assembly_" + name + "_sub" + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + "_extended.tsv"
        tsvs.append(annotation_tsv)

    create_directory_with_soft_links(tsvs, target)
    print("Start pangenome analysis")
    pangenome_analysis(target)
    print("Start creating faa")
    create_fasta_file(tsvs, directory)
    print("Start stat creation")
    make_stat_file(target)

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
    target = directory + "/matrix_tsv"

    for i in files:
        name = os.path.split(i)[1].partition('.')[0]
        print(i, name)
        bakta_annotation(i, name[-24:])
        annotation_tsv = os.path.split(i)[0] + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + ".tsv"
        annotation(annotation_tsv)
        annotation_tsv = os.path.split(i)[0] + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + "_extended.tsv"
        tsvs.append(annotation_tsv)

    create_directory_with_soft_links(tsvs, target)
    print("Start pangenome")
    pangenome_analysis(target)
    print("Start creating faa")
    create_fasta_file(tsvs, directory)
    print("Start stat creation")
    make_stat_file(target)

    asyncio.run(send_smth(text=["Ends annotation"]))


def pipeline_assembly(directory):
    """
    Fortunately, .fasta files could be in one directory.
    """
    files = []
    for filename in os.listdir(directory):
        if filename.endswith('.fasta'):
            file_path = os.path.join(directory, filename)
            files.append(file_path)

    tsvs = []
    target = directory + "/matrix_tsv"

    for i in files:
        name = os.path.split(i)[1].partition('.')[0]
        print(i, name)
        bakta_annotation(i, name[-24:])
        annotation_tsv = os.path.split(i)[0] + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + ".tsv"
        annotation(annotation_tsv)
        annotation_tsv = os.path.split(i)[0] + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + "_extended.tsv"
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
    Add polishing to bakta annotation.
    """
    files = []
    for filename in os.listdir(directory):
        if filename.endswith('.tsv') and not filename.endswith('.hypotheticals.tsv'):
            file_path = os.path.join(directory, filename)
            files.append(file_path)

    for i in files:
        name = os.path.split(i)[1].partition('.')[0]
        print(i, name)
        annotation(i)

    asyncio.run(send_smth(text=["Ends polishing annotation"]))


def pipeline_assembly_bakta_only(directory):
    """
    Annotate all .fasta in directory by bakta with custom db.
    """
    files = []
    for filename in os.listdir(directory):
        if filename.endswith('.fasta'):
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


def polishing_annotation_debug(directory):
    """
    Fortunately, .fasta files could be in one directory.
    """
    files = []
    for filename in os.listdir(directory):
        if filename.endswith('.fasta'):
            file_path = os.path.join(directory, filename)
            files.append(file_path)

    tsvs = []
    target = directory + "/matrix_tsv"

    for i in files:
        name = os.path.split(i)[1].partition('.')[0]
        print(i, name)
        annotation_tsv = os.path.split(i)[0] + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + ".tsv"
        # annotation(annotation_tsv)
        annotation_tsv = os.path.split(i)[0] + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + "_extended.tsv"
        tsvs.append(annotation_tsv)

    create_directory_with_soft_links(tsvs, target)
    print("Start pangenome")
    pangenome_analysis(target)
    print("Start creating faa")
    create_fasta_file(tsvs, directory)
    print("Start stat creation")
    make_stat_file(target)
    asyncio.run(send_smth(text=["Ends annotation"]))




# polishing_annotation_debug("/storage/data1/marmi/annotation_project/fidelity_dataset_v3")


def polishing_annotation_for_cohort(directory):
    """
    Fortunately, .fasta files could be in one directory.
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
        if os.path.exists(os.path.split(i)[0] + "/bakta_annotation"):
            os.rename(os.path.split(i)[0] + "/bakta_annotation", os.path.split(i)[0] + "/bakta_annotation_" + name[-24:])
        annotation_tsv = os.path.split(i)[0] + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + ".tsv"
        # annotation(annotation_tsv)
        annotation_tsv = os.path.split(i)[0] + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + "_extended.tsv"
        tsvs.append(annotation_tsv)

    create_directory_with_soft_links(tsvs, target)
    print("Start pangenome")
    #pangenome_analysis(target)
    print("Start creating faa")
    create_fasta_file(tsvs, directory)
    print("Start stat creation")
    make_stat_file(target)
    asyncio.run(send_smth(text=["Ends annotation"]))
