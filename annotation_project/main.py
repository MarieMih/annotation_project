import os
import sys
import asyncio
from .preparation import bakta_annotation
from .preparation import send_smth
from .preparation import filtering_fastq_pe
from .preparation import assembly_unicycler_pe
from .annotation import annotation
from .pangenome.pangenome_analysis import create_directory_with_soft_links
from .pangenome.pangenome_analysis import pangenome_analysis
from .create_union_protein_fasta_from_gffs import create_fasta_file


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

    asyncio.run(send_smth(text=["Ends annotation"]))
