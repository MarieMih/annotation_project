import os
import asyncio
from annotation import annotation
from converting_to_gtf import convert_gff_to_gtf
from correct_annotation_files import correct_annotation_files
from create_union_protein_fasta_from_gffs import create_fasta_file
from pangenome.pangenome_analysis import pangenome_analysis
from pangenome.pangenome_analysis import create_directory_with_soft_links
from preparation import assembly_unicycler_pe
from preparation import bakta_annotation
from preparation import send_smth
from preparation import filtering_fastq_pe
from make_common_protein_fasta import make_common_protein_fasta
from metrics.stat import make_stat_file


# def pipeline_since_fastq_old(directory):
#     """"
#     Full pipeline with filteration, assembling, annotation
#     and pangenome analysis.
#     """
#     files = []
#     for filename in os.listdir(directory):
#         if filename.endswith('_1.fq.gz'):
#             file_path = os.path.join(directory, filename)
#             files.append(file_path)

#     tsvs = []
#     target = directory + "/matrix_tsv"

#     for i in files:
#         name = os.path.split(i)[1].partition('.')[0]
#         print(i, name)
#         read_1, read_2 = filtering_fastq_pe(i)
#         assembly = assembly_unicycler_pe(read_1, read_2)
#         assembly = os.path.split(i)[0] + "/assembly_" + name + "_sub" + "/assembly.fasta"
#         bakta_annotation(assembly, name[-24:])
#         annotation_tsv = os.path.split(i)[0] + "/assembly_" + name + "_sub" + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + ".tsv"
#         annotation(annotation_tsv)
#         annotation_tsv = os.path.split(i)[0] + "/assembly_" + name + "_sub" + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + "_extended.tsv"
#         tsvs.append(annotation_tsv)

#     create_directory_with_soft_links(tsvs, target)
#     print("Start pangenome analysis")
#     pangenome_analysis(target)
#     print("Start creating faa")
#     create_fasta_file(tsvs, directory)
#     print("Start stat creation")
#     make_stat_file(target)

#     asyncio.run(send_smth(text=["Ends annotation"]))


# def pipeline_assembly_file_old(file):
#     """
#     Unfortunately, .fasta files must be in different directories
#     which lies in the same one as sample.txt.
#     """
#     files = []
#     directory = os.path.split(file)[0]
#     with open(file, "r") as f:
#         for i in f:
#             files.append(i[:-1])

#     tsvs = []
#     target = directory + "/matrix_tsv"

#     for i in files:
#         name = os.path.split(i)[1].partition('.')[0]
#         print(i, name)
#         bakta_annotation(i, name[-24:])
#         annotation_tsv = os.path.split(i)[0] + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + ".tsv"
#         annotation(annotation_tsv)
#         annotation_tsv = os.path.split(i)[0] + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + "_extended.tsv"
#         tsvs.append(annotation_tsv)

#     create_directory_with_soft_links(tsvs, target)
#     print("Start pangenome")
#     pangenome_analysis(target)
#     print("Start creating faa")
#     create_fasta_file(tsvs, directory)
#     print("Start stat creation")
#     make_stat_file(target)

#     asyncio.run(send_smth(text=["Ends annotation"]))

# def pipeline_assembly(directory):
#     """
#     Fortunately, .fasta files could be in one directory.
#     """
#     files = []
#     for filename in os.listdir(directory):
#         if filename.endswith('.fasta') or filename.endswith('.fa'):
#             file_path = os.path.join(directory, filename)
#             files.append(file_path)

#     tsvs = []
#     target = directory + "/matrix_tsv"

#     for i in files:
#         name = os.path.split(i)[1].partition('.')[0]
#         print(i, name)
#         bakta_annotation(i, name[-24:])
#         annotation_tsv = os.path.split(i)[0] + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + ".tsv"
#         annotation(annotation_tsv)
#         annotation_tsv = os.path.split(i)[0] + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + "_extended.tsv"
#         tsvs.append(annotation_tsv)

#     create_directory_with_soft_links(tsvs, target)
#     print("Start pangenome")
#     pangenome_analysis(target)
#     print("Start creating faa")
#     create_fasta_file(tsvs, directory)
#     print("Start stat creation")
#     make_stat_file(target)
#     asyncio.run(send_smth(text=["Ends annotation"]))


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
    target = directory + "/matrix_tsv"

    for i in files:
        name = os.path.split(i)[1].partition('.')[0]
        print(i, name)
        bakta_annotation(i, name[-24:])
        annotation_tsv = os.path.split(i)[0] + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + ".tsv"
        annotation(annotation_tsv)
        annotation_tsv = os.path.split(i)[0] + "/bakta_annotation_" + name[-24:] + "/" + name[-24:] + "_extended.tsv"
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
        if os.path.exists(os.path.split(i)[0] + "/bakta_annotation"):
            os.rename(os.path.split(i)[0] + "/bakta_annotation", os.path.split(i)[0] + "/bakta_annotation_" + name[-24:])
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


# polishing_annotation("/storage/data1/marmi/assembly_project/ecoli_crohn_isolates_annotation")

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


# pipeline_assembly("/storage/data1/marmi/assembly_project/ecoli_crohn_isolates_annotation/assemblies/annotation_of_plasmid")
