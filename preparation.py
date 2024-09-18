"""Module for assembling and annotating fastq."""

import os
import sys
import subprocess
import asyncio
from telegram_send import send

BAKTA_DB = "/storage/data1/marmi/bakta_db/db"
PROTEINS = ""


async def send_smth(text=None, image=None, file=None):
    """
    Function for send to telegram.
    """
    if text:
        await send(messages=text)

    if image:
        for j in image:
            with open(j, "rb") as f:
                await send(images=[f])

    if file:
        for j in file:
            with open(j, "rb") as f:
                await send(files=[f])


def filtering_fastq_pe(fastq,
                       genome_size=460000000,
                       avg_qual=28,
                       length_required=50):
    """
    Filter NGS-data with fastp and subsample with rasusa.
    """
    if fastq.endswith("_1.fq.gz") or fastq.endswith("_2.fq.gz"):
        fastq_1_in = fastq[:-8] + "_1.fq.gz"
        fastq_2_in = fastq[:-8] + "_2.fq.gz"
        fastq_1_out = fastq[:-8] + "_1_filtered.fq.gz"
        fastq_2_out = fastq[:-8] + "_2_filtered.fq.gz"
        fastq_1_ready = fastq[:-8] + "_1_sub.fq.gz"
        fastq_2_ready = fastq[:-8] + "_2_sub.fq.gz"
    elif fastq.endswith("_1.fq") or fastq.endswith("_2.fq"):
        fastq_1_in = fastq[:-5] + "_1.fq"
        fastq_2_in = fastq[:-5] + "_2.fq"
        fastq_1_out = fastq[:-5] + "_1_filtered.fq"
        fastq_2_out = fastq[:-5] + "_2_filtered.fq"
        fastq_1_ready = fastq[:-5] + "_1_sub.fq.gz"
        fastq_2_ready = fastq[:-5] + "_2_sub.fq.gz"
    elif fastq.endswith("_1.fastq") or fastq.endswith("_2.fastq"):
        fastq_1_in = fastq[:-8] + "_1.fq"
        fastq_2_in = fastq[:-8] + "_2.fq"
        fastq_1_out = fastq[:-8] + "_1_filtered.fq"
        fastq_2_out = fastq[:-8] + "_2_filtered.fq"
        fastq_1_ready = fastq[:-8] + "_1_sub.fq.gz"
        fastq_2_ready = fastq[:-8] + "_2_sub.fq.gz"
    elif fastq.endswith("_1.fastq.gz") or fastq.endswith("_2.fastq.gz"):
        fastq_1_in = fastq[:-11] + "_1.fq.gz"
        fastq_2_in = fastq[:-11] + "_2.fq.gz"
        fastq_1_out = fastq[:-11] + "_1_filtered.fq"
        fastq_2_out = fastq[:-11] + "_2_filtered.fq"
        fastq_1_ready = fastq[:-11] + "_1_sub.fq.gz"
        fastq_2_ready = fastq[:-11] + "_2_sub.fq.gz"
    else:
        print("Must be something wrong with extension.")
        sys.exit()

    pref = fastq.rpartition('.')[0]
    f = open(pref + "_fastp.log", "w", encoding="utf-8")
    try:
        subprocess.run(['fastp', '-e', str(avg_qual), '-w', "8",
                        '--length_required', str(length_required),
                        '--in1', fastq_1_in, '--in2', fastq_2_in,
                        '--out1', fastq_1_out, '--out2', fastq_2_out
                        ], check=True, stderr=f)
    except Exception as e:
        print(e, "\n", "Error while running fastp")
        f.close()
        sys.exit()
    f.close()

    f = open(pref + "_rasusa.log", "w", encoding="utf-8")
    try:
        subprocess.run(['rasusa', 'reads', '-O', 'g',
                        '-b', str(genome_size),
                        '-o', fastq_1_ready, '-o', fastq_2_ready,
                        fastq_1_out, fastq_2_out
                        ], check=True, stderr=f)
    except Exception as e:
        f.close()
        print(e, "\n", "Error while running rasusa")
        sys.exit()
    f.close()

    return fastq_1_ready, fastq_2_ready


def filtering_fastq_se(fastq,
                       genome_size=460000000,
                       avg_qual=28,
                       length_required=50):
    """
    Filter NGS-data with fastp and subsample with rasusa.
    """
    if fastq.endswith(".fastq")     \
       or fastq.endswith(".fq")     \
       or fastq.endswith(".fq.gz")  \
       or fastq.endswith(".fastq.gz"):
        fastq_in = fastq
        parts = os.path.split(fastq)
        name = parts[1].partition('.')[0]
        extension = parts[1].partition('.')[1]
        fastq_out = parts[0] + name + "_filtered." + extension
        fastq_ready = parts[0] + name + "_sub." + "fq.gz"
    else:
        print("Must be something wrong with extension.")
        sys.exit()

    pref = fastq.rpartition('.')[0]
    f = open(pref + "_fastp.log", "w", encoding="utf-8")
    try:
        subprocess.run(['fastp', '-e', str(avg_qual), '-w', '8',
                        '--length_required', str(length_required),
                        '--in', fastq_in, '--out', fastq_out
                        ], check=True, stderr=f)
    except Exception as e:
        print(e)
        print("Error while running fastp")
        f.close()
        sys.exit()
    f.close()

    f = open(pref + "_rasusa.log", "w", encoding="utf-8")  
    try:
        subprocess.run(['rasusa', 'reads', '-O', 'g',
                        '-b', str(genome_size),
                        '-o', fastq_ready,
                        fastq_out
                        ], check=True, stderr=f)
    except Exception as e:
        f.close()
        print(e)
        print("Error while running rasusa")
        sys.exit()
    f.close()

    return fastq_ready


def assembly_unicycler_pe(fastq_1, fastq_2):
    """
    Assembly with unicycler.
    """
    parts = os.path.split(fastq_1)
    name = parts[1].partition('.')[0]
    pth = parts[0] + "/assembly_" + name
    try:
        subprocess.run(['unicycler',
                        '-1', fastq_1, '-2', fastq_2,
                        '-o', pth,
                        '--threads', '8'
                        ], check=True)
    except Exception as e:
        print(e, "\n", "Error with unicycler.")
        sys.exit()

    return pth + "/assembly.fasta"


def bakta_annotation(fasta, locus_tag):
    """
    Annotation of bacterial assembly.
    """
    parts = os.path.split(fasta)
    pth = parts[0] + "/bakta_annotation"
    f = open(parts[0] + "/bakta_annotation" + ".log", "w", encoding="utf-8")
    try:
        subprocess.run(['bakta',
                        '--db', BAKTA_DB,
                        # '--proteins', PROTEINS,
                        '--locus-tag', locus_tag,
                        '--prefix', locus_tag,
                        '--threads', '8',
                        '--output', pth,
                        fasta
                        ], check=True, stdout=f)
    except Exception as e:
        print(e, "\n", "Error with bakta.")
        f.close()
        sys.exit()
    f.close()


if __name__ == "__main__":
    directory = sys.argv[1]
    files = []
    for filename in os.listdir(directory):
        if filename.endswith('_1.fq.gz'):
            file_path = os.path.join(directory, filename)
            files.append(file_path)

    for i in files:
        name = os.path.split(i)[1].partition('.')[0]
        print(i, name)
        read_1, read_2 = filtering_fastq_pe(i)
        assembly = assembly_unicycler_pe(read_1, read_2)
        # bakta_annotation(assembly, name)

    asyncio.run(send_smth(text=["Ends filtering"]))
