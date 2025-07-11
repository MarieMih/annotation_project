"""Module for assembling and annotating fastq."""

import os
import sys
import subprocess
from telegram_send import send
sys.path.append(os.path.dirname(__file__))
import common_variables


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
        subprocess.run(['fastp', '-e', str(avg_qual), '-w', common_variables.N_THREADS,
                        '--length_required', str(length_required),
                        '--in1', fastq_1_in, '--in2', fastq_2_in,
                        '--out1', fastq_1_out, '--out2', fastq_2_out,
                        '-j', pref + "_fastp.json",
                        '-h', pref + "_fastp.html"
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
        fastq_out = os.path.join(parts[0], name + "_filtered." + extension)
        fastq_ready = os.path.join(parts[0], name + "_sub." + "fq.gz")
    else:
        print("Must be something wrong with extension.")
        sys.exit()

    pref = fastq.rpartition('.')[0]
    f = open(pref + "_fastp.log", "w", encoding="utf-8")
    try:
        subprocess.run(['fastp', '-e', str(avg_qual), '-w', common_variables.N_THREADS,
                        '--length_required', str(length_required),
                        '--in1', fastq_in, '--out1', fastq_out,
                        '-j', pref + "_fastp.json",
                        '-h', pref + "_fastp.html"
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
    pth = os.path.join(parts[0], "assembly_" + name)
    try:
        subprocess.run(['unicycler',
                        '-1', fastq_1, '-2', fastq_2,
                        '-o', pth,
                        '--threads', common_variables.N_THREADS
                        ], check=True)
    except Exception as e:
        print(e, "\n", "Error with unicycler.")
        sys.exit()

    return os.path.join(pth, "assembly.fasta")


def bakta_annotation(fasta, locus_tag):
    """
    Annotation of bacterial assembly.
    """
    print(common_variables.BAKTA_DB)
    parts = os.path.split(fasta)
    pth = os.path.join(parts[0], "bakta_annotation_" + locus_tag)
    f = open(os.path.join(parts[0], "bakta_annotation_" + locus_tag + ".log"), "w", encoding="utf-8")
    try:
        subprocess.run(['bakta',
                        '--db', common_variables.BAKTA_DB,
                        '--proteins', common_variables.PROTEINS,
                        '--force',
                        '--skip-plot',
                        '--locus-tag', locus_tag,
                        '--prefix', locus_tag,
                        '--threads', common_variables.N_THREADS,
                        '--output', pth,
                        fasta
                        ], check=True, stdout=f)
    except Exception as e:
        print(e, "\n", "Error with bakta.")
        f.close()
        # sys.exit()
    f.close()
