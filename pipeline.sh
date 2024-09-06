#!/bin/bash

# аннотирует фасту с белковыми последовательностями
upimapi -i /storage/data1/marmi/annotation_project/phastest/zvl_predicted_genes.faa -o /storage/data1/marmi/annotation_project/phastest/zvl_predicted_genes_upimapi -db uniprot -t 7

# ищем интегразы
isescan.py --seqfile genomes/ --output --nthread 7


# кластеризация
mmseqs easy-cluster examples/DB.fasta clusterRes tmp --min-seq-id 0.5 -c 0.8 --cov-mode 1
# вытаскивание представителей 
mmseqs createsubdb DB_clu DB DB_clu_rep
mmseqs convert2fasta DB_clu_rep DB_clu_rep.fasta
