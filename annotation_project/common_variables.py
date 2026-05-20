import os

QUERY_ALIGN_LENGTH_UP     = 1.05
QUERY_ALIGN_LENGTH_BOTTOM = 0.95
PIDENT                    = 90
RATIO_LEN_UP              = 1.1
RATIO_LEN_BOTTOM          = 0.9

BAKTA_DB = None
PROTEINS = None
N_THREADS = str(int(os.cpu_count() * 0.75 if (os.cpu_count() is not None) else 1))

TOOL = "MMSEQS2"

TOOL_FOR_RNA = "None"

SEND_NOTIFICATION = None

BAKTA_TSV_HEADER = "Sequence Id,Type,Start,Stop,Strand,Locus Tag,Gene,Product,DbXrefs"
BAKTA_INFERENCE_TSV_HEADER = "Sequence Id,Type,Start,Stop,Strand,Locus Tag,Score,Evalue,Query Cov,Subject Cov,Id,Accession"
PANGENOME_TSV_HEADER = "Sequence Id,Type,Start,Stop,Strand,Gene,Gene Name,Gene synonymes,Product,DbXrefs,Organism,Inference,KEGG,GO,gene_id,transcript_id"