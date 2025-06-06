import sys
import subprocess
import os

QUERY_ALIGN_LENGTH_UP = 1.05
QUERY_ALIGN_LENGTH_BOTTOM = 0.95
PIDENT = 90
RATIO_LEN_UP = 1.1
RATIO_LEN_BOTTOM = 0.9

UPIMAPI_RESOURCES = None
UPIMAPI_DATABASE = "uniprot"

BAKTA_DB = None
PROTEINS = None
N_THREADS = str(int(os.cpu_count() * 0.75 if (os.cpu_count() is not None) else 1))

SEND_NOTIFICATION = None

