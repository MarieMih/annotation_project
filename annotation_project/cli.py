#!/storage/data1/marmi/miniconda/installation/envs/test_env_annotation/bin/python
import argparse
import sys
import os
sys.path.append(os.path.dirname(__file__))
import common_variables
from main import pipeline_assembly
from main import pipeline_since_fastq
from main import polishing_annotation
from main import pipeline_assembly_file
from main import pipeline_stat_all_tsv_in_dir
from main import pipeline_assembly_bakta_only
from main import pipeline_setting
from main import pipeline_headers
from resource_monitor import ResourceMonitor
from helpers import return_str_with_date_and_time


def launch_pipeline_assembly(args):
    data_line = return_str_with_date_and_time()
    report_path = os.path.join(args.directory, "monitor_report_" + data_line + ".txt")
    with ResourceMonitor(interval=args.interval, report_path=report_path, label="assembly_annotation") as monitor:
        pipeline_assembly(args.directory, args.jobs)


def launch_pipeline_fastq(args):
    data_line = return_str_with_date_and_time()
    report_path = os.path.join(args.directory, "monitor_report_" + data_line + ".txt")
    with ResourceMonitor(interval=args.interval, report_path=report_path, label="fastq_annotation") as monitor:
        pipeline_since_fastq(args.directory, args.jobs)


def launch_annotation_polishing(args):
    data_line = return_str_with_date_and_time()
    report_path = os.path.join(args.directory, "monitor_report_" + data_line + ".txt")
    with ResourceMonitor(interval=args.interval, report_path=report_path, label="polishing_annotation") as monitor:
        polishing_annotation(args.directory)


def launch_pipeline_assembly_file(args):
    pipeline_assembly_file(args.file)

def launch_pipeline_headers(args):
    pipeline_headers(args.input, args.output)

def launch_pipeline_stat_all_tsv_in_dir(args):
    pipeline_stat_all_tsv_in_dir(args.directory)


def launch_pipeline_assembly_bakta_only(args):
    data_line = return_str_with_date_and_time()
    report_path = os.path.join(args.directory, "monitor_report_" + data_line + ".txt")
    with ResourceMonitor(interval=args.interval, report_path=report_path, label="assembly_bakta_only") as monitor:
        pipeline_assembly_bakta_only(args.directory, args.jobs)


def launch_pipeline_setting(args):
    pipeline_setting()


def create_parser():
    prog = sys.argv[0]
    if ".py" in prog:
        prog = "python {}".format(prog)
    else:
        prog = None

    parser = argparse.ArgumentParser(prog=prog, description='Launch annotation')
    subparsers = parser.add_subparsers(metavar="command", required=True)

    parser_assembly = subparsers.add_parser("assembly", help="annotate all assembly.fasta in directory")
    parser_assembly.add_argument("-d", "--directory", metavar="DIRECTORY")
    parser_assembly.add_argument("--bakta-db", metavar="BAKTA_DB")
    parser_assembly.add_argument("--user-db", metavar="USER_DB")
    parser_assembly.add_argument("--tool", metavar="CLUSTER_TOOL")
    parser_assembly.add_argument("--send-tg", help="telegram messages", action="store_true")
    parser_assembly.add_argument("-t", "--threads", metavar="THREADS")
    parser_assembly.add_argument("--jobs", type=int, metavar="JOBS", default=1, help="Number of sample-level annotation jobs to run in parallel")
    parser_assembly.add_argument("--interval", type=float, metavar="INTERVAL", default=60.0, help="A snapshot for the monitor_report will be taken at each interval (in sec). Default = 60 sec")
    parser_assembly.set_defaults(func=launch_pipeline_assembly)

    parser_assembly_file = subparsers.add_parser("assembly_file", help="annotate all assembly.fasta written in .txt")
    parser_assembly_file.add_argument("-f", "--file", metavar="TXT_FILE")
    parser_assembly_file.add_argument("--bakta-db", metavar="BAKTA_DB")
    parser_assembly_file.add_argument("--user-db", metavar="USER_DB")
    parser_assembly_file.add_argument("--tool", metavar="CLUSTER_TOOL")
    parser_assembly_file.add_argument("--jobs", type=int, metavar="JOBS", default=1, help="Number of sample-level annotation jobs to run in parallel")
    parser_assembly_file.add_argument("--send-tg", help="telegram messages", action="store_true")
    parser_assembly_file.add_argument("-t", "--threads", metavar="THREADS")
    parser_assembly_file.add_argument("--interval", type=float, metavar="INTERVAL", default=60.0, help="A snapshot for the monitor_report will be taken at each interval (in sec). Default = 60 sec")
    parser_assembly_file.set_defaults(func=launch_pipeline_assembly_file)

    parser_fastq = subparsers.add_parser("fastq", help="filter, assembly and annotate all fastqs in directory")
    parser_fastq.add_argument("-d", "--directory", metavar="DIRECTORY")
    parser_fastq.add_argument("--bakta-db", metavar="BAKTA_DB")
    parser_fastq.add_argument("--user-db", metavar="USER_DB")
    parser_fastq.add_argument("--tool", metavar="CLUSTER_TOOL")
    parser_fastq.add_argument("--send-tg", help="telegram messages", action="store_true")
    parser_fastq.add_argument("-t", "--threads", metavar="THREADS")
    parser_fastq.add_argument("--jobs", type=int, metavar="JOBS", default=1, help="Number of sample-level annotation jobs to run in parallel")
    parser_fastq.add_argument("--interval", type=float, metavar="INTERVAL", default=60.0, help="A snapshot for the monitor_report will be taken at each interval (in sec). Default = 60 sec")
    parser_fastq.set_defaults(func=launch_pipeline_fastq)

    parser_polish = subparsers.add_parser("polish", help="polish bakta annotation in directory")
    parser_polish.add_argument("-d", "--directory", metavar="DIRECTORY")
    parser_polish.add_argument("--tool", metavar="CLUSTER_TOOL")
    parser_polish.add_argument("--send-tg", help="telegram messages", action="store_true")
    parser_polish.add_argument("-t", "--threads", metavar="THREADS")
    parser_polish.add_argument("--interval", type=float, metavar="INTERVAL", default=60.0, help="A snapshot for the monitor_report will be taken at each interval (in sec). Default = 60 sec")
    parser_polish.set_defaults(func=launch_annotation_polishing)

    parser_header = subparsers.add_parser("header", help="change header in faa for compatibility")
    parser_header.add_argument("-i", "--input", metavar="INPUT")
    parser_header.add_argument("-o", "--output", metavar="OUTPUT")
    parser_header.set_defaults(func=launch_pipeline_headers)

    parser_stat = subparsers.add_parser("stat", help="stat all tsv files in directory")
    parser_stat.add_argument("-d", "--directory", metavar="DIRECTORY")
    parser_stat.set_defaults(func=launch_pipeline_stat_all_tsv_in_dir)

    parser_stat_bakta = subparsers.add_parser("bakta", help="annotate all .fasta in directory by only custom bakta")
    parser_stat_bakta.add_argument("-d", "--directory", metavar="DIRECTORY")
    parser_stat_bakta.add_argument("--bakta-db", metavar="BAKTA_DB")
    parser_stat_bakta.add_argument("--user-db", metavar="USER_DB")
    parser_stat_bakta.add_argument("-t", "--threads", metavar="THREADS")
    parser_stat_bakta.add_argument("--jobs", type=int, metavar="JOBS", default=1, help="Number of sample-level annotation jobs to run in parallel")
    parser_stat_bakta.add_argument("--interval", type=float, metavar="INTERVAL", default=60.0, help="A snapshot for the monitor_report will be taken at each interval (in sec). Default = 60 sec")
    parser_stat_bakta.set_defaults(func=launch_pipeline_assembly_bakta_only)

    parser_stat = subparsers.add_parser("setting", help="help to set up all databases")
    parser_stat.set_defaults(func=launch_pipeline_setting)

    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()
    if "func" in args:
        common_variables.BAKTA_DB = args.bakta_db if (hasattr(args, "bakta_db") and (args.bakta_db != "") and (args.bakta_db is not None)) else None
        common_variables.PROTEINS = args.user_db if (hasattr(args, "user_db") and (args.user_db != "") and (args.user_db is not None)) else None
        common_variables.TOOL = args.tool if (hasattr(args, "tool") and (args.tool != "") and (args.tool is not None)) else "MMSEQS2"
        common_variables.N_THREADS = args.threads if (hasattr(args, "threads") and (args.threads is not None)) else str(int(os.cpu_count() * 0.75 if (os.cpu_count() is not None) else 1))
        common_variables.SEND_NOTIFICATION = args.send_tg if hasattr(args, "send_tg") else False
        args.func(args)
    # else:
    #     parser.print_help()


if __name__ == '__main__':
    main()
