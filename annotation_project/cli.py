import argparse
import sys
from main import pipeline_assembly
from main import pipeline_since_fastq
from main import polishing_annotation


def launch_pipeline_assembly(args):
    pipeline_assembly(args.directory)


def launch_pipeline_fastq(args):
    pipeline_since_fastq(args.directory)


def launch_annotation_polishing(args):
    polishing_annotation(args.directory)


def create_parser():
    prog = sys.argv[0]
    if ".py" in prog:
        prog = "python {}".format(prog)
    else:
        prog = None
    parser = argparse.ArgumentParser(prog=prog,
                                     description='Launch annotation')

    subparsers = parser.add_subparsers(metavar="command", required=True)

    parser_assembly = subparsers.add_parser("assembly",
                                            help="annotate all assembly.fasta in directory")
    parser_assembly.add_argument("-d", "--directory", metavar="DIRECTORY")
    parser_assembly.set_defaults(func=launch_pipeline_assembly)

    parser_fastq = subparsers.add_parser("fastq",
                                         help="filter, assembly and annotate all fastqs in directory")
    parser_fastq.add_argument("-d", "--directory", metavar="DIRECTORY")
    parser_fastq.set_defaults(func=launch_pipeline_fastq)

    parser_polish = subparsers.add_parser("polish",
                                          help="polish bakta annotation in directory")
    parser_polish.add_argument("-d", "--directory", metavar="DIRECTORY")
    parser_polish.set_defaults(func=launch_annotation_polishing)

    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()
    if "func" in args:
        args.func(args)
    else:
        parser.print_help()
    return 0


if __name__ == '__main__':
    main()
