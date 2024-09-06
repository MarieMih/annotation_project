#import Bio
import argparse
import sys, os
import subprocess
from multiprocessing import Process
import shutil

def pipe (args):
    try:
        chck=os.path.exists(args.input_dir)
        if (chck == True):
            retcode = subprocess.call(['ls', args.input_dir, args.outdir])
            if retcode.returncode < 0:
                print("Child was terminated by signal", -retcode.returncode, file=sys.stderr)
            else:
                print("Success create_dir. Child returned", retcode.returncode, file=sys.stderr)              
        else:
            print("Input directory doesn't exist.")
    except OSError as e:
        print("Execution failed:", e, file=sys.stderr)

def clustering (dir_to_clust, output_dir, tmp_dir):
    try:
        chck=os.path.exists(dir_to_clust)
        program=shutil.which("mmseqs")
        if (chck == True):
            retcode = subprocess.run([program, 'easy-cluster', dir_to_clust, output_dir, tmp_dir], shell = True)
            if retcode.returncode < 0:
                print("Child was terminated by signal", -retcode.returncode, file=sys.stderr)
            else:
                print("Success create_dir. Child returned", retcode.returncode, file=sys.stderr)              
        else:
            print("Input directory doesn't exist.")
    except OSError as e:
        print("Execution failed:", e, file=sys.stderr)

def create_parser():
    parser = argparse.ArgumentParser(
        prog="proann",
        description="This program is designed to make the most precise annotation of bacteria."
    )

    subparsers = parser.add_subparsers(metavar="command", required=True)
    
    parser_pipe = subparsers.add_parser("pipe", help="all analysis")
    parser_pipe.add_argument("-i", "--input-dir", metavar="INPUTDIR", required=True,
                        help="Input directory with all genomes")
    parser_pipe.add_argument("-o", "--outdir", metavar="OUTDIR", default=os.getcwd(),
                        help="Output directory")
    parser_pipe.set_defaults(func=pipe)

    parser_clustering = subparsers.add_parser("clustering", help="all analysis")
    parser_clustering.add_argument("-i", "--dir_to_clust", metavar="INPUTDIR", required=True,
                        help="Input directory with all genomes")
    parser_clustering.add_argument("-o", "--output_dir", metavar="OUTDIR", default=os.getcwd(),
                        help="Output directory")
    parser_clustering.set_defaults(func=clustering)

    return parser

def main():
    clustering("/storage/data1/marmi/annotation_project/uniprot_db/uniprotkb_escherichia_phage.faa",
               "/storage/data1/marmi/annotation_project/uniprot_db/test_proann_clustering/test",
               "/storage/data1/marmi/annotation_project/uniprot_db/test_proann_clustering/tmp")
    parser = create_parser()
    args = parser.parse_args()
    if "func" in args: # Check subcommand function argument
        args.func(args)
    else:
        parser.print_help()
    return 0

if __name__ == '__main__':
    main()
