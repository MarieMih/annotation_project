import os
import shutil
import sys
from annotation_project.main import pipeline_assembly_file


TEST_FASTA_DIRECTORY = os.path.join(os.path.split(os.path.realpath(sys.argv[0]))[0], "test_fasta_report_dataset")


if __name__ == "__main__":
    cur_dir = os.path.split(os.path.realpath(sys.argv[0]))[0]
    tmp = os.path.join(cur_dir, "tmp")
    try:
        if os.path.exists(tmp):
            shutil.rmtree(tmp)
        shutil.copytree(TEST_FASTA_DIRECTORY, tmp)
    except OSError as e:
        print(e)
        sys.exit()
    TEST_FASTA_FILE = os.path.join(tmp, "long.txt")
    pipeline_assembly_file(TEST_FASTA_FILE)
