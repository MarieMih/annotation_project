import os
import shutil
import sys
from annotation_project.main import pipeline_since_fastq


TEST_FASTQ_FOLDER = os.path.join(os.path.split(os.path.realpath(sys.argv[0]))[0], "test_fastq_dataset")

if __name__ == "__main__":
    cur_dir = os.path.split(os.path.realpath(sys.argv[0]))[0]
    tmp = os.path.join(cur_dir, "tmp")
    try:
        if os.path.exists(tmp):
            shutil.rmtree(tmp)
        shutil.copytree(TEST_FASTQ_FOLDER, tmp)
    except OSError as e:
        print(e)
        sys.exit()
    pipeline_since_fastq(tmp)
