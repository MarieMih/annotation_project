import subprocess


def convert_gff_to_gtf(file):
    """
    Use AGAT.
    """
    new_gtf = file.rpartition(".")[0] + ".gtf"
    subprocess.run(["agat_convert_sp_gff2gtf.pl",
                    "--gff", file,
                    "-o", new_gtf],
                   check=True)

convert_gff_to_gtf("/storage/data1/marmi/annotation_project_dev/annotation_project/crohn/long/MG_glu_ho/bakta_annotation/MG_glu_ho_extended.gff3")
