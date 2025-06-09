![Colinca](https://github.com/MarieMih/annotation_project/blob/rebase/other/colinca.png)



# Colinca

Tool that creates colinear annotations for group of samples of one species.

## Contents
- [Installation](#installation)
- [Usage](#usage)
- [Input](#input)
- [Options](#options)
- [Output](#output)
- [Example data](#example-data)
- [Examples](#examples)
- [FAQ](#faq)
- [Acknowledgements](#acknowledgements)
- [Contact](#contact)
- [License](#license)

## Installation

    git clone -b rebase https://github.com/MarieMih/annotation_project.git
    cd annotation_project
    conda env create -f environment.yml
    conda activate annotation_project

## Usage

Before launching it's needed to download databases for Bakta, UPIMAPI and desired taxids.

Please, run firstly:

    python cli.py setting

While interacting with command line, you could choose what type of Bakta database use (full/light), for what taxons build protein-trusted list and configure telegram-bot for getting messages.


After that, you can run pipeline:

    python cli.py assembly -d </path/to/folder>

## Input

If you want to start from PE reads, use this command:

    python cli.py fastq -d </path/to/folder-with-fastq>
  
If you want to start with assemblies, use these commands:

    python cli.py assembly -d </path/to/folder-with-fasta>

    #or 

    python cli.py assembly_file -f </path/to/file-with-paths-to-assembly>

If you want to start with annotations, use these commands:

    python cli.py polish -d </path/to/file-with-paths-to-annotation>

If you need to add samples to you current group, use this:

    python cli.py bakta -d </path/to/folder-with-fasta>

and after run polish on new group.

#### News
Latest version is 0.91.

## Options
Command line view:
```
usage: python cli.py [-h] command ...

Launch annotation

positional arguments:
  command
    assembly     annotate all assembly.fasta in directory
    assembly_file
                 annotate all assembly.fasta written in .txt
    fastq        filter, assembly and annotate all fastqs in directory
    polish       polish bakta annotation in directory
    stat         stat all tsv files in directory
    bakta        annotate all .fasta in directory by only custom bakta
    setting      help to set up all databases

options:
  -h, --help     show this help message and exit
```

#### The -d parameter [string]
Common parameter for set directory with files.
#### The -t parameter [integer]
Common parameter for set number of threads.
#### The --send-tg parameter [boolean]
Common parameter that switch on telegram messages.
#### The --user-db parameter [string]
Parameter that set path to protein-trusted list file.
#### The --bakta-db parameter [string]
Parameter that set path to bakta database folder.
#### The --upimapi-db parameter [string]
Parameter that set path to UPIMAPI database.


## Output
Principal file is __extended.tsv_.

The results file contains the following columns:

| No. | Column             | Description |
|-----|--------------------|-------------|
| 1   | **Sequence Id**    | Identifier of the sequence, such as contig or chromosome name (e.g., `contig_1`). |
| 2   | **Type**           | Type of annotated feature: `cds` stands for coding sequence. |
| 3   | **Start**          | Start position of the feature on the sequence (1-based coordinate system). |
| 4   | **Stop**           | End position of the feature on the sequence. |
| 5   | **Strand**         | DNA strand the gene is located on: `+` for forward, `-` for reverse. |
| 6   | **Locus Tag**      | Unique gene identifier within the annotation set. |
| 7   | **Gene**           | Standard or provisional gene symbol, if available (e.g., `rspA`). |
| 8   | **Product**        | Name or description of the protein encoded by this gene (e.g., `Starvation-sensing protein RspA`). |
| 9   | **DbXrefs**        | External databases and their identifiers referring to this gene or protein (e.g., `RefSeq`, `UniParc`, `UniRef`, `SO:`). |
| 10  | **Organism**       | Name of the organism where this gene was identified (e.g., `Escherichia coli`). |
| 11  | **Entry UniProtKB**| Unique protein ID from the UniProtKB database (e.g., `P38104`). |
| 12  | **GO**             | Functional annotations from the Gene Ontology database: molecular function, biological process, and cellular component (e.g., `DNA binding [GO:0003677]`). |
| 13  | **KEGG**           | KEGG gene and pathway identifiers (e.g., `ecj:JW1573`). |
| 14  | **UniPathway**     | Metabolic pathway identifiers from UniPathway (e.g., `UPA00211;UPA00250`). |
| 15  | **Pathway**        | Descriptive names of metabolic or biological pathways involving the protein (e.g., `spermidine degradation`). |
| 16  | **Keywords**       | Keywords describing protein properties: enzymatic activity, localization, reference proteome status, etc. |
| 17  | **Transcript_id**  | Unique transcript identifier in the format `cds\|GENE\|UniProtID`. |
| 18  | **Gene_id**        | Unique gene identifier, usually matching `Transcript_id`. |

## Example data
TODO

## Examples
TODO

#### 1. Case of assemblies
TODO

> [!NOTE]
>  If you want to repolish annotations, run in directory with all annotation folders:
>
>     rm -r ./bakta_annotation_*/userprotein* ./bakta_annotation_*/*extended.* ./bakta_annotation_*/*_uniref100* ./bakta_annotation_*/*upimapi_ref2ref ./bakta_annotation_*/*cds_sorf.tsv ./bakta_annotation_*/*_detected.faa  ./bakta_annotation_*/*_semidefined.tsv ./bakta_annotation_*/*_rna.tsv ./bakta_annotation_*/*_unknown.faa union_detected_faa union_unknown_faa upimapi_output matrix_tsv

> [!WARNING]
> Please, understand that this _**is not fully tested**_ .


## FAQ
- **Why use this tool instead of Bakta?**

Particularly, this is most suitable if you want to compare a lot of similar samples.

- **Why Colinca?**

Colinca is an anagram of "Colling" and "annotation". It's also famous dance and tree.

## Acknowledgements
- Work is done as a magister diplome in MIPT.

## Contact
Maria Mikhailycheva (mikhailycheva.mv@phystech.edu)

## License
MIT license.
