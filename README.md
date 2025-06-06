![Colinca](https://github.com/MarieMih/annotation_project/other/colinca.png)



# Colinca

New annotation tool.

## Contents
- [Installation](#installation)
- [Usage](#usage)
- [Input](#input)
- [Output](#output)
- [Options](#options)
- [Example data](#example-data)
- [Examples](#examples)
- [FAQ](#faq)
- [Acknowledgements](#acknowledgements)
- [How to cite](#how-to-cite)
- [Contact](#contact)
- [License](#license)

## Installation

    conda env create -f environment.yml
    git clone https://github.com/MarieMih/annotation_project.git
    python cli.py --help

## Usage

    python cli.py assembly -d </path/to/folder>

## Input
Different input could be used.

#### Some topic
Latest version is 0.91.


## Output
Principal file is __extended.tsv_.

The results file contains the following columns:

| Column name | Explanation |
| ----------- | ----------- |
| Gene | The gene name |
| Annotation | Annotation |

## Options
Command line view:
```
usage: python ../annotation_project/cli.py [-h] command ...

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

#### The -d parameter
TODO

## Example data
TODO

## Examples
TODO

#### 1. Case of assemblies
TODO

> [!NOTE]
>  Use `--help` option.  

> [!WARNING]
> Please, understand that this _**is not fully tested**_ .


## FAQ
- **Why use this tool instead of Bakta?**

Particularly, this is most suitable if you want to compare a lot of similar samples.

- **Why Colinca?**

Colinca is an anagram of "Colling" and "annotation". It's also famous dance and tree.

## Acknowledgements
- Work is done as a magister diplome in MIPT.

## How to cite
If you use Colinca, please consider citing:
**Cool authors** Great annotation. _Nature_. 2060;17:238

The article gonna be Open Access some day.

## Contact
Maria Mikhailycheva (mikhailycheva.mv@phystech.edu)

## License
MIT license.
