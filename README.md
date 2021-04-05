# Leishmania Mundinia Genomes Assembly and Annotation Pipeline (LMGAAP)

[Hatim Almutairi](mailto:hatim.almutiairi@hotmail.com) 2021

This repository was made for the purpose of reproducing de novo Genome assembly and annotation for the following genomes:
 - [Leishmania (Mundinia) martiniquensis strain: LV760 (isolate:LSCM1)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA691531)
 - [Leishmania (Mundinia) orientalis strain: LV768 (isolate:LSCM4)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA691532)
 - [Leishmania (Mundinia) enrietti strain: LV763 (isolate:CUR178)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA691534)
 - [Leishmania sp. Ghana strain: LV757 (isolate:GH5)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA691536)
 - [Leishmania (Mundinia) sp. Namibia strain: LV425 (isolate:253)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA689706)
 - [Porcisia hertigi strain: LV43 (isolate:C119)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA691541)

## References:
- [LMGAAP on Github](https://github.com/hatimalmutairi/LMGAAP.git)

## Pipeline Outline
![Pipeline Outline](https://github.com/hatimalmutairi/LMGAAP/blob/main/Pipline_Outline.png)

## Content
This repository contains An automated pipline for the assembly and annotation of six Leishmania genomes and was writen and excuted using [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html).

### Friendly warning
Running the pipeline requires around 500GB of free space including 140GB for SRAs.

## How to use
 1. Use  ```git clone``` command to clone the repository to your working directory. you can also download it from [Lancaster University Research Directory]() or from [here]()
```sh
$ git clone https://github.com/hatimalmutairi/LMGAAP.git
```
 2. Create a [conda](https://docs.conda.io/en/latest/) environment with [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) installed.
 Anaconda or miniconda must be installed already in youe machine.
```sh
$ conda create -n LMGAAP -c bioconda -c conda-forge snakemake -y
```
 3. Activate the new environment
```
$ conda activate LMGAAP
```
 4. Set the downloaded directory ad the working directory and  run the entire pipline
```sh
$ snakemake --cores 8 --use-conda
```
The option ```--cores``` lets you choose how many CPUs to use.
 
 5. (Optional) Perform a dry run of Snakemake using the option ```--dry-run```
```sh
$ snakemake --dry-run
```
 6. (Optional) Create a rulegraph wth the following command
```sh
$ snakemake --rulegraph | dot -Tpng > rulegraph.png 
```
For further instrucion how to use snakemake, visit [Sankemake](https://snakemake.readthedocs.io/en/stable/index.html) documentation

## Where to find the final output after runnning the pipeline
 - The final assemblies with annotated proteins and transcripts will be in ```results/annotation/{sample}/Final_Annotation/```
