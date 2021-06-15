# LGAAP: *Leishmaniinae* Genome Assembly and Annotation Pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4707445.svg)](https://doi.org/10.5281/zenodo.4707445)

[Hatim Almutairi](mailto:hatim.almutiairi@hotmail.com) 2021

This repository was made for the purpose of reproducing *de novo* Genome assembly and annotation for the following genomes:
 - [*Leishmania (Mundinia) martiniquensis* strain: LV760 (isolate:LSCM1)](https://www.ncbi.nlm.nih.gov/assembly/GCA_017916325.1)
 - [*Leishmania (Mundinia) orientalis* strain: LV768 (isolate:LSCM4)](https://www.ncbi.nlm.nih.gov/assembly/GCA_017916335.1)
 - [*Leishmania (Mundinia) enriettii* strain: LV763 (isolate:CUR178)](https://www.ncbi.nlm.nih.gov/assembly/GCA_017916305.1)
 - [*Leishmania (Mundinia) sp.* Ghana strain: LV757 (isolate:GH5)](https://www.ncbi.nlm.nih.gov/assembly/GCA_017918215.1)
 - [*Leishmania (Mundinia) sp.* Namibia strain: LV425 (isolate:253)](https://www.ncbi.nlm.nih.gov/assembly/GCA_017918225.1)
 - [*Porcisia hertigi* strain: LV43 (isolate:C119)](https://www.ncbi.nlm.nih.gov/assembly/GCA_017918235.1)


## Pipeline Outline
![Pipeline Outline](https://github.com/hatimalmutairi/LGAAP/blob/main/Pipline_Outline.png)

## Content
This repository contains An automated pipline for the assembly and annotation of six *Leishmaniinae* genomes and was writen and excuted using [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html). the pipeline run [314 steps](https://github.com/hatimalmutairi/LGAAP/blob/main/rulegraph.svg) to complete. 

### Friendly warning
Running the pipeline requires around 500GB of free space including 140GB for SRAs.

## How to use
Before you start, make sure that [docker software](https://docs.docker.com/get-docker/) is already installed in your machine. 
 1. Use  ```git clone``` command to clone the repository to your working directory. you can also download it from  [zenodo](https://zenodo.org/record/4707445)
```sh
$ git clone https://github.com/hatimalmutairi/LGAAP.git
```
 2. Create a [conda](https://docs.conda.io/en/latest/) environment with [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) installed.
 Anaconda or miniconda must be installed already in youe machine.
```sh
$ conda create -n LGAAP -c bioconda -c conda-forge snakemake -y
```
 3. Activate the new environment
```
$ conda activate LGAAP
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
```results/annotation/{sample}/Final_Annotation/ ``` contains the final assemblies with annotated proteins and transcripts.

## Estimated Time of run
We attempted to run the pipeline on a machine with 24 CPUs and 120 GB of RAM. The entire operation took approximately 15377 minutes (around 10 days). The round of evidence-based annotation consumed the majority of the running time, with each genome taking approximately 30 hours to complete.


![Runtime](https://github.com/hatimalmutairi/LGAAP/blob/main/Runtime.png)

## Pipeline blueprint
Each genome assembly was created in a sperate exctuion file (all of them can be seen in ```workflow/rules``` directory) 

Click [here](https://github.com/hatimalmutairi/LGAAP/blob/main/rulegraph.svg) to see the rulegraph

Take note that the following command lines were included to demonstrate the main commands; however, attempting to use them directly may fail because the majority of the commands were configured to be executed inside either a conda environment or a docker container. Rather than that, it is recommended to run the snakemake piplene commands, which will take care of all software and its dependencies.

### Assembly using Flye v2.8.2
```sh
flye --nano-raw {input.long} --genome-size 35m --threads {threads} -o results/LSCM1_Assembly_LongRead_by_Flye
```
### Polishing using minimap2 v2.17
```sh
minimap2 -t {threads} -ax sr {input.assembly} {input.sam}  > {output.sam}
```
### Polishing using Samtools v1.11
```sh
samtools faidx {input} > {output}
samtools view -@ {threads} -bt {input.idx} {input.sam} > {output.bam}
samtools sort -@ {threads} {input.bam} -o {output.bam_sort}
samtools index -@ {threads} {input.bam_sort} > {output.bam_idx}
```
### Polishing using BCFtools v1.11
```sh
bcftools mpileup --threads {threads} -Ou -d 1000000 -f {input.assembly} {input.bam_sort}
bcftools call -mv -Oz -o {output.PE}
bcftools norm --threads {threads} -f {input.assembly} {input.PE} -Ob -o {output}
bcftools filter --threads {threads} --IndelGap 5 {input.PE} -Ob -o {output}
bcftools index --threads {threads} {input} -o {output}
bcftools consensus -f {input.assembly} {input.PE} > {output}
```
### Polishing using Pilon v1.23
```sh
pilon --threads {threads} --genome {input.assembly} --bam {input.bam_sort} --outdir {output}
```
### Arrange Assembly using RaGOO v1.1
```sh
ragoo.py -t {threads} -C -T sr -b -g 0 -R {input.reads} {input.assembly} {input.reference}
```
### Clean and sort the Polished Assembly using  funannotate v1.7.1
```sh
funannotate clean --exhaustive --minlen 50 --pident 10 -i {input} -o {output}
```
### Scan contaminants using Blast+ v2.10.1
```sh
blastn -db {input.UniVec} -max_hsps 1 -max_target_seqs 1 -outfmt '6 qseqid qstart qend' -num_threads {threads} -query {input.assembly} -out {output}
```
### Mask contaminants using bedtools v2.30.0
```sh
bedtools maskfasta -fi {input.fa1} -bed {input.bed1} -fo {output}
```
### QC using GAAS v1.2.0
```sh
gaas_fasta_statistics.pl -f {input} -o {output}
```
### Detecting Repeats using RepeatModeler
```sh
BuildDatabase -name {sample}_RepeatDB -engine ncbi {input}
RepeatModeler -pa {threads} -engine ncbi -database {wildcards.sample}_RepeatDB
```
### Mask repeats using RepeatMasker v1.4
```sh
RepeatMasker -a -gff -x -pa {threads} -e ncbi -lib {input.Repeat} {input.assembly}
```
### Classify transpsosable element (TE) by TEclass v2.1.3b
```sh
TEclassTest.pl -r {input} > {output}
```
### Annotation using MAKER v2.31.10
```sh
maker -fix_nucleotides -genome {input.g1} -base {wildcards.sample}_Round1 maker_opts.ctl maker_bopts.ctl maker_exe.ctl
gff3_merge -n -s -d {input} > {output}
fasta_merge -d {input} > {output}
maker_map_ids --prefix {wildcards.sample}_ --justify 5 {input} > {output}
maker_functional_fasta {input.Uniprot} {input} {input.gff} > {output.gff}
maker_functional_fasta {input.Uniprot} {input} {input.proteins} > {output.proteins}
maker_functional_fasta {input.Uniprot} {input} {input.transcripts} > {output.transcripts}
```
### Annotation QC using Genometools v1.2.1
```sh
gt gff3 -sort -tidy {input} | gt stat > {output}
```
### processing annotation using GAAS v1.2.0
```sh
gaas_maker_merge_outputs_from_datastore.pl -i {input} -o {output}
```
### Assign Functional Annotations using interproscan v5.22-61.0 and Blast+ v2.10.1
```sh
blastp -num_threads {threads} -query {input} -db {input} -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out {output}
interproscan.sh -appl pfam -dp -f TSV -goterms -iprlookup -pa -t p -i {input} -o {output}
```
### Keep Longest Isoform using AGAT v0.6.0
```sh
agat_sp_keep_longest_isoform.pl -gff {input} -o {output}
```
### Create the Final Annotation file using MAKER v2.31.10
```sh
ipr_update_gff {input.gff} {input.iprscan} > {output}
```
