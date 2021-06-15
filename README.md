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
 6. (Optional) Create a rulegraph wth the following command similar to [this](https://github.com/hatimalmutairi/LGAAP/blob/main/rulegraph.svg)
```sh
$ snakemake --rulegraph | dot -Tpng > rulegraph.png 
```
For further instrucion how to use snakemake, visit [Sankemake](https://snakemake.readthedocs.io/en/stable/index.html) documentation

## Where to find the final output after runnning the pipeline
```results/annotation/{sample}/Final_Annotation/ ``` contains the final assemblies with annotated proteins and transcripts.

## Estimated Time of run
We used a virtual machine with 24 CPUs and 120 GB of RAM to run this pipeline. It took [15377 minutes to complete the job](https://github.com/hatimalmutairi/LGAAP/blob/main/Runtime.png) (around 10 days). The evidence-based annotation took up the majority of the time, taking about 30 hours to finish for each genome.

## Pipeline blueprint
Each genome assembly was created in a sperate exctuion file (all of them can be seen in ```workflow/rules``` directory) 

Take note that the following command lines were included to demonstrate the main commands; however, attempting to use them directly may fail because the majority of the commands were configured to be executed inside either a conda environment or a docker container. Rather than that, it is recommended to run the snakemake piplene commands, which will take care of all software and its dependencies.

1. Assembly using [Flye](https://github.com/fenderglass/Flye) (v2.8.2)
```sh
flye --nano-raw {input.longreads} --genome-size 35m --threads {threads} -o {output}
```
2. Mapping short reads using [minimap2](https://github.com/lh3/minimap2#map-long-genomic) (v2.17)
```sh
minimap2 -t {threads} -ax sr {input.assembly} {input.sam}  > {output.sam}
```
3. Viewing, sorting, and indexing the mapped reads using [Samtools](https://github.com/samtools/samtools) (v1.11)
```sh
samtools faidx {input} > {output}
samtools view -@ {threads} -bt {input.idx} {input.sam} > {output.bam}
samtools sort -@ {threads} {input.bam} -o {output.bam_sort}
samtools index -@ {threads} {input.bam_sort} > {output.bam_idx}
```
4. Creating consensus sequncing using [BCFtools](https://github.com/samtools/bcftools) (v1.11)
```sh
bcftools mpileup --threads {threads} -Ou -d 1000000 -f {input.assembly} {input.bam_sort}
bcftools call -mv -Oz -o {output.PE}
bcftools norm --threads {threads} -f {input.assembly} {input.PE} -Ob -o {output}
bcftools filter --threads {threads} --IndelGap 5 {input.PE} -Ob -o {output}
bcftools index --threads {threads} {input} -o {output}
bcftools consensus -f {input.assembly} {input.PE} > {output}
```
5. Polishing the assemblies using [Pilon](https://github.com/broadinstitute/pilon) (v1.23)
```sh
pilon --threads {threads} --genome {input.assembly} --bam {input.bam_sort} --outdir {output}
```
6. Arranging sssemblies using [RaGOO](https://github.com/malonge/RaGOO) (v1.1)
```sh
ragoo.py -t {threads} -C -T sr -b -g 0 -R {input.reads} {input.assembly} {input.reference}
```
7. Cleaning the Polished Assembly using [funannotate](https://github.com/nextgenusfs/funannotate) (v1.7.1)
```sh
funannotate clean --exhaustive --minlen 50 --pident 10 -i {input} -o {output}
```
8. Scan contaminants using [Blast+](https://github.com/ncbi/blast_plus_docs)(v2.10.1) usnig [UniVec database (build 10.0)](https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/) 
```sh
blastn -db {input.UniVec} -max_hsps 1 -max_target_seqs 1 -outfmt '6 qseqid qstart qend' -num_threads {threads} -query {input.assembly} -out {output}
```
9. Mask contaminants using [bedtools](https://github.com/arq5x/bedtools2) (v2.30.0)
```sh
bedtools maskfasta -fi {input.fa1} -bed {input.bed1} -fo {output}
```
10. QC using [GAAS](https://github.com/NBISweden/GAAS) (v1.2.0)
```sh
gaas_fasta_statistics.pl -f {input} -o {output}
```
11. Detecting Repeats using [RepeatModeler](https://github.com/Dfam-consortium/TETools)
```sh
BuildDatabase -name {sample}_RepeatDB -engine ncbi {input}
RepeatModeler -pa {threads} -engine ncbi -database {wildcards.sample}_RepeatDB
```
12. Mask repeats using [RepeatMasker](https://github.com/Dfam-consortium/TETools) (v1.4)
```sh
RepeatMasker -a -gff -x -pa {threads} -e ncbi -lib {input.Repeat} {input.assembly}
```
13. Classify transpsosable element (TE) by [TEclass](https://hub.docker.com/repository/docker/hatimalmutairi/teclass-2.1.3b) (v2.1.3b)
```sh
TEclassTest.pl -r {input} > {output}
```
14. Annotation using [MAKER](https://hub.docker.com/r/hatimalmutairi/lmgaap-maker) (v2.31.10)
```sh
maker -fix_nucleotides -genome {input.g1} -base {wildcards.sample}_Round{1} maker_opts.ctl maker_bopts.ctl maker_exe.ctl
gff3_merge -n -s -d {input} > {output}
fasta_merge -d {input} > {output}
maker_map_ids --prefix {wildcards.sample}_ --justify 5 {input} > {output}
maker_functional_fasta {input.Uniprot} {input} {input.gff} > {output.gff}
maker_functional_fasta {input.Uniprot} {input} {input.proteins} > {output.proteins}
maker_functional_fasta {input.Uniprot} {input} {input.transcripts} > {output.transcripts}
```
15. Annotation QC using [Genometools](https://quay.io/repository/biocontainers/genometools?tag=1.2.1--py27_0&tab=tags) (v1.2.1)
```sh
gt gff3 -sort -tidy {input} | gt stat > {output}
```
16. processing annotation using [GAAS](https://github.com/NBISweden/GAAS) (v1.2.0)
```sh
gaas_maker_merge_outputs_from_datastore.pl -i {input} -o {output}
```
17. Assign Functional Annotations using [interproscan](https://github.com/blaxterlab/interproscan-docker) (v5.22-61.0) and [Blast+](https://github.com/ncbi/blast_plus_docs)
```sh
blastp -num_threads {threads} -query {input} -db {input} -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out {output}
interproscan.sh -appl pfam -dp -f TSV -goterms -iprlookup -pa -t p -i {input} -o {output}
```
18. Keep Longest Isoform using [AGAT](https://github.com/NBISweden/AGAT) (v0.6.0)
```sh
agat_sp_keep_longest_isoform.pl -gff {input} -o {output}
```
19. Create the Final Annotation file using [MAKER](https://hub.docker.com/r/hatimalmutairi/lmgaap-maker) (v2.31.10)
```sh
ipr_update_gff {input.gff} {input.iprscan} > {output}
```
