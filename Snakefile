samples=["LSCM1","LSCM4","CUR178","GH5","JKF63","JIQ42"]
report: "Snakefile"

rule all:
	input:
		expand("results/{sample}_QC_Assemblies_by_QUAST/report.html", sample=samples),
		expand("results/annotation/{sample}/Final_Annotation/Final_Annotation.gff", sample=samples),
		expand("results/annotation/{sample}/Final_Annotation/Final_Annotation_Longest_Isoform.gff", sample=samples),
		expand("results/annotation/{sample}/Final_Annotation/annotation.proteins.fasta", sample=samples),
		expand("results/annotation/{sample}/Final_Annotation/annotation.transcripts.fasta", sample=samples),
include: "workflow/rules/get_fastq.smk",
include: "workflow/rules/LSCM1.smk",
include: "workflow/rules/LSCM4.smk",
include: "workflow/rules/CUR178.smk",
include: "workflow/rules/GH5.smk",
include: "workflow/rules/JKF63.smk",
include: "workflow/rules/JIQ42.smk",
include: "workflow/rules/Annotation.smk",


#echo "s,h:m:s,max_rss,max_vms,max_uss,max_pss,io_in,io_out,mean_load,cpu_time" > benchmarks_ALL.csv && grep -v "^s" *.tsv| grep ""| sed 's/:/\t/g' | sed 's/\t/,/g' >> benchmarks_ALL.csv
