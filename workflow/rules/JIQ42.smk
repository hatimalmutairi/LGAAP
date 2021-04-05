###
rule JIQ42_Assembly_LongRead_by_Flye:
	input:
		long_1= "data/JIQ42/SRR13558767.fastq",
	output:
		"results/JIQ42_Assembly_LongRead_by_Flye/assembly.fasta",
	threads: workflow.cores
	group: "Assembly"
	log:
		"results/log/JIQ42_Assembly_LongRead_by_Flye.log"
	benchmark:
		"results/benchmarks/JIQ42_Assembly_LongRead_by_Flye.tsv"
	conda:
		"envs/Flye.yaml"
	shell:
		"flye --nano-raw {input.long_1} --genome-size 35m --threads {threads} -o results/JIQ42_Assembly_LongRead_by_Flye 2>&1 | tee {log}"
###
rule JIQ42_Index_Assembly_by_Samtools:
	input:
		"results/JIQ42_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		"results/JIQ42_Assembly_LongRead_by_Flye/assembly.fasta.fai",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/JIQ42_Index_Assembly_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools faidx {input} > {output}"
###	
rule JIQ42_Map_ShortReads_by_minimap2:
	input:
		PE1_1="data/JIQ42/SRR13558765_1.fastq",
		PE1_2="data/JIQ42/SRR13558765_2.fastq",
		PE2_1="data/JIQ42/SRR13558764_1.fastq",
		PE2_2="data/JIQ42/SRR13558764_2.fastq",
		PE3_1="data/JIQ42/SRR13558766_1.fastq",
		PE3_2="data/JIQ42/SRR13558766_2.fastq",
		assembly="results/JIQ42_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		sam1="results/JIQ42_Assembly_LongRead_by_Flye/PE1.sam",
		sam2="results/JIQ42_Assembly_LongRead_by_Flye/PE2.sam",
		sam3="results/JIQ42_Assembly_LongRead_by_Flye/PE3.sam",
	threads: workflow.cores
	group: "Polishing_1"
	log:
		"results/log/JIQ42_Map_ShortReads_by_minimap2.log"
	benchmark:
		"results/benchmarks/JIQ42_Map_ShortReads_by_minimap2.tsv"
	conda:
		"envs/minimap2.yaml"
	shell:
		"""
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE1_1} {input.PE1_2} > {output.sam1}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE2_1} {input.PE2_2} > {output.sam2}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE3_1} {input.PE3_2} > {output.sam3}
		"""
###
rule JIQ42_View_ShortReads_by_Samtools:
	input:
		sam1="results/JIQ42_Assembly_LongRead_by_Flye/PE1.sam",
		sam2="results/JIQ42_Assembly_LongRead_by_Flye/PE2.sam",
		sam3="results/JIQ42_Assembly_LongRead_by_Flye/PE3.sam",
		assembly="results/JIQ42_Assembly_LongRead_by_Flye/assembly.fasta",
		idx="results/JIQ42_Assembly_LongRead_by_Flye/assembly.fasta.fai",
	output:
		bam1="results/JIQ42_Assembly_LongRead_by_Flye/PE1.bam",
		bam2="results/JIQ42_Assembly_LongRead_by_Flye/PE2.bam",
		bam3="results/JIQ42_Assembly_LongRead_by_Flye/PE3.bam",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/JIQ42_View_ShortReads_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools view -@ {threads} -bt {input.idx} {input.sam1} > {output.bam1}
		samtools view -@ {threads} -bt {input.idx} {input.sam2} > {output.bam2}
		samtools view -@ {threads} -bt {input.idx} {input.sam3} > {output.bam3}
		"""
###
rule JIQ42_Sort_ShortReads_by_Samtools:
	input:
		bam1="results/JIQ42_Assembly_LongRead_by_Flye/PE1.bam",
		bam2="results/JIQ42_Assembly_LongRead_by_Flye/PE2.bam",
		bam3="results/JIQ42_Assembly_LongRead_by_Flye/PE3.bam",
	output:
		bam_sort1="results/JIQ42_Assembly_LongRead_by_Flye/PE_sort1.bam",
		bam_sort2="results/JIQ42_Assembly_LongRead_by_Flye/PE_sort2.bam",
		bam_sort3="results/JIQ42_Assembly_LongRead_by_Flye/PE_sort3.bam",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/JIQ42_Sort_ShortReads_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools sort -@ {threads} {input.bam1} -o {output.bam_sort1}
		samtools sort -@ {threads} {input.bam2} -o {output.bam_sort2}
		samtools sort -@ {threads} {input.bam3} -o {output.bam_sort3}
		"""
###
rule JIQ42_Index_ShortReads_by_Samtools:
	input:
		bam_sort1="results/JIQ42_Assembly_LongRead_by_Flye/PE_sort1.bam",
		bam_sort2="results/JIQ42_Assembly_LongRead_by_Flye/PE_sort2.bam",
		bam_sort3="results/JIQ42_Assembly_LongRead_by_Flye/PE_sort3.bam",
	output:
		bam_idx1="results/JIQ42_Assembly_LongRead_by_Flye/PE_idx1.bai",
		bam_idx2="results/JIQ42_Assembly_LongRead_by_Flye/PE_idx2.bai",
		bam_idx3="results/JIQ42_Assembly_LongRead_by_Flye/PE_idx3.bai",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/JIQ42_Index_ShortReads_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools index -@ {threads} {input.bam_sort1} > {output.bam_idx1}
		samtools index -@ {threads} {input.bam_sort2} > {output.bam_idx2}
		samtools index -@ {threads} {input.bam_sort3} > {output.bam_idx3}
		"""
###
rule JIQ42_Mpileup_Assembly_by_BCFtools:
	input:
		bam_sort1="results/JIQ42_Assembly_LongRead_by_Flye/PE_sort1.bam",
		bam_sort2="results/JIQ42_Assembly_LongRead_by_Flye/PE_sort2.bam",
		bam_sort3="results/JIQ42_Assembly_LongRead_by_Flye/PE_sort3.bam",
		assembly="results/JIQ42_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		PE="results/JIQ42_Assembly_LongRead_by_Flye/PE.vcf.gz",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/JIQ42_Mpileup_Assembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools mpileup --threads {threads} -Ou -d 1000000 -f {input.assembly} {input.bam_sort1} {input.bam_sort2} {input.bam_sort3} | bcftools call -mv -Oz -o {output.PE}"
###
rule JIQ42_Norm_Assembly_by_BCFtools:
	input:
		PE="results/JIQ42_Assembly_LongRead_by_Flye/PE.vcf.gz",
		assembly="results/JIQ42_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		"results/JIQ42_Assembly_LongRead_by_Flye/PE.norm.bcf",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/JIQ42_Norm_Assembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools norm --threads {threads} -f {input.assembly} {input.PE} -Ob -o {output}"
###
rule JIQ42_Filter_Assembly_by_BCFtools:
	input:
		PE="results/JIQ42_Assembly_LongRead_by_Flye/PE.norm.bcf",
	output:
		"results/JIQ42_Assembly_LongRead_by_Flye/PE.norm.flt-indels.bcf",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/JIQ42_Filter_Assembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools filter --threads {threads} --IndelGap 5 {input.PE} -Ob -o {output}"
###
rule JIQ42_Index_Assembly_by_BCFtools:
	input:
		"results/JIQ42_Assembly_LongRead_by_Flye/PE.norm.flt-indels.bcf",
	output:
		"results/JIQ42_Assembly_LongRead_by_Flye/PE.norm.flt-indels.bcf.csi",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/JIQ42_Index_Assembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools index --threads {threads} {input} -o {output}"
###
rule JIQ42_Consensus_Assembly_by_BCFtools:
	input:
		"results/JIQ42_Assembly_LongRead_by_Flye/PE.norm.flt-indels.bcf.csi",
		PE="results/JIQ42_Assembly_LongRead_by_Flye/PE.norm.flt-indels.bcf",
		assembly="results/JIQ42_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		"results/JIQ42_Assembly_LongRead_by_Flye/Flye_consensus.fasta",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/JIQ42_Consensus_Assembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools consensus -f {input.assembly} {input.PE} > {output}"
###
rule JIQ42_Polish_Assembly_by_Pilon:
	input:
		bam_idx1="results/JIQ42_Assembly_LongRead_by_Flye/PE_idx1.bai",
		bam_idx2="results/JIQ42_Assembly_LongRead_by_Flye/PE_idx2.bai",
		bam_idx3="results/JIQ42_Assembly_LongRead_by_Flye/PE_idx3.bai",
		bam_sort1="results/JIQ42_Assembly_LongRead_by_Flye/PE_sort1.bam",
		bam_sort2="results/JIQ42_Assembly_LongRead_by_Flye/PE_sort2.bam",
		bam_sort3="results/JIQ42_Assembly_LongRead_by_Flye/PE_sort3.bam",
		assembly="results/JIQ42_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		"results/JIQ42_Polish_Assembly_by_Pilon/pilon.fasta",
	threads: workflow.cores
	group: "Polishing_2"
	log:
		"results/log/JIQ42_Polish_Assembly_by_Pilon.log"
	benchmark:
		"results/benchmarks/JIQ42_Polish_Assembly_by_Pilon.tsv"
	conda:
		"envs/pilon.yaml"
	shell:	
		"pilon --threads {threads} -Xmx300000M --genome {input.assembly} --bam {input.bam_sort1} --bam {input.bam_sort2} --bam {input.bam_sort3} --outdir results/JIQ42_Polish_Assembly_by_Pilon 2>&1 | tee {log}"
###
rule JIQ42_index_PolishedAssembly_JIQ42: 
	input:
		"results/JIQ42_Polish_Assembly_by_Pilon/pilon.fasta",
	output:
		"results/JIQ42_Polish_Assembly_by_Pilon/pilon.fai",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/JIQ42_index_PolishedAssembly_JIQ42.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools faidx {input} > {output}"
###	
checkpoint JIQ42_Map_PolishedAssembly_by_minimap2:
	input:
		PE1_1="data/JIQ42/SRR13558765_1.fastq",
		PE1_2="data/JIQ42/SRR13558765_2.fastq",
		PE2_1="data/JIQ42/SRR13558764_1.fastq",
		PE2_2="data/JIQ42/SRR13558764_2.fastq",
		PE3_1="data/JIQ42/SRR13558766_1.fastq",
		PE3_2="data/JIQ42/SRR13558766_2.fastq",
		assembly="results/JIQ42_Polish_Assembly_by_Pilon/pilon.fasta",
	output:
		sam1="results/JIQ42_Polish_Assembly_by_Pilon/PE1.sam",
		sam2="results/JIQ42_Polish_Assembly_by_Pilon/PE2.sam",
		sam3="results/JIQ42_Polish_Assembly_by_Pilon/PE3.sam",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/JIQ42_Map_PolishedAssembly_by_minimap2.tsv"
	conda:
		"envs/minimap2.yaml"
	shell:
		"""
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE1_1} {input.PE1_2} > {output.sam1}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE2_1} {input.PE2_2} > {output.sam2}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE3_1} {input.PE3_2} > {output.sam3}
		"""
###
rule JIQ42_View_PolishedAssembly_by_Samtools:
	input:
		sam1="results/JIQ42_Polish_Assembly_by_Pilon/PE1.sam",
		sam2="results/JIQ42_Polish_Assembly_by_Pilon/PE2.sam",
		sam3="results/JIQ42_Polish_Assembly_by_Pilon/PE3.sam",
		assembly="results/JIQ42_Polish_Assembly_by_Pilon/pilon.fasta",
		idx="results/JIQ42_Polish_Assembly_by_Pilon/pilon.fai",
	output:
		bam1="results/JIQ42_Polish_Assembly_by_Pilon/PE1.bam",
		bam2="results/JIQ42_Polish_Assembly_by_Pilon/PE2.bam",
		bam3="results/JIQ42_Polish_Assembly_by_Pilon/PE3.bam",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/JIQ42_View_PolishedAssembly_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools view -@ {threads} -bt {input.idx} {input.sam1} > {output.bam1}
		samtools view -@ {threads} -bt {input.idx} {input.sam2} > {output.bam2}
		samtools view -@ {threads} -bt {input.idx} {input.sam3} > {output.bam3}
		"""
###
rule JIQ42_Sort_PolishedAssembly_by_Samtools:
	input:
		bam1="results/JIQ42_Polish_Assembly_by_Pilon/PE1.bam",
		bam2="results/JIQ42_Polish_Assembly_by_Pilon/PE2.bam",
		bam3="results/JIQ42_Polish_Assembly_by_Pilon/PE3.bam",
	output:
		bam_sort1="results/JIQ42_Polish_Assembly_by_Pilon/PE_sort1.bam",
		bam_sort2="results/JIQ42_Polish_Assembly_by_Pilon/PE_sort2.bam",
		bam_sort3="results/JIQ42_Polish_Assembly_by_Pilon/PE_sort3.bam",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/JIQ42_Sort_PolishedAssembly_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools sort -@ {threads} {input.bam1} -o {output.bam_sort1}
		samtools sort -@ {threads} {input.bam2} -o {output.bam_sort2}
		samtools sort -@ {threads} {input.bam3} -o {output.bam_sort3}
		"""
###
rule JIQ42_Index_PolishedAssembly_by_Samtools:
	input:
		bam_sort1="results/JIQ42_Polish_Assembly_by_Pilon/PE_sort1.bam",
		bam_sort2="results/JIQ42_Polish_Assembly_by_Pilon/PE_sort2.bam",
		bam_sort3="results/JIQ42_Polish_Assembly_by_Pilon/PE_sort3.bam",
	output:
		bam_sort1="results/JIQ42_Polish_Assembly_by_Pilon/PE_sort1.bai",
		bam_sort2="results/JIQ42_Polish_Assembly_by_Pilon/PE_sort2.bai",
		bam_sort3="results/JIQ42_Polish_Assembly_by_Pilon/PE_sort3.bai",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/JIQ42_Index_PolishedAssembly_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools index -@ {threads} {input.bam_sort1} > {output.bam_sort1}
		samtools index -@ {threads} {input.bam_sort2} > {output.bam_sort2}
		samtools index -@ {threads} {input.bam_sort3} > {output.bam_sort3}
		"""
###
rule JIQ42_Mpileup_PolishedAssembly_by_BCFtools:
	input:
		bam_sort1="results/JIQ42_Polish_Assembly_by_Pilon/PE_sort1.bam",
		bam_sort2="results/JIQ42_Polish_Assembly_by_Pilon/PE_sort2.bam",
		bam_sort3="results/JIQ42_Polish_Assembly_by_Pilon/PE_sort3.bam",
		assembly="results/JIQ42_Polish_Assembly_by_Pilon/pilon.fasta",
	output:
		PE="results/JIQ42_Polish_Assembly_by_Pilon/PE.vcf.gz",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/JIQ42_Mpileup_PolishedAssembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools mpileup --threads {threads} -Ou -d 1000000 -f {input.assembly} {input.bam_sort1} {input.bam_sort2} {input.bam_sort3}| bcftools call -mv -Oz -o {output.PE}"
###
rule JIQ42_Norm_PolishedAssembly_by_BCFtools:
	input:
		PE="results/JIQ42_Polish_Assembly_by_Pilon/PE.vcf.gz",
		assembly="results/JIQ42_Polish_Assembly_by_Pilon/pilon.fasta",
	output:
		"results/JIQ42_Polish_Assembly_by_Pilon/PE.norm.bcf",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/JIQ42_Norm_PolishedAssembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools norm --threads {threads} -f {input.assembly} {input.PE} -Ob -o {output}"
###
rule JIQ42_Filter_PolishedAssembly_by_BCFtools:
	input:
		PE="results/JIQ42_Polish_Assembly_by_Pilon/PE.norm.bcf",
	output:
		"results/JIQ42_Polish_Assembly_by_Pilon/PE.norm.flt-indels.bcf",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/JIQ42_Filter_PolishedAssembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools filter --threads {threads} --IndelGap 5 {input.PE} -Ob -o {output}"
###
rule JIQ42_Index_PolishedAssembly_by_BCFtools:
	input:
		"results/JIQ42_Polish_Assembly_by_Pilon/PE.norm.flt-indels.bcf",
	output:
		"results/JIQ42_Polish_Assembly_by_Pilon/PE.norm.flt-indels.bcf.csi",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/JIQ42_Index_PolishedAssembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools index --threads {threads} {input} -o {output}"
###
rule JIQ42_Consensus_PolishedAssembly_by_BCFtools:
	input:
		"results/JIQ42_Polish_Assembly_by_Pilon/PE.norm.flt-indels.bcf.csi",
		PE="results/JIQ42_Polish_Assembly_by_Pilon/PE.norm.flt-indels.bcf",
		assembly="results/JIQ42_Polish_Assembly_by_Pilon/pilon.fasta",
	output:
		"results/JIQ42_Polish_Assembly_by_Pilon/Pilon_consensus.fasta",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/JIQ42_Consensus_PolishedAssembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools consensus -f {input.assembly} {input.PE} > {output}"
###
rule JIQ42_aggregate_reads:
	input:
		PE1_1="data/JIQ42/SRR13558765_1.fastq",
		PE1_2="data/JIQ42/SRR13558765_2.fastq",
		PE2_1="data/JIQ42/SRR13558764_1.fastq",
		PE2_2="data/JIQ42/SRR13558764_2.fastq",
		PE3_1="data/JIQ42/SRR13558766_1.fastq",
		PE3_2="data/JIQ42/SRR13558766_2.fastq",
	output:
		"data/JIQ42/JIQ42_Short_reads.fastq"
	group: "Polishing_3"
	shell:
		"cat {input.PE1_1} {input.PE1_2} {input.PE2_1} {input.PE2_2} {input.PE3_1} {input.PE3_2} > {output}"
###
rule JIQ42_pre_RaGOO_Short_reads:
	input:
		"data/JIQ42/JIQ42_Short_reads.fastq",
	output:
		"JIQ42_Short_reads.fastq",
	shell:
		"cp {input} {output}"
###
###
rule JIQ42_Sort_PolishedAssembly_by_funannotate:
	input:
		"results/JIQ42_Polish_Assembly_by_Pilon/Pilon_consensus.fasta",
	output:
		"results/JIQ42_Sort_PolishedAssembly_by_funannotate/JIQ42_Draft_Sort.fasta",
	threads: workflow.cores	
	group: "Polishing_3"
	log:
		"results/log/JIQ42_Sort_PolishedAssembly_by_funannotate.log"
	benchmark:
		"results/benchmarks/JIQ42_Sort_PolishedAssembly_by_funannotate.tsv"
	container:
		"docker://nextgenusfs/funannotate:latest"
	shell:
		"docker run -u 0 -v $(pwd):/home/data -w /home/data nextgenusfs/funannotate:latest funannotate sort -i {input} -o {output}"
###
rule JIQ42_pre_RaGOO_Polished_Assembly:
	input:
		"results/JIQ42_Sort_PolishedAssembly_by_funannotate/JIQ42_Draft_Sort.fasta",
	output:
		"JIQ42_Polished_Assembly.fasta",
	shell:
		"cp {input} {output}"
###
rule JIQ42_Arrange_Assembly_by_RaGOO:
	input:
		reads= "JIQ42_Short_reads.fastq",
		assembly= "JIQ42_Polished_Assembly.fasta",
		ref= "LmajorFriedlin_Genome.fasta",
	output:
		"results/JIQ42_Arrange_Assembly_by_RaGOO/ragoo_output/ragoo.fasta",
	threads: workflow.cores
	group: "Polishing_3"
	log:
		"results/log/JIQ42_Arrange_Assembly_by_RaGOO.log"
	benchmark:
		"results/benchmarks/JIQ42_Arrange_Assembly_by_RaGOO.tsv"
	conda:
		"envs/ragoo.yaml"
	shell:
		"rm -rf ragoo_output; ragoo.py -t {threads} -C -T sr -b -g 0 -R {input.reads} {input.assembly} {input.ref}; cp -r ragoo_output results/JIQ42_Arrange_Assembly_by_RaGOO; rm -rf JIQ42_Short_reads.fastq; rm -rf JIQ42_Polished_Assembly.fasta; rm -rf ragoo_output"
###
rule JIQ42_Clean_PolishedAssembly_by_funannotate:
	input:
		"results/JIQ42_Arrange_Assembly_by_RaGOO/ragoo_output/ragoo.fasta",
	output:
		"results/JIQ42_Clean_PolishedAssembly_by_funannotate/funannotate_Assembly.fasta",
	threads: workflow.cores	
	group: "Polishing_3"	
	params:
		"--exhaustive --minlen 50 --pident 10",
	log:
		"results/log/Clean_JIQ42_Clean_PolishedAssembly_by_funannotate.log"
	benchmark:
		"results/benchmarks/Clean_JIQ42_Clean_PolishedAssembly_by_funannotate.tsv"
	container: "docker://nextgenusfs/funannotate:latest"
	shell:
		"docker run -u 0 -v $(pwd):/home/data -w /home/data nextgenusfs/funannotate:latest funannotate clean {params} -i {input} -o {output}"
###
rule JIQ42_Final_assembly:
	input:
		"results/JIQ42_Clean_PolishedAssembly_by_funannotate/funannotate_Assembly.fasta",
	output:
		"results/Final_assembly/JIQ42_Genome.fasta",
	shell:
		"cp {input} {output}; sed -i 's/_RaGOO//g' {output}; sed -i 's/LmjF/JIQ42/g' {output}"
###
# Assessing Quality and Completeness:
###
checkpoint JIQ42_QC_Assemblies_by_QUAST:
	input:
		"results/JIQ42_Assembly_LongRead_by_Flye/assembly.fasta",
		"results/JIQ42_Assembly_LongRead_by_Flye/Flye_consensus.fasta",
		"results/JIQ42_Polish_Assembly_by_Pilon/pilon.fasta",
		"results/JIQ42_Polish_Assembly_by_Pilon/Pilon_consensus.fasta",
		"results/JIQ42_Arrange_Assembly_by_RaGOO/ragoo_output/ragoo.fasta",
		"results/JIQ42_Clean_PolishedAssembly_by_funannotate/funannotate_Assembly.fasta",
		"results/Final_assembly/JIQ42_Genome.fasta",
	output:
		"results/JIQ42_QC_Assemblies_by_QUAST/report.html",
	params:
		"--eukaryote --circos",
	threads: workflow.cores
	group: "Quality_Assessment"
	log:
		"results/log/JIQ42_QC_Assemblies_by_QUAST.log"
	benchmark:
		"results/benchmarks/JIQ42_QC_Assemblies_by_QUAST.tsv"
	container:
		"docker://quay.io/biocontainers/quast:5.0.1--py27pl526ha92aebf_0"
	#conda:
		#"envs/quast.yaml"
	shell:
		"docker run -v $(pwd):/home/data -w /home/data quay.io/biocontainers/quast:5.0.1--py27pl526ha92aebf_0 quast.py {params} --threads {threads} --output-dir results/JIQ42_QC_Assemblies_by_QUAST {input} 2>&1 | tee {log}" 
###