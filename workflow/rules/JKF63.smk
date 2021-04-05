###
rule JKF63_Assembly_LongRead_by_Flye:
	input:
		long_1= "data/JKF63/SRR13558757.fastq",
	output:
		"results/JKF63_Assembly_LongRead_by_Flye/assembly.fasta",
	threads: workflow.cores
	group: "Assembly"
	log:
		"results/log/JKF63_Assembly_LongRead_by_Flye.log"
	benchmark:
		"results/benchmarks/JKF63_Assembly_LongRead_by_Flye.tsv"
	conda:
		"envs/Flye.yaml"
	shell:
		"flye --nano-raw {input.long_1} --genome-size 35m --threads {threads} -o results/JKF63_Assembly_LongRead_by_Flye 2>&1 | tee {log}"
###
rule JKF63_Index_Assembly_by_Samtools:
	input:
		"results/JKF63_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		"results/JKF63_Assembly_LongRead_by_Flye/assembly.fasta.fai",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/JKF63_Index_Assembly_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools faidx {input} > {output}"
###	
rule JKF63_Map_ShortReads_by_minimap2:
	input:
		PE1_1="data/JKF63/SRR13558756_1.fastq",
		PE1_2="data/JKF63/SRR13558756_2.fastq",
		PE2_1="data/JKF63/SRR13558754_1.fastq",
		PE2_2="data/JKF63/SRR13558754_2.fastq",
		PE3_1="data/JKF63/SRR13558755_1.fastq",
		PE3_2="data/JKF63/SRR13558755_2.fastq",
		assembly="results/JKF63_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		sam1="results/JKF63_Assembly_LongRead_by_Flye/PE1.sam",
		sam2="results/JKF63_Assembly_LongRead_by_Flye/PE2.sam",
		sam3="results/JKF63_Assembly_LongRead_by_Flye/PE3.sam",
	threads: workflow.cores
	group: "Polishing_1"
	log:
		"results/log/JKF63_Map_ShortReads_by_minimap2.log"
	benchmark:
		"results/benchmarks/JKF63_Map_ShortReads_by_minimap2.tsv"
	conda:
		"envs/minimap2.yaml"
	shell:
		"""
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE1_1} {input.PE1_2} > {output.sam1}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE2_1} {input.PE2_2} > {output.sam2}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE3_1} {input.PE3_2} > {output.sam3}
		"""
###
rule JKF63_View_ShortReads_by_Samtools:
	input:
		sam1="results/JKF63_Assembly_LongRead_by_Flye/PE1.sam",
		sam2="results/JKF63_Assembly_LongRead_by_Flye/PE2.sam",
		sam3="results/JKF63_Assembly_LongRead_by_Flye/PE3.sam",
		assembly="results/JKF63_Assembly_LongRead_by_Flye/assembly.fasta",
		idx="results/JKF63_Assembly_LongRead_by_Flye/assembly.fasta.fai",
	output:
		bam1="results/JKF63_Assembly_LongRead_by_Flye/PE1.bam",
		bam2="results/JKF63_Assembly_LongRead_by_Flye/PE2.bam",
		bam3="results/JKF63_Assembly_LongRead_by_Flye/PE3.bam",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/JKF63_View_ShortReads_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools view -@ {threads} -bt {input.idx} {input.sam1} > {output.bam1}
		samtools view -@ {threads} -bt {input.idx} {input.sam2} > {output.bam2}
		samtools view -@ {threads} -bt {input.idx} {input.sam3} > {output.bam3}
		"""
###
rule JKF63_Sort_ShortReads_by_Samtools:
	input:
		bam1="results/JKF63_Assembly_LongRead_by_Flye/PE1.bam",
		bam2="results/JKF63_Assembly_LongRead_by_Flye/PE2.bam",
		bam3="results/JKF63_Assembly_LongRead_by_Flye/PE3.bam",
	output:
		bam_sort1="results/JKF63_Assembly_LongRead_by_Flye/PE_sort1.bam",
		bam_sort2="results/JKF63_Assembly_LongRead_by_Flye/PE_sort2.bam",
		bam_sort3="results/JKF63_Assembly_LongRead_by_Flye/PE_sort3.bam",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/JKF63_Sort_ShortReads_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools sort -@ {threads} {input.bam1} -o {output.bam_sort1}
		samtools sort -@ {threads} {input.bam2} -o {output.bam_sort2}
		samtools sort -@ {threads} {input.bam3} -o {output.bam_sort3}
		"""
###
rule JKF63_Index_ShortReads_by_Samtools:
	input:
		bam_sort1="results/JKF63_Assembly_LongRead_by_Flye/PE_sort1.bam",
		bam_sort2="results/JKF63_Assembly_LongRead_by_Flye/PE_sort2.bam",
		bam_sort3="results/JKF63_Assembly_LongRead_by_Flye/PE_sort3.bam",
	output:
		bam_idx1="results/JKF63_Assembly_LongRead_by_Flye/PE_idx1.bai",
		bam_idx2="results/JKF63_Assembly_LongRead_by_Flye/PE_idx2.bai",
		bam_idx3="results/JKF63_Assembly_LongRead_by_Flye/PE_idx3.bai",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/JKF63_Index_ShortReads_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools index -@ {threads} {input.bam_sort1} > {output.bam_idx1}
		samtools index -@ {threads} {input.bam_sort2} > {output.bam_idx2}
		samtools index -@ {threads} {input.bam_sort3} > {output.bam_idx3}
		"""
###
rule JKF63_Mpileup_Assembly_by_BCFtools:
	input:
		bam_sort1="results/JKF63_Assembly_LongRead_by_Flye/PE_sort1.bam",
		bam_sort2="results/JKF63_Assembly_LongRead_by_Flye/PE_sort2.bam",
		bam_sort3="results/JKF63_Assembly_LongRead_by_Flye/PE_sort3.bam",
		assembly="results/JKF63_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		PE="results/JKF63_Assembly_LongRead_by_Flye/PE.vcf.gz",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/JKF63_Mpileup_Assembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools mpileup --threads {threads} -Ou -d 1000000 -f {input.assembly} {input.bam_sort1} {input.bam_sort2} {input.bam_sort3} | bcftools call -mv -Oz -o {output.PE}"
###
rule JKF63_Norm_Assembly_by_BCFtools:
	input:
		PE="results/JKF63_Assembly_LongRead_by_Flye/PE.vcf.gz",
		assembly="results/JKF63_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		"results/JKF63_Assembly_LongRead_by_Flye/PE.norm.bcf",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/JKF63_Norm_Assembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools norm --threads {threads} -f {input.assembly} {input.PE} -Ob -o {output}"
###
rule JKF63_Filter_Assembly_by_BCFtools:
	input:
		PE="results/JKF63_Assembly_LongRead_by_Flye/PE.norm.bcf",
	output:
		"results/JKF63_Assembly_LongRead_by_Flye/PE.norm.flt-indels.bcf",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/JKF63_Filter_Assembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools filter --threads {threads} --IndelGap 5 {input.PE} -Ob -o {output}"
###
rule JKF63_Index_Assembly_by_BCFtools:
	input:
		"results/JKF63_Assembly_LongRead_by_Flye/PE.norm.flt-indels.bcf",
	output:
		"results/JKF63_Assembly_LongRead_by_Flye/PE.norm.flt-indels.bcf.csi",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/JKF63_Index_Assembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools index --threads {threads} {input} -o {output}"
###
rule JKF63_Consensus_Assembly_by_BCFtools:
	input:
		"results/JKF63_Assembly_LongRead_by_Flye/PE.norm.flt-indels.bcf.csi",
		PE="results/JKF63_Assembly_LongRead_by_Flye/PE.norm.flt-indels.bcf",
		assembly="results/JKF63_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		"results/JKF63_Assembly_LongRead_by_Flye/Flye_consensus.fasta",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/JKF63_Consensus_Assembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools consensus -f {input.assembly} {input.PE} > {output}"
###
rule JKF63_Polish_Assembly_by_Pilon:
	input:
		bam_idx1="results/JKF63_Assembly_LongRead_by_Flye/PE_idx1.bai",
		bam_idx2="results/JKF63_Assembly_LongRead_by_Flye/PE_idx2.bai",
		bam_idx3="results/JKF63_Assembly_LongRead_by_Flye/PE_idx3.bai",
		bam_sort1="results/JKF63_Assembly_LongRead_by_Flye/PE_sort1.bam",
		bam_sort2="results/JKF63_Assembly_LongRead_by_Flye/PE_sort2.bam",
		bam_sort3="results/JKF63_Assembly_LongRead_by_Flye/PE_sort3.bam",
		assembly="results/JKF63_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		"results/JKF63_Polish_Assembly_by_Pilon/pilon.fasta",
	threads: workflow.cores
	group: "Polishing_2"
	log:
		"results/log/JKF63_Polish_Assembly_by_Pilon.log"
	benchmark:
		"results/benchmarks/JKF63_Polish_Assembly_by_Pilon.tsv"
	conda:
		"envs/pilon.yaml"
	shell:	
		"pilon --threads {threads} -Xmx300000M --genome {input.assembly} --bam {input.bam_sort1} --bam {input.bam_sort2} --bam {input.bam_sort3} --outdir results/JKF63_Polish_Assembly_by_Pilon 2>&1 | tee {log}"
###
rule JKF63_index_PolishedAssembly_JKF63: 
	input:
		"results/JKF63_Polish_Assembly_by_Pilon/pilon.fasta",
	output:
		"results/JKF63_Polish_Assembly_by_Pilon/pilon.fai",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/JKF63_index_PolishedAssembly_JKF63.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools faidx {input} > {output}"
###	
checkpoint JKF63_Map_PolishedAssembly_by_minimap2:
	input:
		PE1_1="data/JKF63/SRR13558756_1.fastq",
		PE1_2="data/JKF63/SRR13558756_2.fastq",
		PE2_1="data/JKF63/SRR13558754_1.fastq",
		PE2_2="data/JKF63/SRR13558754_2.fastq",
		PE3_1="data/JKF63/SRR13558755_1.fastq",
		PE3_2="data/JKF63/SRR13558755_2.fastq",
		assembly="results/JKF63_Polish_Assembly_by_Pilon/pilon.fasta",
	output:
		sam1="results/JKF63_Polish_Assembly_by_Pilon/PE1.sam",
		sam2="results/JKF63_Polish_Assembly_by_Pilon/PE2.sam",
		sam3="results/JKF63_Polish_Assembly_by_Pilon/PE3.sam",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/JKF63_Map_PolishedAssembly_by_minimap2.tsv"
	conda:
		"envs/minimap2.yaml"
	shell:
		"""
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE1_1} {input.PE1_2} > {output.sam1}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE2_1} {input.PE2_2} > {output.sam2}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE3_1} {input.PE3_2} > {output.sam3}
		"""
###
rule JKF63_View_PolishedAssembly_by_Samtools:
	input:
		sam1="results/JKF63_Polish_Assembly_by_Pilon/PE1.sam",
		sam2="results/JKF63_Polish_Assembly_by_Pilon/PE2.sam",
		sam3="results/JKF63_Polish_Assembly_by_Pilon/PE3.sam",
		assembly="results/JKF63_Polish_Assembly_by_Pilon/pilon.fasta",
		idx="results/JKF63_Polish_Assembly_by_Pilon/pilon.fai",
	output:
		bam1="results/JKF63_Polish_Assembly_by_Pilon/PE1.bam",
		bam2="results/JKF63_Polish_Assembly_by_Pilon/PE2.bam",
		bam3="results/JKF63_Polish_Assembly_by_Pilon/PE3.bam",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/JKF63_View_PolishedAssembly_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools view -@ {threads} -bt {input.idx} {input.sam1} > {output.bam1}
		samtools view -@ {threads} -bt {input.idx} {input.sam2} > {output.bam2}
		samtools view -@ {threads} -bt {input.idx} {input.sam3} > {output.bam3}
		"""
###
rule JKF63_Sort_PolishedAssembly_by_Samtools:
	input:
		bam1="results/JKF63_Polish_Assembly_by_Pilon/PE1.bam",
		bam2="results/JKF63_Polish_Assembly_by_Pilon/PE2.bam",
		bam3="results/JKF63_Polish_Assembly_by_Pilon/PE3.bam",
	output:
		bam_sort1="results/JKF63_Polish_Assembly_by_Pilon/PE_sort1.bam",
		bam_sort2="results/JKF63_Polish_Assembly_by_Pilon/PE_sort2.bam",
		bam_sort3="results/JKF63_Polish_Assembly_by_Pilon/PE_sort3.bam",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/JKF63_Sort_PolishedAssembly_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools sort -@ {threads} {input.bam1} -o {output.bam_sort1}
		samtools sort -@ {threads} {input.bam2} -o {output.bam_sort2}
		samtools sort -@ {threads} {input.bam3} -o {output.bam_sort3}
		"""
###
rule JKF63_Index_PolishedAssembly_by_Samtools:
	input:
		bam_sort1="results/JKF63_Polish_Assembly_by_Pilon/PE_sort1.bam",
		bam_sort2="results/JKF63_Polish_Assembly_by_Pilon/PE_sort2.bam",
		bam_sort3="results/JKF63_Polish_Assembly_by_Pilon/PE_sort3.bam",
	output:
		bam_sort1="results/JKF63_Polish_Assembly_by_Pilon/PE_sort1.bai",
		bam_sort2="results/JKF63_Polish_Assembly_by_Pilon/PE_sort2.bai",
		bam_sort3="results/JKF63_Polish_Assembly_by_Pilon/PE_sort3.bai",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/JKF63_Index_PolishedAssembly_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools index -@ {threads} {input.bam_sort1} > {output.bam_sort1}
		samtools index -@ {threads} {input.bam_sort2} > {output.bam_sort2}
		samtools index -@ {threads} {input.bam_sort3} > {output.bam_sort3}
		"""
###
rule JKF63_Mpileup_PolishedAssembly_by_BCFtools:
	input:
		bam_sort1="results/JKF63_Polish_Assembly_by_Pilon/PE_sort1.bam",
		bam_sort2="results/JKF63_Polish_Assembly_by_Pilon/PE_sort2.bam",
		bam_sort3="results/JKF63_Polish_Assembly_by_Pilon/PE_sort3.bam",
		assembly="results/JKF63_Polish_Assembly_by_Pilon/pilon.fasta",
	output:
		PE="results/JKF63_Polish_Assembly_by_Pilon/PE.vcf.gz",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/JKF63_Mpileup_PolishedAssembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools mpileup --threads {threads} -Ou -d 1000000 -f {input.assembly} {input.bam_sort1} {input.bam_sort2} {input.bam_sort3}| bcftools call -mv -Oz -o {output.PE}"
###
rule JKF63_Norm_PolishedAssembly_by_BCFtools:
	input:
		PE="results/JKF63_Polish_Assembly_by_Pilon/PE.vcf.gz",
		assembly="results/JKF63_Polish_Assembly_by_Pilon/pilon.fasta",
	output:
		"results/JKF63_Polish_Assembly_by_Pilon/PE.norm.bcf",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/JKF63_Norm_PolishedAssembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools norm --threads {threads} -f {input.assembly} {input.PE} -Ob -o {output}"
###
rule JKF63_Filter_PolishedAssembly_by_BCFtools:
	input:
		PE="results/JKF63_Polish_Assembly_by_Pilon/PE.norm.bcf",
	output:
		"results/JKF63_Polish_Assembly_by_Pilon/PE.norm.flt-indels.bcf",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/JKF63_Filter_PolishedAssembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools filter --threads {threads} --IndelGap 5 {input.PE} -Ob -o {output}"
###
rule JKF63_Index_PolishedAssembly_by_BCFtools:
	input:
		"results/JKF63_Polish_Assembly_by_Pilon/PE.norm.flt-indels.bcf",
	output:
		"results/JKF63_Polish_Assembly_by_Pilon/PE.norm.flt-indels.bcf.csi",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/JKF63_Index_PolishedAssembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools index --threads {threads} {input} -o {output}"
###
rule JKF63_Consensus_PolishedAssembly_by_BCFtools:
	input:
		"results/JKF63_Polish_Assembly_by_Pilon/PE.norm.flt-indels.bcf.csi",
		PE="results/JKF63_Polish_Assembly_by_Pilon/PE.norm.flt-indels.bcf",
		assembly="results/JKF63_Polish_Assembly_by_Pilon/pilon.fasta",
	output:
		"results/JKF63_Polish_Assembly_by_Pilon/Pilon_consensus.fasta",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/JKF63_Consensus_PolishedAssembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools consensus -f {input.assembly} {input.PE} > {output}"
###
rule JKF63_aggregate_reads:
	input:
		PE1_1="data/JKF63/SRR13558756_1.fastq",
		PE1_2="data/JKF63/SRR13558756_2.fastq",
		PE2_1="data/JKF63/SRR13558754_1.fastq",
		PE2_2="data/JKF63/SRR13558754_2.fastq",
		PE3_1="data/JKF63/SRR13558755_1.fastq",
		PE3_2="data/JKF63/SRR13558755_2.fastq",
	output:
		"data/JKF63/JKF63_Short_reads.fastq"
	group: "Polishing_3"
	shell:
		"cat {input.PE1_1} {input.PE1_2} {input.PE2_1} {input.PE2_2} {input.PE3_1} {input.PE3_2} > {output}"
###
rule JKF63_pre_RaGOO_Short_reads:
	input:
		"data/JKF63/JKF63_Short_reads.fastq",
	output:
		"JKF63_Short_reads.fastq",
	shell:
		"cp {input} {output}"
###
rule JKF63_Sort_PolishedAssembly_by_funannotate:
	input:
		"results/JKF63_Polish_Assembly_by_Pilon/Pilon_consensus.fasta",
	output:
		"results/JKF63_Sort_PolishedAssembly_by_funannotate/JKF63_Draft_Sort.fasta",
	threads: workflow.cores	
	group: "Polishing_3"
	log:
		"results/log/JKF63_Sort_PolishedAssembly_by_funannotate.log"
	benchmark:
		"results/benchmarks/JKF63_Sort_PolishedAssembly_by_funannotate.tsv"
	container:
		"docker://nextgenusfs/funannotate:latest"
	shell:
		"docker run -u 0 -v $(pwd):/home/data -w /home/data nextgenusfs/funannotate:latest funannotate sort -i {input} -o {output}"
###
rule JKF63_pre_RaGOO_Polished_Assembly:
	input:
		"results/JKF63_Sort_PolishedAssembly_by_funannotate/JKF63_Draft_Sort.fasta",
	output:
		"JKF63_Polished_Assembly.fasta",
	shell:
		"cp {input} {output}"
###
rule JKF63_Arrange_Assembly_by_RaGOO:
	input:
		reads= "JKF63_Short_reads.fastq",
		assembly= "JKF63_Polished_Assembly.fasta",
		ref= "EmonterogeiiLV88_Genome.fasta",
	output:
		"results/JKF63_Arrange_Assembly_by_RaGOO/ragoo_output/ragoo.fasta",
	threads: workflow.cores
	group: "Polishing_3"
	log:
		"results/log/JKF63_Arrange_Assembly_by_RaGOO.log"
	benchmark:
		"results/benchmarks/JKF63_Arrange_Assembly_by_RaGOO.tsv"
	conda:
		"envs/ragoo.yaml"
	shell:
		"rm -rf ragoo_output; ragoo.py -t {threads} -C -T sr -b -g 0 -R {input.reads} {input.assembly} {input.ref}; cp -r ragoo_output results/JKF63_Arrange_Assembly_by_RaGOO; rm -rf JKF63_Short_reads.fastq; rm -rf JKF63_Polished_Assembly.fasta; rm -rf ragoo_output"
###
rule JKF63_Clean_PolishedAssembly_by_funannotate:
	input:
		"results/JKF63_Arrange_Assembly_by_RaGOO/ragoo_output/ragoo.fasta",
	output:
		"results/JKF63_Clean_PolishedAssembly_by_funannotate/funannotate_Assembly.fasta",
	threads: workflow.cores	
	group: "Polishing_3"	
	params:
		"--exhaustive --minlen 50 --pident 10",
	log:
		"results/log/Clean_JKF63_Clean_PolishedAssembly_by_funannotate.log"
	benchmark:
		"results/benchmarks/Clean_JKF63_Clean_PolishedAssembly_by_funannotate.tsv"
	container: "docker://nextgenusfs/funannotate:latest"
	shell:
		"docker run -u 0 -v $(pwd):/home/data -w /home/data nextgenusfs/funannotate:latest funannotate clean {params} -i {input} -o {output}"
###
rule JKF63_Final_assembly:
	input:
		"results/JKF63_Clean_PolishedAssembly_by_funannotate/funannotate_Assembly.fasta",
	output:
		"results/Final_assembly/JKF63_Genome.fasta",
	shell:
		"cp {input} {output}; sed -i 's/_RaGOO//g' {output}; sed -i 's/EmoLV88/JKF63/g' {output}"
###
# Assessing Quality and Completeness:
###
checkpoint JKF63_QC_Assemblies_by_QUAST:
	input:
		"results/JKF63_Assembly_LongRead_by_Flye/assembly.fasta",
		"results/JKF63_Assembly_LongRead_by_Flye/Flye_consensus.fasta",
		"results/JKF63_Polish_Assembly_by_Pilon/pilon.fasta",
		"results/JKF63_Polish_Assembly_by_Pilon/Pilon_consensus.fasta",
		"results/JKF63_Arrange_Assembly_by_RaGOO/ragoo_output/ragoo.fasta",
		"results/JKF63_Clean_PolishedAssembly_by_funannotate/funannotate_Assembly.fasta",
		"results/Final_assembly/JKF63_Genome.fasta",
	output:
		"results/JKF63_QC_Assemblies_by_QUAST/report.html",
	params:
		"--eukaryote --circos",
	threads: workflow.cores
	group: "Quality_Assessment"
	log:
		"results/log/JKF63_QC_Assemblies_by_QUAST.log"
	benchmark:
		"results/benchmarks/JKF63_QC_Assemblies_by_QUAST.tsv"
	container:
		"docker://quay.io/biocontainers/quast:5.0.1--py27pl526ha92aebf_0"
	#conda:
		#"envs/quast.yaml"
	shell:
		"docker run -v $(pwd):/home/data -w /home/data quay.io/biocontainers/quast:5.0.1--py27pl526ha92aebf_0 quast.py {params} --threads {threads} --output-dir results/JKF63_QC_Assemblies_by_QUAST {input} 2>&1 | tee {log}" 
###