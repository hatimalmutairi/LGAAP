###
rule GH5_Assembly_LongRead_by_Flye:
	input:
		long_1= "data/GH5/SRR13558805.fastq",
	output:
		"results/GH5_Assembly_LongRead_by_Flye/assembly.fasta",
	threads: workflow.cores
	group: "Assembly"
	log:
		"results/log/GH5_Assembly_LongRead_by_Flye.log"
	benchmark:
		"results/benchmarks/GH5_Assembly_LongRead_by_Flye.tsv"
	conda:
		"envs/Flye.yaml"
	shell:
		"flye --nano-raw {input.long_1} --genome-size 35m --threads {threads} -o results/GH5_Assembly_LongRead_by_Flye 2>&1 | tee {log}"
###
rule GH5_Index_Assembly_by_Samtools:
	input:
		"results/GH5_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		"results/GH5_Assembly_LongRead_by_Flye/assembly.fasta.fai",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/GH5_Index_Assembly_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools faidx {input} > {output}"
###	
rule GH5_Map_ShortReads_by_minimap2:
	input:
		PE1_1="data/GH5/SRR13558800_1.fastq",
		PE1_2="data/GH5/SRR13558800_2.fastq",
		PE2_1="data/GH5/SRR13558802_1.fastq",
		PE2_2="data/GH5/SRR13558802_2.fastq",
		PE3_1="data/GH5/SRR13558803_1.fastq",
		PE3_2="data/GH5/SRR13558803_2.fastq",
		PE4_1="data/GH5/SRR13558804_1.fastq",
		PE4_2="data/GH5/SRR13558804_2.fastq",
		PE5_1="data/GH5/SRR13558801_1.fastq",
		PE5_2="data/GH5/SRR13558801_2.fastq",
		assembly="results/GH5_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		sam1="results/GH5_Assembly_LongRead_by_Flye/PE1.sam",
		sam2="results/GH5_Assembly_LongRead_by_Flye/PE2.sam",
		sam3="results/GH5_Assembly_LongRead_by_Flye/PE3.sam",
		sam4="results/GH5_Assembly_LongRead_by_Flye/PE4.sam",
		sam5="results/GH5_Assembly_LongRead_by_Flye/PE5.sam",
	threads: workflow.cores
	group: "Polishing_1"
	log:
		"results/log/GH5_Map_ShortReads_by_minimap2.log"
	benchmark:
		"results/benchmarks/GH5_Map_ShortReads_by_minimap2.tsv"
	conda:
		"envs/minimap2.yaml"
	shell:
		"""
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE1_1} {input.PE1_2} > {output.sam1}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE2_1} {input.PE2_2} > {output.sam2}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE3_1} {input.PE3_2} > {output.sam3}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE4_1} {input.PE4_2} > {output.sam4}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE5_1} {input.PE5_2} > {output.sam5}
		"""
###
rule GH5_View_ShortReads_by_Samtools:
	input:
		sam1="results/GH5_Assembly_LongRead_by_Flye/PE1.sam",
		sam2="results/GH5_Assembly_LongRead_by_Flye/PE2.sam",
		sam3="results/GH5_Assembly_LongRead_by_Flye/PE3.sam",
		sam4="results/GH5_Assembly_LongRead_by_Flye/PE4.sam",
		sam5="results/GH5_Assembly_LongRead_by_Flye/PE5.sam",
		assembly="results/GH5_Assembly_LongRead_by_Flye/assembly.fasta",
		idx="results/GH5_Assembly_LongRead_by_Flye/assembly.fasta.fai",
	output:
		bam1="results/GH5_Assembly_LongRead_by_Flye/PE1.bam",
		bam2="results/GH5_Assembly_LongRead_by_Flye/PE2.bam",
		bam3="results/GH5_Assembly_LongRead_by_Flye/PE3.bam",
		bam4="results/GH5_Assembly_LongRead_by_Flye/PE4.bam",
		bam5="results/GH5_Assembly_LongRead_by_Flye/PE5.bam",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/GH5_View_ShortReads_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools view -@ {threads} -bt {input.idx} {input.sam1} > {output.bam1}
		samtools view -@ {threads} -bt {input.idx} {input.sam2} > {output.bam2}
		samtools view -@ {threads} -bt {input.idx} {input.sam3} > {output.bam3}
		samtools view -@ {threads} -bt {input.idx} {input.sam4} > {output.bam4}
		samtools view -@ {threads} -bt {input.idx} {input.sam5} > {output.bam5}
		"""
###
rule GH5_Sort_ShortReads_by_Samtools:
	input:
		bam1="results/GH5_Assembly_LongRead_by_Flye/PE1.bam",
		bam2="results/GH5_Assembly_LongRead_by_Flye/PE2.bam",
		bam3="results/GH5_Assembly_LongRead_by_Flye/PE3.bam",
		bam4="results/GH5_Assembly_LongRead_by_Flye/PE4.bam",
		bam5="results/GH5_Assembly_LongRead_by_Flye/PE5.bam",
	output:
		bam_sort1="results/GH5_Assembly_LongRead_by_Flye/PE_sort1.bam",
		bam_sort2="results/GH5_Assembly_LongRead_by_Flye/PE_sort2.bam",
		bam_sort3="results/GH5_Assembly_LongRead_by_Flye/PE_sort3.bam",
		bam_sort4="results/GH5_Assembly_LongRead_by_Flye/PE_sort4.bam",
		bam_sort5="results/GH5_Assembly_LongRead_by_Flye/PE_sort5.bam",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/GH5_Sort_ShortReads_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools sort -@ {threads} {input.bam1} -o {output.bam_sort1}
		samtools sort -@ {threads} {input.bam2} -o {output.bam_sort2}
		samtools sort -@ {threads} {input.bam3} -o {output.bam_sort3}
		samtools sort -@ {threads} {input.bam4} -o {output.bam_sort4}
		samtools sort -@ {threads} {input.bam5} -o {output.bam_sort5}
		"""
###
rule GH5_Index_ShortReads_by_Samtools:
	input:
		bam_sort1="results/GH5_Assembly_LongRead_by_Flye/PE_sort1.bam",
		bam_sort2="results/GH5_Assembly_LongRead_by_Flye/PE_sort2.bam",
		bam_sort3="results/GH5_Assembly_LongRead_by_Flye/PE_sort3.bam",
		bam_sort4="results/GH5_Assembly_LongRead_by_Flye/PE_sort4.bam",
		bam_sort5="results/GH5_Assembly_LongRead_by_Flye/PE_sort5.bam",
	output:
		bam_idx1="results/GH5_Assembly_LongRead_by_Flye/PE_idx1.bai",
		bam_idx2="results/GH5_Assembly_LongRead_by_Flye/PE_idx2.bai",
		bam_idx3="results/GH5_Assembly_LongRead_by_Flye/PE_idx3.bai",
		bam_idx4="results/GH5_Assembly_LongRead_by_Flye/PE_idx4.bai",
		bam_idx5="results/GH5_Assembly_LongRead_by_Flye/PE_idx5.bai",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/GH5_Index_ShortReads_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools index -@ {threads} {input.bam_sort1} > {output.bam_idx1}
		samtools index -@ {threads} {input.bam_sort2} > {output.bam_idx2}
		samtools index -@ {threads} {input.bam_sort3} > {output.bam_idx3}
		samtools index -@ {threads} {input.bam_sort4} > {output.bam_idx4}
		samtools index -@ {threads} {input.bam_sort5} > {output.bam_idx5}
		"""
###
rule GH5_Mpileup_Assembly_by_BCFtools:
	input:
		bam_sort1="results/GH5_Assembly_LongRead_by_Flye/PE_sort1.bam",
		bam_sort2="results/GH5_Assembly_LongRead_by_Flye/PE_sort2.bam",
		bam_sort3="results/GH5_Assembly_LongRead_by_Flye/PE_sort3.bam",
		bam_sort4="results/GH5_Assembly_LongRead_by_Flye/PE_sort4.bam",
		bam_sort5="results/GH5_Assembly_LongRead_by_Flye/PE_sort5.bam",
		assembly="results/GH5_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		PE="results/GH5_Assembly_LongRead_by_Flye/PE.vcf.gz",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/GH5_Mpileup_Assembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools mpileup --threads {threads} -Ou -d 1000000 -f {input.assembly} {input.bam_sort1} {input.bam_sort2} {input.bam_sort3} {input.bam_sort4} {input.bam_sort5} | bcftools call -mv -Oz -o {output.PE}"
###
rule GH5_Norm_Assembly_by_BCFtools:
	input:
		PE="results/GH5_Assembly_LongRead_by_Flye/PE.vcf.gz",
		assembly="results/GH5_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		"results/GH5_Assembly_LongRead_by_Flye/PE.norm.bcf",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/GH5_Norm_Assembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools norm --threads {threads} -f {input.assembly} {input.PE} -Ob -o {output}"
###
rule GH5_Filter_Assembly_by_BCFtools:
	input:
		PE="results/GH5_Assembly_LongRead_by_Flye/PE.norm.bcf",
	output:
		"results/GH5_Assembly_LongRead_by_Flye/PE.norm.flt-indels.bcf",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/GH5_Filter_Assembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools filter --threads {threads} --IndelGap 5 {input.PE} -Ob -o {output}"
###
rule GH5_Index_Assembly_by_BCFtools:
	input:
		"results/GH5_Assembly_LongRead_by_Flye/PE.norm.flt-indels.bcf",
	output:
		"results/GH5_Assembly_LongRead_by_Flye/PE.norm.flt-indels.bcf.csi",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/GH5_Index_Assembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools index --threads {threads} {input} -o {output}"
###
rule GH5_Consensus_Assembly_by_BCFtools:
	input:
		"results/GH5_Assembly_LongRead_by_Flye/PE.norm.flt-indels.bcf.csi",
		PE="results/GH5_Assembly_LongRead_by_Flye/PE.norm.flt-indels.bcf",
		assembly="results/GH5_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		"results/GH5_Assembly_LongRead_by_Flye/Flye_consensus.fasta",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/GH5_Consensus_Assembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools consensus -f {input.assembly} {input.PE} > {output}"
###
rule GH5_Polish_Assembly_by_Pilon:
	input:
		bam_idx1="results/GH5_Assembly_LongRead_by_Flye/PE_idx1.bai",
		bam_idx2="results/GH5_Assembly_LongRead_by_Flye/PE_idx2.bai",
		bam_idx3="results/GH5_Assembly_LongRead_by_Flye/PE_idx3.bai",
		bam_idx4="results/GH5_Assembly_LongRead_by_Flye/PE_idx4.bai",
		bam_idx5="results/GH5_Assembly_LongRead_by_Flye/PE_idx5.bai",
		bam_sort1="results/GH5_Assembly_LongRead_by_Flye/PE_sort1.bam",
		bam_sort2="results/GH5_Assembly_LongRead_by_Flye/PE_sort2.bam",
		bam_sort3="results/GH5_Assembly_LongRead_by_Flye/PE_sort3.bam",
		bam_sort4="results/GH5_Assembly_LongRead_by_Flye/PE_sort4.bam",
		bam_sort5="results/GH5_Assembly_LongRead_by_Flye/PE_sort5.bam",
		assembly="results/GH5_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		"results/GH5_Polish_Assembly_by_Pilon/pilon.fasta",
	threads: workflow.cores
	group: "Polishing_2"
	log:
		"results/log/GH5_Polish_Assembly_by_Pilon.log"
	benchmark:
		"results/benchmarks/GH5_Polish_Assembly_by_Pilon.tsv"
	conda:
		"envs/pilon.yaml"
	shell:	
		"pilon --threads {threads} -Xmx300000M --genome {input.assembly} --bam {input.bam_sort1} --bam {input.bam_sort2} --bam {input.bam_sort3} --bam {input.bam_sort4} --bam {input.bam_sort5} --outdir results/GH5_Polish_Assembly_by_Pilon 2>&1 | tee {log}"
###
rule GH5_index_PolishedAssembly_GH5: 
	input:
		"results/GH5_Polish_Assembly_by_Pilon/pilon.fasta",
	output:
		"results/GH5_Polish_Assembly_by_Pilon/pilon.fai",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/GH5_index_PolishedAssembly_GH5.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools faidx {input} > {output}"
###	
checkpoint GH5_Map_PolishedAssembly_by_minimap2:
	input:
		PE1_1="data/GH5/SRR13558800_1.fastq",
		PE1_2="data/GH5/SRR13558800_2.fastq",
		PE2_1="data/GH5/SRR13558802_1.fastq",
		PE2_2="data/GH5/SRR13558802_2.fastq",
		PE3_1="data/GH5/SRR13558803_1.fastq",
		PE3_2="data/GH5/SRR13558803_2.fastq",
		PE4_1="data/GH5/SRR13558804_1.fastq",
		PE4_2="data/GH5/SRR13558804_2.fastq",
		PE5_1="data/GH5/SRR13558801_1.fastq",
		PE5_2="data/GH5/SRR13558801_2.fastq",
		assembly="results/GH5_Polish_Assembly_by_Pilon/pilon.fasta",
	output:
		sam1="results/GH5_Polish_Assembly_by_Pilon/PE1.sam",
		sam2="results/GH5_Polish_Assembly_by_Pilon/PE2.sam",
		sam3="results/GH5_Polish_Assembly_by_Pilon/PE3.sam",
		sam4="results/GH5_Polish_Assembly_by_Pilon/PE4.sam",
		sam5="results/GH5_Polish_Assembly_by_Pilon/PE5.sam",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/GH5_Map_PolishedAssembly_by_minimap2.tsv"
	conda:
		"envs/minimap2.yaml"
	shell:
		"""
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE1_1} {input.PE1_2} > {output.sam1}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE2_1} {input.PE2_2} > {output.sam2}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE3_1} {input.PE3_2} > {output.sam3}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE4_1} {input.PE4_2} > {output.sam4}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE5_1} {input.PE5_2} > {output.sam5}
		"""
###
rule GH5_View_PolishedAssembly_by_Samtools:
	input:	
		sam1="results/GH5_Polish_Assembly_by_Pilon/PE1.sam",
		sam2="results/GH5_Polish_Assembly_by_Pilon/PE2.sam",
		sam3="results/GH5_Polish_Assembly_by_Pilon/PE3.sam",
		sam4="results/GH5_Polish_Assembly_by_Pilon/PE4.sam",
		sam5="results/GH5_Polish_Assembly_by_Pilon/PE5.sam",
		assembly="results/GH5_Polish_Assembly_by_Pilon/pilon.fasta",
		idx="results/GH5_Polish_Assembly_by_Pilon/pilon.fai",
	output:
		bam1="results/GH5_Polish_Assembly_by_Pilon/PE1.bam",
		bam2="results/GH5_Polish_Assembly_by_Pilon/PE2.bam",
		bam3="results/GH5_Polish_Assembly_by_Pilon/PE3.bam",
		bam4="results/GH5_Polish_Assembly_by_Pilon/PE4.bam",
		bam5="results/GH5_Polish_Assembly_by_Pilon/PE5.bam",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/GH5_View_PolishedAssembly_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools view -@ {threads} -bt {input.idx} {input.sam1} > {output.bam1}
		samtools view -@ {threads} -bt {input.idx} {input.sam2} > {output.bam2}
		samtools view -@ {threads} -bt {input.idx} {input.sam3} > {output.bam3}
		samtools view -@ {threads} -bt {input.idx} {input.sam4} > {output.bam4}
		samtools view -@ {threads} -bt {input.idx} {input.sam5} > {output.bam5}
		"""
###
rule GH5_Sort_PolishedAssembly_by_Samtools:
	input:
		bam1="results/GH5_Polish_Assembly_by_Pilon/PE1.bam",
		bam2="results/GH5_Polish_Assembly_by_Pilon/PE2.bam",
		bam3="results/GH5_Polish_Assembly_by_Pilon/PE3.bam",
		bam4="results/GH5_Polish_Assembly_by_Pilon/PE4.bam",
		bam5="results/GH5_Polish_Assembly_by_Pilon/PE5.bam",
	output:
		bam_sort1="results/GH5_Polish_Assembly_by_Pilon/PE_sort1.bam",
		bam_sort2="results/GH5_Polish_Assembly_by_Pilon/PE_sort2.bam",
		bam_sort3="results/GH5_Polish_Assembly_by_Pilon/PE_sort3.bam",
		bam_sort4="results/GH5_Polish_Assembly_by_Pilon/PE_sort4.bam",
		bam_sort5="results/GH5_Polish_Assembly_by_Pilon/PE_sort5.bam",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/GH5_Sort_PolishedAssembly_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools sort -@ {threads} {input.bam1} -o {output.bam_sort1}
		samtools sort -@ {threads} {input.bam2} -o {output.bam_sort2}
		samtools sort -@ {threads} {input.bam3} -o {output.bam_sort3}
		samtools sort -@ {threads} {input.bam4} -o {output.bam_sort4}
		samtools sort -@ {threads} {input.bam5} -o {output.bam_sort5}
		"""
###
rule GH5_Index_PolishedAssembly_by_Samtools:
	input:
		bam_sort1="results/GH5_Polish_Assembly_by_Pilon/PE_sort1.bam",
		bam_sort2="results/GH5_Polish_Assembly_by_Pilon/PE_sort2.bam",
		bam_sort3="results/GH5_Polish_Assembly_by_Pilon/PE_sort3.bam",
		bam_sort4="results/GH5_Polish_Assembly_by_Pilon/PE_sort4.bam",
		bam_sort5="results/GH5_Polish_Assembly_by_Pilon/PE_sort5.bam",
	output:
		bam_idx1="results/GH5_Polish_Assembly_by_Pilon/PE_sort1.bai",
		bam_idx2="results/GH5_Polish_Assembly_by_Pilon/PE_sort2.bai",
		bam_idx3="results/GH5_Polish_Assembly_by_Pilon/PE_sort3.bai",
		bam_idx4="results/GH5_Polish_Assembly_by_Pilon/PE_sort4.bai",
		bam_idx5="results/GH5_Polish_Assembly_by_Pilon/PE_sort5.bai",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/GH5_Index_PolishedAssembly_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools index -@ {threads} {input.bam_sort1} > {output.bam_idx1}
		samtools index -@ {threads} {input.bam_sort2} > {output.bam_idx2}
		samtools index -@ {threads} {input.bam_sort3} > {output.bam_idx3}
		samtools index -@ {threads} {input.bam_sort4} > {output.bam_idx4}
		samtools index -@ {threads} {input.bam_sort5} > {output.bam_idx5}
		"""
###
rule GH5_Mpileup_PolishedAssembly_by_BCFtools:
	input:
		bam_sort1="results/GH5_Polish_Assembly_by_Pilon/PE_sort1.bam",
		bam_sort2="results/GH5_Polish_Assembly_by_Pilon/PE_sort2.bam",
		bam_sort3="results/GH5_Polish_Assembly_by_Pilon/PE_sort3.bam",
		bam_sort4="results/GH5_Polish_Assembly_by_Pilon/PE_sort4.bam",
		bam_sort5="results/GH5_Polish_Assembly_by_Pilon/PE_sort5.bam",
		assembly="results/GH5_Polish_Assembly_by_Pilon/pilon.fasta",
	output:
		PE="results/GH5_Polish_Assembly_by_Pilon/PE.vcf.gz",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/GH5_Mpileup_PolishedAssembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools mpileup --threads {threads} -Ou -d 1000000 -f {input.assembly} {input.bam_sort1} {input.bam_sort2} {input.bam_sort3} {input.bam_sort4} {input.bam_sort5} | bcftools call -mv -Oz -o {output.PE}"
###
rule GH5_Norm_PolishedAssembly_by_BCFtools:
	input:
		PE="results/GH5_Polish_Assembly_by_Pilon/PE.vcf.gz",
		assembly="results/GH5_Polish_Assembly_by_Pilon/pilon.fasta",
	output:
		"results/GH5_Polish_Assembly_by_Pilon/PE.norm.bcf",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/GH5_Norm_PolishedAssembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools norm --threads {threads} -f {input.assembly} {input.PE} -Ob -o {output}"
###
rule GH5_Filter_PolishedAssembly_by_BCFtools:
	input:
		PE="results/GH5_Polish_Assembly_by_Pilon/PE.norm.bcf",
	output:
		"results/GH5_Polish_Assembly_by_Pilon/PE.norm.flt-indels.bcf",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/GH5_Filter_PolishedAssembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools filter --threads {threads} --IndelGap 5 {input.PE} -Ob -o {output}"
###
rule GH5_Index_PolishedAssembly_by_BCFtools:
	input:
		"results/GH5_Polish_Assembly_by_Pilon/PE.norm.flt-indels.bcf",
	output:
		"results/GH5_Polish_Assembly_by_Pilon/PE.norm.flt-indels.bcf.csi",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/GH5_Index_PolishedAssembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools index --threads {threads} {input} -o {output}"
###
rule GH5_Consensus_PolishedAssembly_by_BCFtools:
	input:
		"results/GH5_Polish_Assembly_by_Pilon/PE.norm.flt-indels.bcf.csi",
		PE="results/GH5_Polish_Assembly_by_Pilon/PE.norm.flt-indels.bcf",
		assembly="results/GH5_Polish_Assembly_by_Pilon/pilon.fasta",
	output:
		"results/GH5_Polish_Assembly_by_Pilon/Pilon_consensus.fasta",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/GH5_Consensus_PolishedAssembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools consensus -f {input.assembly} {input.PE} > {output}"
###
rule GH5_aggregate_reads:
	input:
		PE1_1="data/GH5/SRR13558800_1.fastq",
		PE1_2="data/GH5/SRR13558800_2.fastq",
		PE2_1="data/GH5/SRR13558802_1.fastq",
		PE2_2="data/GH5/SRR13558802_2.fastq",
		PE3_1="data/GH5/SRR13558803_1.fastq",
		PE3_2="data/GH5/SRR13558803_2.fastq",
		PE4_1="data/GH5/SRR13558804_1.fastq",
		PE4_2="data/GH5/SRR13558804_2.fastq",
		PE5_1="data/GH5/SRR13558801_1.fastq",
		PE5_2="data/GH5/SRR13558801_2.fastq",
	output:
		"data/GH5/GH5_Short_reads.fastq"
	group: "Polishing_3"
	shell:
		"cat {input.PE1_1} {input.PE1_2} {input.PE2_1} {input.PE2_2} {input.PE3_1} {input.PE3_2} {input.PE4_1} {input.PE4_2} {input.PE5_1} {input.PE5_2} > {output}"
###
rule GH5_pre_RaGOO_Short_reads:
	input:
		"data/GH5/GH5_Short_reads.fastq",
	output:
		"GH5_Short_reads.fastq",
	shell:
		"cp {input} {output}"
###
rule GH5_Sort_PolishedAssembly_by_funannotate:
	input:
		"results/GH5_Polish_Assembly_by_Pilon/Pilon_consensus.fasta",
	output:
		"results/GH5_Sort_PolishedAssembly_by_funannotate/GH5_Draft_Sort.fasta",
	threads: workflow.cores	
	group: "Polishing_3"
	log:
		"results/log/GH5_Sort_PolishedAssembly_by_funannotate.log"
	benchmark:
		"results/benchmarks/GH5_Sort_PolishedAssembly_by_funannotate.tsv"
	container:
		"docker://nextgenusfs/funannotate:latest"
	shell:
		"docker run -u 0 -v $(pwd):/home/data -w /home/data nextgenusfs/funannotate:latest funannotate sort -i {input} -o {output}"
###
rule GH5_pre_RaGOO_Polished_Assembly:
	input:
		"results/GH5_Sort_PolishedAssembly_by_funannotate/GH5_Draft_Sort.fasta",
	output:
		"GH5_Polished_Assembly.fasta",
	shell:
		"cp {input} {output}"
###
rule GH5_Arrange_Assembly_by_RaGOO:
	input:
		reads= "GH5_Short_reads.fastq",
		assembly= "GH5_Polished_Assembly.fasta",
		ref= "LmajorFriedlin_Genome.fasta",
	output:
		"results/GH5_Arrange_Assembly_by_RaGOO/ragoo_output/ragoo.fasta",
	threads: workflow.cores
	group: "Polishing_3"
	log:
		"results/log/GH5_Arrange_Assembly_by_RaGOO.log"
	benchmark:
		"results/benchmarks/GH5_Arrange_Assembly_by_RaGOO.tsv"
	conda:
		"envs/ragoo.yaml"
	shell:
		"rm -rf ragoo_output; ragoo.py -t {threads} -C -T sr -b -g 0 -R {input.reads} {input.assembly} {input.ref}; cp -r ragoo_output results/GH5_Arrange_Assembly_by_RaGOO; rm -rf GH5_Short_reads.fastq; rm -rf GH5_Polished_Assembly.fasta; rm -rf ragoo_output"

###
rule GH5_Clean_PolishedAssembly_by_funannotate:
	input:
		"results/GH5_Arrange_Assembly_by_RaGOO/ragoo_output/ragoo.fasta",
	output:
		"results/GH5_Clean_PolishedAssembly_by_funannotate/funannotate_Assembly.fasta",
	threads: workflow.cores	
	group: "Polishing_3"	
	params:
		"--exhaustive --minlen 50 --pident 10",
	log:
		"results/log/Clean_GH5_Clean_PolishedAssembly_by_funannotate.log"
	benchmark:
		"results/benchmarks/Clean_GH5_Clean_PolishedAssembly_by_funannotate.tsv"
	container: "docker://nextgenusfs/funannotate:latest"
	shell:
		"docker run -u 0 -v $(pwd):/home/data -w /home/data nextgenusfs/funannotate:latest funannotate clean {params} -i {input} -o {output}"
###
rule GH5_Final_assembly:
	input:
		"results/GH5_Clean_PolishedAssembly_by_funannotate/funannotate_Assembly.fasta",
	output:
		"results/Final_assembly/GH5_Genome.fasta",
	shell:
		"cp {input} {output}; sed -i 's/_RaGOO//g' {output}; sed -i 's/LmjF/GH5/g' {output}"
###
# Assessing Quality and Completeness:
###
checkpoint GH5_QC_Assemblies_by_QUAST:
	input:
		"results/GH5_Assembly_LongRead_by_Flye/assembly.fasta",
		"results/GH5_Assembly_LongRead_by_Flye/Flye_consensus.fasta",
		"results/GH5_Polish_Assembly_by_Pilon/pilon.fasta",
		"results/GH5_Polish_Assembly_by_Pilon/Pilon_consensus.fasta",
		"results/GH5_Arrange_Assembly_by_RaGOO/ragoo_output/ragoo.fasta",
		"results/GH5_Clean_PolishedAssembly_by_funannotate/funannotate_Assembly.fasta",
		"results/Final_assembly/GH5_Genome.fasta",
	output:
		"results/GH5_QC_Assemblies_by_QUAST/report.html",
	params:
		"--eukaryote --circos",
	threads: workflow.cores
	group: "Quality_Assessment"
	log:
		"results/log/GH5_QC_Assemblies_by_QUAST.log"
	benchmark:
		"results/benchmarks/GH5_QC_Assemblies_by_QUAST.tsv"
	container:
		"docker://quay.io/biocontainers/quast:5.0.1--py27pl526ha92aebf_0"
	#conda:
		#"envs/quast.yaml"
	shell:
		"docker run -v $(pwd):/home/data -w /home/data quay.io/biocontainers/quast:5.0.1--py27pl526ha92aebf_0 quast.py {params} --threads {threads} --output-dir results/GH5_QC_Assemblies_by_QUAST {input} 2>&1 | tee {log}" 
###