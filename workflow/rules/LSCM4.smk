###
rule LSCM4_Assembly_LongRead_by_Flye:
	input:
		long_1= "data/LSCM4/SRR13558782.fastq",
	output:
		"results/LSCM4_Assembly_LongRead_by_Flye/assembly.fasta",
	threads: workflow.cores
	group: "Assembly"
	log:
		"results/log/LSCM4_Assembly_LongRead_by_Flye.log"
	benchmark:
		"results/benchmarks/LSCM4_Assembly_LongRead_by_Flye.tsv"
	conda:
		"envs/Flye.yaml"
	shell:
		"flye --nano-raw {input.long_1} --genome-size 35m --threads {threads} -o results/LSCM4_Assembly_LongRead_by_Flye 2>&1 | tee {log}"
###
rule LSCM4_Index_Assembly_by_Samtools:
	input:
		"results/LSCM4_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		"results/LSCM4_Assembly_LongRead_by_Flye/assembly.fasta.fai",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/LSCM4_Index_Assembly_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools faidx {input} > {output}"
###	
rule LSCM4_Map_ShortReads_by_minimap2:
	input:
		PE1_1="data/LSCM4/SRR13558779_1.fastq",
		PE1_2="data/LSCM4/SRR13558779_2.fastq",
		PE2_1="data/LSCM4/SRR13558777_1.fastq",
		PE2_2="data/LSCM4/SRR13558777_2.fastq",
		PE3_1="data/LSCM4/SRR13558775_1.fastq",
		PE3_2="data/LSCM4/SRR13558775_2.fastq",
		PE4_1="data/LSCM4/SRR13558774_1.fastq",
		PE4_2="data/LSCM4/SRR13558774_2.fastq",
		PE5_1="data/LSCM4/SRR13558780_1.fastq",
		PE5_2="data/LSCM4/SRR13558780_2.fastq",
		PE6_1="data/LSCM4/SRR13558776_1.fastq",
		PE6_2="data/LSCM4/SRR13558776_2.fastq",
		PE7_1="data/LSCM4/SRR13558778_1.fastq",
		PE7_2="data/LSCM4/SRR13558778_2.fastq",
		PE8_1="data/LSCM4/SRR13558781_1.fastq",
		PE8_2="data/LSCM4/SRR13558781_2.fastq",
		assembly="results/LSCM4_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		sam1="results/LSCM4_Assembly_LongRead_by_Flye/PE1.sam",
		sam2="results/LSCM4_Assembly_LongRead_by_Flye/PE2.sam",
		sam3="results/LSCM4_Assembly_LongRead_by_Flye/PE3.sam",
		sam4="results/LSCM4_Assembly_LongRead_by_Flye/PE4.sam",
		sam5="results/LSCM4_Assembly_LongRead_by_Flye/PE5.sam",
		sam6="results/LSCM4_Assembly_LongRead_by_Flye/PE6.sam",
		sam7="results/LSCM4_Assembly_LongRead_by_Flye/PE7.sam",
		sam8="results/LSCM4_Assembly_LongRead_by_Flye/PE8.sam",
	threads: workflow.cores
	group: "Polishing_1"
	log:
		"results/log/LSCM4_Map_ShortReads_by_minimap2.log"
	benchmark:
		"results/benchmarks/LSCM4_Map_ShortReads_by_minimap2.tsv"
	conda:
		"envs/minimap2.yaml"
	shell:
		"""
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE1_1} {input.PE1_2} > {output.sam1}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE2_1} {input.PE2_2} > {output.sam2}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE3_1} {input.PE3_2} > {output.sam3}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE4_1} {input.PE4_2} > {output.sam4}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE5_1} {input.PE5_2} > {output.sam5}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE6_1} {input.PE6_2} > {output.sam6}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE7_1} {input.PE7_2} > {output.sam7}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE8_1} {input.PE8_2} > {output.sam8}
		"""
###
rule LSCM4_View_ShortReads_by_Samtools:
	input:
		sam1="results/LSCM4_Assembly_LongRead_by_Flye/PE1.sam",
		sam2="results/LSCM4_Assembly_LongRead_by_Flye/PE2.sam",
		sam3="results/LSCM4_Assembly_LongRead_by_Flye/PE3.sam",
		sam4="results/LSCM4_Assembly_LongRead_by_Flye/PE4.sam",
		sam5="results/LSCM4_Assembly_LongRead_by_Flye/PE5.sam",
		sam6="results/LSCM4_Assembly_LongRead_by_Flye/PE6.sam",
		sam7="results/LSCM4_Assembly_LongRead_by_Flye/PE7.sam",
		sam8="results/LSCM4_Assembly_LongRead_by_Flye/PE8.sam",
		assembly="results/LSCM4_Assembly_LongRead_by_Flye/assembly.fasta",
		idx="results/LSCM4_Assembly_LongRead_by_Flye/assembly.fasta.fai",
	output:
		bam1="results/LSCM4_Assembly_LongRead_by_Flye/PE1.bam",
		bam2="results/LSCM4_Assembly_LongRead_by_Flye/PE2.bam",
		bam3="results/LSCM4_Assembly_LongRead_by_Flye/PE3.bam",
		bam4="results/LSCM4_Assembly_LongRead_by_Flye/PE4.bam",
		bam5="results/LSCM4_Assembly_LongRead_by_Flye/PE5.bam",
		bam6="results/LSCM4_Assembly_LongRead_by_Flye/PE6.bam",
		bam7="results/LSCM4_Assembly_LongRead_by_Flye/PE7.bam",
		bam8="results/LSCM4_Assembly_LongRead_by_Flye/PE8.bam",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/LSCM4_View_ShortReads_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools view -@ {threads} -bt {input.idx} {input.sam1} > {output.bam1}
		samtools view -@ {threads} -bt {input.idx} {input.sam2} > {output.bam2}
		samtools view -@ {threads} -bt {input.idx} {input.sam3} > {output.bam3}
		samtools view -@ {threads} -bt {input.idx} {input.sam4} > {output.bam4}
		samtools view -@ {threads} -bt {input.idx} {input.sam5} > {output.bam5}
		samtools view -@ {threads} -bt {input.idx} {input.sam6} > {output.bam6}
		samtools view -@ {threads} -bt {input.idx} {input.sam7} > {output.bam7}
		samtools view -@ {threads} -bt {input.idx} {input.sam8} > {output.bam8}
		"""
###
rule LSCM4_Sort_ShortReads_by_Samtools:
	input:
		bam1="results/LSCM4_Assembly_LongRead_by_Flye/PE1.bam",
		bam2="results/LSCM4_Assembly_LongRead_by_Flye/PE2.bam",
		bam3="results/LSCM4_Assembly_LongRead_by_Flye/PE3.bam",
		bam4="results/LSCM4_Assembly_LongRead_by_Flye/PE4.bam",
		bam5="results/LSCM4_Assembly_LongRead_by_Flye/PE5.bam",
		bam6="results/LSCM4_Assembly_LongRead_by_Flye/PE6.bam",
		bam7="results/LSCM4_Assembly_LongRead_by_Flye/PE7.bam",
		bam8="results/LSCM4_Assembly_LongRead_by_Flye/PE8.bam",
	output:
		bam_sort1="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort1.bam",
		bam_sort2="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort2.bam",
		bam_sort3="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort3.bam",
		bam_sort4="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort4.bam",
		bam_sort5="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort5.bam",
		bam_sort6="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort6.bam",
		bam_sort7="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort7.bam",
		bam_sort8="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort8.bam",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/LSCM4_Sort_ShortReads_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools sort -@ {threads} {input.bam1} -o {output.bam_sort1}
		samtools sort -@ {threads} {input.bam2} -o {output.bam_sort2}
		samtools sort -@ {threads} {input.bam3} -o {output.bam_sort3}
		samtools sort -@ {threads} {input.bam4} -o {output.bam_sort4}
		samtools sort -@ {threads} {input.bam5} -o {output.bam_sort5}
		samtools sort -@ {threads} {input.bam6} -o {output.bam_sort6}
		samtools sort -@ {threads} {input.bam7} -o {output.bam_sort7}
		samtools sort -@ {threads} {input.bam8} -o {output.bam_sort8}
		"""
###
rule LSCM4_Index_ShortReads_by_Samtools:
	input:
		bam_sort1="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort1.bam",
		bam_sort2="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort2.bam",
		bam_sort3="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort3.bam",
		bam_sort4="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort4.bam",
		bam_sort5="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort5.bam",
		bam_sort6="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort6.bam",
		bam_sort7="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort7.bam",
		bam_sort8="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort8.bam",
	output:
		bam_idx1="results/LSCM4_Assembly_LongRead_by_Flye/PE_idx1.bai",
		bam_idx2="results/LSCM4_Assembly_LongRead_by_Flye/PE_idx2.bai",
		bam_idx3="results/LSCM4_Assembly_LongRead_by_Flye/PE_idx3.bai",
		bam_idx4="results/LSCM4_Assembly_LongRead_by_Flye/PE_idx4.bai",
		bam_idx5="results/LSCM4_Assembly_LongRead_by_Flye/PE_idx5.bai",
		bam_idx6="results/LSCM4_Assembly_LongRead_by_Flye/PE_idx6.bai",
		bam_idx7="results/LSCM4_Assembly_LongRead_by_Flye/PE_idx7.bai",
		bam_idx8="results/LSCM4_Assembly_LongRead_by_Flye/PE_idx8.bai",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/LSCM4_Index_ShortReads_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools index -@ {threads} {input.bam_sort1} > {output.bam_idx1}
		samtools index -@ {threads} {input.bam_sort2} > {output.bam_idx2}
		samtools index -@ {threads} {input.bam_sort3} > {output.bam_idx3}
		samtools index -@ {threads} {input.bam_sort4} > {output.bam_idx4}
		samtools index -@ {threads} {input.bam_sort5} > {output.bam_idx5}
		samtools index -@ {threads} {input.bam_sort6} > {output.bam_idx6}
		samtools index -@ {threads} {input.bam_sort7} > {output.bam_idx7}
		samtools index -@ {threads} {input.bam_sort8} > {output.bam_idx8}
		"""
###
rule LSCM4_Mpileup_Assembly_by_BCFtools:
	input:
		bam_sort1="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort1.bam",
		bam_sort2="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort2.bam",
		bam_sort3="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort3.bam",
		bam_sort4="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort4.bam",
		bam_sort5="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort5.bam",
		bam_sort6="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort6.bam",
		bam_sort7="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort7.bam",
		bam_sort8="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort8.bam",
		assembly="results/LSCM4_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		PE="results/LSCM4_Assembly_LongRead_by_Flye/PE.vcf.gz",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/LSCM4_Mpileup_Assembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools mpileup --threads {threads} -Ou -d 1000000 -f {input.assembly} {input.bam_sort1} {input.bam_sort2} {input.bam_sort3} {input.bam_sort4} {input.bam_sort5} {input.bam_sort6} {input.bam_sort7} {input.bam_sort8} | bcftools call -mv -Oz -o {output.PE}"
###
rule LSCM4_Norm_Assembly_by_BCFtools:
	input:
		PE="results/LSCM4_Assembly_LongRead_by_Flye/PE.vcf.gz",
		assembly="results/LSCM4_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		"results/LSCM4_Assembly_LongRead_by_Flye/PE.norm.bcf",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/LSCM4_Norm_Assembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools norm --threads {threads} -f {input.assembly} {input.PE} -Ob -o {output}"
###
rule LSCM4_Filter_Assembly_by_BCFtools:
	input:
		PE="results/LSCM4_Assembly_LongRead_by_Flye/PE.norm.bcf",
	output:
		"results/LSCM4_Assembly_LongRead_by_Flye/PE.norm.flt-indels.bcf",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/LSCM4_Filter_Assembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools filter --threads {threads} --IndelGap 5 {input.PE} -Ob -o {output}"
###
rule LSCM4_Index_Assembly_by_BCFtools:
	input:
		"results/LSCM4_Assembly_LongRead_by_Flye/PE.norm.flt-indels.bcf",
	output:
		"results/LSCM4_Assembly_LongRead_by_Flye/PE.norm.flt-indels.bcf.csi",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/LSCM4_Index_Assembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools index --threads {threads} {input} -o {output}"
###
rule LSCM4_Consensus_Assembly_by_BCFtools:
	input:
		"results/LSCM4_Assembly_LongRead_by_Flye/PE.norm.flt-indels.bcf.csi",
		PE="results/LSCM4_Assembly_LongRead_by_Flye/PE.norm.flt-indels.bcf",
		assembly="results/LSCM4_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		"results/LSCM4_Assembly_LongRead_by_Flye/Flye_consensus.fasta",
	threads: workflow.cores
	group: "Polishing_1"
	benchmark:
		"results/benchmarks/LSCM4_Consensus_Assembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools consensus -f {input.assembly} {input.PE} > {output}"
###
rule LSCM4_Polish_Assembly_by_Pilon:
	input:
		bam_idx1="results/LSCM4_Assembly_LongRead_by_Flye/PE_idx1.bai",
		bam_idx2="results/LSCM4_Assembly_LongRead_by_Flye/PE_idx2.bai",
		bam_idx3="results/LSCM4_Assembly_LongRead_by_Flye/PE_idx3.bai",
		bam_idx4="results/LSCM4_Assembly_LongRead_by_Flye/PE_idx4.bai",
		bam_idx5="results/LSCM4_Assembly_LongRead_by_Flye/PE_idx5.bai",
		bam_idx6="results/LSCM4_Assembly_LongRead_by_Flye/PE_idx6.bai",
		bam_idx7="results/LSCM4_Assembly_LongRead_by_Flye/PE_idx7.bai",
		bam_idx8="results/LSCM4_Assembly_LongRead_by_Flye/PE_idx8.bai",
		bam_sort1="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort1.bam",
		bam_sort2="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort2.bam",
		bam_sort3="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort3.bam",
		bam_sort4="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort4.bam",
		bam_sort5="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort5.bam",
		bam_sort6="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort6.bam",
		bam_sort7="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort7.bam",
		bam_sort8="results/LSCM4_Assembly_LongRead_by_Flye/PE_sort8.bam",
		assembly="results/LSCM4_Assembly_LongRead_by_Flye/assembly.fasta",
	output:
		"results/LSCM4_Polish_Assembly_by_Pilon/pilon.fasta",
	threads: workflow.cores
	group: "Polishing_2"
	log:
		"results/log/LSCM4_Polish_Assembly_by_Pilon.log"
	benchmark:
		"results/benchmarks/LSCM4_Polish_Assembly_by_Pilon.tsv"
	conda:
		"envs/pilon.yaml"
	shell:	
		"pilon --threads {threads} -Xmx300000M --genome {input.assembly} --bam {input.bam_sort1} --bam {input.bam_sort2} --bam {input.bam_sort3} --bam {input.bam_sort4} --bam {input.bam_sort5} --bam {input.bam_sort6} --bam {input.bam_sort7} --bam {input.bam_sort8} --outdir results/LSCM4_Polish_Assembly_by_Pilon 2>&1 | tee {log}"
###
rule LSCM4_index_PolishedAssembly_LSCM4: 
	input:
		"results/LSCM4_Polish_Assembly_by_Pilon/pilon.fasta",
	output:
		"results/LSCM4_Polish_Assembly_by_Pilon/pilon.fai",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/LSCM4_index_PolishedAssembly_LSCM4.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools faidx {input} > {output}"
###	
checkpoint LSCM4_Map_PolishedAssembly_by_minimap2:
	input:
		PE1_1="data/LSCM4/SRR13558779_1.fastq",
		PE1_2="data/LSCM4/SRR13558779_2.fastq",
		PE2_1="data/LSCM4/SRR13558777_1.fastq",
		PE2_2="data/LSCM4/SRR13558777_2.fastq",
		PE3_1="data/LSCM4/SRR13558775_1.fastq",
		PE3_2="data/LSCM4/SRR13558775_2.fastq",
		PE4_1="data/LSCM4/SRR13558774_1.fastq",
		PE4_2="data/LSCM4/SRR13558774_2.fastq",
		PE5_1="data/LSCM4/SRR13558780_1.fastq",
		PE5_2="data/LSCM4/SRR13558780_2.fastq",
		PE6_1="data/LSCM4/SRR13558776_1.fastq",
		PE6_2="data/LSCM4/SRR13558776_2.fastq",
		PE7_1="data/LSCM4/SRR13558778_1.fastq",
		PE7_2="data/LSCM4/SRR13558778_2.fastq",
		PE8_1="data/LSCM4/SRR13558781_1.fastq",
		PE8_2="data/LSCM4/SRR13558781_2.fastq",
		assembly="results/LSCM4_Polish_Assembly_by_Pilon/pilon.fasta",
	output:
		sam1="results/LSCM4_Polish_Assembly_by_Pilon/PE1.sam",
		sam2="results/LSCM4_Polish_Assembly_by_Pilon/PE2.sam",
		sam3="results/LSCM4_Polish_Assembly_by_Pilon/PE3.sam",
		sam4="results/LSCM4_Polish_Assembly_by_Pilon/PE4.sam",
		sam5="results/LSCM4_Polish_Assembly_by_Pilon/PE5.sam",
		sam6="results/LSCM4_Polish_Assembly_by_Pilon/PE6.sam",
		sam7="results/LSCM4_Polish_Assembly_by_Pilon/PE7.sam",
		sam8="results/LSCM4_Polish_Assembly_by_Pilon/PE8.sam",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/LSCM4_Map_PolishedAssembly_by_minimap2.tsv"
	conda:
		"envs/minimap2.yaml"
	shell:
		"""
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE1_1} {input.PE1_2} > {output.sam1}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE2_1} {input.PE2_2} > {output.sam2}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE3_1} {input.PE3_2} > {output.sam3}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE4_1} {input.PE4_2} > {output.sam4}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE5_1} {input.PE5_2} > {output.sam5}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE6_1} {input.PE6_2} > {output.sam6}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE7_1} {input.PE7_2} > {output.sam7}
		minimap2 -t {threads} -ax sr {input.assembly} {input.PE8_1} {input.PE8_2} > {output.sam8}
		"""
###
rule LSCM4_View_PolishedAssembly_by_Samtools:
	input:	
		sam1="results/LSCM4_Polish_Assembly_by_Pilon/PE1.sam",
		sam2="results/LSCM4_Polish_Assembly_by_Pilon/PE2.sam",
		sam3="results/LSCM4_Polish_Assembly_by_Pilon/PE3.sam",
		sam4="results/LSCM4_Polish_Assembly_by_Pilon/PE4.sam",
		sam5="results/LSCM4_Polish_Assembly_by_Pilon/PE5.sam",
		sam6="results/LSCM4_Polish_Assembly_by_Pilon/PE6.sam",
		sam7="results/LSCM4_Polish_Assembly_by_Pilon/PE7.sam",
		sam8="results/LSCM4_Polish_Assembly_by_Pilon/PE8.sam",
		assembly="results/LSCM4_Polish_Assembly_by_Pilon/pilon.fasta",
		idx="results/LSCM4_Polish_Assembly_by_Pilon/pilon.fai",
	output:
		bam1="results/LSCM4_Polish_Assembly_by_Pilon/PE1.bam",
		bam2="results/LSCM4_Polish_Assembly_by_Pilon/PE2.bam",
		bam3="results/LSCM4_Polish_Assembly_by_Pilon/PE3.bam",
		bam4="results/LSCM4_Polish_Assembly_by_Pilon/PE4.bam",
		bam5="results/LSCM4_Polish_Assembly_by_Pilon/PE5.bam",
		bam6="results/LSCM4_Polish_Assembly_by_Pilon/PE6.bam",
		bam7="results/LSCM4_Polish_Assembly_by_Pilon/PE7.bam",
		bam8="results/LSCM4_Polish_Assembly_by_Pilon/PE8.bam",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/LSCM4_View_PolishedAssembly_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools view -@ {threads} -bt {input.idx} {input.sam1} > {output.bam1}
		samtools view -@ {threads} -bt {input.idx} {input.sam2} > {output.bam2}
		samtools view -@ {threads} -bt {input.idx} {input.sam3} > {output.bam3}
		samtools view -@ {threads} -bt {input.idx} {input.sam4} > {output.bam4}
		samtools view -@ {threads} -bt {input.idx} {input.sam5} > {output.bam5}
		samtools view -@ {threads} -bt {input.idx} {input.sam6} > {output.bam6}
		samtools view -@ {threads} -bt {input.idx} {input.sam7} > {output.bam7}
		samtools view -@ {threads} -bt {input.idx} {input.sam8} > {output.bam8}
		"""
###
rule LSCM4_Sort_PolishedAssembly_by_Samtools:
	input:
		bam1="results/LSCM4_Polish_Assembly_by_Pilon/PE1.bam",
		bam2="results/LSCM4_Polish_Assembly_by_Pilon/PE2.bam",
		bam3="results/LSCM4_Polish_Assembly_by_Pilon/PE3.bam",
		bam4="results/LSCM4_Polish_Assembly_by_Pilon/PE4.bam",
		bam5="results/LSCM4_Polish_Assembly_by_Pilon/PE5.bam",
		bam6="results/LSCM4_Polish_Assembly_by_Pilon/PE6.bam",
		bam7="results/LSCM4_Polish_Assembly_by_Pilon/PE7.bam",
		bam8="results/LSCM4_Polish_Assembly_by_Pilon/PE8.bam",
	output:
		bam_sort1="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort1.bam",
		bam_sort2="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort2.bam",
		bam_sort3="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort3.bam",
		bam_sort4="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort4.bam",
		bam_sort5="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort5.bam",
		bam_sort6="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort6.bam",
		bam_sort7="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort7.bam",
		bam_sort8="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort8.bam",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/LSCM4_Sort_PolishedAssembly_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools sort -@ {threads} {input.bam1} -o {output.bam_sort1}
		samtools sort -@ {threads} {input.bam2} -o {output.bam_sort2}
		samtools sort -@ {threads} {input.bam3} -o {output.bam_sort3}
		samtools sort -@ {threads} {input.bam4} -o {output.bam_sort4}
		samtools sort -@ {threads} {input.bam5} -o {output.bam_sort5}
		samtools sort -@ {threads} {input.bam6} -o {output.bam_sort6}
		samtools sort -@ {threads} {input.bam7} -o {output.bam_sort7}
		samtools sort -@ {threads} {input.bam8} -o {output.bam_sort8}
		"""
###
rule LSCM4_Index_PolishedAssembly_by_Samtools:
	input:
		bam_sort1="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort1.bam",
		bam_sort2="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort2.bam",
		bam_sort3="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort3.bam",
		bam_sort4="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort4.bam",
		bam_sort5="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort5.bam",
		bam_sort6="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort6.bam",
		bam_sort7="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort7.bam",
		bam_sort8="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort8.bam",
	output:
		bam_idx1="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort1.bai",
		bam_idx2="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort2.bai",
		bam_idx3="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort3.bai",
		bam_idx4="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort4.bai",
		bam_idx5="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort5.bai",
		bam_idx6="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort6.bai",
		bam_idx7="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort7.bai",
		bam_idx8="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort8.bai",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/LSCM4_Index_PolishedAssembly_by_Samtools.tsv"
	conda:
		"envs/samtools.yaml"
	shell:
		"""
		samtools index -@ {threads} {input.bam_sort1} > {output.bam_idx1}
		samtools index -@ {threads} {input.bam_sort2} > {output.bam_idx2}
		samtools index -@ {threads} {input.bam_sort3} > {output.bam_idx3}
		samtools index -@ {threads} {input.bam_sort4} > {output.bam_idx4}
		samtools index -@ {threads} {input.bam_sort5} > {output.bam_idx5}
		samtools index -@ {threads} {input.bam_sort6} > {output.bam_idx6}
		samtools index -@ {threads} {input.bam_sort7} > {output.bam_idx7}
		samtools index -@ {threads} {input.bam_sort8} > {output.bam_idx8}
		"""
###
rule LSCM4_Mpileup_PolishedAssembly_by_BCFtools:
	input:
		bam_sort1="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort1.bam",
		bam_sort2="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort2.bam",
		bam_sort3="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort3.bam",
		bam_sort4="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort4.bam",
		bam_sort5="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort5.bam",
		bam_sort6="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort6.bam",
		bam_sort7="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort7.bam",
		bam_sort8="results/LSCM4_Polish_Assembly_by_Pilon/PE_sort8.bam",
		assembly="results/LSCM4_Polish_Assembly_by_Pilon/pilon.fasta",
	output:
		PE="results/LSCM4_Polish_Assembly_by_Pilon/PE.vcf.gz",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/LSCM4_Mpileup_PolishedAssembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools mpileup --threads {threads} -Ou -d 1000000 -f {input.assembly} {input.bam_sort1} {input.bam_sort2} {input.bam_sort3} {input.bam_sort4} {input.bam_sort5} {input.bam_sort6} {input.bam_sort7} {input.bam_sort8} | bcftools call -mv -Oz -o {output.PE}"
###
rule LSCM4_Norm_PolishedAssembly_by_BCFtools:
	input:
		PE="results/LSCM4_Polish_Assembly_by_Pilon/PE.vcf.gz",
		assembly="results/LSCM4_Polish_Assembly_by_Pilon/pilon.fasta",
	output:
		"results/LSCM4_Polish_Assembly_by_Pilon/PE.norm.bcf",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/LSCM4_Norm_PolishedAssembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools norm --threads {threads} -f {input.assembly} {input.PE} -Ob -o {output}"
###
rule LSCM4_Filter_PolishedAssembly_by_BCFtools:
	input:
		PE="results/LSCM4_Polish_Assembly_by_Pilon/PE.norm.bcf",
	output:
		"results/LSCM4_Polish_Assembly_by_Pilon/PE.norm.flt-indels.bcf",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/LSCM4_Filter_PolishedAssembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools filter --threads {threads} --IndelGap 5 {input.PE} -Ob -o {output}"
###
rule LSCM4_Index_PolishedAssembly_by_BCFtools:
	input:
		"results/LSCM4_Polish_Assembly_by_Pilon/PE.norm.flt-indels.bcf",
	output:
		"results/LSCM4_Polish_Assembly_by_Pilon/PE.norm.flt-indels.bcf.csi",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/LSCM4_Index_PolishedAssembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools index --threads {threads} {input} -o {output}"
###
rule LSCM4_Consensus_PolishedAssembly_by_BCFtools:
	input:
		"results/LSCM4_Polish_Assembly_by_Pilon/PE.norm.flt-indels.bcf.csi",
		PE="results/LSCM4_Polish_Assembly_by_Pilon/PE.norm.flt-indels.bcf",
		assembly="results/LSCM4_Polish_Assembly_by_Pilon/pilon.fasta",
	output:
		"results/LSCM4_Polish_Assembly_by_Pilon/Pilon_consensus.fasta",
	threads: workflow.cores
	group: "Polishing_2"
	benchmark:
		"results/benchmarks/LSCM4_Consensus_PolishedAssembly_by_BCFtools.tsv"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools consensus -f {input.assembly} {input.PE} > {output}"
###
rule LSCM4_aggregate_reads:
	input:
		PE1_1="data/LSCM4/SRR13558779_1.fastq",
		PE1_2="data/LSCM4/SRR13558779_2.fastq",
		PE2_1="data/LSCM4/SRR13558777_1.fastq",
		PE2_2="data/LSCM4/SRR13558777_2.fastq",
		PE3_1="data/LSCM4/SRR13558775_1.fastq",
		PE3_2="data/LSCM4/SRR13558775_2.fastq",
		PE4_1="data/LSCM4/SRR13558774_1.fastq",
		PE4_2="data/LSCM4/SRR13558774_2.fastq",
		PE5_1="data/LSCM4/SRR13558780_1.fastq",
		PE5_2="data/LSCM4/SRR13558780_2.fastq",
		PE6_1="data/LSCM4/SRR13558776_1.fastq",
		PE6_2="data/LSCM4/SRR13558776_2.fastq",
		PE7_1="data/LSCM4/SRR13558778_1.fastq",
		PE7_2="data/LSCM4/SRR13558778_2.fastq",
		PE8_1="data/LSCM4/SRR13558781_1.fastq",
		PE8_2="data/LSCM4/SRR13558781_2.fastq",
	output:
		"data/LSCM4/LSCM4_Short_reads.fastq"
	group: "Polishing_3"
	shell:
		"cat {input.PE1_1} {input.PE1_2} {input.PE2_1} {input.PE2_2} {input.PE3_1} {input.PE3_2} {input.PE4_1} {input.PE4_2} {input.PE5_1} {input.PE5_2} {input.PE6_1} {input.PE6_2} {input.PE7_1} {input.PE7_2} {input.PE8_1} {input.PE8_2} > {output}"

###
rule LSCM4_pre_RaGOO_Short_reads:
	input:
		"data/LSCM4/LSCM4_Short_reads.fastq",
	output:
		"LSCM4_Short_reads.fastq",
	shell:
		"cp {input} {output}"
###
rule LSCM4_Sort_PolishedAssembly_by_funannotate:
	input:
		"results/LSCM4_Polish_Assembly_by_Pilon/Pilon_consensus.fasta",
	output:
		"results/LSCM4_Sort_PolishedAssembly_by_funannotate/LSCM4_Draft_Sort.fasta",
	threads: workflow.cores	
	group: "Polishing_3"
	log:
		"results/log/LSCM4_Sort_PolishedAssembly_by_funannotate.log"
	benchmark:
		"results/benchmarks/LSCM4_Sort_PolishedAssembly_by_funannotate.tsv"
	container:
		"docker://nextgenusfs/funannotate:latest"
	shell:
		"docker run -u 0 -v $(pwd):/home/data -w /home/data nextgenusfs/funannotate:latest funannotate sort -i {input} -o {output}"
###
rule LSCM4_pre_RaGOO_Polished_Assembly:
	input:
		"results/LSCM1_Sort_PolishedAssembly_by_funannotate/LSCM1_Draft_Sort.fasta",
	output:
		"LSCM4_Polished_Assembly.fasta",
	shell:
		"cp {input} {output}"
###
rule LSCM4_Arrange_Assembly_by_RaGOO:
	input:
		reads= "LSCM4_Short_reads.fastq",
		assembly= "LSCM4_Polished_Assembly.fasta",
		ref= "LmajorFriedlin_Genome.fasta",
	output:
		"results/LSCM4_Arrange_Assembly_by_RaGOO/ragoo_output/ragoo.fasta",
	threads: workflow.cores
	group: "Polishing_3"
	log:
		"results/log/LSCM4_Arrange_Assembly_by_RaGOO.log"
	benchmark:
		"results/benchmarks/LSCM4_Arrange_Assembly_by_RaGOO.tsv"
	conda:
		"envs/ragoo.yaml"
	shell:
		"rm -rf ragoo_output; ragoo.py -t {threads} -C -T sr -b -g 0 -R {input.reads} {input.assembly} {input.ref}; cp -r ragoo_output results/LSCM4_Arrange_Assembly_by_RaGOO; rm -rf LSCM4_Short_reads.fastq; rm -rf LSCM4_Polished_Assembly.fasta; rm -rf ragoo_output"
###
rule LSCM4_Clean_PolishedAssembly_by_funannotate:
	input:
		"results/LSCM4_Arrange_Assembly_by_RaGOO/ragoo_output/ragoo.fasta",
	output:
		"results/LSCM4_Clean_PolishedAssembly_by_funannotate/funannotate_Assembly.fasta",
	threads: workflow.cores	
	group: "Polishing_3"	
	params:
		"--exhaustive --minlen 50 --pident 10",
	log:
		"results/log/Clean_LSCM4_Clean_PolishedAssembly_by_funannotate.log"
	benchmark:
		"results/benchmarks/Clean_LSCM4_Clean_PolishedAssembly_by_funannotate.tsv"
	container: "docker://nextgenusfs/funannotate:latest"
	shell:
		"docker run -u 0 -v $(pwd):/home/data -w /home/data nextgenusfs/funannotate:latest funannotate clean {params} -i {input} -o {output}"
###
rule LSCM4_Final_assembly:
	input:
		"results/LSCM4_Clean_PolishedAssembly_by_funannotate/funannotate_Assembly.fasta",
	output:
		"results/Final_assembly/LSCM4_Genome.fasta",
	shell:
		"cp {input} {output}; sed -i 's/_RaGOO//g' {output}; sed -i 's/LmjF/LSCM4/g' {output}"
###
# Assessing Quality and Completeness:
###
checkpoint LSCM4_QC_Assemblies_by_QUAST:
	input:
		"results/LSCM4_Assembly_LongRead_by_Flye/assembly.fasta",
		"results/LSCM4_Assembly_LongRead_by_Flye/Flye_consensus.fasta",
		"results/LSCM4_Polish_Assembly_by_Pilon/pilon.fasta",
		"results/LSCM4_Polish_Assembly_by_Pilon/Pilon_consensus.fasta",
		"results/LSCM4_Arrange_Assembly_by_RaGOO/ragoo_output/ragoo.fasta",
		"results/LSCM4_Clean_PolishedAssembly_by_funannotate/funannotate_Assembly.fasta",
		"results/Final_assembly/LSCM4_Genome.fasta",
	output:
		"results/LSCM4_QC_Assemblies_by_QUAST/report.html",
	params:
		"--eukaryote --circos",
	threads: workflow.cores
	group: "Quality_Assessment"
	log:
		"results/log/LSCM4_QC_Assemblies_by_QUAST.log"
	benchmark:
		"results/benchmarks/LSCM4_QC_Assemblies_by_QUAST.tsv"
	container:
		"docker://quay.io/biocontainers/quast:5.0.1--py27pl526ha92aebf_0"
	#conda:
		#"envs/quast.yaml"
	shell:
		"docker run -v $(pwd):/home/data -w /home/data quay.io/biocontainers/quast:5.0.1--py27pl526ha92aebf_0 quast.py {params} --threads {threads} --output-dir results/LSCM4_QC_Assemblies_by_QUAST {input} 2>&1 | tee {log}" 
###