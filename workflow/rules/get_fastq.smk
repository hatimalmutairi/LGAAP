###
rule LSCM1_Download_Short_Read:
	output:
		"data/LSCM1/SRR13558784_1.fastq",
		"data/LSCM1/SRR13558784_2.fastq",
		"data/LSCM1/SRR13558785_1.fastq",
		"data/LSCM1/SRR13558785_2.fastq",
		"data/LSCM1/SRR13558792_1.fastq",
		"data/LSCM1/SRR13558792_2.fastq",
	threads: workflow.cores
	benchmark:
		"results/benchmarks/LSCM1_Download_Short_Read.tsv"
	conda:
		"envs/sra-tools.yaml"
	shell:
		"fasterq-dump -t data/temp --split-files --outdir data/LSCM1/ SRR13558784 SRR13558785 SRR13558792"
###
rule LSCM4_Download_Short_Read:
	output:
		"data/LSCM4/SRR13558774_1.fastq",
		"data/LSCM4/SRR13558774_2.fastq",
		"data/LSCM4/SRR13558775_1.fastq",
		"data/LSCM4/SRR13558775_2.fastq",
		"data/LSCM4/SRR13558776_1.fastq",
		"data/LSCM4/SRR13558776_2.fastq",
		"data/LSCM4/SRR13558777_1.fastq",
		"data/LSCM4/SRR13558777_2.fastq",
		"data/LSCM4/SRR13558778_1.fastq",
		"data/LSCM4/SRR13558778_2.fastq",
		"data/LSCM4/SRR13558779_1.fastq",
		"data/LSCM4/SRR13558779_2.fastq",
		"data/LSCM4/SRR13558780_1.fastq",
		"data/LSCM4/SRR13558780_2.fastq",
		"data/LSCM4/SRR13558781_1.fastq",
		"data/LSCM4/SRR13558781_2.fastq",
	threads: workflow.cores
	benchmark:
		"results/benchmarks/LSCM4_Download_Short_Read.tsv"
	conda:
		"envs/sra-tools.yaml"
	shell:
		"fasterq-dump -t data/temp --split-files --outdir data/LSCM4/ SRR13558774 SRR13558775 SRR13558776 SRR13558777 SRR13558778 SRR13558779 SRR13558780 SRR13558781"
###
rule CUR178_Download_Short_Read:
	output:
		"data/CUR178/SRR13558795_1.fastq",
		"data/CUR178/SRR13558795_2.fastq",
		"data/CUR178/SRR13558796_1.fastq",
		"data/CUR178/SRR13558796_2.fastq",
		"data/CUR178/SRR13558797_1.fastq",
		"data/CUR178/SRR13558797_2.fastq",
	threads: workflow.cores
	benchmark:
		"results/benchmarks/CUR178_Download_Short_Read.tsv"
	conda:
		"envs/sra-tools.yaml"
	shell:
		"fasterq-dump -t data/temp --split-files --outdir data/CUR178/ SRR13558795 SRR13558796 SRR13558797"
###
rule GH5_Download_Short_Read:
	output:
		"data/GH5/SRR13558800_1.fastq",
		"data/GH5/SRR13558800_2.fastq",
		"data/GH5/SRR13558801_1.fastq",
		"data/GH5/SRR13558801_2.fastq",
		"data/GH5/SRR13558802_1.fastq",
		"data/GH5/SRR13558802_2.fastq",
		"data/GH5/SRR13558803_1.fastq",
		"data/GH5/SRR13558803_2.fastq",
		"data/GH5/SRR13558804_1.fastq",
		"data/GH5/SRR13558804_2.fastq",
	threads: workflow.cores
	benchmark:
		"results/benchmarks/GH5_Download_Short_Read.tsv"
	conda:
		"envs/sra-tools.yaml"
	shell:
		"fasterq-dump -t data/temp --split-files --outdir data/GH5/ SRR13558800 SRR13558801 SRR13558802 SRR13558803 SRR13558804"
###
rule JKF63_Download_Short_Read:
	output:
		"data/JKF63/SRR13558754_1.fastq",
		"data/JKF63/SRR13558754_2.fastq",
		"data/JKF63/SRR13558755_1.fastq",
		"data/JKF63/SRR13558755_2.fastq",
		"data/JKF63/SRR13558756_1.fastq",
		"data/JKF63/SRR13558756_2.fastq",
	threads: workflow.cores
	benchmark:
		"results/benchmarks/JKF63_Download_Short_Read.tsv"
	conda:
		"envs/sra-tools.yaml"
	shell:
		"fasterq-dump -t data/temp --split-files --outdir data/JKF63/ SRR13558754 SRR13558755 SRR13558756"
###
rule JIQ42_Download_Short_Read:
	output:
		"data/JIQ42/SRR13558764_1.fastq",
		"data/JIQ42/SRR13558764_2.fastq",
		"data/JIQ42/SRR13558765_1.fastq",
		"data/JIQ42/SRR13558765_2.fastq",
		"data/JIQ42/SRR13558766_1.fastq",
		"data/JIQ42/SRR13558766_2.fastq",
	threads: workflow.cores
	benchmark:
		"results/benchmarks/JIQ42_Download_Short_Read.tsv"
	conda:
		"envs/sra-tools.yaml"
	shell:
		"fasterq-dump -t data/temp --split-files --outdir data/JIQ42/ SRR13558764 SRR13558765 SRR13558766"
###
rule LSCM1_Download_Long_Read:
	output:
		"data/LSCM1/SRR13558786.fastq",
		"data/LSCM1/SRR13558788.fastq",
		"data/LSCM1/SRR13558790.fastq",
		"data/LSCM1/SRR13558793.fastq",
	threads: workflow.cores
	benchmark:
		"results/benchmarks/LSCM1_Download_Long_Read.tsv"
	conda:
		"envs/sra-tools.yaml"
	shell:
		"fasterq-dump -t data/temp --outdir data/LSCM1/ SRR13558786 SRR13558788 SRR13558790 SRR13558793"
###
rule LSCM4_Download_Long_Read:
	output:
		"data/LSCM4/SRR13558782.fastq",
	threads: workflow.cores
	benchmark:
		"results/benchmarks/LSCM4_Download_Long_Read.tsv"
	conda:
		"envs/sra-tools.yaml"
	shell:
		"fasterq-dump -t data/temp --outdir data/LSCM4/ SRR13558782"
###
rule CUR178_Download_Long_Read:
	output:
		"data/CUR178/SRR13558798.fastq",
	threads: workflow.cores
	benchmark:
		"results/benchmarks/CUR178_Download_Long_Read.tsv"
	conda:
		"envs/sra-tools.yaml"
	shell:
		"fasterq-dump -t data/temp --outdir data/CUR178/ SRR13558798"
###
rule GH5_Download_Long_Read:
	output:
		"data/GH5/SRR13558805.fastq",
	threads: workflow.cores
	benchmark:
		"results/benchmarks/GH5_Download_Long_Read.tsv"
	conda:
		"envs/sra-tools.yaml"
	shell:
		"fasterq-dump -t data/temp --outdir data/GH5/ SRR13558805"
###
rule JKF63_Download_Long_Read:
	output:
		"data/JKF63/SRR13558757.fastq",
	threads: workflow.cores
	benchmark:
		"results/benchmarks/JKF63_Download_Long_Read.tsv"
	conda:
		"envs/sra-tools.yaml"
	shell:
		"fasterq-dump -t data/temp --outdir data/JKF63/ SRR13558757"
###
rule JIQ42_Download_Long_Read:
	output:
		"data/JIQ42/SRR13558767.fastq",
	threads: workflow.cores
	benchmark:
		"results/benchmarks/JIQ42_Download_Long_Read.tsv"
	conda:
		"envs/sra-tools.yaml"
	shell:
		"fasterq-dump -t data/temp --outdir data/JIQ42/ SRR13558767"
###
###
rule download_Evidance1:
	output:
		"data/ref/LmajorFriedlin_Genome.fasta",
	group: "Download"
	benchmark:
		"results/benchmarks/download_Evidance1.tsv"
	shell:
		"wget https://tritrypdb.org/common/downloads/release-47/LmajorFriedlin/fasta/data/TriTrypDB-47_LmajorFriedlin_Genome.fasta; cp TriTrypDB-47_LmajorFriedlin_Genome.fasta {output}"
rule pre_RaGOO_ref:
	input:
		"data/ref/LmajorFriedlin_Genome.fasta",
	output:
		"LmajorFriedlin_Genome.fasta",
	benchmark:
		"results/benchmarks/pre_RaGOO_ref.tsv"
	shell:
		"cp {input} {output}"
###
rule download_Evidance2:
	output:
		"data/ref/EmonterogeiiLV88_Genome.fasta",
	group: "Download"
	benchmark:
		"results/benchmarks/download_Evidance2.tsv"
	shell:
		"wget https://tritrypdb.org/common/downloads/release-47/EmonterogeiiLV88/fasta/data/TriTrypDB-47_EmonterogeiiLV88_Genome.fasta; cp TriTrypDB-47_EmonterogeiiLV88_Genome.fasta {output}"
###
rule pre_RaGOO_ref2:
	input:
		"data/ref/EmonterogeiiLV88_Genome.fasta",
	output:
		"EmonterogeiiLV88_Genome.fasta",
	benchmark:
		"results/benchmarks/pre_RaGOO_ref2.tsv"
	shell:
		"cp {input} {output}"
###
rule download_Annotation_Evidance:
	output:
		n1="data/Evidance/TriTrypDB-47_Linfantum_ESTs.fasta",
		n2="data/Evidance/TriTrypDB-47_Lamazonensis_ESTs.fasta",
		n3="data/Evidance/TriTrypDB-47_LpanamensisMHOMCOL81L13_AnnotatedTranscripts.fasta",
		n4="data/Evidance/TriTrypDB-47_LpanamensisMHOMCOL81L13_AnnotatedProteins.fasta",
		n5="data/Evidance/TriTrypDB-47_LpanamensisMHOMCOL81L13.gff",
		n6="data/Evidance/TriTrypDB-47_LaethiopicaL147_AnnotatedProteins.fasta",
		n7="data/Evidance/TriTrypDB-47_LaethiopicaL147_AnnotatedTranscripts.fasta",
		n8="data/Evidance/TriTrypDB-47_LaethiopicaL147.gff",
		n9="data/Evidance/TriTrypDB-47_Lmexicana_ESTs.fasta",
		n10="data/Evidance/TriTrypDB-47_LamazonensisMHOMBR71973M2269_AnnotatedTranscripts.fasta",
		n11="data/Evidance/TriTrypDB-47_LamazonensisMHOMBR71973M2269_AnnotatedProteins.fasta",
		n12="data/Evidance/TriTrypDB-47_LamazonensisMHOMBR71973M2269.gff",
		n13="data/Evidance/TriTrypDB-47_LdonovaniCL-SL_AnnotatedProteins.fasta",
		n14="data/Evidance/TriTrypDB-47_LdonovaniCL-SL_AnnotatedTranscripts.fasta",
		n15="data/Evidance/TriTrypDB-47_LdonovaniCL-SL.gff",
		n16="data/Evidance/TriTrypDB-47_LmexicanaMHOMGT2001U1103_AnnotatedProteins.fasta",
		n17="data/Evidance/TriTrypDB-47_LmexicanaMHOMGT2001U1103_AnnotatedTranscripts.fasta",
		n18="data/Evidance/TriTrypDB-47_LmexicanaMHOMGT2001U1103.gff",
		n19="data/Evidance/TriTrypDB-47_LturanicaLEM423_AnnotatedProteins.fasta",
		n20="data/Evidance/TriTrypDB-47_LturanicaLEM423_AnnotatedTranscripts.fasta",
		n21="data/Evidance/TriTrypDB-47_LturanicaLEM423.gff",
		n22="data/Evidance/TriTrypDB-47_LinfantumJPCM5_AnnotatedTranscripts.fasta",
		n23="data/Evidance/TriTrypDB-47_LinfantumJPCM5_AnnotatedProteins.fasta",
		n24="data/Evidance/TriTrypDB-47_LinfantumJPCM5.gff",
		n25="data/Evidance/TriTrypDB-47_LbraziliensisMHOMBR75M2903_AnnotatedTranscripts.fasta",
		n26="data/Evidance/TriTrypDB-47_LbraziliensisMHOMBR75M2903_AnnotatedProteins.fasta",
		n27="data/Evidance/TriTrypDB-47_LbraziliensisMHOMBR75M2903.gff",
		n28="data/Evidance/TriTrypDB-47_LmajorSD75.1_AnnotatedTranscripts.fasta",
		n29="data/Evidance/TriTrypDB-47_LmajorSD75.1_AnnotatedProteins.fasta",
		n30="data/Evidance/TriTrypDB-47_LmajorSD75.1.gff",
		n31="data/Evidance/TriTrypDB-47_LbraziliensisMHOMBR75M2904_AnnotatedTranscripts.fasta",
		n32="data/Evidance/TriTrypDB-47_LbraziliensisMHOMBR75M2904_AnnotatedProteins.fasta",
		n33="data/Evidance/TriTrypDB-47_LbraziliensisMHOMBR75M2904.gff",
		n34="data/Evidance/TriTrypDB-47_LtarentolaeParrotTarII_AnnotatedProteins.fasta",
		n35="data/Evidance/TriTrypDB-47_LtarentolaeParrotTarII_AnnotatedTranscripts.fasta",
		n36="data/Evidance/TriTrypDB-47_LtarentolaeParrotTarII.gff",
		n37="data/Evidance/TriTrypDB-47_LtropicaL590_AnnotatedTranscripts.fasta",
		n38="data/Evidance/TriTrypDB-47_LtropicaL590_AnnotatedProteins.fasta",
		n39="data/Evidance/TriTrypDB-47_LtropicaL590.gff",
		n40="data/Evidance/TriTrypDB-47_LbraziliensisMHOMBR75M2904_2019_AnnotatedProteins.fasta",
		n41="data/Evidance/TriTrypDB-47_LbraziliensisMHOMBR75M2904_2019_AnnotatedTranscripts.fasta",
		n42="data/Evidance/TriTrypDB-47_LbraziliensisMHOMBR75M2904_2019.gff",
		n43="data/Evidance/TriTrypDB-47_LgerbilliLEM452_AnnotatedTranscripts.fasta",
		n44="data/Evidance/TriTrypDB-47_LgerbilliLEM452_AnnotatedProteins.fasta",
		n45="data/Evidance/TriTrypDB-47_LgerbilliLEM452.gff",
		n46="data/Evidance/TriTrypDB-47_LdonovaniBPK282A1_AnnotatedTranscripts.fasta",
		n47="data/Evidance/TriTrypDB-47_LdonovaniBPK282A1_AnnotatedProteins.fasta",
		n48="data/Evidance/TriTrypDB-47_LdonovaniBPK282A1.gff",
		n49="data/Evidance/TriTrypDB-47_Lmajor_ESTs.fasta",
		n50="data/Evidance/TriTrypDB-47_LmajorLV39c5_AnnotatedTranscripts.fasta",
		n51="data/Evidance/TriTrypDB-47_LmajorLV39c5_AnnotatedProteins.fasta",
		n52="data/Evidance/TriTrypDB-47_LmajorLV39c5.gff",
		n53="data/Evidance/TriTrypDB-47_Lbraziliensis_ESTs.fasta",
		n54="data/Evidance/TriTrypDB-47_LdonovaniLV9_AnnotatedTranscripts.fasta",
		n55="data/Evidance/TriTrypDB-47_LdonovaniLV9_AnnotatedProteins.fasta",
		n56="data/Evidance/TriTrypDB-47_LdonovaniLV9.gff",
		n57="data/Evidance/TriTrypDB-47_LmajorFriedlin_AnnotatedTranscripts.fasta",
		n58="data/Evidance/TriTrypDB-47_LmajorFriedlin_AnnotatedProteins.fasta",
		n59="data/Evidance/TriTrypDB-47_LmajorFriedlin.gff",
		n60="data/Evidance/TriTrypDB-47_LenriettiiLEM3045_AnnotatedProteins.fasta",
		n61="data/Evidance/TriTrypDB-47_LenriettiiLEM3045_AnnotatedTranscripts.fasta",
		n62="data/Evidance/TriTrypDB-47_LenriettiiLEM3045.gff",
		n63="data/Evidance/TriTrypDB-47_LpanamensisMHOMPA94PSC1_AnnotatedTranscripts.fasta",
		n64="data/Evidance/TriTrypDB-47_LpanamensisMHOMPA94PSC1_AnnotatedProteins.fasta",
		n65="data/Evidance/TriTrypDB-47_LpanamensisMHOMPA94PSC1.gff",
		n66="data/Evidance/TriTrypDB-47_LarabicaLEM1108_AnnotatedProteins.fasta",
		n67="data/Evidance/TriTrypDB-47_LarabicaLEM1108_AnnotatedTranscripts.fasta",
		n68="data/Evidance/TriTrypDB-47_LarabicaLEM1108.gff",
	group: "Download"
	benchmark:
		"results/benchmarks/download_Annotation_Evidance.tsv"
	shell:
		"""
		wget https://tritrypdb.org/common/downloads/release-47/Linfantum/fasta/TriTrypDB-47_Linfantum_ESTs.fasta -O {output.n1}
		wget https://tritrypdb.org/common/downloads/release-47/Lamazonensis/fasta/TriTrypDB-47_Lamazonensis_ESTs.fasta -O {output.n2}
		wget https://tritrypdb.org/common/downloads/release-47/LpanamensisMHOMCOL81L13/fasta/data/TriTrypDB-47_LpanamensisMHOMCOL81L13_AnnotatedTranscripts.fasta -O {output.n3}
		wget https://tritrypdb.org/common/downloads/release-47/LpanamensisMHOMCOL81L13/fasta/data/TriTrypDB-47_LpanamensisMHOMCOL81L13_AnnotatedProteins.fasta -O {output.n4}
		wget https://tritrypdb.org/common/downloads/release-47/LpanamensisMHOMCOL81L13/gff/data/TriTrypDB-47_LpanamensisMHOMCOL81L13.gff -O {output.n5}
		wget https://tritrypdb.org/common/downloads/release-47/LaethiopicaL147/fasta/data/TriTrypDB-47_LaethiopicaL147_AnnotatedProteins.fasta -O {output.n6}
		wget https://tritrypdb.org/common/downloads/release-47/LaethiopicaL147/fasta/data/TriTrypDB-47_LaethiopicaL147_AnnotatedTranscripts.fasta -O {output.n7}
		wget https://tritrypdb.org/common/downloads/release-47/LaethiopicaL147/gff/data/TriTrypDB-47_LaethiopicaL147.gff -O {output.n8}
		wget https://tritrypdb.org/common/downloads/release-47/Lmexicana/fasta/TriTrypDB-47_Lmexicana_ESTs.fasta -O {output.n9}
		wget https://tritrypdb.org/common/downloads/release-47/LamazonensisMHOMBR71973M2269/fasta/data/TriTrypDB-47_LamazonensisMHOMBR71973M2269_AnnotatedTranscripts.fasta -O {output.n10}
		wget https://tritrypdb.org/common/downloads/release-47/LamazonensisMHOMBR71973M2269/fasta/data/TriTrypDB-47_LamazonensisMHOMBR71973M2269_AnnotatedProteins.fasta -O {output.n11}
		wget https://tritrypdb.org/common/downloads/release-47/LamazonensisMHOMBR71973M2269/gff/data/TriTrypDB-47_LamazonensisMHOMBR71973M2269.gff -O {output.n12}
		wget https://tritrypdb.org/common/downloads/release-47/LdonovaniCL-SL/fasta/data/TriTrypDB-47_LdonovaniCL-SL_AnnotatedProteins.fasta -O {output.n13}
		wget https://tritrypdb.org/common/downloads/release-47/LdonovaniCL-SL/fasta/data/TriTrypDB-47_LdonovaniCL-SL_AnnotatedTranscripts.fasta -O {output.n14}
		wget https://tritrypdb.org/common/downloads/release-47/LdonovaniCL-SL/gff/data/TriTrypDB-47_LdonovaniCL-SL.gff -O {output.n15}
		wget https://tritrypdb.org/common/downloads/release-47/LmexicanaMHOMGT2001U1103/fasta/data/TriTrypDB-47_LmexicanaMHOMGT2001U1103_AnnotatedProteins.fasta -O {output.n16}
		wget https://tritrypdb.org/common/downloads/release-47/LmexicanaMHOMGT2001U1103/fasta/data/TriTrypDB-47_LmexicanaMHOMGT2001U1103_AnnotatedTranscripts.fasta -O {output.n17}
		wget https://tritrypdb.org/common/downloads/release-47/LmexicanaMHOMGT2001U1103/gff/data/TriTrypDB-47_LmexicanaMHOMGT2001U1103.gff -O {output.n18}
		wget https://tritrypdb.org/common/downloads/release-47/LturanicaLEM423/fasta/data/TriTrypDB-47_LturanicaLEM423_AnnotatedProteins.fasta -O {output.n19}
		wget https://tritrypdb.org/common/downloads/release-47/LturanicaLEM423/fasta/data/TriTrypDB-47_LturanicaLEM423_AnnotatedTranscripts.fasta -O {output.n20}
		wget https://tritrypdb.org/common/downloads/release-47/LturanicaLEM423/gff/data/TriTrypDB-47_LturanicaLEM423.gff -O {output.n21}
		wget https://tritrypdb.org/common/downloads/release-47/LinfantumJPCM5/fasta/data/TriTrypDB-47_LinfantumJPCM5_AnnotatedTranscripts.fasta -O {output.n22}
		wget https://tritrypdb.org/common/downloads/release-47/LinfantumJPCM5/fasta/data/TriTrypDB-47_LinfantumJPCM5_AnnotatedProteins.fasta -O {output.n23}
		wget https://tritrypdb.org/common/downloads/release-47/LinfantumJPCM5/gff/data/TriTrypDB-47_LinfantumJPCM5.gff -O {output.n24}
		wget https://tritrypdb.org/common/downloads/release-47/LbraziliensisMHOMBR75M2903/fasta/data/TriTrypDB-47_LbraziliensisMHOMBR75M2903_AnnotatedTranscripts.fasta -O {output.n25}
		wget https://tritrypdb.org/common/downloads/release-47/LbraziliensisMHOMBR75M2903/fasta/data/TriTrypDB-47_LbraziliensisMHOMBR75M2903_AnnotatedProteins.fasta -O {output.n26}
		wget https://tritrypdb.org/common/downloads/release-47/LbraziliensisMHOMBR75M2903/gff/data/TriTrypDB-47_LbraziliensisMHOMBR75M2903.gff -O {output.n27}
		wget https://tritrypdb.org/common/downloads/release-47/LmajorSD75.1/fasta/data/TriTrypDB-47_LmajorSD75.1_AnnotatedTranscripts.fasta -O {output.n28}
		wget https://tritrypdb.org/common/downloads/release-47/LmajorSD75.1/fasta/data/TriTrypDB-47_LmajorSD75.1_AnnotatedProteins.fasta -O {output.n29}
		wget https://tritrypdb.org/common/downloads/release-47/LmajorSD75.1/gff/data/TriTrypDB-47_LmajorSD75.1.gff -O {output.n30}
		wget https://tritrypdb.org/common/downloads/release-47/LbraziliensisMHOMBR75M2904/fasta/data/TriTrypDB-47_LbraziliensisMHOMBR75M2904_AnnotatedTranscripts.fasta -O {output.n31}
		wget https://tritrypdb.org/common/downloads/release-47/LbraziliensisMHOMBR75M2904/fasta/data/TriTrypDB-47_LbraziliensisMHOMBR75M2904_AnnotatedProteins.fasta -O {output.n32}
		wget https://tritrypdb.org/common/downloads/release-47/LbraziliensisMHOMBR75M2904/gff/data/TriTrypDB-47_LbraziliensisMHOMBR75M2904.gff -O {output.n33}
		wget https://tritrypdb.org/common/downloads/release-47/LtarentolaeParrotTarII/fasta/data/TriTrypDB-47_LtarentolaeParrotTarII_AnnotatedProteins.fasta -O {output.n34}
		wget https://tritrypdb.org/common/downloads/release-47/LtarentolaeParrotTarII/fasta/data/TriTrypDB-47_LtarentolaeParrotTarII_AnnotatedTranscripts.fasta -O {output.n35}
		wget https://tritrypdb.org/common/downloads/release-47/LtarentolaeParrotTarII/gff/data/TriTrypDB-47_LtarentolaeParrotTarII.gff -O {output.n36}
		wget https://tritrypdb.org/common/downloads/release-47/LtropicaL590/fasta/data/TriTrypDB-47_LtropicaL590_AnnotatedTranscripts.fasta -O {output.n37}
		wget https://tritrypdb.org/common/downloads/release-47/LtropicaL590/fasta/data/TriTrypDB-47_LtropicaL590_AnnotatedProteins.fasta -O {output.n38}
		wget https://tritrypdb.org/common/downloads/release-47/LtropicaL590/gff/data/TriTrypDB-47_LtropicaL590.gff -O {output.n39}
		wget https://tritrypdb.org/common/downloads/release-47/LbraziliensisMHOMBR75M2904_2019/fasta/data/TriTrypDB-47_LbraziliensisMHOMBR75M2904_2019_AnnotatedProteins.fasta -O {output.n40}
		wget https://tritrypdb.org/common/downloads/release-47/LbraziliensisMHOMBR75M2904_2019/fasta/data/TriTrypDB-47_LbraziliensisMHOMBR75M2904_2019_AnnotatedTranscripts.fasta -O {output.n41}
		wget https://tritrypdb.org/common/downloads/release-47/LbraziliensisMHOMBR75M2904_2019/gff/data/TriTrypDB-47_LbraziliensisMHOMBR75M2904_2019.gff -O {output.n42}
		wget https://tritrypdb.org/common/downloads/release-47/LgerbilliLEM452/fasta/data/TriTrypDB-47_LgerbilliLEM452_AnnotatedTranscripts.fasta -O {output.n43}
		wget https://tritrypdb.org/common/downloads/release-47/LgerbilliLEM452/fasta/data/TriTrypDB-47_LgerbilliLEM452_AnnotatedProteins.fasta -O {output.n44}
		wget https://tritrypdb.org/common/downloads/release-47/LgerbilliLEM452/gff/data/TriTrypDB-47_LgerbilliLEM452.gff -O {output.n45}
		wget https://tritrypdb.org/common/downloads/release-47/LdonovaniBPK282A1/fasta/data/TriTrypDB-47_LdonovaniBPK282A1_AnnotatedTranscripts.fasta -O {output.n46}
		wget https://tritrypdb.org/common/downloads/release-47/LdonovaniBPK282A1/fasta/data/TriTrypDB-47_LdonovaniBPK282A1_AnnotatedProteins.fasta -O {output.n47}
		wget https://tritrypdb.org/common/downloads/release-47/LdonovaniBPK282A1/gff/data/TriTrypDB-47_LdonovaniBPK282A1.gff -O {output.n48}
		wget https://tritrypdb.org/common/downloads/release-47/Lmajor/fasta/TriTrypDB-47_Lmajor_ESTs.fasta -O {output.n49}
		wget https://tritrypdb.org/common/downloads/release-47/LmajorLV39c5/fasta/data/TriTrypDB-47_LmajorLV39c5_AnnotatedTranscripts.fasta -O {output.n50}
		wget https://tritrypdb.org/common/downloads/release-47/LmajorLV39c5/fasta/data/TriTrypDB-47_LmajorLV39c5_AnnotatedProteins.fasta -O {output.n51}
		wget https://tritrypdb.org/common/downloads/release-47/LmajorLV39c5/gff/data/TriTrypDB-47_LmajorLV39c5.gff -O {output.n52}
		wget https://tritrypdb.org/common/downloads/release-47/Lbraziliensis/fasta/TriTrypDB-47_Lbraziliensis_ESTs.fasta -O {output.n53}
		wget https://tritrypdb.org/common/downloads/release-47/LdonovaniLV9/fasta/data/TriTrypDB-47_LdonovaniLV9_AnnotatedTranscripts.fasta -O {output.n54}
		wget https://tritrypdb.org/common/downloads/release-47/LdonovaniLV9/fasta/data/TriTrypDB-47_LdonovaniLV9_AnnotatedProteins.fasta -O {output.n55}
		wget https://tritrypdb.org/common/downloads/release-47/LdonovaniLV9/gff/data/TriTrypDB-47_LdonovaniLV9.gff -O {output.n56}
		wget https://tritrypdb.org/common/downloads/release-47/LmajorFriedlin/fasta/data/TriTrypDB-47_LmajorFriedlin_AnnotatedTranscripts.fasta -O {output.n57}
		wget https://tritrypdb.org/common/downloads/release-47/LmajorFriedlin/fasta/data/TriTrypDB-47_LmajorFriedlin_AnnotatedProteins.fasta -O {output.n58}
		wget https://tritrypdb.org/common/downloads/release-47/LmajorFriedlin/gff/data/TriTrypDB-47_LmajorFriedlin.gff -O {output.n59}
		wget https://tritrypdb.org/common/downloads/release-47/LenriettiiLEM3045/fasta/data/TriTrypDB-47_LenriettiiLEM3045_AnnotatedProteins.fasta -O {output.n60}
		wget https://tritrypdb.org/common/downloads/release-47/LenriettiiLEM3045/fasta/data/TriTrypDB-47_LenriettiiLEM3045_AnnotatedTranscripts.fasta -O {output.n61}
		wget https://tritrypdb.org/common/downloads/release-47/LenriettiiLEM3045/gff/data/TriTrypDB-47_LenriettiiLEM3045.gff -O {output.n62}
		wget https://tritrypdb.org/common/downloads/release-47/LpanamensisMHOMPA94PSC1/fasta/data/TriTrypDB-47_LpanamensisMHOMPA94PSC1_AnnotatedTranscripts.fasta -O {output.n63}
		wget https://tritrypdb.org/common/downloads/release-47/LpanamensisMHOMPA94PSC1/fasta/data/TriTrypDB-47_LpanamensisMHOMPA94PSC1_AnnotatedProteins.fasta -O {output.n64}
		wget https://tritrypdb.org/common/downloads/release-47/LpanamensisMHOMPA94PSC1/gff/data/TriTrypDB-47_LpanamensisMHOMPA94PSC1.gff -O {output.n65}
		wget https://tritrypdb.org/common/downloads/release-47/LarabicaLEM1108/fasta/data/TriTrypDB-47_LarabicaLEM1108_AnnotatedProteins.fasta -O {output.n66}
		wget https://tritrypdb.org/common/downloads/release-47/LarabicaLEM1108/fasta/data/TriTrypDB-47_LarabicaLEM1108_AnnotatedTranscripts.fasta -O {output.n67}
		wget https://tritrypdb.org/common/downloads/release-47/LarabicaLEM1108/gff/data/TriTrypDB-47_LarabicaLEM1108.gff -O {output.n68}
		"""
###