###
rule QC1_by_GAAS:
	input:
		expand("results/Final_assembly/{sample}_Genome.fasta", sample=samples),
	output:
		expand("results/annotation/{sample}/QC1.stats.txt", sample=samples),
	conda:
		"envs/GAAS.yaml"
	shell: 
		"gaas_fasta_statistics.pl -f {input} -o {output}"
###
rule Download_Contamination_Database:
	output:
		"UniVec"
	shell:
		"wget https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec"
###
ext1=[".nhr",".nin",".nsq"]
ext2=["",".nhr",".nin",".nsq"]
ext3=[".pdb",".phr",".pin",".pot",".psq",".ptf",".pto"]
###
rule Build_Contamination_Database_by_Blast:
	input:
		"UniVec",
	output:
		expand("UniVec{ext}", ext=ext1),
	threads: workflow.cores
	conda:
		"envs/blast.yaml"		
	shell:
		"makeblastdb -in UniVec -dbtype nucl -out UniVec"
###
rule Scan_Contamination_by_blast:
	input:
		expand("UniVec{ext}", ext=ext2),
		g1="results/Final_assembly/{sample}_Genome.fasta",
		UniVec="UniVec",
	output:
		"results/annotation/{sample}/contamination.bed",
	params:
		"-max_hsps 1 -max_target_seqs 1 -outfmt '6 qseqid qstart qend'",
	threads: workflow.cores
	conda:
		"envs/blast.yaml"
	shell:
		"blastn -db {input.UniVec} {params} -num_threads {threads} -query {input.g1} -out {output}"
###
rule Mask_contaminants_by_bedtools:
	input:
		bed1="results/annotation/{sample}/contamination.bed",
		fa1="results/Final_assembly/{sample}_Genome.fasta",
	output:
		"{sample}_Clean_Genome.fasta",
	threads: workflow.cores
	conda:
		"envs/bedtools.yaml"
	shell:
		"bedtools maskfasta -fi {input.fa1} -bed {input.bed1} -fo {output}"
###
rule QC2_by_GAAS:
	input:
		"{sample}_Clean_Genome.fasta",
	output:
		"results/annotation/{sample}/QC2.stats.txt",
	threads: workflow.cores
	conda:
		"envs/GAAS.yaml"
	shell: 
		"gaas_fasta_statistics.pl -f {input} -o {output}"
###
rule Build_Database_by_RepeatModeler:
	input:
		"{sample}_Clean_Genome.fasta",
	output:
		"{sample}_RepeatDB.nhr",
		"{sample}_RepeatDB.nin",
		"{sample}_RepeatDB.nnd",
		"{sample}_RepeatDB.nni",
		"{sample}_RepeatDB.nog",
		"{sample}_RepeatDB.nsq",
		"{sample}_RepeatDB.translation",
	threads: workflow.cores
	container:
		"docker://dfam/tetools:latest"
	shell:
		"docker run --rm -v $(pwd):/data -w /data/ --mount type=bind,source=/usr/local/bin/trf,target=/opt/trf,ro dfam/tetools:latest BuildDatabase -name {wildcards.sample}_RepeatDB -engine ncbi {input}"
###
rule Detecting_Repeats_by_RepeatModeler:
	input:
		"{sample}_RepeatDB.nhr",
		"{sample}_RepeatDB.nin",
		"{sample}_RepeatDB.nnd",
		"{sample}_RepeatDB.nni",
		"{sample}_RepeatDB.nog",
		"{sample}_RepeatDB.nsq",
		"{sample}_RepeatDB.translation",
	output:
		"{sample}_RepeatDB-families.fa",
		"{sample}_RepeatDB-families.stk",
	threads: workflow.cores
	container:
		"docker://dfam/tetools:latest"
	shell:
		"docker run --rm -v $(pwd):/data -w /data/ --mount type=bind,source=/usr/local/bin/trf,target=/opt/trf,ro dfam/tetools:latest RepeatModeler -pa {threads} -engine ncbi -database {wildcards.sample}_RepeatDB"
###
rule Classify_TE_by_TEclass:
	input:
		"{sample}_RepeatDB-families.fa",
	output:
		"{sample}_TEclass-families/{sample}_RepeatDB-families.fa.lib",
	threads: workflow.cores
	container:
		"docker://hatimalmutairi/teclass-2.1.3b:latest"
	shell:
		"docker run --rm -v $(pwd):/home/data hatimalmutairi/teclass-2.1.3b:latest TEclassTest.pl -r {input}; mv {wildcards.sample}_RepeatDB-families.fa_*/* {wildcards.sample}_TEclass-families/"
###
rule Reformat_TEclass:
	input:
		"{sample}_TEclass-families/{sample}_RepeatDB-families.fa.lib",
	output:
		"{sample}_Reformat_Repeat.fa",
	shell:
		"sed -E 's/^>(\w*.*)#(Unknown)(.*)\|TEclass result: (\w+)\|.*/>\1#\4 \3/;s/unclear/Unknown/g' {input} > {output}"
###
rule Mask_repeats_by_RepeatMasker:
	input:
		a1="{sample}_Reformat_Repeat.fa",
		a2="{sample}_Clean_Genome.fasta",
	output:
		"{sample}_Clean_Genome.fasta.out",
		"{sample}_Clean_Genome.fasta.cat.gz",
		"{sample}_Clean_Genome.fasta.tbl",
		"{sample}_Clean_Genome.fasta.align",
		"{sample}_Clean_Genome.fasta.masked",
		"{sample}_Clean_Genome.fasta.out.gff",
	threads: workflow.cores
	container:
		"docker://dfam/tetools:latest"
	shell:
		"docker run --rm -v $(pwd):/data -w /data/ --mount type=bind,source=/usr/local/bin/trf,target=/opt/trf,ro dfam/tetools:latest RepeatMasker -a -gff -x -pa {threads} -e ncbi -lib {input.a1} {input.a2}"
###
rule Mask_repeats2_by_RepeatMasker:
	input:
		"{sample}_Clean_Genome.fasta.cat.gz",
	output:
		"{sample}_Clean_Genome.fasta.cat",
	shell:
		"gunzip -k {input} > {output}"
###
rule Mask_repeats3_by_RepeatMasker:
	input:
		"{sample}_Clean_Genome.fasta.out",
	output:
		"{sample}_Clean_Genome_Repeats.gff3",
	threads: workflow.cores
	container:
		"docker://dfam/tetools:latest"
	shell:
		"docker run --rm -v $(pwd):/data -w /data/ --mount type=bind,source=/usr/local/bin/trf,target=/opt/trf,ro dfam/tetools:latest rmOutToGFF3.pl {input} > {output}"
###
rule Repeats_reformat1:
	input:
		"{sample}_Clean_Genome_Repeats.gff3",
	output:
		"{sample}_Clean_Genome_Repeats2.gff3",
	threads: workflow.cores
	shell:
		"grep -v -e 'Satellite' -e ')n' -e '-rich' {input} > {output}"
###
rule Repeats_reformat2:
	input:
		"{sample}_Clean_Genome_Repeats2.gff3",
	output:
		"{sample}_Clean_Genome_Repeats3.gff3",
	threads: workflow.cores
	shell:
		"python scripts/reformatGff.py {input} {output}"
###
evidance= ["LarabicaLEM1108_AnnotatedProteins.fasta","LpanamensisMHOMCOL81L13_AnnotatedProteins.fasta","LbraziliensisMHOMBR75M2904_AnnotatedProteins.fasta","LtarentolaeParrotTarII_AnnotatedProteins.fasta","LturanicaLEM423_AnnotatedProteins.fasta","LdonovaniLV9_AnnotatedProteins.fasta","LmajorSD75.1_AnnotatedProteins.fasta","LgerbilliLEM452_AnnotatedProteins.fasta","LbraziliensisMHOMBR75M2904_2019_AnnotatedProteins.fasta","LtropicaL590_AnnotatedProteins.fasta","LdonovaniBPK282A1_AnnotatedProteins.fasta","LpanamensisMHOMPA94PSC1_AnnotatedProteins.fasta","LenriettiiLEM3045_AnnotatedProteins.fasta","LmajorLV39c5_AnnotatedProteins.fasta","LdonovaniCL-SL_AnnotatedProteins.fasta","LamazonensisMHOMBR71973M2269_AnnotatedProteins.fasta","LmexicanaMHOMGT2001U1103_AnnotatedProteins.fasta","LaethiopicaL147_AnnotatedProteins.fasta","LinfantumJPCM5_AnnotatedProteins.fasta","LmajorFriedlin_AnnotatedProteins.fasta","LenriettiiLEM3045.gff","LdonovaniLV9_AnnotatedTranscripts.fasta","LdonovaniCL-SL_AnnotatedTranscripts.fasta","LinfantumJPCM5_AnnotatedTranscripts.fasta","LarabicaLEM1108_AnnotatedTranscripts.fasta","LmajorLV39c5_AnnotatedTranscripts.fasta","LbraziliensisMHOMBR75M2904_AnnotatedTranscripts.fasta","LdonovaniBPK282A1_AnnotatedTranscripts.fasta","LpanamensisMHOMCOL81L13_AnnotatedTranscripts.fasta","LgerbilliLEM452_AnnotatedTranscripts.fasta","LmajorFriedlin_AnnotatedTranscripts.fasta","LpanamensisMHOMPA94PSC1_AnnotatedTranscripts.fasta","LmajorSD75.1_AnnotatedTranscripts.fasta","LturanicaLEM423_AnnotatedTranscripts.fasta","LenriettiiLEM3045_AnnotatedTranscripts.fasta","LaethiopicaL147_AnnotatedTranscripts.fasta","LmexicanaMHOMGT2001U1103_AnnotatedTranscripts.fasta","LtarentolaeParrotTarII_AnnotatedTranscripts.fasta","LbraziliensisMHOMBR75M2904_2019_AnnotatedTranscripts.fasta","LamazonensisMHOMBR71973M2269_AnnotatedTranscripts.fasta","LtropicaL590_AnnotatedTranscripts.fasta","LbraziliensisMHOMBR75M2903_AnnotatedTranscripts.fasta","Lbraziliensis_ESTs.fasta","Linfantum_ESTs.fasta","Lmexicana_ESTs.fasta","Lmajor_ESTs.fasta","Lamazonensis_ESTs.fasta"]

rule Round1_Annotation_by_MAKER:
	input:
		expand("data/Evidance/TriTrypDB-47_{evidance}",evidance=evidance),
		"{sample}_Clean_Genome_Repeats3.gff3",
		g1="{sample}_Clean_Genome.fasta",
		c1="config/maker_opts.{sample}.R1.ctl",
		c2="config/maker_bopts.docker.ctl",
		c3="config/maker_exe.docker.ctl",
	output:
		#directory("{sample}_Round1.maker.output"),
		"{sample}_Round1.maker.output/{sample}_Round1_master_datastore_index.log",
	threads: workflow.cores
	container:
		"docker://hatimalmutairi/lmgaap-maker:latest"
	shell:
		"docker run --rm -v $(pwd):/home/data -w /home/data hatimalmutairi/lmgaap-maker mpirun --allow-run-as-root -n {threads} maker -fix_nucleotides -genome {input.g1} -base {wildcards.sample}_Round1 {input.c1} {input.c2} {input.c3}"
### 
checkpoint Post_Round1_Annotation_by_MAKER:
	input:
		"{sample}_Round1.maker.output/{sample}_Round1_master_datastore_index.log",
	output:
		out1="{sample}_Round1.maker.output/Round1_all.gff",
		out2="{sample}_Round1.maker.output/Round1_Noseq.gff",
	threads: workflow.cores
	container:
		"docker://hatimalmutairi/lmgaap-maker:latest"
	shell:
		"""
		docker run --rm -v $(pwd):/home/data -w /home/data hatimalmutairi/lmgaap-maker:latest gff3_merge -s -d {input} > {output.out1}
		docker run --rm -v $(pwd):/home/data -w /home/data hatimalmutairi/lmgaap-maker:latest gff3_merge -n -s -d {input} > {output.out2}
		docker run --rm -v $(pwd):/home/data -w /home/data hatimalmutairi/lmgaap-maker:latest fasta_merge -d {input}
		"""
###
rule QC_Round1_by_Genometools:
	input:
		"{sample}_Round1.maker.output/Round1_all.gff",
	output:
		"{sample}_Round1_QC_report_1.txt",
	container:
		"docker://quay.io/biocontainers/genometools:1.2.1--py27_0"
	shell:
		"docker run --rm -v $(pwd):/home/data -w /home/data quay.io/biocontainers/genometools:1.2.1--py27_0 gt gff3 -sort -tidy {input} | gt stat > {output}"
###
rule QC_Round1_by_GAAS:
	input:
		"{sample}_Round1.maker.output/Round1_all.gff",
	output:
		"{sample}_Round1_QC_report_2.txt",
	threads: workflow.cores
	conda:
		"envs/GAAS.yaml"
	shell: 
		"agat_sp_statistics.pl -d --gff {input} -o {output}"
###
checkpoint Round2_Annotation_by_MAKER:
	input:
		"{sample}_Round1.maker.output/{sample}_Round1_master_datastore_index.log",
		"{sample}_Round1.maker.output/Round1_all.gff",
		"{sample}_Clean_Genome_Repeats3.gff3",
		g1="{sample}_Clean_Genome.fasta",
		c1="config/maker_opts.{sample}.R2.ctl",
		c2="config/maker_bopts.docker.ctl",
		c3="config/maker_exe.docker.ctl",
	output:
		"{sample}_Round2.maker.output/{sample}_Round2_master_datastore_index.log",
		directory("{sample}_Round2.maker.output"),
	threads: workflow.cores
	container:
		"docker://hatimalmutairi/lmgaap-maker:latest"
	shell:
		"docker run --rm -v $(pwd):/home/data -w /home/data hatimalmutairi/lmgaap-maker mpirun --allow-run-as-root -n {threads} maker -genome {input.g1} -base {wildcards.sample}_Round2 {input.c1} {input.c2} {input.c3}"
###
rule Post_Round2_Annotation_by_MAKER:
	input:		
		"{sample}_Round2.maker.output/Round2_master_datastore_index.log",
	output:
		out1="{sample}_Round2.maker.output/Round2_all.gff",
		out2="{sample}_Round2.maker.output/Round2_Noseq.gff",
	threads: workflow.cores
	container:
		"docker://hatimalmutairi/lmgaap-maker:latest"
	shell:
		"""
		docker run --rm -v $(pwd):/home/data -w /home/data hatimalmutairi/lmgaap-maker:latest gff3_merge -s -d {input} > {output.out1}
		docker run --rm -v $(pwd):/home/data -w /home/data hatimalmutairi/lmgaap-maker:latest gff3_merge -n -s -d {input} > {output.out2}
		docker run --rm -v $(pwd):/home/data -w /home/data hatimalmutairi/lmgaap-maker:latest fasta_merge -d {input}
		"""
###
rule QC_Round2_by_Genometools:
	input:
		"{sample}_Round2.maker.output/Round2_all.gff",
	output:
		"results/annotation/{sample}/Round2_QC_report_1.txt",
	threads: workflow.cores
	container:
		"docker://quay.io/biocontainers/genometools:1.2.1--py27_0"
	shell:
		"docker run --rm -v $(pwd):/home/data -w /home/data quay.io/biocontainers/genometools:1.2.1--py27_0 gt gff3 -sort -tidy {input} | gt stat > {output}"
###
rule QC_Round2_by_GAAS:
	input:
		"{sample}_Round2.maker.output/Round2_all.gff",
	output:
		"results/annotation/{sample}/Round2_QC_report_2.txt",
	threads: workflow.cores
	conda:
		"envs/GAAS.yaml"
	shell: 
		"agat_sp_statistics.pl -d --gff {input} -o {output}"
###
checkpoint Proccess_Annotation_Rounds:
	input:
		"{sample}_Round2.maker.output",
	output:
		"maker_output_processed_{sample}_Round2/maker_annotation.gff",
		"maker_output_processed_{sample}_Round2/maker_annotation.proteins.fasta",
		"maker_output_processed_{sample}_Round2/maker_annotation.transcripts.fasta",
		d1=directory("maker_output_processed_{sample}_Round2")
	threads: workflow.cores
	conda:
		"envs/GAAS.yaml"
	shell: 
		"gaas_maker_merge_outputs_from_datastore.pl -i {input} -o {output.d1}"
###
checkpoint create_Final_Annotation:
	input:
		a1="maker_output_processed_{sample}_Round2/maker_annotation.gff",
		a2="maker_output_processed_{sample}_Round2/maker_annotation.proteins.fasta",
		a3="maker_output_processed_{sample}_Round2/maker_annotation.transcripts.fasta",
		a4="{sample}_Clean_Genome.fasta",
	output:
		a1="results/annotation/{sample}/Final_Annotation/annotation.gff",
		a2="results/annotation/{sample}/Final_Annotation/annotation.proteins.fasta",
		a3="results/annotation/{sample}/Final_Annotation/annotation.transcripts.fasta",
		a4="results/annotation/{sample}/Final_Annotation/Genome.fasta",
	threads: workflow.cores
	shell:
		"""
		cp {input.a1} {output.a1}
		cp {input.a2} {output.a2}
		cp {input.a3} {output.a3}
		cp {input.a4} {output.a4}
		"""
###
rule Download_Uniprot_Database:
	output:
		"Uniprot",
	shell:
		"wget -O - ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz | gunzip -c > {output}"
###
rule Build_Uniprot_Database_by_Blast:
	input:
		"Uniprot",
	output:
		expand("Uniprot{ext}", ext=ext3),
	threads: workflow.cores
	conda:
		"envs/blast.yaml"		
	shell:
		"makeblastdb -in {input} -dbtype prot -out Uniprot"
###
rule Search_for_Annotation_by_blast:
	input:
		u0="results/annotation/{sample}/Final_Annotation/annotation.proteins.fasta",
		u1=expand("Uniprot{ext}", ext=ext3),
		u2="Uniprot"
	output:
		"results/annotation/{sample}/Final_Annotation/annotation.proteins.blastp",
	threads: workflow.cores
	conda:
		"envs/blast.yaml"
	shell:
		"blastp -num_threads {threads} -query {input.u0} -db {input.u2} -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out {output}"
###
rule Search_for_Annotation_by_interproscan:
	input:
		"results/annotation/{sample}/Final_Annotation/annotation.proteins.fasta",
	output:
		"results/annotation/{sample}/Final_Annotation/annotation.proteins.iprscan",
	threads: workflow.cores
	container:
		"docker://blaxterlab/interproscan:latest"
	shell:
		"docker run --rm -v $(pwd):/in -w /in blaxterlab/interproscan:latest interproscan.sh -appl pfam -dp -f TSV -goterms -iprlookup -pa -t p -i {input} -o {output}"
###
rule Assign_Annotations_Map:
	input:
		"results/annotation/{sample}/Final_Annotation/annotation.gff",
	output:
		"results/annotation/{sample}/Final_Annotation/annotation.map",
	threads: workflow.cores
	container:
		"docker://hatimalmutairi/lmgaap-maker:latest"
	shell:
		"docker run --rm -v $(pwd):/home/data -w /home/data  hatimalmutairi/lmgaap-maker:latest maker_map_ids --prefix {wildcards.sample}_ --justify 5 {input} > {output}"
###
rule Assign_new_ID_to_Annotations:
	input:
		map1="results/annotation/{sample}/Final_Annotation/annotation.map",
		a1="results/annotation/{sample}/Final_Annotation/annotation.gff",
		a2="results/annotation/{sample}/Final_Annotation/annotation.proteins.fasta",
		a3="results/annotation/{sample}/Final_Annotation/annotation.transcripts.fasta",
		a4="results/annotation/{sample}/Final_Annotation/annotation.proteins.blastp",
		a5="results/annotation/{sample}/Final_Annotation/annotation.proteins.iprscan",
	threads: workflow.cores
	container:
		"docker://hatimalmutairi/lmgaap-maker:latest"
	shell:
		"""
		docker run --rm -v $(pwd):/home/data -w /home/data  hatimalmutairi/lmgaap-maker:latest map_gff_ids {input.map1} {input.a1}
		docker run --rm -v $(pwd):/home/data -w /home/data  hatimalmutairi/lmgaap-maker:latest map_fasta_ids {input.map1} {input.a2}
		docker run --rm -v $(pwd):/home/data -w /home/data  hatimalmutairi/lmgaap-maker:latest map_fasta_ids {input.map1} {input.a3}
		docker run --rm -v $(pwd):/home/data -w /home/data  hatimalmutairi/lmgaap-maker:latest map_data_ids {input.map1} {input.a4}
		docker run --rm -v $(pwd):/home/data -w /home/data  hatimalmutairi/lmgaap-maker:latest map_data_ids {input.map1} {input.a5}
		"""
###
rule Assign_Functional_Annotations:
	input:
		gff="results/annotation/{sample}/Final_Annotation/annotation.gff",
		proteins="results/annotation/{sample}/Final_Annotation/annotation.proteins.fasta",
		transcripts="results/annotation/{sample}/Final_Annotation/annotation.transcripts.fasta",
		blastp="results/annotation/{sample}/Final_Annotation/annotation.proteins.blastp",
		Uniprot="Uniprot"
	output:
		gff="results/annotation/{sample}/Final_Annotation/putative_function.gff",
		proteins="results/annotation/{sample}/Final_Annotation/putative_function.proteins.fasta",
		transcripts="results/annotation/{sample}/Final_Annotation/putative_function.transcripts.fasta",
	threads: workflow.cores
	container:
		"docker://hatimalmutairi/lmgaap-maker:latest"
	shell:
		"""
		docker run --rm -v $(pwd):/home/data -w /home/data  hatimalmutairi/lmgaap-maker:latest maker_functional_fasta {input.Uniprot} {input.blastp} {input.gff} > {output.gff}
		docker run --rm -v $(pwd):/home/data -w /home/data  hatimalmutairi/lmgaap-maker:latest maker_functional_fasta {input.Uniprot} {input.blastp} {input.proteins} > {output.proteins}
		docker run --rm -v $(pwd):/home/data -w /home/data  hatimalmutairi/lmgaap-maker:latest maker_functional_fasta {input.Uniprot} {input.blastp} {input.transcripts} > {output.transcripts}
		"""
###
rule Assign_Final_GFF:
	input:
		gff="results/annotation/{sample}/Final_Annotation/putative_function.gff",
		iprscan="results/annotation/{sample}/Final_Annotation/annotation.proteins.iprscan",
	output:
		"results/annotation/{sample}/Final_Annotation/Final_Annotation.gff",
	threads: workflow.cores
	container:
		"docker://hatimalmutairi/lmgaap-maker:latest"
	shell:
		"docker run --rm -v $(pwd):/home/data -w /home/data hatimalmutairi/lmgaap-maker:latest ipr_update_gff {input.gff} {input.iprscan} > {output}"
###
rule Keep_Longest_Isoform_by_AGAT:
	input:
		"results/annotation/{sample}/Final_Annotation/Final_Annotation.gff",
	output:
		"results/annotation/{sample}/Final_Annotation/Final_Annotation_Longest_Isoform.gff",
	threads: workflow.cores
	conda:
		"envs/AGAT.yaml"
	shell: 
		"agat_sp_keep_longest_isoform.pl -gff {input} -o {output}"
###
