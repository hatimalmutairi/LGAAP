#-----Genome (these are always required)
genome=
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=1 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=1 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=1 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=1 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=1 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=1 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=1 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est=data/Evidance/TriTrypDB-47_Lmexicana_ESTs.fasta:Lmexicana_ESTs,data/Evidance/TriTrypDB-47_Linfantum_ESTs.fasta:Linfantum_ESTs,data/Evidance/TriTrypDB-47_Lmajor_ESTs.fasta:Lmajor_ESTs,data/Evidance/TriTrypDB-47_Lbraziliensis_ESTs.fasta:Lbraziliensis_ESTs,data/Evidance/TriTrypDB-47_Lamazonensis_ESTs.fasta:Lamazonensis_ESTs
altest=data/Evidance/TriTrypDB-47_LbraziliensisMHOMBR75M2904_2019_AnnotatedTranscripts.fasta:LbraziliensisMHOMBR75M2904_2019,data/Evidance/TriTrypDB-47_LdonovaniLV9_AnnotatedTranscripts.fasta:LdonovaniLV9,data/Evidance/TriTrypDB-47_LdonovaniCL-SL_AnnotatedTranscripts.fasta:LdonovaniCL-SL,data/Evidance/TriTrypDB-47_LpanamensisMHOMPA94PSC1_AnnotatedTranscripts.fasta:LpanamensisMHOMPA94PSC1,data/Evidance/TriTrypDB-47_LturanicaLEM423_AnnotatedTranscripts.fasta:LturanicaLEM423,data/Evidance/TriTrypDB-47_LtropicaL590_AnnotatedTranscripts.fasta:LtropicaL590,data/Evidance/TriTrypDB-47_LtarentolaeParrotTarII_AnnotatedTranscripts.fasta:LtarentolaeParrotTarII,data/Evidance/TriTrypDB-47_LpanamensisMHOMCOL81L13_AnnotatedTranscripts.fasta:LpanamensisMHOMCOL81L13,data/Evidance/TriTrypDB-47_LmexicanaMHOMGT2001U1103_AnnotatedTranscripts.fasta:LmexicanaMHOMGT2001U1103,data/Evidance/TriTrypDB-47_LmajorSD75.1_AnnotatedTranscripts.fasta:LmajorSD75.1,data/Evidance/TriTrypDB-47_LmajorLV39c5_AnnotatedTranscripts.fasta:LmajorLV39c5,data/Evidance/TriTrypDB-47_LmajorFriedlin_AnnotatedTranscripts.fasta:LmajorFriedlin,data/Evidance/TriTrypDB-47_LinfantumJPCM5_AnnotatedTranscripts.fasta:LinfantumJPCM5,data/Evidance/TriTrypDB-47_LgerbilliLEM452_AnnotatedTranscripts.fasta:LgerbilliLEM452,data/Evidance/TriTrypDB-47_LenriettiiLEM3045_AnnotatedTranscripts.fasta:LenriettiiLEM3045,data/Evidance/TriTrypDB-47_LdonovaniBPK282A1_AnnotatedTranscripts.fasta:LdonovaniBPK282A1,data/Evidance/TriTrypDB-47_LbraziliensisMHOMBR75M2904_AnnotatedTranscripts.fasta:LbraziliensisMHOMBR75M2904,data/Evidance/TriTrypDB-47_LbraziliensisMHOMBR75M2903_AnnotatedTranscripts.fasta:LbraziliensisMHOMBR75M2903,data/Evidance/TriTrypDB-47_LarabicaLEM1108_AnnotatedTranscripts.fasta:LarabicaLEM1108,data/Evidance/TriTrypDB-47_LamazonensisMHOMBR71973M2269_AnnotatedTranscripts.fasta:LamazonensisMHOMBR71973M2269,data/Evidance/TriTrypDB-47_LaethiopicaL147_AnnotatedTranscripts.fasta:LaethiopicaL147
est_gff=data/Evidance/TriTrypDB-47_LenriettiiLEM3045.gff #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff=data/Evidance/TriTrypDB-47_LbraziliensisMHOMBR75M2904_2019.gff:LbraziliensisMHOMBR75M2904_2019,data/Evidance/TriTrypDB-47_LdonovaniLV9.gff:LdonovaniLV9,data/Evidance/TriTrypDB-47_LdonovaniCL-SL.gff:LdonovaniCL-SL,data/Evidance/TriTrypDB-47_LpanamensisMHOMPA94PSC1.gff:LpanamensisMHOMPA94PSC1,data/Evidance/TriTrypDB-47_LturanicaLEM423.gff:LturanicaLEM423,data/Evidance/TriTrypDB-47_LtropicaL590.gff:LtropicaL590,data/Evidance/TriTrypDB-47_LtarentolaeParrotTarII.gff:LtarentolaeParrotTarII,data/Evidance/TriTrypDB-47_LpanamensisMHOMCOL81L13.gff:LpanamensisMHOMCOL81L13,data/Evidance/TriTrypDB-47_LmexicanaMHOMGT2001U1103.gff:LmexicanaMHOMGT2001U1103,data/Evidance/TriTrypDB-47_LmajorSD75.1.gff:LmajorSD75.1,data/Evidance/TriTrypDB-47_LmajorLV39c5.gff:LmajorLV39c5,data/Evidance/TriTrypDB-47_LmajorFriedlin.gff:LmajorFriedlin,data/Evidance/TriTrypDB-47_LgerbilliLEM452.gff:LgerbilliLEM452,data/Evidance/TriTrypDB-47_LinfantumJPCM5.gff:LinfantumJPCM5,data/Evidance/TriTrypDB-47_LenriettiiLEM3045.gff:LenriettiiLEM3045,data/Evidance/TriTrypDB-47_LdonovaniBPK282A1.gff:LdonovaniBPK282A1,data/Evidance/TriTrypDB-47_LbraziliensisMHOMBR75M2904.gff:LbraziliensisMHOMBR75M2904,data/Evidance/TriTrypDB-47_LbraziliensisMHOMBR75M2903.gff:LbraziliensisMHOMBR75M2903,data/Evidance/TriTrypDB-47_LarabicaLEM1108.gff:LarabicaLEM1108,data/Evidance/TriTrypDB-47_LamazonensisMHOMBR71973M2269.gff:LamazonensisMHOMBR71973M2269,data/Evidance/TriTrypDB-47_LaethiopicaL147.gff:LaethiopicaL147

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=data/Evidance/TriTrypDB-47_LbraziliensisMHOMBR75M2904_2019_AnnotatedProteins.fasta:LbraziliensisMHOMBR75M2904_2019,data/Evidance/TriTrypDB-47_LdonovaniLV9_AnnotatedProteins.fasta:LdonovaniLV9,data/Evidance/TriTrypDB-47_LdonovaniCL-SL_AnnotatedProteins.fasta:LdonovaniCL-SL,data/Evidance/TriTrypDB-47_LpanamensisMHOMPA94PSC1_AnnotatedProteins.fasta:LpanamensisMHOMPA94PSC1,data/Evidance/TriTrypDB-47_LturanicaLEM423_AnnotatedProteins.fasta:LturanicaLEM423,data/Evidance/TriTrypDB-47_LtropicaL590_AnnotatedProteins.fasta:LtropicaL590,data/Evidance/TriTrypDB-47_LtarentolaeParrotTarII_AnnotatedProteins.fasta:LtarentolaeParrotTarII,data/Evidance/TriTrypDB-47_LpanamensisMHOMCOL81L13_AnnotatedProteins.fasta:LpanamensisMHOMCOL81L13,data/Evidance/TriTrypDB-47_LmexicanaMHOMGT2001U1103_AnnotatedProteins.fasta:LmexicanaMHOMGT2001U1103,data/Evidance/TriTrypDB-47_LmajorSD75.1_AnnotatedProteins.fasta:LmajorSD75.1,data/Evidance/TriTrypDB-47_LmajorLV39c5_AnnotatedProteins.fasta:LmajorLV39c5,data/Evidance/TriTrypDB-47_LmajorFriedlin_AnnotatedProteins.fasta:LmajorFriedlin,data/Evidance/TriTrypDB-47_LinfantumJPCM5_AnnotatedProteins.fasta:LinfantumJPCM5,data/Evidance/TriTrypDB-47_LgerbilliLEM452_AnnotatedProteins.fasta:LgerbilliLEM452,data/Evidance/TriTrypDB-47_LdonovaniBPK282A1_AnnotatedProteins.fasta:LdonovaniBPK282A1,data/Evidance/TriTrypDB-47_LenriettiiLEM3045_AnnotatedProteins.fasta:LenriettiiLEM3045,data/Evidance/TriTrypDB-47_LbraziliensisMHOMBR75M2904_AnnotatedProteins.fasta:LbraziliensisMHOMBR75M2904,data/Evidance/TriTrypDB-47_LbraziliensisMHOMBR75M2903_AnnotatedProteins.fasta:LbraziliensisMHOMBR75M2903,data/Evidance/TriTrypDB-47_LarabicaLEM1108_AnnotatedProteins.fasta:LarabicaLEM1108,data/Evidance/TriTrypDB-47_LamazonensisMHOMBR71973M2269_AnnotatedProteins.fasta:LamazonensisMHOMBR71973M2269,data/Evidance/TriTrypDB-47_LaethiopicaL147_AnnotatedProteins.fasta:LaethiopicaL147
protein_gff=  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org=all #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff=JIQ42_Clean_Genome_Repeats3.gff3 #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm= #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species= #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=1 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
