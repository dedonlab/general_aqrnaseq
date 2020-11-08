# general_aqrnaseq
a general AQRNAseq pipeline that can process both eukaryote and prokaryote RNAseq samples on unix clusters using slurm workload manager

# Section 1: sequence level analyses
Step1: count the number of sequencing reads in original fastq files

1.1 run the sequence counting submission shell script submit_seqcount.sh

./submit_seqcount.sh path_to_the_directory_containing_original_fastq_files/*fastq

example: ./submit_seqcount.sh /home/fastqbin/*fastq

1.2 concatenate the counting results

head *count >raw_seq_count

1.3 format the counting results to tablulate summary

./format.pl raw_seq_count raw_seqcount_formatted

Step2: Pear assembly and linker stripping

2.1 submit pear assembly jobs

./submit_pear.sh path_to_the_directory_containing_original_sequencing_fastq_files/*1_sequence.fastq

example: ./submit_pear.sh /home/fastqbin/*1_sequence.fastq

2.2 extract assembled read counts from pear statistics output files and format them to tabular summary

./submit_pearCount.sh

cat *assembledByPear >assembledByPearCounts

Step3: trim 2 random nucleotides added at 3 prime of the inserts

3.1 submit trimming scripts

./submit_trim3p2N.sh

3.2 count the number of sequencing reads in 3p2N stripped fq files

./submit_seqcount_trim3p2N.sh

Step4: select sqeunces with the correct length. For miRNA, inserts with length 20-30bp are selected allowing nuclutide additions at both ends

4.1 run  slength selection scripts

./submit_len20_30.sh

4.2 Count the number of seqeunces with length 20bp-30bp and format the counts to tabular summary

./submit_seqcount20_30bp.sh

head *count20_30bp >seqcount20_30bp

./format20_30bp.pl seqcount20_30bp seqcount20_30bp_formatted

Step5: convert fastq files to fasta files

./submit_fq2fa.sh

Step6: count the abundance of each unique read

6.1 count unique reads in each sample

./submit_fastxcollapser.sh

6.2 formatting the counts

./submit_format.sh

Step7: merge counts to tabular summary each sequence per row and each sample per column

sbatch merge_count.sh

# Section2: analyses of specific types of RNA

Step8: Prep a fasta files containing unique seqeunces in all samples for blast

cut -f1 merged_count >seq

./fa.pl seq seq.fa

Step9: Blast

sbatch blast.sh blastDabase

example: sbatch blast.sh /home/blastdb/mature_HomoSapiens_miRNA.fa

Step10: select blast hits 17bp or longer

./blastlen17.pl seq.blast seq17bp.blast

Step11: count the number of non redundant sequences hitting RNA species of interest

cut -f1 seq17b0.blast|sort -u|wc -l

Step12: Count the number of RNA genes got a hit in blast

cut -f2 seq17b0.blast|sort -u|wc -l

Step13: conclude counting results

XXXX total sequences. Among them, XXXX unique sequences hit XXX RNA species of interest

Step14: subset results to sequences belong to RNA species of interest only

./faToTable.pl seq.fa seq.tab

./overlaphash.pl seq17bp.blast seq.tab seq17bpWithSeq.blast

cat header seq17bpWithSeq.blast >seq17bpWithSeq.blasth

./overlaphash2.pl seq17bpWithSeq.blasth merged_count merged_miRNA_17bpPerfect_count

# Section3: Genome wide analyses

Step15: mapping unique reads to reference genome

sbatch bwa.sh seq.fa seq path_to_bwa_reference_file

example: sbatch bwa.sh seq.fa seq /home/Genomes/bwa_indexes/hg38simple.fa

Step16: Samtools reverted reverse strand sequences. So incorprated original reads

./seqAndHeader.pl seq seqTable

./origSeq.pl seq.sam seqTable origSeq.sam

./formatSamFile.pl origSeq.sam samMapped

Step17:Annotate sequences

17.1 annotate sequences within genes.  Grch38Gene.bed should be substitued to the gene boundary bed file of the correpsponding species

sbatch bedtools.sh

17.2 annotate sequences in gene desert regions

./geneDesert.pl samMapped samGeneAnnot geneDesert

17.3 annotate not mapped reads

./notMappedGenes.pl origSeq.sam samMapped notMapped

17.4 put all annotation together

cat samGeneAnnot geneDesert notMapped >allSeq

Step18 calculate row sum of the raw counts

sbatch rowSum.sh

Step19 merge sequence mapping information with counts in each sample

cat headerAllSeq allSeq >allSeqH

./overlapSeqInfoWithCounts.pl allSeqH mergedRawCount.txt mergedRawCountWithAnnotationGenomewise.txt

Step20 sort mergedRawCountWithAnnotationGenomewise.txt by row sum

sbatch sortMergedCount.sh

Step21 Count the number of non redundant sequences in genome wide analyses including reads in gene regions, gene desert regions, and notMapped reads

wc -l samGeneAnnot geneDesert notMapped

You should see something like the following:

1044543 samGeneAnnot

915263 geneDesert

7014408 notMapped
