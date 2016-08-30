## Example of code used to go from fastq.gz files to read count matrix.

genome_path=/path/to/MacaM7_genome.fasta
gtf_path=/path/to/MacaM_Rhesus_Genome_Annotation_v7.6.8.gtf
adapter="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
STAR_1Pass=/path/to/newdir_STAR1
STAR_2Pass=/path/to/newdir_STAR2
STAR_exec=/path/to/STAR_executable
fastq_path=/path/to/fastqs

#Example of file name
sample="stem_of_filename.fastq.gz"
j=${sample/\.fastq\.gz/}

#1) Read trimming of adaptors
#Trim Galore version : 0.2.7 needs to be installed and in your default PATH
#cutadapt version : 1.2.1 needs to be installed and in your default PATH

mkdir $j.trimQ20 ; trim_galore -q 20 -o $j.trimQ20 --phred33 -a $adapter $fastq_path/$j.fastq.gz  ## note that this will place the trimmed fastq file into a new directory in the current directory 


#2) STAR 1st pass mapping :
#a) generate the index 
$STAR_exec --runMode genomeGenerate --genomeDir $STAR_1Pass --genomeFastaFiles $genome_path --runThreadN 12 --sjdbOverhang 99 --sjdbGTFfile $gtf_path


#b) mapping
$STAR_exec --genomeDir $STAR_1Pass --readFilesIn $j.trimQ20/${j}_trimmed.fq.gz --runThreadN 12 --readFilesCommand zcat --outSAMstrandField intronMotif --outFileNamePrefix $j.MacaM.trimQ20.2passALL.MacaM_v7.6.8.10MM.  --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterMismatchNmax 10 --outSAMattributes NH nM MD --outSAMattrRGline  ID:$j PU:Illumina PL:Illumina LB:Library SM:$j CN:CHUSJ


#3) SJDB file generation
#Use all the sjdb files generated from the previous step (i.e., for all samples)
#Concat, sort them and make a unique list (remove duplicated lines).

#4) STAR 2nd pass mapping
#a) index generation (NOTE: This uses a lot of memory, so launch it from a node with a large amount of RAM; e.g., 256GB)
$STAR_exec --runMode genomeGenerate --genomeDir $STAR_2PASS --genomeFastaFiles $genome_path --runThreadN 12 --sjdbGTFfile $gtf_path --sjdbFileChrStartEnd SGE_allCombined.2phases.SJdb.strand.out.tab --sjdbOverhang 99 --limitGenomeGenerateRAM 80000000000


#b) mapping

$STAR_exec --genomeDir $STAR_2PASS --readFilesIn $f --runThreadN 12 --readFilesCommand zcat --outSAMstrandField intronMotif --outFileNamePrefix $j.MacaM.trimQ20.2passALL.MacaM_v7.6.8.10MM.  --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterMismatchNmax 10 --outSAMattributes NH nM MD --outSAMattrRGline  ID:$j PU:Illumina PL:Illumina LB:Library SM:$j CN:CHUSJ --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx


#5) HTseq: get per gene read counts for each sample

#HTseq version : 0.6.1p1

samtools sort -@ 3 -n $j.MacaM.trimQ20.2passALL.MacaM_v7.6.8.10MM.Aligned.sortedByCoord.out.bam $LSCRATCH/$j.MacaM.trimQ20.2passALL.MacaM_v7.6.8.10MM.Aligned.sortedByCoord.out.sort ; python -m HTSeq.scripts.count -s no -f bam -m intersection-nonempty $LSCRATCH/$j.MacaM.trimQ20.2passALL.MacaM_v7.6.8.10MM.Aligned.sortedByCoord.out.sort.bam $gtf_path > Count.genes.HTseq.$j.MacaM.trimQ20.2passALL.MacaM_v7.6.8.10MM.Aligned.sortedByCoord.out.txt

