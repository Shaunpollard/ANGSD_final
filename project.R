

(angsd) bash-4.4$ nano accession_list.txt



(angsd) bash-4.4$ nano sra_download.sh
accession_num=$(cat accession_list.txt)


# make new directory to hold fastq files
#mkdir fastq
#cd fastq

# Read accession numbers from the file and iterate over them
while read -r accession_num; do
echo "$accession_num"  # Print each accession number before downloading
fasterq-dump "$accession_num" -F fastq -O fastq --skip-technical --split-3 --threads 3
done < accession_list.txt

echo "All downloads complete :)"




(angsd) bash-4.4$ bash sra_download.sh 
SRR1804235
spots read      : 39,324,028
reads read      : 78,648,056
reads written   : 78,648,056
SRR1804236
spots read      : 50,710,879
reads read      : 101,421,758
reads written   : 101,421,758
SRR1804237
spots read      : 72,207,041
reads read      : 144,414,082
reads written   : 144,414,082
SRR1804238
spots read      : 43,927,084
reads read      : 87,854,168
reads written   : 87,854,168
SRR1804239
spots read      : 33,887,041
reads read      : 67,774,082
reads written   : 67,774,082
SRR1804240
spots read      : 55,361,657
reads read      : 110,723,314
reads written   : 110,723,314
SRR1804241
spots read      : 28,727,301
reads read      : 57,454,602
reads written   : 57,454,602
SRR1804242
spots read      : 42,629,245
reads read      : 85,258,490
reads written   : 85,258,490

Cannot process ''.
All downloads complete :)


(angsd) bash-4.4$ for file in *_*.fastq; do fastqc "$file"; done 


(base) shaunp@MAC306097 Classwork % scp -i ~/.ssh/Private_Key.txt shp4022@cayuga-login1.cac.cornell.edu:/athena/angsd/scratch/shp4022/project/fastq/SRR1804238_2_trimmed_fastqc.html  /Users/shaunp/Desktop/Weill_Cornell_Graduate/Grad_School/SPRING_2024/Analysis_Next-Gen_Sequencing_Data/Project/Fastqc/SRR1804238_2_trimmed_fastqc.html
SRR1804238_2_trimmed_fastqc.html                                                                    100%  616KB   2.0MB/s   00:00    



(multiqc) [shp4022@cayuga-login1 fastq]$ multiqc *_fastqc.zip 


(multiqc) [shp4022@cayuga-login1 fastq]$ mamba activate trim-galore
for file in *.fastq; do trim_galore --fastqc --stringency 3 $file; done

(angsd) [shp4022@cayuga-login1 genome]$ nano STAR_genomeGenerate.sh 
(angsd) [shp4022@cayuga-login1 genome]$ chmod u+x STAR_genomeGenerate.sh 

STAR --runMode genomeGenerate \
--runThreadN 5 \
--genomeDir galGal6_STAR_Index \
`#--genomeFastaFiles galGal6.fa` \
`#--sjdbGTFfile galGal6.ncbiRefSeq.gtf \
--genomeFastaFiles <(zcat galGal6.fa.gz) \
--sjdbGTFfile <(zcat galGal6.ncbiRefSeq.gtf.gz) \
--sjdbOverhang 99 \
--genomeSAindexNbases 10 


srun -n1 --pty --partition=angsd_class --mem=24G bash -i



# Align reads for E8_rep1
cd Alignment

STAR --runMode alignReads \
--outFilterMultimapNmax 20 \
--runThreadN 8 \
--outFileNamePrefix E8_ \
--genomeDir ../genome/galGal6_STAR_Index \ 
--readFilesIn ../fastq/SRR1804237_1_trimmed.fq.gz ../fastq/SRR1804237_2_trimmed.fq.gz \
--readFilesCommand zcat \
--outFilterType BySJout \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI AS nM MD


(angsd) [shp4022@cayuga-login1 Alignment]$ srun -n1 --pty --partition=angsd_class --mem=24G bash STAR_align.sh

(angsd) [shp4022@cayuga-login1 Alignment]$ srun -n1 --pty --partition=angsd_class --mem=24G bash STAR_align.sh
srun: job 115083 queued and waiting for resources
srun: job 115083 has been allocated resources
/athena/angsd/scratch/mef3005/share/envs/angsd/bin/STAR-avx2 --runMode alignReads --outFilterMultimapNmax 20 --runThreadN 8 --outFileNamePrefix E8_ --genomeDir ../genome/galGal6_STAR_Index --readFilesIn ../fastq/SRR1804237_1_trimmed.fq.gz ../fastq/SRR1804237_2_trimmed.fq.gz --readFilesCommand zcat --outFilterType BySJout --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM MD
STAR version: 2.7.11a   compiled: 2023-09-15T02:58:53+0000 :/opt/conda/conda-bld/star_1694746407721/work/source
Feb 23 10:31:42 ..... started STAR run
Feb 23 10:31:42 ..... loading genome
Feb 23 10:31:50 ..... started mapping

EXITING because of FATAL ERROR in reads input: quality string length is not equal to sequence length
@SRR1804237.72003183
+
  .)<<AFFF..FFF7F)FF7F<FF
SOLUTION: fix your fastq file

Feb 23 11:21:55 ...... FATAL ERROR, exiting
srun: error: c0010: task 0: Exited with exit code 104

# people online are saying that trimming may have impacted the format so I will try this to align again on the untrimmed data


(angsd) [shp4022@cayuga-login1 Alignment]$ samtools index E8_Aligned.sortedByCoord.out.bam

(base) shaunp@MAC306097 Classwork % scp -i ~/.ssh/Private_Key.txt shp4022@cayuga-login1.cac.cornell.edu:/athena/angsd/scratch/shp4022/project/Alignment/E8_Aligned.sortedByCoord.out.bam  /Users/shaunp/Desktop/Weill_Cornell_Graduate/Grad_School/SPRING_2024/Analysis_Next-Gen_Sequencing_Data/Project/E8_Aligned.sortedByCoord.out.bam
E8_Aligned.sortedByCoord.out.bam                                        100% 6149MB   5.9MB/s   17:18    
(base) shaunp@MAC306097 Classwork % 
(base) shaunp@MAC306097 Classwork % 
(base) shaunp@MAC306097 Classwork % scp -i ~/.ssh/Private_Key.txt shp4022@cayuga-login1.cac.cornell.edu:/athena/angsd/scratch/shp4022/project/Alignment/E8_Aligned.sortedByCoord.out.bam.bai  /Users/shaunp/Desktop/Weill_Cornell_Graduate/Grad_School/SPRING_2024/Analysis_Next-Gen_Sequencing_Data/Project/E8_Aligned.sortedByCoord.out.bam.bai
E8_Aligned.sortedByCoord.out.bam.bai                                    100% 1538KB   2.5MB/s   00:00    







#STAR --runMode alignReads \
#     --outFilterMultimapNmax 20 \
#     --runThreadN 8 \
#     --outFileNamePrefix E8_ \
#     --genomeDir ../genome/galGal6_STAR_Index \
#     --readFilesIn ../fastq/SRR1804237_1.fastq.gz ../fastq/SRR1804237_2.fastq.gz \
#     --readFilesCommand zcat \
#     --outFilterType BySJout \
#     --outSAMtype BAM SortedByCoordinate \
#     --outSAMattributes NH HI AS nM MD


# E8 retina rep 2 Align reads

STAR --runMode alignReads \
--outFilterMultimapNmax 20 \
--runThreadN 8 \
--outFileNamePrefix E8_retina_rep2_ \
--genomeDir ../genome/galGal6_STAR_Index \
--readFilesIn ../fastq/SRR1804238_1.fastq.gz ../fastq/SRR1804238_2.fastq.gz \
--readFilesCommand zcat \
--outFilterType BySJout \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI AS nM MD

# E16 retina rep 1 Align reads

STAR --runMode alignReads \
--outFilterMultimapNmax 20 \
--runThreadN 8 \
--outFileNamePrefix E16_retina_rep1_ \
--genomeDir ../genome/galGal6_STAR_Index \
--readFilesIn ../fastq/SRR1804235_1.fastq.gz ../fastq/SRR1804235_2.fastq.gz \
--readFilesCommand zcat \
--outFilterType BySJout \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI AS nM MD
# E16 retina rep 2 Align reads

STAR --runMode alignReads \
--outFilterMultimapNmax 20 \
--runThreadN 8 \
--outFileNamePrefix E16_retina_rep2_ \
--genomeDir ../genome/galGal6_STAR_Index \
--readFilesIn ../fastq/SRR1804236_1.fastq.gz ../fastq/SRR1804236_2.fastq.gz \
--readFilesCommand zcat \
--outFilterType BySJout \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI AS nM MD
# E18 retina rep 1 Align reads

STAR --runMode alignReads \
--outFilterMultimapNmax 20 \
--runThreadN 8 \
--outFileNamePrefix E18_retina_rep1_ \
--genomeDir ../genome/galGal6_STAR_Index \
--readFilesIn ../fastq/SRR1804239_1.fastq.gz ../fastq/SRR1804239_2.fastq.gz \
--readFilesCommand zcat \
--outFilterType BySJout \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI AS nM MD


# E18 retina rep 2 Align reads

STAR --runMode alignReads \
--outFilterMultimapNmax 20 \
--runThreadN 8 \
--outFileNamePrefix E18_retina_rep2_ \
--genomeDir ../genome/galGal6_STAR_Index \
--readFilesIn ../fastq/SRR1804240_1.fastq.gz ../fastq/SRR1804240_2.fastq.gz \
--readFilesCommand zcat \
--outFilterType BySJout \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI AS nM MD


# E18 cornea rep 1 Align reads

STAR --runMode alignReads \
--outFilterMultimapNmax 20 \
--runThreadN 8 \
--outFileNamePrefix E18_cornea_rep1_ \
--genomeDir ../genome/galGal6_STAR_Index \
--readFilesIn ../fastq/SRR1804241_1.fastq.gz ../fastq/SRR1804241_2.fastq.gz \
--readFilesCommand zcat \
--outFilterType BySJout \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI AS nM MD

# E18 cornea rep 2 Align reads

STAR --runMode alignReads \
--outFilterMultimapNmax 20 \
--runThreadN 8 \
--outFileNamePrefix E18_cornea_rep2_ \
--genomeDir ../genome/galGal6_STAR_Index \
--readFilesIn ../fastq/SRR1804242_1.fastq.gz ../fastq/SRR1804242_2.fastq.gz \
--readFilesCommand zcat \
--outFilterType BySJout \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI AS nM MD


qualimap bamqc -bam *_Aligned.sortedByCoord.out.bam -outformat pdf -outfile chI_BWA_ERR458878_report.pdf -outdir /athena/angsd/scratch/shp4022/angsd_hw/Homework5/

srun -n1 --pty --partition=angsd_class --cpus-per-task 2 --mem=24G bash -i


