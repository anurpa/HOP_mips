#merge de-multiplexed, UMI trimmed reads using PEAR

/src/PEAR/pear-0.9.6/pear -y 10G -j 4 -f sample.R1.fq -r sample.R2.fq -o sample.merge.fq

#align merged reads using bwa mem plus add read group info, remove unmapped reads, and sort resulting bam

/src/bwa/bwa-0.7.12/bwa mem -M -t 12 -R '@RG\tID:Control_DNA\tSM:sample.id\tPL:illumina\tLB:MIPs\tPU:NextSeq500' /refs/GRCh37/Homo_sapiens_assembly19.fasta_0.7.12 sample.merge.fq.assembled.fastq  | /src/samtools/samtools-1.1/samtools view -F 4 -b - | /src/samtools/samtools-1.1/samtools sort - sample.merged.initial

#convert sam to bam to make file readable with our current script
#you may not need to do this step if using a different tool to collapse
/src/samtools/samtools-1.1/samtools view sample.initial.bam > sample.initial.sam

#collapse to unique reads
python /scripts/mip_extract_unique_tags.py sample.initial.sam > sample.unique.sam

#trim arms
#note, need to have an arm info file. This is a tab-delimited file with chromosome, start position (start of the ligation if negative strand MIP, start of extension arm if positive strand mip), strand (0 if negative strand mip, 16 if positive strand mip. These are flags in Sam file), length ligation arm, length extension arm. 
#i is the insert size. For these mips 248 is the correct insert size
python /scripts/mip_trim_arms_se.py -s sample.unique.sam -i 248 -q 1 Cancer_Mips_Arm_Info.txt > sample.trimmed.unique.sam

#convert the unique, trimmed sam back to bam 
#first create a sequence dictionary from the human reference genome
java -Xmx10g -jar /src/picard/picard-tools-1.110/CreateSequenceDictionary.jar OUTPUT=dict.sam R=/refs/GRCh37/Homo_sapiens_assembly19.fasta

#re-add readgroups
at dict.sam sample.trimmed.unique.sam | java -Xmx10g -jar /src/picard/picard-tools-1.110/SortSam.jar I=/dev/stdin O=/dev/stdout SO=coordinate VALIDATION_STRINGENCY=SILENT | java -Xmx10g -jar /src/picard/picard-tools-1.110/AddOrReplaceReadGroups.jar I=/dev/stdin O=sample.named.trimmed.unique.400k.sam RGID=Test_DNA RGLB=Cancer_MIPs RGPL=PE300Mid RGPU=NextSeq500 RGSM=$name VALIDATION_STRINGENCY=SILENT

#convert sam to bam
/src/samtools/samtools-1.1/samtools view -bS -t \
/refs/GRCh37/Homo_sapiens_assembly19.fasta.fai \
sample.named.trimmed.unique.sam | /src/samtools/samtools-1.1/samtools sort \
- sample.trimmed.unique

#build a bam index
/src/samtools/samtools-1.1/samtools index sample.trimmed.unique.bam

#go ahead and delete intermediary files
#sample.unique.sam, sample.trimmed.unique.sam, sample.named.trimmed.unique.sam
#proceed to variant calling
