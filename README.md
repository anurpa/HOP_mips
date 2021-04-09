# HOP_mips

This pipeline was developed to validate library prep using Molecular Inversion Probes(MIPs). 

# Workflow
1. Merge and trim UMI of de-multiplexed reads using PEAR.
2. Align merged reads using bwa mem plus add read group info, remove unmapped reads, and sort resulting bam.
3. Convert Bam to Sam file.
4. Collapse to unique reads.
5. Trim MIPs trim arms.
6. Convert the unique, trimmed sam back to bam 
7. Build a bam index
8. Germline mutation callling using GATK.
