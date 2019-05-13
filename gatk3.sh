#Snakefile to run gatk haplotype caller 

SAMPLES=["NA12878","NA02633","NA11254","NA13705","NA13715","NA13961","NA14090","NA14094","NA14170",\
"NA14622","NA14623","NA14624","NA14626","NA14634","NA14636","NA14637","NA14638","NA14639","NA14684",\
"NA16658","NA16660","NA21849"]

rule all:
  input: expand("/mips/analysis/{sample}_genotype.vcf",sample=SAMPLES)
  
rule gatk:
  input:"/mips/analysis/{sample}.final.bam"
  output:"/mips/analysis/{sample}.g.vcf"
  params: name="gatk_{sample}"
  
  run:
    shell("""java -Xmx20G -jar /tools/GenomeAnalysisTK.jar -T HaplotypeCaller -mbq 20  \
    -R ../GRch37/Homo_sapiens_assembly19.fasta -I {input} -o {output} -ERC GVCF --disable_auto_index_creation_and_locking_when_reading_rods \
    --dbsnp /mips/GRch37/DBSNP_150_uw.vcf """)
    
rule genotype:
  input:"/mips/analysis/{sample}.g.vcf"
  output:"/mips/analysis/{sample}_genotype.vcf"
  params: name="gatk_{sample}"

  run:
    shell("""java -Xmx20G -jar /tools/GenomeAnalysisTK.jar -T GenotypeGVCFs   \
    -R ../GRch37/Homo_sapiens_assembly19.fasta -I {input} -o {output} \
    --dbsnp /mips/GRch37/DBSNP_150_uw.vcf """)    
