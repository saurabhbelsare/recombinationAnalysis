numChr=1
chromosomes=[str(i) for i in range(1,numChr+1)]
#regions=["AFR", "AMR", "EUR", "EAS", "SAS"]
regions=["SAS"]

rule all:
    input:
        expand("output/chr{chrNum}/{region}/{region}.chr{chrNum}_GRCh38.genotypes.20170504.cleanedup_removeDuplicates.vcf", chrNum=chromosomes, region=regions)

#rule make_directories:
#    output:
#        "output/chr{chrNum}/{region}"
#    shell:
#        "mkdir {output}"

rule get_vcf:
    output:
        "output/chr{chromosomes}/{regions}/{regions}.chr{chromosomes}_GRCh38.genotypes.20170504.cleanedup_removeDuplicates.vcf"
    shell:
        """
        vcftools --vcf /media/walllab/saurabh/data/1000genomes/vcf/phase3_GRCh38_cleanedup/ALL.chr{chromosomes}_GRCh38.genotypes.20170504.cleanedup.recode.vcf --keep individualData/{regions}_1KG.list --remove individualData/{regions}_10X.list --chr {chromosomes} --recode --out {output}
        """
