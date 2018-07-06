numChr=1
chromosomes=[str(i) for i in range(1,numChr+1)]
#regions=["AFR", "AMR", "EUR", "EAS", "SAS"]
regions=["SAS"]

rule all:
    input:
        expand("output/chr{chrNum}/{region}/freqs.txt", chrNum=chromosomes, region=regions)

#rule make_directories:
#    output:
#        "output/chr{chrNum}/{region}"
#    shell:
#        "mkdir {output}"

rule get_vcf:
    output:
        "output/chr{chromosomes}/{regions}/{regions}.chr{chromosomes}_GRCh38.genotypes.20170504.cleanedup_removeDuplicates.recode.vcf"
    params:
        "output/chr{chromosomes}/{regions}/{regions}.chr{chromosomes}_GRCh38.genotypes.20170504.cleanedup_removeDuplicates"
    shell:
        """
        vcftools --vcf /media/walllab/saurabh/data/1000genomes/vcf/phase3_GRCh38_cleanedup_removeDuplicates/ALL.chr{chromosomes}_GRCh38.genotypes.20170504.cleanedup_removeDuplicates.vcf --keep individualData/{regions}_1KG.list --remove individualData/{regions}_10X.list --chr {chromosomes} --recode --out {params}
        """

rule convert_to_ldhat_format:
    input:
        "output/chr{chromosomes}/{regions}/{regions}.chr{chromosomes}_GRCh38.genotypes.20170504.cleanedup_removeDuplicates.recode.vcf"
    output:
        "output/chr{chromosomes}/{regions}/{regions}.chr{chromosomes}_GRCh38_cleanedup_removeDuplicates_first1869.ldhat.sites",
        "output/chr{chromosomes}/{regions}/{regions}.chr{chromosomes}_GRCh38_cleanedup_removeDuplicates_first1869.ldhat.locs"
    params:
        "output/chr{chromosomes}/{regions}/{regions}.chr{chromosomes}_GRCh38_cleanedup_removeDuplicates_first1869"
    shell:
        """
        vcftools --vcf {input} --chr {chromosomes} --to-bp 634115 --ldhat --phased --out {params}
        """

rule run_ldhat_pairwise:
    input:
        sitesFile="output/chr{chromosomes}/{regions}/{regions}.chr{chromosomes}_GRCh38_cleanedup_removeDuplicates_first1869.ldhat.sites",
        locsFile="output/chr{chromosomes}/{regions}/{regions}.chr{chromosomes}_GRCh38_cleanedup_removeDuplicates_first1869.ldhat.locs"
    output:
        freqsOut="output/chr{chromosomes}/{regions}/freqs.txt",
        tableOut="output/chr{chromosomes}/{regions}/type_table.txt"
    params:
        activeDir="output/chr{chromosomes}/{regions}/"
    shell:
        """
        pairwise -seq {input.sitesFile} -loc {input.locsFile}
        """





