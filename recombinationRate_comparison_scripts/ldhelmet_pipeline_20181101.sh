##### New LDhelmet pipeline Nov. 1, 2018

### Extract variants 
# Include 24 olive baboon founders and filter out monomorphic sites and singletons 
# (including singleton REF), and thin so that no two SNPs are closer than 10 bp
# Execute: ./run_vcftools_VarsForBeagle_20181029.sh

cd /media/walllab/jacqueline/baboon/beagle

VCFTOOLS=/media/walllab/jacqueline/programs/vcftools_0.1.13/bin/vcftools
FILE=Joint_Calling_High_and_Locoverage_eff_AddAnnot_HighCov_Masked_Filtered_chrALL_PASS_rehead.vcf.gz
KEEP=PANfounder_24pure.list
OUT=24PANfounders_beagle

for CHR in {1..20}; do
	${VCFTOOLS} --gzvcf ${FILE} --keep ${KEEP} --chr ${CHR} --out ${OUT}_chr${CHR} \
	--recode --remove-indels --remove-filtered-all --non-ref-ac-any 1 \
	--min-alleles 2 --max-alleles 2 --mac 2 --thin 10 &
done


### Phase with Beagle
# Set ne=30000, consistent with SMC++ recent population size for olive baboons
# Execute: for i in {1..20}; do ./run_beagle_20181101.sh ${i} & done

BEAGLE=/media/walllab/jacqueline/programs/beagle5.0/beagle.03Jul18.40b.jar
cd /media/walllab/jacqueline/baboon/beagle
FILE=${1}

java -jar -Xmx16g ${BEAGLE} \
impute=false \
ne=30000 \
nthreads=1 \
gt=${FILE} \
out=${FILE%.vcf.gz}_phased &>${FILE%.vcf.gz}_phased.log


### Index phased VCF files

for i in *phased.vcf.gz; do tabix -p vcf ${i} & done


### Convert phased VCF files to fasta files

VCF2FASTA=/media/walllab/jacqueline/programs/vcflib/bin/vcf2fasta
REF=/media/walllab/jacqueline/baboon/PapAnu2.0/Papio_anubis.PapAnu2.0.dna.toplevel.fa

for i in *phased.vcf.gz; do ${VCF2FASTA} --reference ${REF} --prefix ${i%.vcf.gz}_ --default-ploidy 2 ${i} & done


### Combine fasta files to make one file per chromosome

OUT=/media/walllab/jacqueline/baboon/ldhelmet/24_baboons/0_fastafiles
for chr in {1..20}; do cat 24PANfounders_beagle_chr${chr}.recode_Annot_phased_*_${chr}:*.fa > ${OUT}/24PANfounders_chr${chr}.fa & done


### Run LDhelmet find_confs
# Execute: ./run_ldhelmet_1_findconfs_20181101.sh

LDHELMET=/media/walllab/jacqueline/programs/LDhelmet_v1.10/ldhelmet
FASTADIR=/media/walllab/jacqueline/baboon/ldhelmet/24_baboons/0_fastafiles
cd /media/walllab/jacqueline/baboon/ldhelmet/24_baboons

for chr in {1..20}; do ${LDHELMET} find_confs --num_threads 2 -w 50 -o 1_confs/chr${chr}.conf ${FASTADIR}/24PANfounders_chr${chr}.fa & done


### Run LDhelmet table_gen
# Set theta (-t) to 0.0008, approximation of Watterson's theta calculated from data
# Execute: for i in {1..20}; do ./run_ldhelmet_2_table_gen_20181101.sh & done

LDHELMET=/media/walllab/jacqueline/programs/LDhelmet_v1.10/ldhelmet
CHR=${1}
cd /media/walllab/jacqueline/baboon/ldhelmet/24_baboons
CONF=1_confs/24PANfounders_chr${CHR}.conf
LK=2_tables/24PANfounders_chr${CHR}.lk
LOG=2_tables/24PANfounders_chr${CHR}.lk.log

date > ${LOG}
${LDHELMET} table_gen --num_threads 2 -t 0.0008 -r 0.0 0.1 10.0 1.0 100.0 -c ${CONF} -o ${LK} &>> ${LOG}
date >> ${LOG}


### Run LDhelmet pade
# Set theta (-t) as above
# Execute: for i in {1..20}; do ./run_ldhelmet_3_pade_20181101.sh ${i} & done

LDHELMET=/media/walllab/jacqueline/programs/LDhelmet_v1.10/ldhelmet
CHR=${1}
cd /media/walllab/jacqueline/baboon/ldhelmet/24_baboons
CONF=1_confs/24PANfounders_chr${CHR}.conf
PADE=3_pade/24PANfounders_chr${CHR}.pade
LOG=3_pade/24PANfounders_chr${CHR}.pade.log

date > ${LOG}
${LDHELMET} pade --num_threads 3 -t 0.0008 -x 11 -c ${CONF} -o ${PADE} &>> ${LOG}
date >> ${LOG}


### Run LDhelmet rjmcmc
# Use multiple block penalties
# Execute: for i in {1..20}; do ./run_ldhelmet_4_rjmcmc_20181101.sh ${1} 5.0 & done
# Execute: for i in {1..20}; do ./run_ldhelmet_4_rjmcmc_20181101.sh ${1} 25.0 & done
# Execute: for i in {1..20}; do ./run_ldhelmet_4_rjmcmc_20181101.sh ${1} 50.0 & done

LDHELMET=/media/walllab/jacqueline/programs/LDhelmet_v1.10/ldhelmet
CHR=${1}
BLOCK=${2}
cd /media/walllab/jacqueline/baboon/ldhelmet/24_baboons
FASTA=0_fastafiles/24PANfounders_chr${CHR}.fa
LK=2_tables/24PANfounders_chr${CHR}.lk
PADE=3_pade/24PANfounders_chr${CHR}.pade
POST=4_rjmcmc/24PANfounders_chr${CHR}_bp${BLOCK}.post
LOG=4_rjmcmc/24PANfounders_chr${CHR}_bp${BLOCK}.post.log

THR=3
WIN=50
BURN=100000
ITER=1000000

date > ${LOG}
${LDHELMET} rjmcmc --num_threads ${THR} -w ${WIN} -b ${BLOCK} --burn_in ${BURN} -n ${ITER} -s ${FASTA} -l ${LK} -p ${PADE} -o ${POST} &>> ${LOG}
date >> ${LOG}


### Run LDhelmet post_to_text

LDHELMET=/media/walllab/jacqueline/programs/LDhelmet_v1.10/ldhelmet
cd /media/walllab/jacqueline/baboon/ldhelmet/24_baboons/4_rjmcmc

for i in *.post; do ${LDHELMET} post_to_text -m -p 0.025 -p 0.50 -p 0.975 -o ${i}.txt ${i} ; done
for i in *.post.txt; do gzip ${i} ; done



