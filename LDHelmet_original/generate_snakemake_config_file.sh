printf "threads: 6\n"
printf "\n"
printf "ldhelmet_params:\n"
printf "\twindow: 50\n"
printf "\trho_grid: 0.0 0.1 10.0 1.0 100.0\n"
printf "\tnum_pade_coeff: 11\n"
printf "\tblock_penalty: 50.0\n"
printf "\tburn_in: 100000\n"
printf "\ttotal_iter: 1000000\n"
printf "\tlow_percentile: 0.025\n"
printf "\tmid_percentile: 0.5\n"
printf "\thigh_percentile: 0.975\n"
printf "\n"
printf "theta_values:\n"

chromosomeLengthArray=(248956422 242193529 198295559 190214555 181538259 170805979 159345973 145138636 138394717 133797422 135086622 133275309 114364328 107043718 101991189 90338345 83257441 80373285 58617616 64444167 46709983 50818468)

for region in "SAS"
do
    for ((i=1;i<23;i++))
    do
        numSegregatingSites="$(grep -v "#" ~/data/1000genomes/phase3/data/phase3_GRCh38_cleanedup_removeDuplicates/${region}/${region}.chr${i}_GRCh38.genotypes.20170504.cleanedup_removeDuplicates.recode.vcf | wc -l)"
        numIndividuals="$(awk '{if ($1 == "#CHROM"){print NF-9; exit}}' ~/data/1000genomes/phase3/data/phase3_GRCh38_cleanedup_removeDuplicates/${region}/${region}.chr${i}_GRCh38.genotypes.20170504.cleanedup_removeDuplicates.recode.vcf)"
        theta=$(~/software/anaconda3/bin/python wattersonTheta.py ${numSegregatingSites} ${numIndividuals} ${chromosomeLengthArray[$i-1]})
        printf "\t${i}: %.4f\n" "${theta}"
    done
done
