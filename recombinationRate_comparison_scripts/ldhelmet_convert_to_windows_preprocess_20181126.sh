# Converting LDhelmet results to window-based estimates of rho

# Input data files are gzipped human-readable versions of ldhelmet output
# Format:
# Col. 1: left_snp
# Col. 2: right_snp
# Col. 3: mean rho/bp
# Col.4+: other stats


# Convert to bed files
# New format will be:
# Col. 1: chromosome (scaffold) name
# Col. 2: left_snp
# Col. 3: right_snp
# Col. 4: mean rho/bp
for i in {1..20}; do zcat 24PANfounders_chr${i}_bp5.post.txt.gz | tail -n+4 | awk -v var=${i} '{printf "%s\t%s\t%s\t%s\n",var,$1,$2,$3}' | bgzip > 24PANfounders_chr${i}_bp5.post.bed.gz; done


# Index bed files with tabix
for i in *bed.gz; do tabix -p bed ${i} ; done


# Run script to convert results to window-based estimates
SCRIPT=/media/walllab/jacqueline/scripts/ldhelmet/ldhelmet_convert_to_windows_20181126.py
WINSIZE=100000
STEPSIZE=100000
for i in {1..20}; do python ${SCRIPT} 24PANfounders_chr${i}_bp5.post.bed.gz ${WINSIZE} ${STEPSIZE} ${i}; done
