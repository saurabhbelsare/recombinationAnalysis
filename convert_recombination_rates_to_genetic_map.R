### This script converts recombination rates (r), say generated from pyrho, 
### to a genetic map.

library(readr)

chr_num <- 20

recombination_rates_file <- paste0('/Users/saurabhbelsare/Downloads/baboon_data/recombination_rates/baboon_pyrho_chr',chr_num,'.rmap.bed')
recombination_rates <- read_table(recombination_rates_file)
genetic_map_file <- paste0('/Users/saurabhbelsare/Downloads/baboon_data/recombination_rates/baboon_pyrho_chr',chr_num,'.txt')

#genetic_map_file <- 'Downloads/hg38/YRI/YRI_recombination_map_hapmap_format_hg38_chr_1.txt'
#genetic_map <- read_delim(genetic_map_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE)

num_sites <- length(recombination_rates$Start) + 1

Chromosome <- rep(chr_num, num_sites)
`Position(bp)` <- c(recombination_rates$Start, recombination_rates$End[length(recombination_rates$End)])
`Rate(cM/Mb)` <- c(recombination_rates$reco_rate_per_base_per_generation * 10^8, 0.00)
`Map(cM)` <- rep(0.00, length(recombination_rates$End) + 1)
for (i in 1:length(recombination_rates$Start)){
  `Map(cM)`[i + 1] <- recombination_rates$reco_rate_per_base_per_generation[i] * (recombination_rates$End[i] - recombination_rates$Start[i]) * 100 + `Map(cM)`[i] 
}

genetic_map <- data.frame(Chromosome, `Position(bp)`, `Rate(cM/Mb)`, `Map(cM)`)

write.table(genetic_map, genetic_map_file, row.names = FALSE, sep = '\t')


rec_rates <- c()
for (chr_num in 1:20){
  recombination_rates_file <- paste0('/Users/saurabhbelsare/Downloads/baboon_data/recombination_rates/baboon_pyrho_chr',chr_num,'.rmap.bed')
  recombination_rates <- read_table(recombination_rates_file)
  rec_rates <- c(rec_rates, print(mean(recombination_rates$reco_rate_per_base_per_generation)))
}
print(rec_rates)

mean(rec_rates)