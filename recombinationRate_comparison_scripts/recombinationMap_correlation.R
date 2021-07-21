library(readr)

## BEB

pop_1KG_firstHalf_chr22_GRCh38_cleanedup_removeDuplicates <- read_table2("/home/user/workfiles/recombinationAnalysis/human/MAA/30Indv/unzippedFiles/MAA_gSAsia_30Indv.chr22_GRCh38_cleanedup_removeDuplicates.txt")
pop_1KG_secondHalf_chr22_GRCh38_cleanedup_removeDuplicates <- read_table2("/home/user/workfiles/recombinationAnalysis/human/MAA/30Indv7/unzippedFiles/MAA_gSAsia_30Indv7.chr22_GRCh38_cleanedup_removeDuplicates.txt")

segmentLength <- 100000
chrLength <- 50818468
numSegments <- ceiling(chrLength/segmentLength)

snpCounts1 <- rep(0,numSegments)
rhoSum1 <- rep(0,numSegments)
rhoMean1 <- rep(0,numSegments)

for(i in 1:length(pop_1KG_firstHalf_chr22_GRCh38_cleanedup_removeDuplicates$left_snp)){
  idx <- ceiling(pop_1KG_firstHalf_chr22_GRCh38_cleanedup_removeDuplicates$left_snp[i]/segmentLength)
  snpCounts1[idx] <- snpCounts1[idx] + 1
  rhoSum1[idx] <- rhoSum1[idx] + pop_1KG_firstHalf_chr22_GRCh38_cleanedup_removeDuplicates$mean[i]
}

rhoMean1 <- rhoSum1/snpCounts1

snpCounts2 <- rep(0,numSegments)
rhoSum2 <- rep(0,numSegments)
rhoMean2 <- rep(0,numSegments)

for(i in 1:length(pop_1KG_secondHalf_chr22_GRCh38_cleanedup_removeDuplicates$left_snp)){
  idx <- ceiling(pop_1KG_secondHalf_chr22_GRCh38_cleanedup_removeDuplicates$left_snp[i]/segmentLength)
  snpCounts2[idx] <- snpCounts2[idx] + 1
  rhoSum2[idx] <- rhoSum2[idx] + pop_1KG_secondHalf_chr22_GRCh38_cleanedup_removeDuplicates$mean[i]
}

rhoMean2 <- rhoSum2/snpCounts2

#cor(rhoMean1[749:(length(rhoMean1)-1)],rhoMean2[749:(length(rhoMean2)-1)], method = "spearman") ### 25kb
#cor(rhoMean1[367:(length(rhoMean1)-1)],rhoMean2[367:(length(rhoMean2)-1)], method = "spearman") ### 50kb
#cor(rhoMean1[164:(length(rhoMean1)-1)],rhoMean2[164:(length(rhoMean2)-1)], method = "spearman") ### 100kb
#cor(rhoMean1[66:length(rhoMean1)],rhoMean2[66:length(rhoMean2)], method = "spearman") ### 250kb
#or(rhoMean1[31:length(rhoMean1)],rhoMean2[31:length(rhoMean2)], method = "spearman") ### 500kb
#cor(rhoMean1[21:length(rhoMean1)],rhoMean2[21:length(rhoMean2)], method = "spearman") ### 750kb
#cor(rhoMean1[16:length(rhoMean1)],rhoMean2[16:length(rhoMean2)], method = "spearman") ### 1Mb
cor(rhoMean1[7:length(rhoMean1)],rhoMean2[7:length(rhoMean2)], method = "spearman") ### 2.5Mb
#cor(rhoMean1[3:length(rhoMean1)],rhoMean2[3:length(rhoMean2)], method = "spearman") ### 5Mb

## GIH

GIH_1KG_firstHalf_chr22_GRCh38_cleanedup_removeDuplicates <- read_table2("Desktop/maf_0_03/GIH_1KG_firstHalf.chr22_GRCh38_cleanedup_removeDuplicates.txt")
GIH_1KG_secondHalf_chr22_GRCh38_cleanedup_removeDuplicates <- read_table2("Desktop/maf_0_03/GIH_1KG_secondHalf.chr22_GRCh38_cleanedup_removeDuplicates.txt")

segmentLength <- 5000000
chrLength <- 50818468
numSegments <- ceiling(chrLength/segmentLength)

snpCounts1 <- rep(0,numSegments)
rhoSum1 <- rep(0,numSegments)
rhoMean1 <- rep(0,numSegments)

for(i in 1:length(GIH_1KG_firstHalf_chr22_GRCh38_cleanedup_removeDuplicates$left_snp)){
  idx <- ceiling(GIH_1KG_firstHalf_chr22_GRCh38_cleanedup_removeDuplicates$left_snp[i]/segmentLength)
  snpCounts1[idx] <- snpCounts1[idx] + 1
  rhoSum1[idx] <- rhoSum1[idx] + GIH_1KG_firstHalf_chr22_GRCh38_cleanedup_removeDuplicates$mean[i]
}

rhoMean1 <- rhoSum1/snpCounts1

snpCounts2 <- rep(0,numSegments)
rhoSum2 <- rep(0,numSegments)
rhoMean2 <- rep(0,numSegments)

for(i in 1:length(GIH_1KG_secondHalf_chr22_GRCh38_cleanedup_removeDuplicates$left_snp)){
  idx <- ceiling(GIH_1KG_secondHalf_chr22_GRCh38_cleanedup_removeDuplicates$left_snp[i]/segmentLength)
  snpCounts2[idx] <- snpCounts2[idx] + 1
  rhoSum2[idx] <- rhoSum2[idx] + GIH_1KG_secondHalf_chr22_GRCh38_cleanedup_removeDuplicates$mean[i]
}

rhoMean2 <- rhoSum2/snpCounts2

#cor(rhoMean1[749:(length(rhoMean1)-1)],rhoMean2[749:(length(rhoMean2)-1)], method = "spearman") ### 25kb
#cor(rhoMean1[367:(length(rhoMean1)-1)],rhoMean2[367:(length(rhoMean2)-1)], method = "spearman") ### 50kb
#cor(rhoMean1[164:(length(rhoMean1)-1)],rhoMean2[164:(length(rhoMean2)-1)], method = "spearman") ### 100kb
#cor(rhoMean1[66:length(rhoMean1)],rhoMean2[66:length(rhoMean2)], method = "spearman") ### 250kb
#cor(rhoMean1[31:length(rhoMean1)],rhoMean2[31:length(rhoMean2)], method = "spearman") ### 500kb
#cor(rhoMean1[21:length(rhoMean1)],rhoMean2[21:length(rhoMean2)], method = "spearman") ### 750kb
#cor(rhoMean1[16:length(rhoMean1)],rhoMean2[16:length(rhoMean2)], method = "spearman") ### 1Mb
#cor(rhoMean1[7:length(rhoMean1)],rhoMean2[7:length(rhoMean2)], method = "spearman") ### 2.5Mb
cor(rhoMean1[3:length(rhoMean1)],rhoMean2[3:length(rhoMean2)]) ### 5Mb























lengthsArray <- c(25000, 50000, 100000, 250000, 500000, 750000, 1000000, 2500000, 5000000)
corrsArray_BEB_chr22_pearson <- c(0.7559602, 0.7944265, 0.8399403, 0.8168536, 0.8939317, 0.87728, 0.9061186, 0.9563531, 0.9850219)
corrsArray_BEB_chr22_spearman <- c(0.9277741, 0.9334574, 0.9242976, 0.9102149, 0.9304135, 0.9181502, 0.9250965, 0.9357143, 0.9666667) # Spearman Coefficent
corrsArray_GIH_chr22_spearman <- c(0.9277846, 0.9395474, 0.9363626, 0.9388087, 0.9408322, 0.9117455, 0.9518662, 0.9571429, 0.9833333) # Spearman Coefficent
corrsArray_GIH_maf_0_03_chr22_spearman <- c(0.9416578, 0.9352084, 0.9213813, 0.9199875, 0.8950736, 0.9016500, 0.8970399, 0.8857143, 0.9333333) # Spearman Coefficent
corrsArray_GIH_maf_0_03_chr22_pearson <- c(0.6986246, 0.7421385,  0.8121344, 0.8802898, 0.8843928, 0.8731264, 0.9271284, 0.9271977, 0.9492897) # Pearson Coefficent


crossCorrArray_BEB_GIH_chr22_spearman <- c(0.9263735, 0.9359126, 0.9323334, 0.9045266, 0.9083542, 0.906535, 0.9454311, 0.8892857, 0.9833333)
crossCorrArray_BEB_GIH_chr22_pearson <- c(0.766283, 0.7971859, 0.8453062, 0.8957853, 0.9135784, 0.9212979, 0.9533464, 0.9598709, 0.9757279)

plot(lengthsArray,corrsArray_BEB_chr22_spearman, log = "x", col="green", ylim = c(0.88,0.99))
points(lengthsArray,corrsArray_GIH_chr22_spearman, col="red")
points(lengthsArray,crossCorrArray_BEB_GIH_chr22_spearman, col="blue")
legend(x="topleft", legend=c("BEB", "GIH", "BEB-GIH"), col=c("green", "red", "blue"), lty=1:3, cex=0.8)

plot(lengthsArray, corrsArray_GIH_maf_0_03_chr22_pearson, log = "x", col="red", ylim = c(0.65,0.99))
legend(x="topleft", legend=c("GIH maf 0.03"), col=c("red"), lty = 1, cex=0.8)



#write(rhoMean1, file = "~/Desktop/BEB_firstHalf_chr22_mean50kb.txt")
#write(rhoMean2, file = "~/Desktop/BEB_secondHalf_chr22_mean50kb.txt")














#########################################################################
# Averaging with Jacqueline's script
########################################################################

Indv1_chr1 <- read_delim("workfiles/recombinationAnalysis/human/MAA/30Indv/MAA_gSAsia_30Indv.chr1_GRCh38_cleanedup_removeDuplicates.bed.gz_100000win_100000step.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
Indv2_chr1 <- read_delim("workfiles/recombinationAnalysis/human/MAA/30Indv2/MAA_gSAsia_30Indv2.chr1_GRCh38_cleanedup_removeDuplicates.bed.gz_100000win_100000step.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
Indv3_chr1 <- read_delim("workfiles/recombinationAnalysis/human/MAA/30Indv2/MAA_gSAsia_30Indv2.chr1_GRCh38_cleanedup_removeDuplicates.bed.gz_100000win_100000step.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
Indv4_chr1 <- read_delim("workfiles/recombinationAnalysis/human/MAA/30Indv2/MAA_gSAsia_30Indv2.chr1_GRCh38_cleanedup_removeDuplicates.bed.gz_100000win_100000step.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
Indv5_chr1 <- read_delim("workfiles/recombinationAnalysis/human/MAA/30Indv2/MAA_gSAsia_30Indv2.chr1_GRCh38_cleanedup_removeDuplicates.bed.gz_100000win_100000step.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
Indv6_chr1 <- read_delim("workfiles/recombinationAnalysis/human/MAA/30Indv2/MAA_gSAsia_30Indv2.chr1_GRCh38_cleanedup_removeDuplicates.bed.gz_100000win_100000step.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
Indv7_chr1 <- read_delim("workfiles/recombinationAnalysis/human/MAA/30Indv2/MAA_gSAsia_30Indv2.chr1_GRCh38_cleanedup_removeDuplicates.bed.gz_100000win_100000step.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
Indv8_chr1 <- read_delim("workfiles/recombinationAnalysis/human/MAA/30Indv2/MAA_gSAsia_30Indv2.chr1_GRCh38_cleanedup_removeDuplicates.bed.gz_100000win_100000step.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
Indv9_chr1 <- read_delim("workfiles/recombinationAnalysis/human/MAA/30Indv2/MAA_gSAsia_30Indv2.chr1_GRCh38_cleanedup_removeDuplicates.bed.gz_100000win_100000step.txt", "\t", escape_double = FALSE, trim_ws = TRUE)




plot(Indv1_chr1$win_rho[1:2481]/sum(Indv1_chr1$win_rho[1:2481]))
plot(Indv1_chr1$win_rho[1:2481]/sum(Indv1_chr1$win_rho[1:2481]), Indv2_chr1$win_rho/sum(Indv2_chr1$win_rho))
plot(Indv1_chr1$win_rho[2:2482]/sum(Indv1_chr1$win_rho[2:2482]), Indv2_chr1$win_rho/sum(Indv2_chr1$win_rho))


cor(Indv1_chr1$win_rho[1:2481],Indv2_chr1$win_rho)

mydf <- data.frame(1:2481,Indv1_chr1$win_rho[1:2481],Indv2_chr1[1:2481])
colnames(mydf) <- c("index","Indv1","Indv2")







MAA_gSAsia_win25000_chr1 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr1.txt")
MAA_gSAsia_win50000_chr1 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr1.txt")
MAA_gSAsia_win75000_chr1 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr1.txt")
MAA_gSAsia_win100000_chr1 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr1.txt")
MAA_gSAsia_win250000_chr1 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr1.txt")
MAA_gSAsia_win500000_chr1 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr1.txt")
MAA_gSAsia_win750000_chr1 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr1.txt")
MAA_gSAsia_win1000000_chr1 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr1.txt")
MAA_gSAsia_win2500000_chr1 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr1.txt")
MAA_gSAsia_win5000000_chr1 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr1.txt")
colnames(MAA_gSAsia_win25000_chr1) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr1) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr1) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr1) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr1) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr1) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr1) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr1) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr1) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr1) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

MAA_gSAsia_win25000_chr2 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr2.txt")
MAA_gSAsia_win50000_chr2 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr2.txt")
MAA_gSAsia_win75000_chr2 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr2.txt")
MAA_gSAsia_win100000_chr2 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr2.txt")
MAA_gSAsia_win250000_chr2 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr2.txt")
MAA_gSAsia_win500000_chr2 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr2.txt")
MAA_gSAsia_win750000_chr2 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr2.txt")
MAA_gSAsia_win1000000_chr2 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr2.txt")
MAA_gSAsia_win2500000_chr2 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr2.txt")
MAA_gSAsia_win5000000_chr2 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr2.txt")
colnames(MAA_gSAsia_win25000_chr2) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr2) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr2) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr2) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr2) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr2) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr2) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr2) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr2) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr2) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

MAA_gSAsia_win25000_chr3 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr3.txt")
MAA_gSAsia_win50000_chr3 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr3.txt")
MAA_gSAsia_win75000_chr3 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr3.txt")
MAA_gSAsia_win100000_chr3 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr3.txt")
MAA_gSAsia_win250000_chr3 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr3.txt")
MAA_gSAsia_win500000_chr3 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr3.txt")
MAA_gSAsia_win750000_chr3 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr3.txt")
MAA_gSAsia_win1000000_chr3 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr3.txt")
MAA_gSAsia_win2500000_chr3 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr3.txt")
MAA_gSAsia_win5000000_chr3 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr3.txt")
colnames(MAA_gSAsia_win25000_chr3) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr3) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr3) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr3) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr3) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr3) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr3) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr3) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr3) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr3) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

MAA_gSAsia_win25000_chr4 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr4.txt")
MAA_gSAsia_win50000_chr4 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr4.txt")
MAA_gSAsia_win75000_chr4 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr4.txt")
MAA_gSAsia_win100000_chr4 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr4.txt")
MAA_gSAsia_win250000_chr4 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr4.txt")
MAA_gSAsia_win500000_chr4 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr4.txt")
MAA_gSAsia_win750000_chr4 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr4.txt")
MAA_gSAsia_win1000000_chr4 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr4.txt")
MAA_gSAsia_win2500000_chr4 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr4.txt")
MAA_gSAsia_win5000000_chr4 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr4.txt")
colnames(MAA_gSAsia_win25000_chr4) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr4) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr4) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr4) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr4) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr4) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr4) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr4) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr4) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr4) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

MAA_gSAsia_win25000_chr5 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr5.txt")
MAA_gSAsia_win50000_chr5 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr5.txt")
MAA_gSAsia_win75000_chr5 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr5.txt")
MAA_gSAsia_win100000_chr5 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr5.txt")
MAA_gSAsia_win250000_chr5 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr5.txt")
MAA_gSAsia_win500000_chr5 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr5.txt")
MAA_gSAsia_win750000_chr5 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr5.txt")
MAA_gSAsia_win1000000_chr5 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr5.txt")
MAA_gSAsia_win2500000_chr5 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr5.txt")
MAA_gSAsia_win5000000_chr5 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr5.txt")
colnames(MAA_gSAsia_win25000_chr5) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr5) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr5) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr5) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr5) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr5) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr5) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr5) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr5) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr5) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

MAA_gSAsia_win25000_chr6 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr6.txt")
MAA_gSAsia_win50000_chr6 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr6.txt")
MAA_gSAsia_win75000_chr6 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr6.txt")
MAA_gSAsia_win100000_chr6 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr6.txt")
MAA_gSAsia_win250000_chr6 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr6.txt")
MAA_gSAsia_win500000_chr6 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr6.txt")
MAA_gSAsia_win750000_chr6 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr6.txt")
MAA_gSAsia_win1000000_chr6 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr6.txt")
MAA_gSAsia_win2500000_chr6 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr6.txt")
MAA_gSAsia_win5000000_chr6 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr6.txt")
colnames(MAA_gSAsia_win25000_chr6) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr6) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr6) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr6) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr6) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr6) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr6) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr6) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr6) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr6) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

MAA_gSAsia_win25000_chr7 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr7.txt")
MAA_gSAsia_win50000_chr7 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr7.txt")
MAA_gSAsia_win75000_chr7 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr7.txt")
MAA_gSAsia_win100000_chr7 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr7.txt")
MAA_gSAsia_win250000_chr7 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr7.txt")
MAA_gSAsia_win500000_chr7 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr7.txt")
MAA_gSAsia_win750000_chr7 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr7.txt")
MAA_gSAsia_win1000000_chr7 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr7.txt")
MAA_gSAsia_win2500000_chr7 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr7.txt")
MAA_gSAsia_win5000000_chr7 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr7.txt")
colnames(MAA_gSAsia_win25000_chr7) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr7) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr7) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr7) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr7) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr7) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr7) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr7) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr7) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr7) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

MAA_gSAsia_win25000_chr8 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr8.txt")
MAA_gSAsia_win50000_chr8 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr8.txt")
MAA_gSAsia_win75000_chr8 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr8.txt")
MAA_gSAsia_win100000_chr8 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr8.txt")
MAA_gSAsia_win250000_chr8 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr8.txt")
MAA_gSAsia_win500000_chr8 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr8.txt")
MAA_gSAsia_win750000_chr8 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr8.txt")
MAA_gSAsia_win1000000_chr8 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr8.txt")
MAA_gSAsia_win2500000_chr8 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr8.txt")
MAA_gSAsia_win5000000_chr8 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr8.txt")
colnames(MAA_gSAsia_win25000_chr8) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr8) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr8) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr8) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr8) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr8) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr8) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr8) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr8) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr8) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

MAA_gSAsia_win25000_chr9 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr9.txt")
MAA_gSAsia_win50000_chr9 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr9.txt")
MAA_gSAsia_win75000_chr9 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr9.txt")
MAA_gSAsia_win100000_chr9 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr9.txt")
MAA_gSAsia_win250000_chr9 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr9.txt")
MAA_gSAsia_win500000_chr9 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr9.txt")
MAA_gSAsia_win750000_chr9 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr9.txt")
MAA_gSAsia_win1000000_chr9 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr9.txt")
MAA_gSAsia_win2500000_chr9 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr9.txt")
MAA_gSAsia_win5000000_chr9 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr9.txt")
colnames(MAA_gSAsia_win25000_chr9) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr9) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr9) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr9) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr9) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr9) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr9) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr9) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr9) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr9) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

MAA_gSAsia_win25000_chr10 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr10.txt")
MAA_gSAsia_win50000_chr10 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr10.txt")
MAA_gSAsia_win75000_chr10 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr10.txt")
MAA_gSAsia_win100000_chr10 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr10.txt")
MAA_gSAsia_win250000_chr10 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr10.txt")
MAA_gSAsia_win500000_chr10 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr10.txt")
MAA_gSAsia_win750000_chr10 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr10.txt")
MAA_gSAsia_win1000000_chr10 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr10.txt")
MAA_gSAsia_win2500000_chr10 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr10.txt")
MAA_gSAsia_win5000000_chr10 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr10.txt")
colnames(MAA_gSAsia_win25000_chr10) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr10) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr10) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr10) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr10) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr10) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr10) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr10) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr10) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr10) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

MAA_gSAsia_win25000_chr11 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr11.txt")
MAA_gSAsia_win50000_chr11 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr11.txt")
MAA_gSAsia_win75000_chr11 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr11.txt")
MAA_gSAsia_win100000_chr11 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr11.txt")
MAA_gSAsia_win250000_chr11 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr11.txt")
MAA_gSAsia_win500000_chr11 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr11.txt")
MAA_gSAsia_win750000_chr11 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr11.txt")
MAA_gSAsia_win1000000_chr11 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr11.txt")
MAA_gSAsia_win2500000_chr11 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr11.txt")
MAA_gSAsia_win5000000_chr11 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr11.txt")
colnames(MAA_gSAsia_win25000_chr11) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr11) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr11) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr11) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr11) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr11) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr11) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr11) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr11) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr11) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

MAA_gSAsia_win25000_chr12 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr12.txt")
MAA_gSAsia_win50000_chr12 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr12.txt")
MAA_gSAsia_win75000_chr12 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr12.txt")
MAA_gSAsia_win100000_chr12 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr12.txt")
MAA_gSAsia_win250000_chr12 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr12.txt")
MAA_gSAsia_win500000_chr12 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr12.txt")
MAA_gSAsia_win750000_chr12 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr12.txt")
MAA_gSAsia_win1000000_chr12 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr12.txt")
MAA_gSAsia_win2500000_chr12 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr12.txt")
MAA_gSAsia_win5000000_chr12 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr12.txt")
colnames(MAA_gSAsia_win25000_chr12) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr12) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr12) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr12) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr12) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr12) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr12) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr12) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr12) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr12) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

MAA_gSAsia_win25000_chr13 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr13.txt")
MAA_gSAsia_win50000_chr13 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr13.txt")
MAA_gSAsia_win75000_chr13 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr13.txt")
MAA_gSAsia_win100000_chr13 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr13.txt")
MAA_gSAsia_win250000_chr13 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr13.txt")
MAA_gSAsia_win500000_chr13 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr13.txt")
MAA_gSAsia_win750000_chr13 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr13.txt")
MAA_gSAsia_win1000000_chr13 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr13.txt")
MAA_gSAsia_win2500000_chr13 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr13.txt")
MAA_gSAsia_win5000000_chr13 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr13.txt")
colnames(MAA_gSAsia_win25000_chr13) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr13) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr13) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr13) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr13) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr13) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr13) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr13) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr13) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr13) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")


MAA_gSAsia_win25000_chr14 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr14.txt")
MAA_gSAsia_win50000_chr14 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr14.txt")
MAA_gSAsia_win75000_chr14 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr14.txt")
MAA_gSAsia_win100000_chr14 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr14.txt")
MAA_gSAsia_win250000_chr14 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr14.txt")
MAA_gSAsia_win500000_chr14 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr14.txt")
MAA_gSAsia_win750000_chr14 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr14.txt")
MAA_gSAsia_win1000000_chr14 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr14.txt")
MAA_gSAsia_win2500000_chr14 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr14.txt")
MAA_gSAsia_win5000000_chr14 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr14.txt")
colnames(MAA_gSAsia_win25000_chr14) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr14) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr14) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr14) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr14) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr14) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr14) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr14) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr14) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr14) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

MAA_gSAsia_win25000_chr15 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr15.txt")
MAA_gSAsia_win50000_chr15 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr15.txt")
MAA_gSAsia_win75000_chr15 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr15.txt")
MAA_gSAsia_win100000_chr15 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr15.txt")
MAA_gSAsia_win250000_chr15 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr15.txt")
MAA_gSAsia_win500000_chr15 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr15.txt")
MAA_gSAsia_win750000_chr15 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr15.txt")
MAA_gSAsia_win1000000_chr15 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr15.txt")
MAA_gSAsia_win2500000_chr15 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr15.txt")
MAA_gSAsia_win5000000_chr15 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr15.txt")
colnames(MAA_gSAsia_win25000_chr15) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr15) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr15) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr15) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr15) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr15) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr15) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr15) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr15) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr15) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

MAA_gSAsia_win25000_chr16 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr16.txt")
MAA_gSAsia_win50000_chr16 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr16.txt")
MAA_gSAsia_win75000_chr16 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr16.txt")
MAA_gSAsia_win100000_chr16 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr16.txt")
MAA_gSAsia_win250000_chr16 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr16.txt")
MAA_gSAsia_win500000_chr16 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr16.txt")
MAA_gSAsia_win750000_chr16 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr16.txt")
MAA_gSAsia_win1000000_chr16 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr16.txt")
MAA_gSAsia_win2500000_chr16 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr16.txt")
MAA_gSAsia_win5000000_chr16 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr16.txt")
colnames(MAA_gSAsia_win25000_chr16) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr16) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr16) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr16) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr16) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr16) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr16) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr16) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr16) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr16) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

MAA_gSAsia_win25000_chr17 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr17.txt")
MAA_gSAsia_win50000_chr17 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr17.txt")
MAA_gSAsia_win75000_chr17 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr17.txt")
MAA_gSAsia_win100000_chr17 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr17.txt")
MAA_gSAsia_win250000_chr17 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr17.txt")
MAA_gSAsia_win500000_chr17 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr17.txt")
MAA_gSAsia_win750000_chr17 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr17.txt")
MAA_gSAsia_win1000000_chr17 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr17.txt")
MAA_gSAsia_win2500000_chr17 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr17.txt")
MAA_gSAsia_win5000000_chr17 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr17.txt")
colnames(MAA_gSAsia_win25000_chr17) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr17) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr17) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr17) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr17) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr17) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr17) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr17) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr17) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr17) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

MAA_gSAsia_win25000_chr18 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr18.txt")
MAA_gSAsia_win50000_chr18 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr18.txt")
MAA_gSAsia_win75000_chr18 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr18.txt")
MAA_gSAsia_win100000_chr18 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr18.txt")
MAA_gSAsia_win250000_chr18 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr18.txt")
MAA_gSAsia_win500000_chr18 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr18.txt")
MAA_gSAsia_win750000_chr18 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr18.txt")
MAA_gSAsia_win1000000_chr18 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr18.txt")
MAA_gSAsia_win2500000_chr18 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr18.txt")
MAA_gSAsia_win5000000_chr18 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr18.txt")
colnames(MAA_gSAsia_win25000_chr18) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr18) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr18) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr18) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr18) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr18) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr18) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr18) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr18) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr18) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

MAA_gSAsia_win25000_chr19 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr19.txt")
MAA_gSAsia_win50000_chr19 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr19.txt")
MAA_gSAsia_win75000_chr19 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr19.txt")
MAA_gSAsia_win100000_chr19 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr19.txt")
MAA_gSAsia_win250000_chr19 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr19.txt")
MAA_gSAsia_win500000_chr19 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr19.txt")
MAA_gSAsia_win750000_chr19 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr19.txt")
MAA_gSAsia_win1000000_chr19 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr19.txt")
MAA_gSAsia_win2500000_chr19 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr19.txt")
MAA_gSAsia_win5000000_chr19 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr19.txt")
colnames(MAA_gSAsia_win25000_chr19) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr19) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr19) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr19) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr19) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr19) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr19) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr19) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr19) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr19) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

MAA_gSAsia_win25000_chr20 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr20.txt")
MAA_gSAsia_win50000_chr20 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr20.txt")
MAA_gSAsia_win75000_chr20 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr20.txt")
MAA_gSAsia_win100000_chr20 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr20.txt")
MAA_gSAsia_win250000_chr20 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr20.txt")
MAA_gSAsia_win500000_chr20 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr20.txt")
MAA_gSAsia_win750000_chr20 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr20.txt")
MAA_gSAsia_win1000000_chr20 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr20.txt")
MAA_gSAsia_win2500000_chr20 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr20.txt")
MAA_gSAsia_win5000000_chr20 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr20.txt")
colnames(MAA_gSAsia_win25000_chr20) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr20) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr20) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr20) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr20) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr20) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr20) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr20) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr20) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr20) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

MAA_gSAsia_win25000_chr21 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr21.txt")
MAA_gSAsia_win50000_chr21 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr21.txt")
MAA_gSAsia_win75000_chr21 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr21.txt")
MAA_gSAsia_win100000_chr21 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr21.txt")
MAA_gSAsia_win250000_chr21 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr21.txt")
MAA_gSAsia_win500000_chr21 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr21.txt")
MAA_gSAsia_win750000_chr21 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr21.txt")
MAA_gSAsia_win1000000_chr21 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr21.txt")
MAA_gSAsia_win2500000_chr21 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr21.txt")
MAA_gSAsia_win5000000_chr21 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr21.txt")
colnames(MAA_gSAsia_win25000_chr21) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr21) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr21) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr21) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr21) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr21) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr21) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr21) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr21) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr21) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

MAA_gSAsia_win25000_chr22 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_25000/MAA_gSAsia_win25000_chr22.txt")
MAA_gSAsia_win50000_chr22 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_50000/MAA_gSAsia_win50000_chr22.txt")
MAA_gSAsia_win75000_chr22 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_75000/MAA_gSAsia_win75000_chr22.txt")
MAA_gSAsia_win100000_chr22 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_100000/MAA_gSAsia_win100000_chr22.txt")
MAA_gSAsia_win250000_chr22 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_250000/MAA_gSAsia_win250000_chr22.txt")
MAA_gSAsia_win500000_chr22 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_500000/MAA_gSAsia_win500000_chr22.txt")
MAA_gSAsia_win750000_chr22 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_750000/MAA_gSAsia_win750000_chr22.txt")
MAA_gSAsia_win1000000_chr22 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_1000000/MAA_gSAsia_win1000000_chr22.txt")
MAA_gSAsia_win2500000_chr22 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_2500000/MAA_gSAsia_win2500000_chr22.txt")
MAA_gSAsia_win5000000_chr22 <- read_csv("workfiles/recombinationAnalysis/human/MAA/combinedFiles/winLength_5000000/MAA_gSAsia_win5000000_chr22.txt")
colnames(MAA_gSAsia_win25000_chr22) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win50000_chr22) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win75000_chr22) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win100000_chr22) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win250000_chr22) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win500000_chr22) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win750000_chr22) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win1000000_chr22) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win2500000_chr22) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")
colnames(MAA_gSAsia_win5000000_chr22) <- c("position","Indv1","Indv2","Indv3","Indv4","Indv5","Indv6","Indv7","Indv8","Indv9")

######### Modify chromosomes 1 (Indv2: 25k, 50k, 75k, 100k), 12 (Indv2,5,7: 25k, 50k, 75k), and 13(Indv4,8: all) so that they are all the same size. 
MAA_gSAsia_win25000_chr1 <- MAA_gSAsia_win25000_chr1[1:length(MAA_gSAsia_win25000_chr1$position)-1,]
MAA_gSAsia_win50000_chr1 <- MAA_gSAsia_win50000_chr1[1:length(MAA_gSAsia_win50000_chr1$position)-1,]
MAA_gSAsia_win75000_chr1 <- MAA_gSAsia_win75000_chr1[1:length(MAA_gSAsia_win75000_chr1$position)-1,]
MAA_gSAsia_win100000_chr1 <- MAA_gSAsia_win100000_chr1[1:length(MAA_gSAsia_win100000_chr1$position)-1,]

MAA_gSAsia_win25000_chr12 <- MAA_gSAsia_win25000_chr12[1:length(MAA_gSAsia_win25000_chr1$position)-2,]
MAA_gSAsia_win50000_chr12 <- MAA_gSAsia_win50000_chr12[1:length(MAA_gSAsia_win50000_chr1$position)-2,]
MAA_gSAsia_win75000_chr12 <- MAA_gSAsia_win75000_chr12[1:length(MAA_gSAsia_win75000_chr1$position)-1,]

MAA_gSAsia_win25000_chr14 <- MAA_gSAsia_win25000_chr1[1:length(MAA_gSAsia_win25000_chr1$position)-150,]
MAA_gSAsia_win50000_chr14 <- MAA_gSAsia_win50000_chr1[1:length(MAA_gSAsia_win50000_chr1$position)-75,]
MAA_gSAsia_win75000_chr14 <- MAA_gSAsia_win75000_chr1[1:length(MAA_gSAsia_win75000_chr1$position)-49,]
MAA_gSAsia_win100000_chr14 <- MAA_gSAsia_win100000_chr1[1:length(MAA_gSAsia_win100000_chr1$position)-38,]
MAA_gSAsia_win250000_chr14 <- MAA_gSAsia_win250000_chr1[1:length(MAA_gSAsia_win250000_chr1$position)-15,]
MAA_gSAsia_win500000_chr14 <- MAA_gSAsia_win500000_chr1[1:length(MAA_gSAsia_win500000_chr1$position)-7,]
MAA_gSAsia_win750000_chr14 <- MAA_gSAsia_win750000_chr1[1:length(MAA_gSAsia_win750000_chr1$position)-5,]
MAA_gSAsia_win1000000_chr14 <- MAA_gSAsia_win1000000_chr1[1:length(MAA_gSAsia_win1000000_chr1$position)-3,]
MAA_gSAsia_win2500000_chr14 <- MAA_gSAsia_win2500000_chr1[1:length(MAA_gSAsia_win2500000_chr1$position)-2,]
MAA_gSAsia_win5000000_chr14 <- MAA_gSAsia_win5000000_chr1[1:length(MAA_gSAsia_win5000000_chr1$position)-1,]

######### Concatenate all chromosomes for each of the window sizes. 

MAA_gSAsia_win25000 <- rbind(MAA_gSAsia_win25000_chr1,MAA_gSAsia_win25000_chr2,MAA_gSAsia_win25000_chr3,MAA_gSAsia_win25000_chr4,MAA_gSAsia_win25000_chr5,MAA_gSAsia_win25000_chr6,
                             MAA_gSAsia_win25000_chr7,MAA_gSAsia_win25000_chr8,MAA_gSAsia_win25000_chr9,MAA_gSAsia_win25000_chr10,MAA_gSAsia_win25000_chr11,MAA_gSAsia_win25000_chr12,
                             MAA_gSAsia_win25000_chr13,MAA_gSAsia_win25000_chr14,MAA_gSAsia_win25000_chr15,MAA_gSAsia_win25000_chr16,MAA_gSAsia_win25000_chr17,MAA_gSAsia_win25000_chr18,
                             MAA_gSAsia_win25000_chr19,MAA_gSAsia_win25000_chr20,MAA_gSAsia_win25000_chr21,MAA_gSAsia_win25000_chr22)

MAA_gSAsia_win50000 <- rbind(MAA_gSAsia_win50000_chr1,MAA_gSAsia_win50000_chr2,MAA_gSAsia_win50000_chr3,MAA_gSAsia_win50000_chr4,MAA_gSAsia_win50000_chr5,MAA_gSAsia_win50000_chr6,
                             MAA_gSAsia_win50000_chr7,MAA_gSAsia_win50000_chr8,MAA_gSAsia_win50000_chr9,MAA_gSAsia_win50000_chr10,MAA_gSAsia_win50000_chr11,MAA_gSAsia_win50000_chr12,
                             MAA_gSAsia_win50000_chr13,MAA_gSAsia_win50000_chr14,MAA_gSAsia_win50000_chr15,MAA_gSAsia_win50000_chr16,MAA_gSAsia_win50000_chr17,MAA_gSAsia_win50000_chr18,
                             MAA_gSAsia_win50000_chr19,MAA_gSAsia_win50000_chr20,MAA_gSAsia_win50000_chr21,MAA_gSAsia_win50000_chr22)

MAA_gSAsia_win75000 <- rbind(MAA_gSAsia_win75000_chr1,MAA_gSAsia_win75000_chr2,MAA_gSAsia_win75000_chr3,MAA_gSAsia_win75000_chr4,MAA_gSAsia_win75000_chr5,MAA_gSAsia_win75000_chr6,
                             MAA_gSAsia_win75000_chr7,MAA_gSAsia_win75000_chr8,MAA_gSAsia_win75000_chr9,MAA_gSAsia_win75000_chr10,MAA_gSAsia_win75000_chr11,MAA_gSAsia_win75000_chr12,
                             MAA_gSAsia_win75000_chr13,MAA_gSAsia_win75000_chr14,MAA_gSAsia_win75000_chr15,MAA_gSAsia_win75000_chr16,MAA_gSAsia_win75000_chr17,MAA_gSAsia_win75000_chr18,
                             MAA_gSAsia_win75000_chr19,MAA_gSAsia_win75000_chr20,MAA_gSAsia_win75000_chr21,MAA_gSAsia_win75000_chr22)

MAA_gSAsia_win100000 <- rbind(MAA_gSAsia_win100000_chr1,MAA_gSAsia_win100000_chr2,MAA_gSAsia_win100000_chr3,MAA_gSAsia_win100000_chr4,MAA_gSAsia_win100000_chr5,MAA_gSAsia_win100000_chr6,
                             MAA_gSAsia_win100000_chr7,MAA_gSAsia_win100000_chr8,MAA_gSAsia_win100000_chr9,MAA_gSAsia_win100000_chr10,MAA_gSAsia_win100000_chr11,MAA_gSAsia_win100000_chr12,
                             MAA_gSAsia_win100000_chr13,MAA_gSAsia_win100000_chr14,MAA_gSAsia_win100000_chr15,MAA_gSAsia_win100000_chr16,MAA_gSAsia_win100000_chr17,MAA_gSAsia_win100000_chr18,
                             MAA_gSAsia_win100000_chr19,MAA_gSAsia_win100000_chr20,MAA_gSAsia_win100000_chr21,MAA_gSAsia_win100000_chr22)

MAA_gSAsia_win250000 <- rbind(MAA_gSAsia_win250000_chr1,MAA_gSAsia_win250000_chr2,MAA_gSAsia_win250000_chr3,MAA_gSAsia_win250000_chr4,MAA_gSAsia_win250000_chr5,MAA_gSAsia_win250000_chr6,
                             MAA_gSAsia_win250000_chr7,MAA_gSAsia_win250000_chr8,MAA_gSAsia_win250000_chr9,MAA_gSAsia_win250000_chr10,MAA_gSAsia_win250000_chr11,MAA_gSAsia_win250000_chr12,
                             MAA_gSAsia_win250000_chr13,MAA_gSAsia_win250000_chr14,MAA_gSAsia_win250000_chr15,MAA_gSAsia_win250000_chr16,MAA_gSAsia_win250000_chr17,MAA_gSAsia_win250000_chr18,
                             MAA_gSAsia_win250000_chr19,MAA_gSAsia_win250000_chr20,MAA_gSAsia_win250000_chr21,MAA_gSAsia_win250000_chr22)

MAA_gSAsia_win500000 <- rbind(MAA_gSAsia_win500000_chr1,MAA_gSAsia_win500000_chr2,MAA_gSAsia_win500000_chr3,MAA_gSAsia_win500000_chr4,MAA_gSAsia_win500000_chr5,MAA_gSAsia_win500000_chr6,
                             MAA_gSAsia_win500000_chr7,MAA_gSAsia_win500000_chr8,MAA_gSAsia_win500000_chr9,MAA_gSAsia_win500000_chr10,MAA_gSAsia_win500000_chr11,MAA_gSAsia_win500000_chr12,
                             MAA_gSAsia_win500000_chr13,MAA_gSAsia_win500000_chr14,MAA_gSAsia_win500000_chr15,MAA_gSAsia_win500000_chr16,MAA_gSAsia_win500000_chr17,MAA_gSAsia_win500000_chr18,
                             MAA_gSAsia_win500000_chr19,MAA_gSAsia_win500000_chr20,MAA_gSAsia_win500000_chr21,MAA_gSAsia_win500000_chr22)

MAA_gSAsia_win750000 <- rbind(MAA_gSAsia_win750000_chr1,MAA_gSAsia_win750000_chr2,MAA_gSAsia_win750000_chr3,MAA_gSAsia_win750000_chr4,MAA_gSAsia_win750000_chr5,MAA_gSAsia_win750000_chr6,
                             MAA_gSAsia_win750000_chr7,MAA_gSAsia_win750000_chr8,MAA_gSAsia_win750000_chr9,MAA_gSAsia_win750000_chr10,MAA_gSAsia_win750000_chr11,MAA_gSAsia_win750000_chr12,
                             MAA_gSAsia_win750000_chr13,MAA_gSAsia_win750000_chr14,MAA_gSAsia_win750000_chr15,MAA_gSAsia_win750000_chr16,MAA_gSAsia_win750000_chr17,MAA_gSAsia_win750000_chr18,
                             MAA_gSAsia_win750000_chr19,MAA_gSAsia_win750000_chr20,MAA_gSAsia_win750000_chr21,MAA_gSAsia_win750000_chr22)

MAA_gSAsia_win1000000 <- rbind(MAA_gSAsia_win1000000_chr1,MAA_gSAsia_win1000000_chr2,MAA_gSAsia_win1000000_chr3,MAA_gSAsia_win1000000_chr4,MAA_gSAsia_win1000000_chr5,MAA_gSAsia_win1000000_chr6,
                             MAA_gSAsia_win1000000_chr7,MAA_gSAsia_win1000000_chr8,MAA_gSAsia_win1000000_chr9,MAA_gSAsia_win1000000_chr10,MAA_gSAsia_win1000000_chr11,MAA_gSAsia_win1000000_chr12,
                             MAA_gSAsia_win1000000_chr13,MAA_gSAsia_win1000000_chr14,MAA_gSAsia_win1000000_chr15,MAA_gSAsia_win1000000_chr16,MAA_gSAsia_win1000000_chr17,MAA_gSAsia_win1000000_chr18,
                             MAA_gSAsia_win1000000_chr19,MAA_gSAsia_win1000000_chr20,MAA_gSAsia_win1000000_chr21,MAA_gSAsia_win1000000_chr22)

MAA_gSAsia_win2500000 <- rbind(MAA_gSAsia_win2500000_chr1,MAA_gSAsia_win2500000_chr2,MAA_gSAsia_win2500000_chr3,MAA_gSAsia_win2500000_chr4,MAA_gSAsia_win2500000_chr5,MAA_gSAsia_win2500000_chr6,
                             MAA_gSAsia_win2500000_chr7,MAA_gSAsia_win2500000_chr8,MAA_gSAsia_win2500000_chr9,MAA_gSAsia_win2500000_chr10,MAA_gSAsia_win2500000_chr11,MAA_gSAsia_win2500000_chr12,
                             MAA_gSAsia_win2500000_chr13,MAA_gSAsia_win2500000_chr14,MAA_gSAsia_win2500000_chr15,MAA_gSAsia_win2500000_chr16,MAA_gSAsia_win2500000_chr17,MAA_gSAsia_win2500000_chr18,
                             MAA_gSAsia_win2500000_chr19,MAA_gSAsia_win2500000_chr20,MAA_gSAsia_win2500000_chr21,MAA_gSAsia_win2500000_chr22)

MAA_gSAsia_win5000000 <- rbind(MAA_gSAsia_win5000000_chr1,MAA_gSAsia_win5000000_chr2,MAA_gSAsia_win5000000_chr3,MAA_gSAsia_win5000000_chr4,MAA_gSAsia_win5000000_chr5,MAA_gSAsia_win5000000_chr6,
                             MAA_gSAsia_win5000000_chr7,MAA_gSAsia_win5000000_chr8,MAA_gSAsia_win5000000_chr9,MAA_gSAsia_win5000000_chr10,MAA_gSAsia_win5000000_chr11,MAA_gSAsia_win5000000_chr12,
                             MAA_gSAsia_win5000000_chr13,MAA_gSAsia_win5000000_chr14,MAA_gSAsia_win5000000_chr15,MAA_gSAsia_win5000000_chr16,MAA_gSAsia_win5000000_chr17,MAA_gSAsia_win5000000_chr18,
                             MAA_gSAsia_win5000000_chr19,MAA_gSAsia_win5000000_chr20,MAA_gSAsia_win5000000_chr21,MAA_gSAsia_win5000000_chr22)




MAA_gSAsia_win25000_allChr <- c(MAA_gSAsia_win25000_chr1,)














































#corVec_win25000_chr1 <- c(cor(MAA_gSAsia_win25000_chr1$Indv1, MAA_gSAsia_win25000_chr1$Indv2),cor(MAA_gSAsia_win25000_chr1$Indv1, MAA_gSAsia_win25000_chr1$Indv3),cor(MAA_gSAsia_win25000_chr1$Indv1, MAA_gSAsia_win25000_chr1$Indv4), 
#                          cor(MAA_gSAsia_win25000_chr1$Indv1, MAA_gSAsia_win25000_chr1$Indv5),cor(MAA_gSAsia_win25000_chr1$Indv1, MAA_gSAsia_win25000_chr1$Indv6),cor(MAA_gSAsia_win25000_chr1$Indv1, MAA_gSAsia_win25000_chr1$Indv7),
#                          cor(MAA_gSAsia_win25000_chr1$Indv1, MAA_gSAsia_win25000_chr1$Indv8),cor(MAA_gSAsia_win25000_chr1$Indv1, MAA_gSAsia_win25000_chr1$Indv9),
#                          cor(MAA_gSAsia_win25000_chr1$Indv2, MAA_gSAsia_win25000_chr1$Indv3),cor(MAA_gSAsia_win25000_chr1$Indv2, MAA_gSAsia_win25000_chr1$Indv4),cor(MAA_gSAsia_win25000_chr1$Indv2, MAA_gSAsia_win25000_chr1$Indv5),
#                          cor(MAA_gSAsia_win25000_chr1$Indv2, MAA_gSAsia_win25000_chr1$Indv6),cor(MAA_gSAsia_win25000_chr1$Indv2, MAA_gSAsia_win25000_chr1$Indv7),cor(MAA_gSAsia_win25000_chr1$Indv2, MAA_gSAsia_win25000_chr1$Indv8),
#                          cor(MAA_gSAsia_win25000_chr1$Indv2, MAA_gSAsia_win25000_chr1$Indv9),
#                          cor(MAA_gSAsia_win25000_chr1$Indv3, MAA_gSAsia_win25000_chr1$Indv4),cor(MAA_gSAsia_win25000_chr1$Indv3, MAA_gSAsia_win25000_chr1$Indv5),cor(MAA_gSAsia_win25000_chr1$Indv3, MAA_gSAsia_win25000_chr1$Indv6),
#                          cor(MAA_gSAsia_win25000_chr1$Indv3, MAA_gSAsia_win25000_chr1$Indv7),cor(MAA_gSAsia_win25000_chr1$Indv3, MAA_gSAsia_win25000_chr1$Indv8),cor(MAA_gSAsia_win25000_chr1$Indv3, MAA_gSAsia_win25000_chr1$Indv9),
#                          cor(MAA_gSAsia_win25000_chr1$Indv4, MAA_gSAsia_win25000_chr1$Indv5),cor(MAA_gSAsia_win25000_chr1$Indv4, MAA_gSAsia_win25000_chr1$Indv6),cor(MAA_gSAsia_win25000_chr1$Indv4, MAA_gSAsia_win25000_chr1$Indv7),
#                          cor(MAA_gSAsia_win25000_chr1$Indv4, MAA_gSAsia_win25000_chr1$Indv8),cor(MAA_gSAsia_win25000_chr1$Indv4, MAA_gSAsia_win25000_chr1$Indv9),
#                          cor(MAA_gSAsia_win25000_chr1$Indv5, MAA_gSAsia_win25000_chr1$Indv6),cor(MAA_gSAsia_win25000_chr1$Indv5, MAA_gSAsia_win25000_chr1$Indv7),cor(MAA_gSAsia_win25000_chr1$Indv5, MAA_gSAsia_win25000_chr1$Indv8),
#                          cor(MAA_gSAsia_win25000_chr1$Indv5, MAA_gSAsia_win25000_chr1$Indv9),
#                          cor(MAA_gSAsia_win25000_chr1$Indv6, MAA_gSAsia_win25000_chr1$Indv7),cor(MAA_gSAsia_win25000_chr1$Indv6, MAA_gSAsia_win25000_chr1$Indv8),cor(MAA_gSAsia_win25000_chr1$Indv6, MAA_gSAsia_win25000_chr1$Indv9),
#                          cor(MAA_gSAsia_win25000_chr1$Indv7, MAA_gSAsia_win25000_chr1$Indv8),cor(MAA_gSAsia_win25000_chr1$Indv7, MAA_gSAsia_win25000_chr1$Indv9),
#                          cor(MAA_gSAsia_win25000_chr1$Indv8, MAA_gSAsia_win25000_chr1$Indv9))
#
#corVec_win25000_chr2 <- c(cor(MAA_gSAsia_win25000_chr2$Indv1, MAA_gSAsia_win25000_chr2$Indv2),cor(MAA_gSAsia_win25000_chr2$Indv1, MAA_gSAsia_win25000_chr2$Indv3),cor(MAA_gSAsia_win25000_chr2$Indv1, MAA_gSAsia_win25000_chr2$Indv4), 
#                          cor(MAA_gSAsia_win25000_chr2$Indv1, MAA_gSAsia_win25000_chr2$Indv5),cor(MAA_gSAsia_win25000_chr2$Indv1, MAA_gSAsia_win25000_chr2$Indv6),cor(MAA_gSAsia_win25000_chr2$Indv1, MAA_gSAsia_win25000_chr2$Indv7),
#                          cor(MAA_gSAsia_win25000_chr2$Indv1, MAA_gSAsia_win25000_chr2$Indv8),cor(MAA_gSAsia_win25000_chr2$Indv1, MAA_gSAsia_win25000_chr2$Indv9),
#                          cor(MAA_gSAsia_win25000_chr2$Indv2, MAA_gSAsia_win25000_chr2$Indv3),cor(MAA_gSAsia_win25000_chr2$Indv2, MAA_gSAsia_win25000_chr2$Indv4),cor(MAA_gSAsia_win25000_chr2$Indv2, MAA_gSAsia_win25000_chr2$Indv5),
#                          cor(MAA_gSAsia_win25000_chr2$Indv2, MAA_gSAsia_win25000_chr2$Indv6),cor(MAA_gSAsia_win25000_chr2$Indv2, MAA_gSAsia_win25000_chr2$Indv7),cor(MAA_gSAsia_win25000_chr2$Indv2, MAA_gSAsia_win25000_chr2$Indv8),
#                          cor(MAA_gSAsia_win25000_chr2$Indv2, MAA_gSAsia_win25000_chr2$Indv9),
#                          cor(MAA_gSAsia_win25000_chr2$Indv3, MAA_gSAsia_win25000_chr2$Indv4),cor(MAA_gSAsia_win25000_chr2$Indv3, MAA_gSAsia_win25000_chr2$Indv5),cor(MAA_gSAsia_win25000_chr2$Indv3, MAA_gSAsia_win25000_chr2$Indv6),
#                          cor(MAA_gSAsia_win25000_chr2$Indv3, MAA_gSAsia_win25000_chr2$Indv7),cor(MAA_gSAsia_win25000_chr2$Indv3, MAA_gSAsia_win25000_chr2$Indv8),cor(MAA_gSAsia_win25000_chr2$Indv3, MAA_gSAsia_win25000_chr2$Indv9),
#                          cor(MAA_gSAsia_win25000_chr2$Indv4, MAA_gSAsia_win25000_chr2$Indv5),cor(MAA_gSAsia_win25000_chr2$Indv4, MAA_gSAsia_win25000_chr2$Indv6),cor(MAA_gSAsia_win25000_chr2$Indv4, MAA_gSAsia_win25000_chr2$Indv7),
#                          cor(MAA_gSAsia_win25000_chr2$Indv4, MAA_gSAsia_win25000_chr2$Indv8),cor(MAA_gSAsia_win25000_chr2$Indv4, MAA_gSAsia_win25000_chr2$Indv9),
#                          cor(MAA_gSAsia_win25000_chr2$Indv5, MAA_gSAsia_win25000_chr2$Indv6),cor(MAA_gSAsia_win25000_chr2$Indv5, MAA_gSAsia_win25000_chr2$Indv7),cor(MAA_gSAsia_win25000_chr2$Indv5, MAA_gSAsia_win25000_chr2$Indv8),
#                          cor(MAA_gSAsia_win25000_chr2$Indv5, MAA_gSAsia_win25000_chr2$Indv9),
#                          cor(MAA_gSAsia_win25000_chr2$Indv6, MAA_gSAsia_win25000_chr2$Indv7),cor(MAA_gSAsia_win25000_chr2$Indv6, MAA_gSAsia_win25000_chr2$Indv8),cor(MAA_gSAsia_win25000_chr2$Indv6, MAA_gSAsia_win25000_chr2$Indv9),
#                          cor(MAA_gSAsia_win25000_chr2$Indv7, MAA_gSAsia_win25000_chr2$Indv8),cor(MAA_gSAsia_win25000_chr2$Indv7, MAA_gSAsia_win25000_chr2$Indv9),
#                          cor(MAA_gSAsia_win25000_chr2$Indv8, MAA_gSAsia_win25000_chr2$Indv9))
#
#corVec_win25000_chr3 <- c(cor(MAA_gSAsia_win25000_chr3$Indv1, MAA_gSAsia_win25000_chr3$Indv2),cor(MAA_gSAsia_win25000_chr3$Indv1, MAA_gSAsia_win25000_chr3$Indv3),cor(MAA_gSAsia_win25000_chr3$Indv1, MAA_gSAsia_win25000_chr3$Indv4), 
#                          cor(MAA_gSAsia_win25000_chr3$Indv1, MAA_gSAsia_win25000_chr3$Indv5),cor(MAA_gSAsia_win25000_chr3$Indv1, MAA_gSAsia_win25000_chr3$Indv6),cor(MAA_gSAsia_win25000_chr3$Indv1, MAA_gSAsia_win25000_chr3$Indv7),
#                          cor(MAA_gSAsia_win25000_chr3$Indv1, MAA_gSAsia_win25000_chr3$Indv8),cor(MAA_gSAsia_win25000_chr3$Indv1, MAA_gSAsia_win25000_chr3$Indv9),
#                          cor(MAA_gSAsia_win25000_chr3$Indv2, MAA_gSAsia_win25000_chr3$Indv3),cor(MAA_gSAsia_win25000_chr3$Indv2, MAA_gSAsia_win25000_chr3$Indv4),cor(MAA_gSAsia_win25000_chr3$Indv2, MAA_gSAsia_win25000_chr3$Indv5),
#                          cor(MAA_gSAsia_win25000_chr3$Indv2, MAA_gSAsia_win25000_chr3$Indv6),cor(MAA_gSAsia_win25000_chr3$Indv2, MAA_gSAsia_win25000_chr3$Indv7),cor(MAA_gSAsia_win25000_chr3$Indv2, MAA_gSAsia_win25000_chr3$Indv8),
#                          cor(MAA_gSAsia_win25000_chr3$Indv2, MAA_gSAsia_win25000_chr3$Indv9),
#                          cor(MAA_gSAsia_win25000_chr3$Indv3, MAA_gSAsia_win25000_chr3$Indv4),cor(MAA_gSAsia_win25000_chr3$Indv3, MAA_gSAsia_win25000_chr3$Indv5),cor(MAA_gSAsia_win25000_chr3$Indv3, MAA_gSAsia_win25000_chr3$Indv6),
#                          cor(MAA_gSAsia_win25000_chr3$Indv3, MAA_gSAsia_win25000_chr3$Indv7),cor(MAA_gSAsia_win25000_chr3$Indv3, MAA_gSAsia_win25000_chr3$Indv8),cor(MAA_gSAsia_win25000_chr3$Indv3, MAA_gSAsia_win25000_chr3$Indv9),
#                          cor(MAA_gSAsia_win25000_chr3$Indv4, MAA_gSAsia_win25000_chr3$Indv5),cor(MAA_gSAsia_win25000_chr3$Indv4, MAA_gSAsia_win25000_chr3$Indv6),cor(MAA_gSAsia_win25000_chr3$Indv4, MAA_gSAsia_win25000_chr3$Indv7),
#                          cor(MAA_gSAsia_win25000_chr3$Indv4, MAA_gSAsia_win25000_chr3$Indv8),cor(MAA_gSAsia_win25000_chr3$Indv4, MAA_gSAsia_win25000_chr3$Indv9),
#                          cor(MAA_gSAsia_win25000_chr3$Indv5, MAA_gSAsia_win25000_chr3$Indv6),cor(MAA_gSAsia_win25000_chr3$Indv5, MAA_gSAsia_win25000_chr3$Indv7),cor(MAA_gSAsia_win25000_chr3$Indv5, MAA_gSAsia_win25000_chr3$Indv8),
#                          cor(MAA_gSAsia_win25000_chr3$Indv5, MAA_gSAsia_win25000_chr3$Indv9),
#                          cor(MAA_gSAsia_win25000_chr3$Indv6, MAA_gSAsia_win25000_chr3$Indv7),cor(MAA_gSAsia_win25000_chr3$Indv6, MAA_gSAsia_win25000_chr3$Indv8),cor(MAA_gSAsia_win25000_chr3$Indv6, MAA_gSAsia_win25000_chr3$Indv9),
#                          cor(MAA_gSAsia_win25000_chr3$Indv7, MAA_gSAsia_win25000_chr3$Indv8),cor(MAA_gSAsia_win25000_chr3$Indv7, MAA_gSAsia_win25000_chr3$Indv9),
#                          cor(MAA_gSAsia_win25000_chr3$Indv8, MAA_gSAsia_win25000_chr3$Indv9))
#
#corVec_win25000_chr4 <- c(cor(MAA_gSAsia_win25000_chr4$Indv1, MAA_gSAsia_win25000_chr4$Indv2),cor(MAA_gSAsia_win25000_chr4$Indv1, MAA_gSAsia_win25000_chr4$Indv3),cor(MAA_gSAsia_win25000_chr4$Indv1, MAA_gSAsia_win25000_chr4$Indv4), 
#                          cor(MAA_gSAsia_win25000_chr4$Indv1, MAA_gSAsia_win25000_chr4$Indv5),cor(MAA_gSAsia_win25000_chr4$Indv1, MAA_gSAsia_win25000_chr4$Indv6),cor(MAA_gSAsia_win25000_chr4$Indv1, MAA_gSAsia_win25000_chr4$Indv7),
#                          cor(MAA_gSAsia_win25000_chr4$Indv1, MAA_gSAsia_win25000_chr4$Indv8),cor(MAA_gSAsia_win25000_chr4$Indv1, MAA_gSAsia_win25000_chr4$Indv9),
#                          cor(MAA_gSAsia_win25000_chr4$Indv2, MAA_gSAsia_win25000_chr4$Indv3),cor(MAA_gSAsia_win25000_chr4$Indv2, MAA_gSAsia_win25000_chr4$Indv4),cor(MAA_gSAsia_win25000_chr4$Indv2, MAA_gSAsia_win25000_chr4$Indv5),
#                          cor(MAA_gSAsia_win25000_chr4$Indv2, MAA_gSAsia_win25000_chr4$Indv6),cor(MAA_gSAsia_win25000_chr4$Indv2, MAA_gSAsia_win25000_chr4$Indv7),cor(MAA_gSAsia_win25000_chr4$Indv2, MAA_gSAsia_win25000_chr4$Indv8),
#                          cor(MAA_gSAsia_win25000_chr4$Indv2, MAA_gSAsia_win25000_chr4$Indv9),
#                          cor(MAA_gSAsia_win25000_chr4$Indv3, MAA_gSAsia_win25000_chr4$Indv4),cor(MAA_gSAsia_win25000_chr4$Indv3, MAA_gSAsia_win25000_chr4$Indv5),cor(MAA_gSAsia_win25000_chr4$Indv3, MAA_gSAsia_win25000_chr4$Indv6),
#                          cor(MAA_gSAsia_win25000_chr4$Indv3, MAA_gSAsia_win25000_chr4$Indv7),cor(MAA_gSAsia_win25000_chr4$Indv3, MAA_gSAsia_win25000_chr4$Indv8),cor(MAA_gSAsia_win25000_chr4$Indv3, MAA_gSAsia_win25000_chr4$Indv9),
#                          cor(MAA_gSAsia_win25000_chr4$Indv4, MAA_gSAsia_win25000_chr4$Indv5),cor(MAA_gSAsia_win25000_chr4$Indv4, MAA_gSAsia_win25000_chr4$Indv6),cor(MAA_gSAsia_win25000_chr4$Indv4, MAA_gSAsia_win25000_chr4$Indv7),
#                          cor(MAA_gSAsia_win25000_chr4$Indv4, MAA_gSAsia_win25000_chr4$Indv8),cor(MAA_gSAsia_win25000_chr4$Indv4, MAA_gSAsia_win25000_chr4$Indv9),
#                          cor(MAA_gSAsia_win25000_chr4$Indv5, MAA_gSAsia_win25000_chr4$Indv6),cor(MAA_gSAsia_win25000_chr4$Indv5, MAA_gSAsia_win25000_chr4$Indv7),cor(MAA_gSAsia_win25000_chr4$Indv5, MAA_gSAsia_win25000_chr4$Indv8),
#                          cor(MAA_gSAsia_win25000_chr4$Indv5, MAA_gSAsia_win25000_chr4$Indv9),
#                          cor(MAA_gSAsia_win25000_chr4$Indv6, MAA_gSAsia_win25000_chr4$Indv7),cor(MAA_gSAsia_win25000_chr4$Indv6, MAA_gSAsia_win25000_chr4$Indv8),cor(MAA_gSAsia_win25000_chr4$Indv6, MAA_gSAsia_win25000_chr4$Indv9),
#                          cor(MAA_gSAsia_win25000_chr4$Indv7, MAA_gSAsia_win25000_chr4$Indv8),cor(MAA_gSAsia_win25000_chr4$Indv7, MAA_gSAsia_win25000_chr4$Indv9),
#                          cor(MAA_gSAsia_win25000_chr4$Indv8, MAA_gSAsia_win25000_chr4$Indv9))
#
#corVec_win25000_chr5 <- c(cor(MAA_gSAsia_win25000_chr5$Indv1, MAA_gSAsia_win25000_chr5$Indv2),cor(MAA_gSAsia_win25000_chr5$Indv1, MAA_gSAsia_win25000_chr5$Indv3),cor(MAA_gSAsia_win25000_chr5$Indv1, MAA_gSAsia_win25000_chr5$Indv4), 
#                          cor(MAA_gSAsia_win25000_chr5$Indv1, MAA_gSAsia_win25000_chr5$Indv5),cor(MAA_gSAsia_win25000_chr5$Indv1, MAA_gSAsia_win25000_chr5$Indv6),cor(MAA_gSAsia_win25000_chr5$Indv1, MAA_gSAsia_win25000_chr5$Indv7),
#                          cor(MAA_gSAsia_win25000_chr5$Indv1, MAA_gSAsia_win25000_chr5$Indv8),cor(MAA_gSAsia_win25000_chr5$Indv1, MAA_gSAsia_win25000_chr5$Indv9),
#                          cor(MAA_gSAsia_win25000_chr5$Indv2, MAA_gSAsia_win25000_chr5$Indv3),cor(MAA_gSAsia_win25000_chr5$Indv2, MAA_gSAsia_win25000_chr5$Indv4),cor(MAA_gSAsia_win25000_chr5$Indv2, MAA_gSAsia_win25000_chr5$Indv5),
#                          cor(MAA_gSAsia_win25000_chr5$Indv2, MAA_gSAsia_win25000_chr5$Indv6),cor(MAA_gSAsia_win25000_chr5$Indv2, MAA_gSAsia_win25000_chr5$Indv7),cor(MAA_gSAsia_win25000_chr5$Indv2, MAA_gSAsia_win25000_chr5$Indv8),
#                          cor(MAA_gSAsia_win25000_chr5$Indv2, MAA_gSAsia_win25000_chr5$Indv9),
#                          cor(MAA_gSAsia_win25000_chr5$Indv3, MAA_gSAsia_win25000_chr5$Indv4),cor(MAA_gSAsia_win25000_chr5$Indv3, MAA_gSAsia_win25000_chr5$Indv5),cor(MAA_gSAsia_win25000_chr5$Indv3, MAA_gSAsia_win25000_chr5$Indv6),
#                          cor(MAA_gSAsia_win25000_chr5$Indv3, MAA_gSAsia_win25000_chr5$Indv7),cor(MAA_gSAsia_win25000_chr5$Indv3, MAA_gSAsia_win25000_chr5$Indv8),cor(MAA_gSAsia_win25000_chr5$Indv3, MAA_gSAsia_win25000_chr5$Indv9),
#                          cor(MAA_gSAsia_win25000_chr5$Indv4, MAA_gSAsia_win25000_chr5$Indv5),cor(MAA_gSAsia_win25000_chr5$Indv4, MAA_gSAsia_win25000_chr5$Indv6),cor(MAA_gSAsia_win25000_chr5$Indv4, MAA_gSAsia_win25000_chr5$Indv7),
#                          cor(MAA_gSAsia_win25000_chr5$Indv4, MAA_gSAsia_win25000_chr5$Indv8),cor(MAA_gSAsia_win25000_chr5$Indv4, MAA_gSAsia_win25000_chr5$Indv9),
#                          cor(MAA_gSAsia_win25000_chr5$Indv5, MAA_gSAsia_win25000_chr5$Indv6),cor(MAA_gSAsia_win25000_chr5$Indv5, MAA_gSAsia_win25000_chr5$Indv7),cor(MAA_gSAsia_win25000_chr5$Indv5, MAA_gSAsia_win25000_chr5$Indv8),
#                          cor(MAA_gSAsia_win25000_chr5$Indv5, MAA_gSAsia_win25000_chr5$Indv9),
#                          cor(MAA_gSAsia_win25000_chr5$Indv6, MAA_gSAsia_win25000_chr5$Indv7),cor(MAA_gSAsia_win25000_chr5$Indv6, MAA_gSAsia_win25000_chr5$Indv8),cor(MAA_gSAsia_win25000_chr5$Indv6, MAA_gSAsia_win25000_chr5$Indv9),
#                          cor(MAA_gSAsia_win25000_chr5$Indv7, MAA_gSAsia_win25000_chr5$Indv8),cor(MAA_gSAsia_win25000_chr5$Indv7, MAA_gSAsia_win25000_chr5$Indv9),
#                          cor(MAA_gSAsia_win25000_chr5$Indv8, MAA_gSAsia_win25000_chr5$Indv9))
#
#corVec_win25000_chr6 <- c(cor(MAA_gSAsia_win25000_chr6$Indv1, MAA_gSAsia_win25000_chr6$Indv2),cor(MAA_gSAsia_win25000_chr6$Indv1, MAA_gSAsia_win25000_chr6$Indv3),cor(MAA_gSAsia_win25000_chr6$Indv1, MAA_gSAsia_win25000_chr6$Indv4), 
#                          cor(MAA_gSAsia_win25000_chr6$Indv1, MAA_gSAsia_win25000_chr6$Indv5),cor(MAA_gSAsia_win25000_chr6$Indv1, MAA_gSAsia_win25000_chr6$Indv6),cor(MAA_gSAsia_win25000_chr6$Indv1, MAA_gSAsia_win25000_chr6$Indv7),
#                          cor(MAA_gSAsia_win25000_chr6$Indv1, MAA_gSAsia_win25000_chr6$Indv8),cor(MAA_gSAsia_win25000_chr6$Indv1, MAA_gSAsia_win25000_chr6$Indv9),
#                          cor(MAA_gSAsia_win25000_chr6$Indv2, MAA_gSAsia_win25000_chr6$Indv3),cor(MAA_gSAsia_win25000_chr6$Indv2, MAA_gSAsia_win25000_chr6$Indv4),cor(MAA_gSAsia_win25000_chr6$Indv2, MAA_gSAsia_win25000_chr6$Indv5),
#                          cor(MAA_gSAsia_win25000_chr6$Indv2, MAA_gSAsia_win25000_chr6$Indv6),cor(MAA_gSAsia_win25000_chr6$Indv2, MAA_gSAsia_win25000_chr6$Indv7),cor(MAA_gSAsia_win25000_chr6$Indv2, MAA_gSAsia_win25000_chr6$Indv8),
#                          cor(MAA_gSAsia_win25000_chr6$Indv2, MAA_gSAsia_win25000_chr6$Indv9),
#                          cor(MAA_gSAsia_win25000_chr6$Indv3, MAA_gSAsia_win25000_chr6$Indv4),cor(MAA_gSAsia_win25000_chr6$Indv3, MAA_gSAsia_win25000_chr6$Indv5),cor(MAA_gSAsia_win25000_chr6$Indv3, MAA_gSAsia_win25000_chr6$Indv6),
#                          cor(MAA_gSAsia_win25000_chr6$Indv3, MAA_gSAsia_win25000_chr6$Indv7),cor(MAA_gSAsia_win25000_chr6$Indv3, MAA_gSAsia_win25000_chr6$Indv8),cor(MAA_gSAsia_win25000_chr6$Indv3, MAA_gSAsia_win25000_chr6$Indv9),
#                          cor(MAA_gSAsia_win25000_chr6$Indv4, MAA_gSAsia_win25000_chr6$Indv5),cor(MAA_gSAsia_win25000_chr6$Indv4, MAA_gSAsia_win25000_chr6$Indv6),cor(MAA_gSAsia_win25000_chr6$Indv4, MAA_gSAsia_win25000_chr6$Indv7),
#                          cor(MAA_gSAsia_win25000_chr6$Indv4, MAA_gSAsia_win25000_chr6$Indv8),cor(MAA_gSAsia_win25000_chr6$Indv4, MAA_gSAsia_win25000_chr6$Indv9),
#                          cor(MAA_gSAsia_win25000_chr6$Indv5, MAA_gSAsia_win25000_chr6$Indv6),cor(MAA_gSAsia_win25000_chr6$Indv5, MAA_gSAsia_win25000_chr6$Indv7),cor(MAA_gSAsia_win25000_chr6$Indv5, MAA_gSAsia_win25000_chr6$Indv8),
#                          cor(MAA_gSAsia_win25000_chr6$Indv5, MAA_gSAsia_win25000_chr6$Indv9),
#                          cor(MAA_gSAsia_win25000_chr6$Indv6, MAA_gSAsia_win25000_chr6$Indv7),cor(MAA_gSAsia_win25000_chr6$Indv6, MAA_gSAsia_win25000_chr6$Indv8),cor(MAA_gSAsia_win25000_chr6$Indv6, MAA_gSAsia_win25000_chr6$Indv9),
#                          cor(MAA_gSAsia_win25000_chr6$Indv7, MAA_gSAsia_win25000_chr6$Indv8),cor(MAA_gSAsia_win25000_chr6$Indv7, MAA_gSAsia_win25000_chr6$Indv9),
#                          cor(MAA_gSAsia_win25000_chr6$Indv8, MAA_gSAsia_win25000_chr6$Indv9))
#
#corVec_win25000_chr1 <- c(cor(MAA_gSAsia_win25000_chr7$Indv1, MAA_gSAsia_win25000_chr7$Indv2),cor(MAA_gSAsia_win25000_chr7$Indv1, MAA_gSAsia_win25000_chr7$Indv3),cor(MAA_gSAsia_win25000_chr7$Indv1, MAA_gSAsia_win25000_chr7$Indv4), 
#                          cor(MAA_gSAsia_win25000_chr7$Indv1, MAA_gSAsia_win25000_chr7$Indv5),cor(MAA_gSAsia_win25000_chr7$Indv1, MAA_gSAsia_win25000_chr7$Indv6),cor(MAA_gSAsia_win25000_chr7$Indv1, MAA_gSAsia_win25000_chr7$Indv7),
#                          cor(MAA_gSAsia_win25000_chr7$Indv1, MAA_gSAsia_win25000_chr7$Indv8),cor(MAA_gSAsia_win25000_chr7$Indv1, MAA_gSAsia_win25000_chr7$Indv9),
#                          cor(MAA_gSAsia_win25000_chr7$Indv2, MAA_gSAsia_win25000_chr7$Indv3),cor(MAA_gSAsia_win25000_chr7$Indv2, MAA_gSAsia_win25000_chr7$Indv4),cor(MAA_gSAsia_win25000_chr7$Indv2, MAA_gSAsia_win25000_chr7$Indv5),
#                          cor(MAA_gSAsia_win25000_chr7$Indv2, MAA_gSAsia_win25000_chr7$Indv6),cor(MAA_gSAsia_win25000_chr7$Indv2, MAA_gSAsia_win25000_chr7$Indv7),cor(MAA_gSAsia_win25000_chr7$Indv2, MAA_gSAsia_win25000_chr7$Indv8),
#                          cor(MAA_gSAsia_win25000_chr7$Indv2, MAA_gSAsia_win25000_chr7$Indv9),
#                          cor(MAA_gSAsia_win25000_chr7$Indv3, MAA_gSAsia_win25000_chr7$Indv4),cor(MAA_gSAsia_win25000_chr7$Indv3, MAA_gSAsia_win25000_chr7$Indv5),cor(MAA_gSAsia_win25000_chr7$Indv3, MAA_gSAsia_win25000_chr7$Indv6),
#                          cor(MAA_gSAsia_win25000_chr7$Indv3, MAA_gSAsia_win25000_chr7$Indv7),cor(MAA_gSAsia_win25000_chr7$Indv3, MAA_gSAsia_win25000_chr7$Indv8),cor(MAA_gSAsia_win25000_chr7$Indv3, MAA_gSAsia_win25000_chr7$Indv9),
#                          cor(MAA_gSAsia_win25000_chr7$Indv4, MAA_gSAsia_win25000_chr7$Indv5),cor(MAA_gSAsia_win25000_chr7$Indv4, MAA_gSAsia_win25000_chr7$Indv6),cor(MAA_gSAsia_win25000_chr7$Indv4, MAA_gSAsia_win25000_chr7$Indv7),
#                          cor(MAA_gSAsia_win25000_chr7$Indv4, MAA_gSAsia_win25000_chr7$Indv8),cor(MAA_gSAsia_win25000_chr7$Indv4, MAA_gSAsia_win25000_chr7$Indv9),
#                          cor(MAA_gSAsia_win25000_chr7$Indv5, MAA_gSAsia_win25000_chr7$Indv6),cor(MAA_gSAsia_win25000_chr7$Indv5, MAA_gSAsia_win25000_chr7$Indv7),cor(MAA_gSAsia_win25000_chr7$Indv5, MAA_gSAsia_win25000_chr7$Indv8),
#                          cor(MAA_gSAsia_win25000_chr7$Indv5, MAA_gSAsia_win25000_chr7$Indv9),
#                          cor(MAA_gSAsia_win25000_chr7$Indv6, MAA_gSAsia_win25000_chr7$Indv7),cor(MAA_gSAsia_win25000_chr7$Indv6, MAA_gSAsia_win25000_chr7$Indv8),cor(MAA_gSAsia_win25000_chr7$Indv6, MAA_gSAsia_win25000_chr7$Indv9),
#                          cor(MAA_gSAsia_win25000_chr7$Indv7, MAA_gSAsia_win25000_chr7$Indv8),cor(MAA_gSAsia_win25000_chr7$Indv7, MAA_gSAsia_win25000_chr7$Indv9),
#                          cor(MAA_gSAsia_win25000_chr7$Indv8, MAA_gSAsia_win25000_chr7$Indv9))
#



#cor(MAA_gSAsia_win25000_chr1$Indv1, MAA_gSAsia_win25000_chr1$Indv2)
#cor(MAA_gSAsia_win25000_chr1$Indv1, MAA_gSAsia_win25000_chr1$Indv3)
#cor(MAA_gSAsia_win25000_chr1$Indv1, MAA_gSAsia_win25000_chr1$Indv4)
#cor(MAA_gSAsia_win25000_chr1$Indv1, MAA_gSAsia_win25000_chr1$Indv5)
#cor(MAA_gSAsia_win25000_chr1$Indv1, MAA_gSAsia_win25000_chr1$Indv6)
#cor(MAA_gSAsia_win25000_chr1$Indv1, MAA_gSAsia_win25000_chr1$Indv7)
#cor(MAA_gSAsia_win25000_chr1$Indv1, MAA_gSAsia_win25000_chr1$Indv8)
#cor(MAA_gSAsia_win25000_chr1$Indv1, MAA_gSAsia_win25000_chr1$Indv9)
#
#cor(MAA_gSAsia_win25000_chr1$Indv2, MAA_gSAsia_win25000_chr1$Indv3)
#cor(MAA_gSAsia_win25000_chr1$Indv2, MAA_gSAsia_win25000_chr1$Indv4)
#cor(MAA_gSAsia_win25000_chr1$Indv2, MAA_gSAsia_win25000_chr1$Indv5)
#cor(MAA_gSAsia_win25000_chr1$Indv2, MAA_gSAsia_win25000_chr1$Indv6)
#cor(MAA_gSAsia_win25000_chr1$Indv2, MAA_gSAsia_win25000_chr1$Indv7)
#cor(MAA_gSAsia_win25000_chr1$Indv2, MAA_gSAsia_win25000_chr1$Indv8)
#cor(MAA_gSAsia_win25000_chr1$Indv2, MAA_gSAsia_win25000_chr1$Indv9)
#
#cor(MAA_gSAsia_win25000_chr1$Indv3, MAA_gSAsia_win25000_chr1$Indv4)
#cor(MAA_gSAsia_win25000_chr1$Indv3, MAA_gSAsia_win25000_chr1$Indv5)
#cor(MAA_gSAsia_win25000_chr1$Indv3, MAA_gSAsia_win25000_chr1$Indv6)
#cor(MAA_gSAsia_win25000_chr1$Indv3, MAA_gSAsia_win25000_chr1$Indv7)
#cor(MAA_gSAsia_win25000_chr1$Indv3, MAA_gSAsia_win25000_chr1$Indv8)
#cor(MAA_gSAsia_win25000_chr1$Indv3, MAA_gSAsia_win25000_chr1$Indv9)
#
#cor(MAA_gSAsia_win25000_chr1$Indv4, MAA_gSAsia_win25000_chr1$Indv5)
#cor(MAA_gSAsia_win25000_chr1$Indv4, MAA_gSAsia_win25000_chr1$Indv6)
#cor(MAA_gSAsia_win25000_chr1$Indv4, MAA_gSAsia_win25000_chr1$Indv7)
#cor(MAA_gSAsia_win25000_chr1$Indv4, MAA_gSAsia_win25000_chr1$Indv8)
#cor(MAA_gSAsia_win25000_chr1$Indv4, MAA_gSAsia_win25000_chr1$Indv9)
#
#cor(MAA_gSAsia_win25000_chr1$Indv5, MAA_gSAsia_win25000_chr1$Indv6)
#cor(MAA_gSAsia_win25000_chr1$Indv5, MAA_gSAsia_win25000_chr1$Indv7)
#cor(MAA_gSAsia_win25000_chr1$Indv5, MAA_gSAsia_win25000_chr1$Indv8)
#cor(MAA_gSAsia_win25000_chr1$Indv5, MAA_gSAsia_win25000_chr1$Indv9)
#
#cor(MAA_gSAsia_win25000_chr1$Indv6, MAA_gSAsia_win25000_chr1$Indv7)
#cor(MAA_gSAsia_win25000_chr1$Indv6, MAA_gSAsia_win25000_chr1$Indv8)
#cor(MAA_gSAsia_win25000_chr1$Indv6, MAA_gSAsia_win25000_chr1$Indv9)
#
#cor(MAA_gSAsia_win25000_chr1$Indv7, MAA_gSAsia_win25000_chr1$Indv8)
#cor(MAA_gSAsia_win25000_chr1$Indv7, MAA_gSAsia_win25000_chr1$Indv9)
#
#cor(MAA_gSAsia_win25000_chr1$Indv8, MAA_gSAsia_win25000_chr1$Indv9)


