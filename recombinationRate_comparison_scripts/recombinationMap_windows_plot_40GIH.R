library(readr)

######Human
directoryName <- "~/workfiles/recombinationAnalysis/human/GIH/40Indv/"
prefixSeg <- "GIH_1KG_firstHalf.chr"
suffixSeg <- "_GRCh38_cleanedup_removeDuplicates.bed.gz_100000win_100000step.txt"

selChrWindowRho_chr1 <- read_delim(paste(directoryName,prefixSeg,"1",suffixSeg,sep=""), 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr2 <- read_delim(paste(directoryName,prefixSeg,"2",suffixSeg,sep=""), 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr3 <- read_delim(paste(directoryName,prefixSeg,"3",suffixSeg,sep=""), 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr4 <- read_delim(paste(directoryName,prefixSeg,"4",suffixSeg,sep=""), 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr5 <- read_delim(paste(directoryName,prefixSeg,"5",suffixSeg,sep=""), 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr6 <- read_delim(paste(directoryName,prefixSeg,"6",suffixSeg,sep=""), 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr7 <- read_delim(paste(directoryName,prefixSeg,"7",suffixSeg,sep=""), 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr8 <- read_delim(paste(directoryName,prefixSeg,"8",suffixSeg,sep=""), 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr9 <- read_delim(paste(directoryName,prefixSeg,"9",suffixSeg,sep=""), 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr10 <- read_delim(paste(directoryName,prefixSeg,"10",suffixSeg,sep=""), 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr11 <- read_delim(paste(directoryName,prefixSeg,"11",suffixSeg,sep=""), 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr12 <- read_delim(paste(directoryName,prefixSeg,"12",suffixSeg,sep=""), 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr13 <- read_delim(paste(directoryName,prefixSeg,"13",suffixSeg,sep=""), 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr14 <- read_delim(paste(directoryName,prefixSeg,"14",suffixSeg,sep=""), 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr15 <- read_delim(paste(directoryName,prefixSeg,"15",suffixSeg,sep=""), 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr16 <- read_delim(paste(directoryName,prefixSeg,"16",suffixSeg,sep=""), 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr17 <- read_delim(paste(directoryName,prefixSeg,"17",suffixSeg,sep=""), 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr18 <- read_delim(paste(directoryName,prefixSeg,"18",suffixSeg,sep=""), 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr19 <- read_delim(paste(directoryName,prefixSeg,"19",suffixSeg,sep=""), 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr20 <- read_delim(paste(directoryName,prefixSeg,"20",suffixSeg,sep=""), 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr21 <- read_delim(paste(directoryName,prefixSeg,"21",suffixSeg,sep=""), 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
selChrWindowRho_chr22 <- read_delim(paste(directoryName,prefixSeg,"22",suffixSeg,sep=""), 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)




allChrRhoVec <- c(selChrWindowRho_chr1$win_rho,selChrWindowRho_chr2$win_rho,selChrWindowRho_chr3$win_rho,selChrWindowRho_chr4$win_rho,
                  selChrWindowRho_chr5$win_rho,selChrWindowRho_chr6$win_rho,selChrWindowRho_chr7$win_rho,selChrWindowRho_chr8$win_rho,
                  selChrWindowRho_chr9$win_rho,selChrWindowRho_chr10$win_rho,selChrWindowRho_chr11$win_rho,selChrWindowRho_chr12$win_rho,
                  selChrWindowRho_chr13$win_rho,selChrWindowRho_chr14$win_rho,selChrWindowRho_chr15$win_rho,selChrWindowRho_chr16$win_rho,
                  selChrWindowRho_chr17$win_rho,selChrWindowRho_chr18$win_rho,selChrWindowRho_chr19$win_rho,selChrWindowRho_chr20$win_rho,
                  selChrWindowRho_chr21$win_rho,selChrWindowRho_chr22$win_rho)

meanRho <- mean(allChrRhoVec)


#plot(selChrWindowRho_chr1$window_start,selChrWindowRho_chr1$win_rho/meanRho,type="l")
#plot(selChrWindowRho_chr2$window_start,selChrWindowRho_chr2$win_rho/meanRho,type="l")
#plot(selChrWindowRho_chr7$window_start,selChrWindowRho_chr7$win_rho/meanRho,type="l")
#plot(selChrWindowRho_chr8$window_start,selChrWindowRho_chr8$win_rho/meanRho,type="l")
#plot(selChrWindowRho_chr12$window_start,selChrWindowRho_chr12$win_rho/meanRho,type="l")
#plot(selChrWindowRho_chr13$window_start,selChrWindowRho_chr13$win_rho/meanRho,type="l")
#
#xVec <- c(rep(1,length(selChrWindowRho_chr1$window_start)),rep(2,length(selChrWindowRho_chr2$window_start)),
#          rep(3,length(selChrWindowRho_chr3$window_start)),rep(4,length(selChrWindowRho_chr4$window_start)),
#          rep(5,length(selChrWindowRho_chr5$window_start)),rep(6,length(selChrWindowRho_chr6$window_start)),
#          rep(7,length(selChrWindowRho_chr7$window_start)),rep(8,length(selChrWindowRho_chr8$window_start)),
#          rep(9,length(selChrWindowRho_chr9$window_start)),rep(10,length(selChrWindowRho_chr10$window_start)),
#          rep(11,length(selChrWindowRho_chr11$window_start)),rep(12,length(selChrWindowRho_chr12$window_start)),
#          rep(13,length(selChrWindowRho_chr13$window_start)),rep(14,length(selChrWindowRho_chr14$window_start)),
#          rep(15,length(selChrWindowRho_chr15$window_start)),rep(16,length(selChrWindowRho_chr16$window_start)),
#          rep(17,length(selChrWindowRho_chr17$window_start)),rep(18,length(selChrWindowRho_chr18$window_start)),
#          rep(19,length(selChrWindowRho_chr19$window_start)),rep(20,length(selChrWindowRho_chr20$window_start)))

allChrRhoVecNorm <- c(selChrWindowRho_chr1$win_rho/meanRho,selChrWindowRho_chr2$win_rho/meanRho,selChrWindowRho_chr3$win_rho/meanRho,selChrWindowRho_chr4$win_rho/meanRho,
                      selChrWindowRho_chr5$win_rho/meanRho,selChrWindowRho_chr6$win_rho/meanRho,selChrWindowRho_chr7$win_rho/meanRho,selChrWindowRho_chr8$win_rho/meanRho,
                      selChrWindowRho_chr9$win_rho/meanRho,selChrWindowRho_chr10$win_rho/meanRho,selChrWindowRho_chr11$win_rho/meanRho,selChrWindowRho_chr12$win_rho/meanRho,
                      selChrWindowRho_chr13$win_rho/meanRho,selChrWindowRho_chr14$win_rho/meanRho,selChrWindowRho_chr15$win_rho/meanRho,selChrWindowRho_chr16$win_rho/meanRho,
                      selChrWindowRho_chr17$win_rho/meanRho,selChrWindowRho_chr18$win_rho/meanRho,selChrWindowRho_chr19$win_rho/meanRho,selChrWindowRho_chr20$win_rho/meanRho,
                      selChrWindowRho_chr21$win_rho/meanRho,selChrWindowRho_chr22$win_rho/meanRho)

length1 <- length(selChrWindowRho_chr1$window_start)
length2 <- length(selChrWindowRho_chr2$window_start)
length3 <- length(selChrWindowRho_chr3$window_start)
length4 <- length(selChrWindowRho_chr4$window_start)
length5 <- length(selChrWindowRho_chr5$window_start)
length6 <- length(selChrWindowRho_chr6$window_start)
length7 <- length(selChrWindowRho_chr7$window_start)
length8 <- length(selChrWindowRho_chr8$window_start)
length9 <- length(selChrWindowRho_chr9$window_start)
length10 <- length(selChrWindowRho_chr10$window_start)
length11 <- length(selChrWindowRho_chr11$window_start)
length12 <- length(selChrWindowRho_chr12$window_start)
length13 <- length(selChrWindowRho_chr13$window_start)
length14 <- length(selChrWindowRho_chr14$window_start)
length15 <- length(selChrWindowRho_chr15$window_start)
length16 <- length(selChrWindowRho_chr16$window_start)
length17 <- length(selChrWindowRho_chr17$window_start)
length18 <- length(selChrWindowRho_chr18$window_start)
length19 <- length(selChrWindowRho_chr19$window_start)
length20 <- length(selChrWindowRho_chr20$window_start)
length21 <- length(selChrWindowRho_chr21$window_start)
length22 <- length(selChrWindowRho_chr22$window_start)


tick0 <- 1
tick1 <- length1
tick2 <- tick1 + length2
tick3 <- tick2 + length3
tick4 <- tick3 + length4
tick5 <- tick4 + length5
tick6 <- tick5 + length6
tick7 <- tick6 + length7
tick8 <- tick7 + length8
tick9 <- tick8 + length9
tick10 <- tick9 + length10
tick11 <- tick10 + length11
tick12 <- tick11 + length12
tick13 <- tick12 + length13
tick14 <- tick13 + length14
tick15 <- tick14 + length15
tick16 <- tick15 + length16
tick17 <- tick16 + length17
tick18 <- tick17 + length18
tick19 <- tick18 + length19
tick20 <- tick19 + length20
tick21 <- tick20 + length21
tick22 <- tick21 + length22








plot(allChrRhoVecNorm,xaxt='n',type = "l",xlab="Chromosome",ylab="Normalized rho/bp",title="GIH 40 Individuals",ylim = c(0,400))
axis(1, at = c(tick1,tick2,tick3,tick4,tick5,tick6,tick7,tick8,tick9,tick10,tick11,tick12,tick13,tick14,tick15,tick16,tick17,tick18,tick19,tick20,tick21,tick22) - 0.5*c(length1,length2,length3,length4,length5,length6,length7,length8,length9,length10,length11,length12,length13,length14,length15,length16,length17,length18,length19,length20,length21,length22), labels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"), tick = FALSE)
axis(1, at = c(tick0,tick1,tick2,tick3,tick4,tick5,tick6,tick7,tick8,tick9,tick10,tick11,tick12,tick13,tick14,tick15,tick16,tick17,tick18,tick19,tick20,tick21,tick22), labels = FALSE)
#axis(side=1,at=c(tick0,tick1,tick2,tick3,tick4,tick5,tick6,tick7,tick8,tick9,tick10,tick11,tick12,tick13,tick14,tick15,tick16,tick17,tick18,tick19,tick20,tick21,tick22),labels=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"))




GIH_recombination_map_hg38_chr_1 <- read_delim("Downloads/hg38/GIH/GIH_recombination_map_hg38_chr_1.bed", 
                                               "\t", escape_double = FALSE, trim_ws = TRUE)
plot(GIH_recombination_map_hg38_chr_1$Start,GIH_recombination_map_hg38_chr_1$reco_rate_per_base_per_generation, type = "l")





