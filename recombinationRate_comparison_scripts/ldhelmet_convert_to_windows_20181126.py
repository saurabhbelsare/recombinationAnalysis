# Script to convert ldhelmet output to window-based estimates of rho
#
# Usage: python ./ldhelmet_convert_to_windows_20181126.py [filename] [window_size] [step_size] [chromosome]
#
# Output: Tab-delimited file with region start coordinates (chrom, window_start), 
# number of sites within window (num_sites), and summed rho for that window (win_rho). 
# Convert to rho/bp by dividing win_rho by num_sites.


import sys
import gzip
import os
import pysam


# set variables
window_size = int(sys.argv[2])
step_size = int(sys.argv[3])
chrom = sys.argv[4]


# open input file (gzipped VCF file), make sure the VCF file is indexed (if not, create index)
filename = sys.argv[1]


# Make sure tabix index is present
if not os.path.exists("%s.tbi" % filename):
    pysam.tabix_index(filename, preset="bed")
parsebed = pysam.Tabixfile(filename)


# Get first and last coordinates
chrom_size={
'1':220367699,
'2':187378091,
'3':180432695,
'4':165299245,
'5':178775436,
'6':174439528,
'7':162156779,
'8':139646187,
'9':125196307,
'10':90941950,
'11':132286798,
'12':104110932,
'13':128036923,
'14':123829720,
'15':107442819,
'16':74645514,
'17':90913898,
'18':72186199,
'19':51301725,
'20':71807805}

BED=gzip.open(filename, 'r')
line=BED.readline()
start_pos = int(line.strip().split()[1])
end_pos = chrom_size[chrom]
BED.close()


# create output file
output = open(filename + '_%swin_%sstep.txt' % (window_size, step_size), 'w')
output.write('chrom\twindow_start\tnum_sites\twin_rho\n')


# Function to calculate rho
def rho_cal(chrom, window_start, window_end):
    window_size=window_end-window_start
    rows = tuple(parsebed.fetch(region="%s:%s-%s" % (chrom, window_start, window_end), parser=pysam.asTuple()))
    rho=[]
    rangestart=int(rows[0][1])
    xshift=window_start-rangestart
    for line in rows:
        startpos=int(line[1])
        endpos=int(line[2])
        rhoperbp=float(line[3])
        interval=endpos-startpos
        rates=[rhoperbp]*interval
        rho=rho+rates
    rho2=rho[xshift:(xshift+window_size)]
    lenrho=len(rho2)
    sumrho=sum(rho2)
    output.write("%s\t%s\t%s\t%s\n" % (chrom, window_start, lenrho, sumrho))
		

# initialize window start and end coordinates
window_start = start_pos
window_end = start_pos + window_size


# calculate rho for window, update window start and end positions, repeat to end of chromosome
while window_end <= end_pos:	
			
	if window_end < end_pos:
		rho_cal(chrom,window_start,window_end)

		window_start = window_start + step_size
		window_end = window_start + window_size

	else:
		rho_cal(chrom,window_start,window_end)
		break	
		
else:
	window_end = end_pos
	rho_cal(chrom,window_start,window_end)


# close files and exit
output.close()
exit()

