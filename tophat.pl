#!/usr/bin/perl -w
##This perl script used to generate tophad command for pair-end fastq files.
#input is the file list of fastq file names, such like
#SRR1842828_1.fastq
#SRR1842828_2.fastq
#etc
#you could get the file list simply by "ls *.fastq > file_list" of the files you need.
while (<>){
	chomp;
	if(/(.+)\_1.fastq/){
	print "tophat2 -p 7 -g 1 -o /export/work/blair/coral/Tophat/$1 /export/work/blair/coral/refs/acropora /export/work/blair/coral/data/$_ /export/work/blair/coral/data/$1\_2\.fastq 1>$1\.log 2>&1 &\n";
	}elsif(/(.+)\_2.fastq/){
}
}
