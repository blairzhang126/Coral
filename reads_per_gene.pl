#!/usr/bin/perl -w
#This perl script used to generate samtools and perl command
#input is the sra names of the files
while (<>){
	chomp;
#write samtools view command, which transform bam file into sam file
	print "samtools view -h -o /export/work/blair/coral/Tophat/$_/accepted_hits\.sam /export/work/blair/coral/Tophat/$_/accepted_hits\.bam &\n";

#write perl command to generate raw reads align to each gene
	print "perl /export/work/blair/coral/scripts/map_sam_reads_to_known_genes_rmdup\.pl /export/work/blair/coral/scripts/output\.psl /export/work/blair/coral/Tophat/$_/accepted_hits\.sam /export/work/blair/coral/align_reads/$_\.rmdup\.txt 1>$_.log 2>&1 &\n";
}
