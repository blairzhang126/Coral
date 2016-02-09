#!/usr/bin/perl -w
use strict;

### Description:

# This perl script parse the output of TopHat [given in the Sequence Alignment/Map (SAM) format (http://samtools.sourceforge.net/)]
# and convert it into the raw number of reads aligned to each Acropora gene.

# The input is (1) the output of TopHat in SAM format, and (2) the alignment of the Acropora genes to the Acropora genome in a PSL format.
# The output is a text file with the raw number of reads mapped to each Acropora gene.  
# To avoid PCR duplicates, only paired-end reads that have unique start position in the genome in both pairs were used (as described in Levin et al., Nature Methods, 2010). 

### Usage:
# perl map_sam_reads_to_known_genes_rmdup.pl GENE_INFORMATION_PSL_FORMAT SEQUENCE_ALIGNMENT_SAM_FORMAT OUTPUT_FILE_TEXT_FORMAT
# Usage example:
# perl map_sam_reads_to_known_genes_rmdup.pl AcroporaGenes.psl ConditionX/accepted_hits.sam all_reads_ConditionX.txt

### Creator:
# Shahar Alon, PhD,
# MIT Media lab
# shaharal@mit.edu

### Start of input

my $start_of_read_in_sam="DCV4KXP1"; # to identify the sequencing reads in the SAM file the prefix string of the read name should be mentioned
# the prefix string of the read name can be easily identified by using the following linux command:
# tail -3 ConditionX/accepted_hits.sam
# which can give, for example:
# DCV4KXP1:270:C3R0AACXX:6:1211:14108:8058	73	NODE_999969_length_32266_cov_11.417375	31698	50	74M2I24M	*	0	0	AGCTAGTGTAGGCTGCACTAAATTAATTCTCACACATCACTATTTAAACTTTACAGCAACAAGAGAACCATACACTGGGAGACTTCAGCTTGTCTGATCA	CCCFFFFFHHHHHJJIJJJJJJJJJJJJJJJJJJJJJJJJJIJJJJJJJJJJJIIGIGHJIIFIIEIJJGEGIFIHEHHED8BCEEEDEEDDDCDDDDDD	AS:i:-16	XN:i:0	XM:i:1	XO:i:1	XG:i:2	NM:i:3	MD:Z:88G9	YT:Z:UU	NH:i:1
# DCV4KXP1:270:C3R0AACXX:6:1309:15271:91063	147	NODE_999969_length_32266_cov_11.417375	31709	50	63M2I35M	=	31689	-118	GCTGCACTAAATTAATTCTCACACATCACTATTTAAACTTTACAGCAACAAGAGAACCATACACTGGGAGACTTCAGCTTGTCTGATCACTGAAGAAAAA	CDCECCFFFFFFEGFGGJIHEJHGHF@AIIJIIIHGGIIGIGHHCGJIGFHEGCHEFDDDD@CHGGHGEGHJIGHIIIIIJJIJIGBFHHFHFFFFFCC@	AS:i:-22	XN:i:0	XM:i:2	XO:i:1	XG:i:2	NM:i:4	MD:Z:77G14A5	YT:Z:UU	NH:i:1
# DCV4KXP1:270:C3R0AACXX:6:2105:7768:18421	137	NODE_999969_length_32266_cov_11.417375	31719	50	53M2I28M2I15M	*	0	0	ATTAATTCTCACACATCACTATTTAAACTTTACAGCAACAAGAGAACCATACACTGGGAGACTTCAGCTTGTCTGATCACTGAAGAAAAAAAAATAGAAA	CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJHIGIJJJJJJJJJJJJJJJJJJJHHHHHHFFFFFFEEEDEDDDDBDCDDD	AS:i:-28	XN:i:0	XM:i:1	XO:i:2	XG:i:4	NM:i:5	MD:Z:67G28	YT:Z:UU	NH:i:1
# in this case, the prefix of the sequencing read name is: DCV4KXP1

my $file_gtf=$ARGV[0];
#example: "AcroporaGenes.psl"; # PSL file format (https://genome.ucsc.edu/FAQ/FAQformat.html#format2) that contains the alignment of the Acropora genes to the Acropora genome 
my $file_sam=$ARGV[1];
#"example: ConditionX/accepted_hits.sam"; # the output of TopHat [given in the Sequence Alignment/Map (SAM) format]
my $file_out=$ARGV[2];
#"example: all_reads_ConditionX.txt" # the output file, which contains two columns (tab-delimited):
# (1) the Acropora gene name, and (2) the number of reads mapped to each gene.
# Example of the output file:
# JR970418.1	1
# JR970422.1	0
# JR970424.1	0
# JR970425.1	979
# JR970427.1	8
# JR970431.1	126
# JR970432.1	61
# JR970433.1	4
# JR970434.1	75
# JR970438.1	2

### End of input

# Variables declaration
my $total_number_of_gene_models;
my $total_number_of_scaffolds;
my $string;
my %hash_gene_models;
my %hash_scaffolds;
my $current_scaffold;
my @all_fields;
my $gene_id;
my $location;
my $i;
my %hash_gene_id_in_each_location;
my $total_number_of_sam_lines;
my $start_location;
my $CIGAR;
my $length_of_read;
my $perfect_match_string;
my $inner_i;
my @all_locations;
my $flag_ok;
my $flag_continue;
my %hash_reads_for_each_gene_model;
my $paired_reads_in_genes;
my $total_reads_in_genes;
my %hash_pair_end;
my $read_name;
my @temp_array;
my %hash_rmdup;
my $rmdup_filtered;
my $before_rmdup;
my $field_block_size;
my @all_block_size;
my $field_location_start;
my @all_location_start;
my $j;

# read the number of gene models and the number of genomic scaffolds from the PSL file that contains the alignment of the Acropora genes to the Acropora genome
open FILE_in_models, "<", $file_gtf or die $!;
$total_number_of_gene_models=0;
$total_number_of_scaffolds=0;
while ($string=<FILE_in_models>)
{
    if ($string =~ /(J\w+\.1).+?(NODE.+?)\t/)
    {
        if (!(exists $hash_gene_models{$1}))
        {
            $hash_gene_models{$1}=1;
            $total_number_of_gene_models++;
        }
        if (!(exists $hash_scaffolds{$2}))
        {
            $hash_scaffolds{$2}=1;
            $total_number_of_scaffolds++;
        }        
    }
}
close FILE_in_models;
print "total number of gene models:\t",$total_number_of_gene_models,"\n";
print "total number of scaffolds:\t",$total_number_of_scaffolds,"\n";

# read the total number of SAM lines (reads aligned against the genome)
open FILE_in_tophat, "<", $file_sam or die $!;
$total_number_of_sam_lines=0;
while ($string=<FILE_in_tophat>)
{
    if ($string =~ /$start_of_read_in_sam/)
    {
        $total_number_of_sam_lines++;
    }
}
close FILE_in_tophat;
print "total number of sam lines:\t",$total_number_of_sam_lines,"\n";

# read the genomic locations of the genes (as they appear in the alignment of the Acropora genes to the Acropora genome) into a hash
open FILE_in_models, "<", $file_gtf or die $!;
while ($string=<FILE_in_models>)
{
    if ($string =~ /(J\w+\.1)/)
    {
        $gene_id=$1;
        @all_fields=split ("\t",$string);           
        $current_scaffold=$all_fields[13];
        $field_block_size=$all_fields[18];
        @all_block_size=split (",",$field_block_size);
        $field_location_start=$all_fields[20];
        @all_location_start=split (",",$field_location_start);
        for ($i=0;$i<$#all_location_start;$i++)
        {
            for ($j=0;$j<$all_block_size[$i];$j++)
            {
                $location=$all_location_start[$i]+$j;
                $hash_gene_id_in_each_location{$current_scaffold}{$location}=$gene_id;
            }               
        }
    }
}
close FILE_in_models;

# Check if the genomic locations of each read (as they appear in the SAM file) overlap the genomic location of a gene model,
# and if so count the read as mapped to the gene model
$paired_reads_in_genes=0;
$total_reads_in_genes=0;
$rmdup_filtered=0;
$before_rmdup=0;
open FILE_in_tophat, "<", $file_sam or die $!;
# read each line of the SAM file (each line represent an alignment of a sequencing read against the genome)
while ($string=<FILE_in_tophat>)
{
    if ($string =~ /$start_of_read_in_sam/)
    {
        # parse each SAM line 
        @all_fields=split("\t",$string);
        $read_name=$all_fields[0];
        $current_scaffold=$all_fields[2];
        $start_location=$all_fields[3];
        $CIGAR=$all_fields[5];
        $length_of_read=length($all_fields[9]);
        $perfect_match_string=$length_of_read."M";
        @all_locations=();
        # Only process the read if its genomic start location overlap a genomic location which map to a gene model
        # Get all the genomic locations to which the read align against and place them in a array
        if (exists $hash_gene_id_in_each_location{$current_scaffold}{$start_location})
        {
            if ($CIGAR =~ /([0-9]+)M([0-9]+)N([0-9]+)M/)
            {
                for ($inner_i=$start_location;$inner_i<($start_location+$1);$inner_i++)
                {
                    push(@all_locations,$inner_i);
                }
                for ($inner_i=($start_location+$1+$2);$inner_i<($start_location+$1+$2+$3);$inner_i++)
                {
                    push(@all_locations,$inner_i);
                }
            }
            elsif ($CIGAR =~ /([0-9]+)M([0-9]+)D([0-9]+)M/)
            {
                for ($inner_i=$start_location;$inner_i<($start_location+$1);$inner_i++)
                {
                    push(@all_locations,$inner_i);
                }
                for ($inner_i=($start_location+$1+$2);$inner_i<($start_location+$1+$2+$3);$inner_i++)
                {
                    push(@all_locations,$inner_i);
                }
            }
            elsif ($CIGAR =~ /([0-9]+)M([0-9]+)I([0-9]+)M/)
            {
                for ($inner_i=$start_location;$inner_i<($start_location+$1);$inner_i++)
                {
                    push(@all_locations,$inner_i);
                }
                for ($inner_i=($start_location+$1);$inner_i<($start_location+$1+$3);$inner_i++)
                {
                    push(@all_locations,$inner_i);
                }
            }
            elsif ($CIGAR =~ /$perfect_match_string/)
            {
                for ($inner_i=$start_location;$inner_i<($start_location+$length_of_read);$inner_i++)
                {
                    push(@all_locations,$inner_i);
                }
            }
            else
            {
                $all_locations[0]="no_where";
            }
            
            # test if all the genomic locations of the processed read correspond to only one gene model
            $flag_ok=1;
            $flag_continue=1;
            if (exists $hash_gene_id_in_each_location{$current_scaffold}{$all_locations[0]})
            {
                $gene_id=$hash_gene_id_in_each_location{$current_scaffold}{$all_locations[0]};
                $i=1;
                while (($flag_ok==1)&&($flag_continue==1))
                {
                    if (exists $hash_gene_id_in_each_location{$current_scaffold}{$all_locations[$i]})
                    {
                        if ($hash_gene_id_in_each_location{$current_scaffold}{$all_locations[$i]} ne $gene_id)
                        {
                            $flag_ok=0;
                        }
                    }
                    else
                    {
                        $flag_ok=0;
                    }
                    if ($i<$#all_locations)
                    {
                        $i++;
                    }
                    else
                    {
                        $flag_continue=0;
                    }                   
                }
            }
            else
            {
                $flag_ok=0;
            }               
        
            # if all these genomic locations appear in only one gene model store the read in a hash which gives the number of raw reads per gene 
            if ($flag_ok==1)
            {
                $total_reads_in_genes=$total_reads_in_genes+1;
                if (exists $hash_pair_end{$read_name})
                {
                    if ($hash_pair_end{$read_name}[0] eq $gene_id)
                    {
                        $before_rmdup=$before_rmdup+1;
                        @temp_array=();
                        push(@temp_array,$start_location);
                        push(@temp_array,$hash_pair_end{$read_name}[1]);
                        
                        # To avoid PCR duplicates, only paired-end reads that have unique start position in the genome in both pairs were used
                        # (as described in Levin et al., Nature Methods, 2010). 
                        if (!(exists $hash_rmdup{$gene_id."_".min(@temp_array)."_".max(@temp_array)}))
                        {
                            #print "YES\t",$gene_id,"\t",min(@temp_array),"\t",max(@temp_array),"\n";
                            $hash_rmdup{$gene_id."_".min(@temp_array)."_".max(@temp_array)}=1;
                            $paired_reads_in_genes=$paired_reads_in_genes+1;
                            if (exists $hash_reads_for_each_gene_model{$gene_id})
                            {
                                $hash_reads_for_each_gene_model{$gene_id}=$hash_reads_for_each_gene_model{$gene_id}+1;
                            }
                            else
                            {
                                $hash_reads_for_each_gene_model{$gene_id}=1;
                            }
                        }
                        else
                        {
                            #print "NO\t",$gene_id,"\t",min(@temp_array),"\t",max(@temp_array),"\n";
                            $rmdup_filtered=$rmdup_filtered+1;
                        }
                        delete $hash_pair_end{$read_name};
                    }
                    else
                    {
                        delete $hash_pair_end{$read_name};
                    }
                }
                else
                {
                    $hash_pair_end{$read_name}[0]=$gene_id;
                    $hash_pair_end{$read_name}[1]=$start_location;
                }
            }
        }
    }
}
close FILE_in_tophat;
# print summary of the analysis
print "total reads aligned to genes:\t",$total_reads_in_genes,"\n";
print "paired reads aligned to genes before rmdup:\t",$before_rmdup,"\n";
print "removed due to PCR duplicates:\t",$rmdup_filtered,"\n";
print "paired reads aligned to genes:\t",$paired_reads_in_genes,"\n";

# print the output file, which contains two columns (tab-delimited):
# (1) the Acropora gene name, and (2) the number of reads mapped to each gene.
open FILE_out, ">", $file_out or die $!;
foreach (sort keys %hash_gene_models)
{
    $gene_id=$_;
    if (exists $hash_reads_for_each_gene_model{$gene_id})
    {
        print FILE_out $gene_id,"\t",$hash_reads_for_each_gene_model{$gene_id},"\n";
    }
    else
    {
        print FILE_out $gene_id,"\t0\n";
    }
}
close FILE_out;

############################################################################################
sub max {
    my @numbers = @_;
    my $max;
    my $running=-1;
    for($running=0;$running<=$#numbers;$running++)
    {
        if (defined($numbers[$running]))
        {
            $max = $numbers[$running];
            last;
        }
    } 
    
    for($running=0;$running<=$#numbers;$running++)
    {
        if (defined($numbers[$running]))
        {
            if ($numbers[$running] > $max)
            {
                $max = $numbers[$running];
            }
        }
    }
    return $max;
}

############################################################################################
sub min {
    my @numbers = @_;
    my $min;
    my $running=-1;
    for($running=0;$running<=$#numbers;$running++)
    {
        if (defined($numbers[$running]))
        {
            $min = $numbers[$running];
            last;
        }
    } 
    
    for($running=0;$running<=$#numbers;$running++)
    {
        if (defined($numbers[$running]))
        {
            if ($numbers[$running] < $min)
            {
                $min = $numbers[$running];
            }
        }
    }
    return $min;
}
