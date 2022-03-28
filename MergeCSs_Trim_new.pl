use strict;
use warnings;

# This takes in input from quantify3seq.pl, i.e. bam_R1_aligned_notrim_sorted_polyA_CS_positions_3UTRfb_quantify.txt
# It outputs CSs associated with CDSs that have number of reads >= $readcutoff and also
# are a local maximum, the locality is determined by $mergebase, for example if $mergebase = 100
# the local region is 100 bases

# example input:
# perl MergeCSs_Trim.pl bam_R1_aligned_notrim_sorted_polyA_CS_positions_3UTRfb_quantify.txt 10 

my $readcutoff = $ARGV[1];
my $mergebase = 50;


open (INPUT, $ARGV[0]);
my @name = split(/\./,$ARGV[0]);
open (OUTPUT, "> $name[0]"."_trim$readcutoff"."_merge$mergebase".".txt");
while (<INPUT>){
    chomp $_;
    #print $_,"\n";
    my @data = split(/\t/,$_);
    my @CSs = split(/\|/,$data[2]);
    my $currentpos = 0;
    my $currentmaxreads = 0;
    my $currentmaxpos = 0;
    my $count = 0;
    my @positions;
    foreach my $CSs (@CSs){
        my @reads = split(/\*/,$CSs);
	if ($reads[1] >= $readcutoff){
	    if ($reads[0] - $currentpos > $mergebase){
	        if ($count != 0){
		    push (@positions,$currentmaxpos."*".$currentmaxreads);
		    #print $data[0]."\t".$data[1]."\t".$currentmaxpos."|".$currentmaxreads,"\n";
		}
		$currentmaxpos = $reads[0];
		$currentmaxreads = $reads[1];
		$count = 1;
	    }
	    else{
	        if ($reads[1] - $currentmaxreads > 0){
		    $currentmaxreads = $reads[1];
		    $currentmaxpos = $reads[0];
		}
	    }
	    $currentpos = $reads[0];
	}
    }
    if ($currentmaxreads != 0){
	push (@positions,$currentmaxpos."*".$currentmaxreads);
	#print $data[0]."\t".$data[1]."\t".$currentmaxpos."|".$currentmaxreads,"\n";
    }
    if (defined($positions[0])){
        print OUTPUT $data[0]."\t".$data[1]."\t";
	print OUTPUT join("|",@positions);
	print OUTPUT "\n";
    }
}
