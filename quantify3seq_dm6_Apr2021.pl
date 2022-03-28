use strict;
use warnings;


# This is going to output a file that connects all of the CSs to their nearest CDS ends (in the correct orientation)

# It outputs the location of the CDS end, the name of the transcripts associated with the cds end,
# a list of the CSs associated with the CDS end and the total number of reads associated with the cds end,
# all separated by tabs

# Example input:
# perl quantify3seq.pl bam_R1_aligned_notrim_sorted_polyA_CS_positions_3UTRfb.txt 


#input, i.e.: OrderedlastCDSexonends
#input, i.e.: upd_R1_alignedCSpos_Reduc10_3UTRfbext500


my %CDSends;
my %CDSendsNames;
open (INPUT, "/labs/mtfuller/cameron/QUANTSEQ_hsTC_dm6_2020/OrderedlastCDSexonends_Apr2021_r636_new.txt");
while (<INPUT>){
    chomp $_;
    #print $_,"\n";
    my @input = split(/\t/,$_);
    my @loc = split(/\|/,$input[0]);
    chomp $input[1];
    push (@{$CDSends{$loc[0]."|".$loc[1]}},$loc[2]);
    $CDSendsNames{$input[0]} = $input[1];
    #print $input[1],"\n";
}

my %CSsites;
my %CSsitesReads;
open (INPUT2, $ARGV[0]);
while (<INPUT2>){
    chomp $_;
    my @input = split(/\t/,$_);
    my @chr = split(/chr/,$input[0]);
    push (@{$CSsites{$chr[1]."|".$input[5]}},$input[1]);
    $CSsitesReads{$chr[1]."|".$input[5]."|".$input[1]} = $input[3];
}

foreach my $CSsites (keys %CSsites){
     @{$CSsites{$CSsites}} = sort {$a <=> $b} @{$CSsites{$CSsites}};
}

my %CDStoCSsite;
foreach my $CSsiteschr (keys %CSsites){
    my $i = 0;
    my $j = 0;
    if ($CSsiteschr =~ /\+/){ 
        while (($i < scalar @{$CSsites{$CSsiteschr}}) && ($j < scalar @{$CDSends{$CSsiteschr}})){
	    #print $CSsites{$CSsiteschr}[$i]."|".$CDSends{$CSsiteschr}[$j],"\n";
	    if ($CSsites{$CSsiteschr}[$i] < $CDSends{$CSsiteschr}[$j]){
		$i=$i+1;
	    }
	    elsif ($CSsites{$CSsiteschr}[$i] >= $CDSends{$CSsiteschr}[$j]){
		if (defined($CDSends{$CSsiteschr}[$j+1])){
		    if ($CDSends{$CSsiteschr}[$j+1] < $CSsites{$CSsiteschr}[$i]){
		        $j=$j+1;
		    }
		    else{
		        push (@{$CDStoCSsite{$CSsiteschr."|".$CDSends{$CSsiteschr}[$j]}},$CSsites{$CSsiteschr}[$i]);
			#print $CSsites{$CSsiteschr}[$i]."|".$CDSends{$CSsiteschr}[$j],"\n"; 
			$i=$i+1;
		    }
		}
		else{
		    $i = $i + 1;
		}
	    }
	}
    }
    else{
        while (($i < scalar @{$CSsites{$CSsiteschr}}) && ($j < scalar @{$CDSends{$CSsiteschr}})){
	    if ($CSsites{$CSsiteschr}[$i] <= $CDSends{$CSsiteschr}[$j]){
	        push (@{$CDStoCSsite{$CSsiteschr."|".$CDSends{$CSsiteschr}[$j]}},$CSsites{$CSsiteschr}[$i]);
		#print $CSsites{$CSsiteschr}[$i]."|".$CDSends{$CSsiteschr}[$j],"\n"; 
		$i=$i+1;
	    }
	    elsif ($CSsites{$CSsiteschr}[$i] > $CDSends{$CSsiteschr}[$j]){
	        $j = $j + 1;
	    }
	}
    }
}

my @name = split(/\./,$ARGV[0]);

open (OUTPUT, "> $name[0]"."_quantify.txt");
foreach my $keys (keys %CDStoCSsite){
    chop $CDSendsNames{$keys};
    print OUTPUT $keys."\t".$CDSendsNames{$keys}."\t";
    my @info = split(/\|/,$keys);
    my $total = 0;
    my $i = 0;
    foreach my $entry (@{$CDStoCSsite{$keys}}){
        print OUTPUT $entry."*".$CSsitesReads{$info[0]."|".$info[1]."|".$entry};
	$i = $i + 1;
	if ($i < scalar (@{$CDStoCSsite{$keys}})){
	    print OUTPUT "|";
        }
	$total = $total + $CSsitesReads{$info[0]."|".$info[1]."|".$entry};
    }
    print OUTPUT "\t";
    print OUTPUT $total;
    print OUTPUT "\n";
}





# open (CDSloc, "/srv/gsfs0/projects/fuller/gohrilla/ForGonzalo_New/Datasets/dmel-all-CDS-r5.55.fasta");
# my %UTRlocations;
# while (<CDSloc>){
#     chomp $_;
#     if ($_ =~ ">"){
#         my @info = split(/\;/,$_);
# 	my @prename = split(/\-P[A-Z]/,$info[0]);
# 	my @name = split(/\>/,$prename[0]);
# 	my @loc = split(/\=/,$info[1]);
# 	my @chr = split(/\:/,$loc[1]);
# 	print $name[1],"\n";
# 	if ($info[1] =~ /\,/){
# 	    my @start = split(/\(/,$chr[1]);
# 	    my @truestart = split(/\.\./,$start[1]);
# 	    my $start = $truestart[0];
# 	    my @end = split(/\)/,$start[1]);
# 	    my @trueend = split(/\.\./,$end[0]);
# 	    my $end = $trueend[-1];
# 	    if ($info[2] =~ "complement"){
# 	        $start = $start - $extension;
# 		for (my $i = $start; $i <= $end; $i++){
# 		    if (defined($UTRlocations{$chr[0]."|-|".$i})){
# 		        if (!($i - $start < $extension)){
# 			    $UTRlocations{$chr[0]."|-|".$i} = $name[1];
# 			}
# 		    }
# 		    else{
# 		        $UTRlocations{$chr[0]."|-|".$i} = $name[1];
# 		    }
# 		}
# 	    }
# 	    else{
# 	        $end = $end + $extension;
# 		for (my $i = $start; $i <= $end; $i++){
# 		    if (defined($UTRlocations{$chr[0]."|+|".$i})){
# 		        if (!($end - $i < $extension)){
# 			    $UTRlocations{$chr[0]."|+|".$i} = $name[1];
# 			}
# 		    }
# 		    else{
# 		        $UTRlocations{$chr[0]."|+|".$i} = $name[1];
# 		    }
# 		}
# 	    }
# 	}
#         elsif ($info[1] =~ "complement"){
# 	    my @start = split(/\(/,$chr[1]);
# 	    my @actstart = split(/\.\./,$start[1]);
# 	    my @end = split(/\)/,$actstart[1]);
# 	    my $start = $actstart[0];
# 	    my $end = $end[0];
# 	    $start = $start - $extension;
# 	    for (my $i = $start; $i <= $end; $i++){
# 	        if (defined($UTRlocations{$chr[0]."|-|".$i})){
# 		    if (!($i - $start < $extension)){
# 		        $UTRlocations{$chr[0]."|-|".$i} = $name[1];
# 		    }
# 		}
# 		else{
# 		    $UTRlocations{$chr[0]."|-|".$i} = $name[1];
# 		}
# 	    }
#         }
# 	else{
# 	    my @start = split (/\.\./,$chr[1]);
# 	    my @end = split(/\;/,$start[1]);
# 	    my $start = $start[0];
# 	    my $end = $end[0];
# 	    $end = $end + $extension;
# 	    for (my $i = $start; $i <= $end; $i++){
# 	        if (defined($UTRlocations{$chr[0]."|+|".$i})){
# 		    if (!($end - $i < $extension)){
# 		        $UTRlocations{$chr[0]."|+|".$i} = $name[1];
# 		    }
# 		}
# 		else{
# 		    $UTRlocations{$chr[0]."|+|".$i} = $name[1];
# 		}
# 	    }
# 	}
#     }
# }

# my %RNAquantification;
# open (
