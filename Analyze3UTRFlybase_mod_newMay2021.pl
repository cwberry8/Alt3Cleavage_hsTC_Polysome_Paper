use strict;
use warnings;

# This script takes as input the flybase file containing a list of all three prime UTRs (dmel-all-three_prime_UTR-r5.55.fasta)
# It also takes in a set of genomic locations (CSs_aly_Goh2_min10red50.bed)
# Outputs the percentage of genomic locations that overlap with three prime UTR locations
# Can also indicate an 'extension' length from known UTRs

# Check the input variables for the CS data set(s), might need to change depending on the input file

# Incorporating a 5'UTR (plus ncRNAs because they can be close to the end of the 3UTR of neighboring genes) mask to prevent extending 3'UTR into adjacent genes! -April+May 2021

open (FIVEUTR, "/labs/mtfuller/cameron/QUANTSEQ_hsTC_dm6_2020/dmel-all-five_prime_UTR-ncRNA-r6.36.fasta");
my %FIVEUTR;
while (<FIVEUTR>){
  chomp $_;
  if ($_ =~ /\>/){
    my @info = split(/\s/,$_);
    my @loc = split(/\=/,$info[2]);

    if ($info[2] =~ "complement"){
      if ($info[2] =~ /\,/){
	my @chr = split(/\:complement\(join\(/,$loc[1]);
	my @chr2 = split(/\)\)\;/,$chr[1]);
	my @startend = split(/\,/,$chr2[0]);

	foreach my $startend (@startend){
	  my @coord = split(/\.\./,$startend);
	  for (my $z = $coord[0]; $z <= $coord[1]; $z++){
	    $FIVEUTR{$chr[0]."|-|".$z} = 1;
	  }	
	}

      }
      else{
	my @chr = split(/\:complement\(/,$loc[1]);
	my @chr2 = split(/\)\;/,$chr[1]);
	my @startend = split(/\.\./,$chr2[0]);
	
	for (my $z = $startend[0]; $z <= $startend[1] ;$z++){
	  $FIVEUTR{$chr[0]."|-|".$z} = 1;
	}
	
      }
    }
    else{
      if ($info[2] =~ /\,/){
	my @chr = split(/\:join\(/,$loc[1]);
	my @chr2 = split(/\);/,$chr[1]);
	my @startend = split(/\,/,$chr2[0]);
	foreach my $startend (@startend){
	  my @coord = split(/\.\./,$startend);
	  for (my $z = $coord[0]; $z <= $coord[1]; $z++){
	    $FIVEUTR{$chr[0]."|+|".$z} = 1;
	  }			    
	}
      }
      else{
	my @chr = split(/\:/,$loc[1]);
	my @chr2 = split(/\;/,$chr[1]);
	my @startend = split(/\.\./,$chr2[0]);
	
	for (my $z = $startend[0]; $z <= $startend[1] ;$z++){
	  $FIVEUTR{$chr[0]."|+|".$z} = 1;
	  my $test = $chr[0]."|+|".$z;
	  #if ($test eq "X|+|7329125"){
	  #  print "It recognized 5UTR of correctly","\n";
	 # }
	}
      }
    }
  }
}


open (INPUT, "/labs/mtfuller/cameron/Polysome_hsTC_3SEQ/dmel-all-three_prime_UTR-r6.36.fasta");
my $extension = 500;
my %UTRlocations;
while (<INPUT>){
    chomp $_;
    if ($_ =~ /\>/){
        my @info = split(/\s/,$_);
	my @loc = split(/\=/,$info[2]);
	my @chr = split(/\:/,$loc[1]);
	my $chr = $chr[0];
	my $dir = "+";
	
	if ($info[2] =~ /\,/){
	    my @start = split(/\(/,$chr[1]);
	    my @end = split(/\)/,$start[1]);
	    my @UTRs = split(/\,/,$end[0]);
	    my $start;
	    my $end;
	    if (($info[2] =~ "complement")){
	        $dir = "-";
		my @start2 = split(/\.\./,$start[2]);
		my $startextension = $start2[0] - $extension;

# check if the extension puts the 3'UTR in the 5'UTR of an adjacent gene. If it is, find the 3' most nucleotide in the extension that does not map to an adjacent 5'UTR
		
		my $done = 0;

		for (my $y = $start2[0]; $y > $startextension; $y--){
		  if ($done == 0){
		    if (defined($FIVEUTR{$chr."|".$dir."|".$y})){
		      $done = 1;
		      $start = $y + 1;
		    }
		  }
		}	       
		if ($done == 0){
		  $start = $start2[0] - $extension;
		}
	      
		my @end2 = split(/\)\)/,$start[2]);
		my @end3 = split(/\.\./,$end2[0]);
		$end = $end3[-1];
	    }
	    else{
		my @start2 = split(/\.\./,$start[1]);
		$start = $start2[0];
		my @end2 = split(/\.\./,$end[0]);
		my $endextension = $end2[-1] + $extension;

# check if the extension puts the 3'UTR in the 5'UTR of an adjacent gene. If it is, find the 3' most nucleotide in the extension that does not map to an adjacent 5'UTR

		my $done = 0;

		for (my $y = $end2[-1]; $y < $endextension; $y++){
		  if ($done == 0){
		    if (defined($FIVEUTR{$chr."|".$dir."|".$y})){
		      $done = 1;
		      $end = $y - 1;
		    }
		  }
		}
		if ($done == 0){
		  $end = $end2[-1] + $extension;
		}
		
	    }
            for (my $i = $start; $i <= $end; $i++){
                $UTRlocations{$chr."|".$i."|".$dir} = 1;
            }
	}
        elsif ($info[2] =~ "complement"){
	  $dir = "-";
	  my @start = split(/\(/,$chr[1]);
	  my @actstart = split(/\.\./,$start[1]);
	  my @end = split(/\)/,$actstart[1]);
	  my $start = $actstart[0];
	  my $end = $end[0];
	  my $startextension = $start - $extension;

	  my $done = 0;

	  for (my $y = $start; $y > $startextension ;$y--){
	    if ($done == 0){
	      if (defined($FIVEUTR{$chr."|".$dir."|".$y})){
		$done = 1;
		$start = $y + 1;
	      }
	    }
	  }	       
	  if ($done == 0){
	    $start = $start - $extension;
	  }
	       
	  for (my $i = $start; $i <= $end; $i++){
	    $UTRlocations{$chr."|".$i."|".$dir} = 1;
	  }
        }
	else{
	    my @start = split (/\.\./,$chr[1]);
	    my @end = split(/\;/,$start[1]);
	    my $start = $start[0];
	    my $end = $end[0];
	    my $endextension = $end + $extension;
	    my $done = 0;
	    
	    for (my $y = $end; $y < $endextension; $y++){
	      if ($done == 0){
		if (defined($FIVEUTR{$chr."|".$dir."|".$y})){
		  $done = 1;
		  $end = $y - 1;
		}
	      }
	    }
	    if ($done == 0){
	      $end = $end + $extension;
	    }


	    for (my $i = $start; $i <= $end; $i++){
	      $UTRlocations{$chr."|".$i."|".$dir} = 1;
	    }
	}
    }
}


my $total = 0;
my $CSinUTR = 0;

open (INPUT2, $ARGV[0]);
my @name = split(/\./,$ARGV[0]);
my %CSsites;
my %readcount;

open (OUTPUT, "> $name[0]"."_3UTRfb_May2021.txt");
#open (OUTPUT, "> ore_R1_alignedCSpos_Reduc10_notin3UTRfb.txt");
while (<INPUT2>){
    chomp $_;
    my @info = split(/\t/,$_);
    my $chr = substr($info[0],3,2);
    $total = $total + 1;
    #print $chr."|".$info[1]."|".$info[3],"\n";
    if (defined($UTRlocations{$chr."|".$info[1]."|".$info[-1]})){
        push (@{$CSsites{$chr."|".$info[-1]}},$info[1]);
	$readcount{$chr."|".$info[-1]."|".$info[1]} = $info[2];
        #print OUTPUT join ("\t",$info[0],$info[1],$info[1],$info[2],$info[2],$info[3]);
	#print OUTPUT "\n";
	$CSinUTR = $CSinUTR + 1;
    }
}

foreach my $keys (keys %CSsites){
    @{$CSsites{$keys}} = sort {$a <=> $b} @{$CSsites{$keys}};
    foreach my $pos (@{$CSsites{$keys}}){
        my @info = split(/\|/,$keys);
	  #if (!($info[0] eq "U" || $info[0] eq "Uextra" || $info[0] eq "M" || $info[0] eq "XHet" ||$info[0] eq "YHet")){
	      print OUTPUT join("\t","chr".$info[0],$pos,$pos,$readcount{$keys."|".$pos},$readcount{$keys."|".$pos},$info[1]);
	      print OUTPUT "\n";
	  #}
    }
}


print $total."|".$CSinUTR,"\n";


# 8436|4195
# in bam, 4195 of 8436 CSs lie in known UTRs
