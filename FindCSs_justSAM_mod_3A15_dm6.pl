use warnings;
use strict;

if (@ARGV != 2) {
        die "need to provide 2 inputs:SAM file and outputname\n";
}
my ($samfile, $outputfile) = ($ARGV[0], $ARGV[1]);
my $Alength = 15;
my $Amismatch = 3;

open (my $INPUT , "<", $samfile) or die "error opening inputfile: $!\n";
open (my $OUTPUT, ">", $outputfile);

my $count = 0;
print "$samfile\n";
while (<$INPUT>) {
        $count++;
        print "$count\n" unless ($count % 500);
        chomp;
        my $line = $_;
        my @bamfields = split;
        my ($alignment, $readstart, $sequence, $chr) = ($bamfields[1], $bamfields[3], $bamfields[9], $bamfields[2]);
        my $length = length $sequence;
        if ($alignment == 0) {
                my $match = '';
                for (my $j = 0; $j < $Alength; $j++) {
                        $match = join '', $match, 'A';
                }
                if ($sequence =~ m/$match/i) {
                        STRINGLOOP: for (my $i = 0; $i < ($length - ($Alength - 1)); $i++) {
                                my $fragment = substr $sequence, $i, $Alength;
                                if ($fragment eq $match) {
                                        my $start = $readstart + $i - 1;
                                        my $end = $start + $Alength;
                                        my $seq = `echo $chr\"\t\"$start\"\t\"$end\"\t\"tmp | /srv/gsfs0/software/bedtools/2.18.0/bin/fastaFromBed -fi /labs/mtfuller/cameron/QUANTSEQ_hsTC_dm6_2020/dm6/$chr.fa -bed stdin -fo stdout`;
                                        chomp $seq;
                                        my @seqs = split(/\n/,$seq);
                                        my $RefA = ( $seqs[1] =~ tr/A/X/ );
                                        if ($RefA <= ($Alength - $Amismatch)) {
                                                print $OUTPUT "$line\n";
                                        }
                                        last STRINGLOOP;
                                }
                        }
    		}
        } elsif ($alignment == 16) {
                my $match = '';
                for (my $j = 0; $j < $Alength; $j++) {
                        $match = join '', $match, 'T';
                }
                if ($sequence =~ m/$match/i) {
                        STRINGLOOP: for (my $i = ($length - $Alength); $i >= 0; $i--) {
                                my $fragment = substr $sequence, $i, $Alength;
                                if ($fragment eq $match) {
                                        my $start = $readstart + $i - 1;
                                        my $end = $start + $Alength;
                                        my $seq = `echo $chr\"\t\"$start\"\t\"$end\"\t\"tmp | /srv/gsfs0/software/bedtools/2.18.0/bin/fastaFromBed -fi /labs/mtfuller/cameron/QUANTSEQ_hsTC_dm6_2020/dm6/$chr.fa -bed stdin -fo stdout`;
                                        chomp $seq;
                                        my @seqs = split(/\n/,$seq);
                                        my $RefA = ( $seqs[1] =~ tr/T/X/ );
                                        if ($RefA <= ($Alength - $Amismatch)) {
                                                print $OUTPUT "$line\n";
                                        }
                                        last STRINGLOOP;
                                }
                        }
                }
        }
}
close $INPUT;
close $OUTPUT;
