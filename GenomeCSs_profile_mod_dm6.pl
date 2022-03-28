use warnings;
use strict;

if (@ARGV != 2) {
        die "need to provide 2 inputs:SAM file and outputname\n";
}
my ($samfile, $outputfile) = ($ARGV[0], $ARGV[1]);
my $Alength = 8;
my $Amismatch = 5;

open (my $INPUT , "<", $samfile) or die "error opening inputfile: $!\n";
open (my $OUTPUT, ">", $outputfile);

my $count = 0;
print "$samfile\n";
my %csshashplus;
my %csshashminus;
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
                                        my $CSSPOS = $readstart + $i - 1;
                                        my $CSSkey = join "\t", $chr, $CSSPOS;
                                        if ($csshashplus{$CSSkey}) {
                                                $csshashplus{$CSSkey} = $csshashplus{$CSSkey} + 1;
                                        } else {
                                                $csshashplus{$CSSkey} = 1;
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
                                        my $CSSPOS = $readstart + $i + $Alength;
                                        my $CSSkey = join "\t", $chr, $CSSPOS;
                                        if ($csshashminus{$CSSkey}) {
                                                $csshashminus{$CSSkey} = $csshashminus{$CSSkey} + 1;
                                        } else {
                                                $csshashminus{$CSSkey} = 1;
                                        }
                                        last STRINGLOOP;
                                }
                        }
                }
        }
}
close $INPUT;
print $OUTPUT "$_\t$csshashplus{$_}\t+\n" foreach (keys %csshashplus);
print $OUTPUT "$_\t$csshashminus{$_}\t-\n" foreach (keys %csshashminus);
close $OUTPUT;
