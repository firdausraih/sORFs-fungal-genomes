#!/usr/bin/perl -w
use strict;
use warnings;
use threads;

################################################################################
# igr_gff.pl                                                                   #
#                                                                              #
# This script calculates the IGRs based on feature.gff and output the IGRs in  #
# gff format.                                                                  #
#                                                                              #
#                                                                              #
# IMPORTANT!                                                                   #
# edit 10 mac2014 using get_IGR.pl                                             #
# lastedit 20may2015                                    		       #
#                                                                              #
################################################################################

my $usage = "perl $0 <input.gff3> <genome.len>\n";

my $key = "gene"; # key for gene features, usually "gene" or "CDS".

my $line;
my @line = ();
my @ctglen = ();

open SEQLEN, $ARGV[1] || die "$usage";
while ($line = <SEQLEN>) {
	chomp $line;
	@line = split(/\t/,$line);
	push @ctglen, [@line];
}

my $c;

for $c (0..$#ctglen) {
	my @gff_fwd = ();
	my @gff_rvs = ();
	
	open GFF, $ARGV[0] || die "$usage";
	while ($line = <GFF>) {
        chomp $line;
        @line = split(/\t/,$line);
        push @gff_fwd, [@line] if ($line[0] eq $ctglen[$c][0] && $line[2] eq $key && $line[6] eq "+");
        push @gff_rvs, [@line] if ($line[0] eq $ctglen[$c][0] && $line[2] eq $key && $line[6] eq "-");
	}
	
	my $i;
	my $start;
	my $end;
	
	my $count = 1;
	
	my @gff = @gff_fwd;
	for $i (0..$#gff) {
		if ( $i == 0 ) {
			$start = 1;
			$end = $gff[$i][3]-1;
		} else {
			$start = $gff[$i-1][4]+1;
			$end = $gff[$i][3]-1;
		}
		print "$ctglen[$c][0]\tIGR\tIGR\t$start\t$end\t.\t+\t.\ttranscript_id \"IGR_$count\"; gene_id \"IGR_$count\";\n";
		$count++;
		
		if ( $i == $#gff ) {
			$start = $gff[$i][4]+1;
			$end = $ctglen[$c][1];
			print "$ctglen[$c][0]\tIGR\tIGR\t$start\t$end\t.\t+\t.\ttranscript_id \"IGR_$count\"; gene_id \"IGR_$count\";\n";
			$count++;
		}
	}
	
	@gff = @gff_rvs;
	for $i (0..$#gff) {
		if ( $i == 0 ) {
			$start = 1;
			$end = $gff[$i][3]-1;
		} else {
			$start = $gff[$i-1][4]+1;
			$end = $gff[$i][3]-1;
		}
		print "$ctglen[$c][0]\tIGR\tIGR\t$start\t$end\t.\t-\t.\ttranscript_id \"IGR_$count\"; gene_id \"IGR_$count\";\n";
		$count++;
		
		if ( $i == $#gff ) {
			$start = $gff[$i][4]+1;
			$end = $ctglen[$c][1];
			print "$ctglen[$c][0]\tIGR\tIGR\t$start\t$end\t.\t-\t.\ttranscript_id \"IGR_$count\"; gene_id \"IGR_$count\";\n";
			$count++;
		}
	}
}
