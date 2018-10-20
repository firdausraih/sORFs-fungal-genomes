#!/usr/bin/perl
use warnings;

my $file = shift;
my $line;
my $cells;

if (!defined ($file)) {
    die "USAGE: perl <script> <input_file>\n";
}

open (GFF, "$file") or die "Cannot open input file$!\n";

while ($line = <GFF>)  {
	chomp $line;
	my @cells = split /\t/, $line;
	#change any type for example here is exon
	if ($cells[2] eq "gene") {
		print "$line \n";
	}
}
