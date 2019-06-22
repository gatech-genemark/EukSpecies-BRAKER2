#!/usr/bin/perl
#---------------------------------
# Alex Lomsadze
# GaTech 2019
#
# move GFF3 to GTF with stops included
# --------------------------------

use strict;
use warnings;

my $gff3 = shift;
my $gtf = shift;
my $v = 1;

if ( !$gff3 ) { die "error, input file name is missing"; }
if ( !$gtf ) { die "error, output file name is missing"; }

my %trans_to_gene = ();

LoadIds($gff3, \%trans_to_gene );

open( my $IN, $gff3) or die "error on open file: $gff3\n";
open( my $OUT, ">", $gtf) or die "error on open file: $gtf\n";
while(<$IN>)
{
	if ( /\tCDS\t/ )
	{
		my $trans = '';

		chomp;

		if ( /Parent=([^;]+)/ )
		{
			$trans = $1;
		}

		if ( !$trans )
			{ die "error, unexpected line format found: $_"; }

		my @arr = split( ',' , $trans );

		foreach my $label (@arr)
		{
			if ( ! exists $trans_to_gene{$label} )
				{ die "error, Parent is missing for: $label\n";}

			my $line = $_;
			$line =~ s/\tParent=.*$/\t/;
			$line =~ s/\tID=.*$/\t/;

			print $OUT $line;

			my $attr = "gene_id \"$trans_to_gene{$label}\"\;";
			$attr .= " transcript_id \"$label\"\;\n";

			print $OUT $attr;
		}
	}
}
close $OUT;
close $IN;

# -------------------
sub LoadIds
{
	my $name = shift;
	my $ref = shift;

	open( my $IN, $name) or die "error on open file: $name\n";
	while(<$IN>)
	{
		if ( /^#/ ) {next;}

		if ( /\tmRNA\t/ )
		{
#			if ($v) { print $_; }

			my $trans = '';
			my $gene = '';

			if ( /ID=([^;]+)/ )
			{
				$trans = $1;
			}

			if ( /Parent=([^;]+)/ )
			{
				$gene = $1;
			}

			if( !$trans or !$gene )
				{ die "error, unexpect format found: $_"; }

			$ref->{$trans} = $gene;

#			if ($v) { print "$trans $gene\n"; }
		}
	}
	close $IN;
}
# -------------------


