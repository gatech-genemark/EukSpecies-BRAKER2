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
	if ( /\tCDS\t/ or /\t[Ss]tart_codon\t/ or /\t[Ss]top_codon\t/ or /\t[Ii]ntron\t/ or /\tgap\t/ )
	{
		chomp;

		my $trans = '';

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

			my @gff = split( '\t' , $_ );

			my $line = '';
			for( my $i = 0; $i < 8; $i += 1 )
			{
				$line .= $gff[$i] ."\t";
			}

			$line .= "gene_id \"$trans_to_gene{$label}\"\; ";
			$line .= "transcript_id \"$label\"\;";

			if( $gff[8] =~ /cds_type=(\S+?);/ )
			{
				$line .= " cds_type \"$1\"\;";
			}
			if( $gff[8] =~ /count=(\S+?);/ )
			{
				$line .= " count \"$1\"\;";
			}
			if( $gff[8] =~ /site_seq=(\S+?);/ )
			{
				$line .= " site_seq \"$1\"\;";;
			}

			$line .= "\n";

			print $OUT $line;
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

		if ( /\tmRNA\t/ or /\t\S_gene_segment\t/ )
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


