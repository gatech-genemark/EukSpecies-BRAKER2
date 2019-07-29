#!/usr/bin/perl
# -------------------
# Alex Lomsadze
# GaTech
# 2019
#
# Seelect records assosiated with pseudogenes from nice GFF3 formatted file
# -------------------

use warnings;
use strict;

my $in = shift;
my $out = shift;

my $text = LoadPseudoRecords($in);

open( my $OUT, ">", $out ) or die "error on open file $out: $!\n";
print $OUT $text;
close $OUT;

# -------------------
sub LoadPseudoRecords
{
	my $name = shift;

	my $pseudo = '';

	open( my $IN, $name ) or die "error on open file $name: $!\n";
	while( <$IN> )
	{
		if ( (/^#/ or /^##/) and !/^###/ )
		{
			$pseudo .= $_;
			next;
		}
		
		next if ( /^\s*$/ );

		my $current = $_;

		while(<$IN>)
		{
			$current .= $_;
			last if ( /^###/ );
		}

		if ( $current =~ /pseudogene/ )  # as in Arabidopsis_thaliana
		{
			$pseudo .= $current;
		}
		elsif ( $current =~ /so_term_name=pseudogene/ ) # as in Caenorhabditis_elegans
		{
			$pseudo .= $current;
		}
	}
	close $IN;

	return $pseudo;
}
# -------------------

