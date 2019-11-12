#!/usr/bin/perl
# ==============================================================
# Alex Lomsadze
# Georgia Institute of Technology, Atlanta, Georgia, US
# Last update 2019
# 
# This script takes as input two files: 
# (1) gene annotation in "enriched" GTF format and 
# (2) intron annotation in GFF format, 
# and outputs genes which have introns matching intron input.
# ==============================================================

use strict;
use warnings;

use Getopt::Long qw( GetOptions );
use Data::Dumper;

my $VERSION = "v1.2";

# ------------------------------------------------
my $v = '';
my $debug = 0;

my $in_gtf = '';
my $out_gtf = '';
my $in_introns = '';
my $no_phase = 0;
my $reverse = 0;
# ------------------------------------------------

Usage() if ( @ARGV < 1 );
ParseCMD();
CheckBeforeRun();

my %h = ();  # introns are here in unique keys
ParseGFF( $in_introns, \%h );

my %transc = ();
ParseTransc( $in_gtf, \%transc );

my $count_no_intron = 0;
my $count_with_intron = 0;
my $count_out = 0;

open( my $OUT, ">" ,$out_gtf ) or die "error on open file $out_gtf: $!\n";
foreach my $key (keys %transc)
{
	if($debug)
	{
		print "# $key\n";
	}

	if ( HasIntron( $transc{$key} ) )
	{
		$count_with_intron += 1;

		if ( CompareIntrons( $transc{$key}, \%h ) )
		{
			if ( ! $reverse )
			{
				print $OUT $transc{$key};
				$count_out += 1;
			}
		}
		else
		{
			if ( $reverse )
			{
				print $OUT $transc{$key};
				$count_out += 1;
			}
		}
	}
	else
	{
		$count_no_intron += 1;
	}
}
close $OUT;

if($v)
{
	print "# genes_no_introns     $count_no_intron\n";
	print "# genes_with_introns   $count_with_intron\n";
	print "# genes_out            $count_out\n";
}

exit 0;

# ================= subs =========================
sub HasIntron
{
	my $txt = shift;

	my $result = 0;

	if ( $txt =~ /\t[Ii]ntron\t/ )
	{
		$result = 1;
	}

	return $result;
}
# ------------------------------------------------
sub CompareIntrons
{
	my $txt = shift;
	my $ref = shift;

	my $is_good = 0;
	my $intron_in = 0;
	my $introns_match = 0;

	my @arr = split( '\n', $txt );

	foreach my $line (@arr)
	{
		if ( $line =~ /\t[Ii]ntron\t/ )
		{
			if( $line =~ /^(\S+)\t\S+\t(\S+)\t(\d+)\t(\d+)\t\S+\t([-+.])\t(\S+)(\t(.*)|\s*)$/ )
			{
				my $id     = $1;
				my $type   = $2;
				my $start  = $3;
				my $end    = $4;
				my $strand = $5;
				my $ph     = $6;
				my $attr   = $7;

				$intron_in += 1;

				my $key = $id ."_". $start ."_". $end ."_". $strand;

				$key .= "_". $ph  if ( ! $no_phase );

				if ( exists $ref->{$key} )
				{
					$introns_match += 1;
				}
			}
			else { die "error, unexpected line format: $line"; }
		}
	}

	if ( $introns_match and ($intron_in == $introns_match ) )
	{
		$is_good = 1;
	}

	if ($debug)
	{
		print "# $intron_in $introns_match\n";
	}

	return $is_good;
}
# ------------------------------------------------
sub ParseTransc
{
	my ($name, $ref) = @_;

	open( my $IN, $name ) or die "error on open file $name: $!\n";
	while( my $line = <$IN> )
	{
		next if( $line =~ /^\s*#/ );
		next if( $line =~ /^\s*$/ );
		next if( $line !~ /\t/ );

		if( $line =~ /^(\S+)\t\S+\t(\S+)\t(\d+)\t(\d+)\t\S+\t([-+.])\t(\S+)(\t(.*)|\s*)$/ )
		{
			my $id     = $1;
			my $type   = $2;
			my $start  = $3;
			my $end    = $4;
			my $strand = $5;
			my $ph     = $6;
			my $attr   = $7;

			if ( $attr =~ /transcript_id \"(\S+)\";/ )
			{
				$ref->{$1} .= $line;
			}
			else { die "error, in parsing transcrip id:\n$line\n"; } 
		}
	}
	close $IN;

	if ($v)
	{
		print "# Transcripts in file $name: ". (scalar keys %$ref) ."\n";
	}
}
# ------------------------------------------------
sub ParseGFF
{
	my ($name, $ref) = @_;

	open( my $IN, $name ) or die "error on open file $name: $!\n";
	while( my $line = <$IN> )
	{
		next if( $line =~ /^\s*#/ );
		next if( $line =~ /^\s*$/ );
		next if( $line !~ /\t/ );
		
		if( $line =~ /^(\S+)\t\S+\t(\S+)\t(\d+)\t(\d+)\t\S+\t([-+.])\t(\S+)(\t(.*)|\s*)$/ )
		{
			my $id     = $1;
			my $type   = $2;
			my $start  = $3;
			my $end    = $4;
			my $strand = $5;
			my $ph     = $6;
			my $attr   = $7;

			$attr = '' if ( !defined $attr );
			$attr =~ s/^\t//;

			if ( $type =~ /^[Ii]ntron/ )
			{;}
			else {next;}

			if ( $v )
			{
				print "warning, strand is not defined: $line" if ($strand eq ".");
			}

			# start of unique key for most of the comparisons
			# 
			my $key = $id ."_". $start ."_". $end ."_". $strand;

			if ( ! $no_phase )
			{
				$key .= "_". $ph;
			}

			if ( exists $ref->{$key} and $v )
			{
				print "warning, more than one record with the same key was detected: $line";
			}

			$ref->{$key} .= $line;
		}
		else
		{
			print "warning, unexpected line format found in gff: $line\n";
		}
	}	
	close $IN;
	
	if ($v)
	{
		print "# Introns in file $name: ". (scalar keys %$ref) ."\n";
	}

#	print Dumper($ref);
}
# ------------------------------------------------
sub CheckBeforeRun
{
	die "error, file not found: option --in_gtf $in_gtf\n" if( ! -e $in_gtf );
	die "error, file not found: option --in_introns $in_introns\n" if( ! -e $in_introns );

	if ( $out_gtf )
	{
		die "error, output file name matches input file\n" if (( $out_gtf eq $in_gtf ) or ( $out_gtf eq $in_introns ));
	}
}
# ------------------------------------------------
sub ParseCMD
{
	my $opt_results = GetOptions
	(
		'in_gtf=s'     => \$in_gtf,
		'out_gtf=s'    => \$out_gtf,
		'in_introns=s' => \$in_introns,
		'no_phase'  => \$no_phase,
		'reverse'   => \$reverse,
		'verbose'   => \$v,
		'debug'     => \$debug,
	);

	die "error on command line\n" if( !$opt_results );
	die "error, unexpected argument found on command line\n" if( @ARGV > 0 );

	$v = 1 if $debug;
}
# ------------------------------------------------
sub Usage
{
	print qq(
Usage:
$0  --in_gtf [name]  --in_introns [name]  --out_gtf [name]

Select genes based on introns
 
Optional:
   --no_phase   exclude intron phase from comparison
   --reverse    reverse selection criteria: print only multi exons genes
                which have at least one mismatch with input introns
General:
   --verbose
   --debug

Version: $VERSION

);
	exit 1;
}
# ------------------------------------------------
 
