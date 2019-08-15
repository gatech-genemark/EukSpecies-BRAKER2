#!/usr/bin/perl
# --------------------------------------------
# Alex Lomsadze
# GaTech
# January 2019
# --------------------------------------------

use strict;
use warnings;
use Getopt::Long;

# --------------------------------------------

if ( $#ARGV == -1 ) { print PrintUsage(); exit 1; }

my $in = '';
my $out = '';
my $list = '';
my $v = '';
my $col = 1;
my $swap = 0;

my $result = GetOptions
(
  'in=s'   => \$in,
  'out=s'  => \$out,
  'list=s' => \$list,
  'verbose'=> \$v,
  'col=i'  => \$col,
  'swap'   => \$swap,
);

if ( !$result ) { print "error on cmd"; print PrintUsage(); exit 1; }
if ( !$in or !$out ) { die "error on cmd: out or in is missing\n"; }
if (( $in eq $out ) or ( $list eq $out)) { die "error on cmd: out equals in or list\n"; }
if ( !$list ) { die "error on cmd: list is missing \n"; }

# --------------------------------------------

my %h;
LoadList( $list, \%h );

# --------------------------------------------

my $count_in = 0;
my $count_out = 0;

my %allowed =
(
	"gene" => 1,
	"mRNA" => 1,
	"CDS" => 1,
	"exon" => 1,

	"Repeat" => 1,
	"similarity" => 1,
	"match" => 1,

#	"five_prime_UTR" => 1,
#	"three_prime_UTR" => 1,
	"intron" => 1,
	"transposable_element" => 1,

	"pseudogenic_transcript" => 1,
	"piRNA" => 1,
	"lincRNA" => 1,
	"miRNA_primary_transcript" => 1,
	"nc_primary_transcript" => 1,
	"pseudogenic_rRNA" => 1,
	"pseudogenic_tRNA" => 1,
	"scRNA" => 1,

	"lnc_RNA" => 1,
	"antisense_lncRNA" => 1,
	"transcript_region" => 1,
	"antisense_RNA" => 1,
	"transposable_element_gene" => 1,

	"ncRNA" => 1,
	"snoRNA" => 1,
	"pseudogene" => 1,
	"snRNA" => 1,
	"tRNA" => 1,
	"rRNA" => 1,
	"pre_miRNA" => 1,
	"miRNA" => 1,

	"ncRNA_gene" => 1,
	"unconfirmed_transcript" => 1,
	"C_gene_segment" => 1,
	"D_gene_segment" => 1,
	"J_gene_segment" => 1,
	"V_gene_segment" => 1,

);

my %found;

open( my $IN, $in ) or die "Can't open $in: $!\n";
open( OUT, ">", $out ) or die "Can't open $out: $!\n";

while ( my $record = <$IN> )
{
	next if ( $record =~ /^\s*$/ );

	++$count_in;

	if ( $record =~ /^##gff-version\s+3/ )
	{
		print OUT $record;
		next;
	}

	if ( $record =~/^(\S+)\t\S+\t(\S+)\t/ )
	{
		my $id = $1;
		my $type = $2;

		next if ( ! exists $allowed{$type} );

		if ( exists( $h{$id} ) )
		{
			$record =~ s/^$id/$h{$id}/;

			print OUT $record;
			++$count_out;
			$found{$id} +=1 ;
		}
		else
		{
#			print "$id $type\n";
		}
	}
	else
	{
#		print $record;
	}
}

close $IN;
close OUT;

if ($list)
{
	CheckAllFound( \%found );
}

print "$count_in in $in file\n";
print "$count_out in $out file\n";

# --------------------------------------------
sub CheckAllFound
{
	my $ref = shift;

	for my $key (keys %$ref )
	{
		if ( ! exists $ref->{$key} )
		{
			print "warning, record was not found: $ref->{$key}\n";
		}
	}	
}
# --------------------------------------------
sub LoadList
{
	my $name = shift;
	my $ref = shift;
	
	if (!$name )
	{
		print "error, file name is empty\n";
		exit 1;
	} 
	
	open( my $F, $name ) or die "error on  open $name: $!\n";
	while( my $line = <$F> )
	{
		next if ( $line =~ /^\s*$/ );
		next if ( $line =~ /^#/ );
		
		if ( $line =~ /^\s*(\S+)\s+(\S+)\s*/ )
		{
			my $col_a = $1;
			my $col_b = $2;

			if ( $col == 2 )
			{
				$col_a = $2;
				$col_b = $1;
			}


			if ( exists $ref->{$col_a} )
			{
				print "error, duplicated entry found in the list: $name\n";
				exit 1;
			}

			if ( $swap )
			{
				$ref->{$col_a} = $col_b;
			}
			else
			{
				$ref->{$col_a} = $col_a;
			}
		}
		else
		{
			print "warning, unexpected format was found in $name: $line";
		}
	}
	close $F;

	if ( scalar keys %$ref < 1 )
	{
		print "error, list is empty: $name\n";
		exit 1;
	}

	if ( $v )
	{
		for my $key (keys %$ref)
		{
			print $key ." ". $ref->{$key} ."\n";
		}
	}
}
# --------------------------------------------
sub PrintUsage
{
	my $txt = "Usage: $0  --in <file name>  --out <file name>  --list <file name>

  This program takes as input GFF formatted anotation and
  outputs GFF records only for sequence ID's listed on the --list file.

  List file may have two columns
  --list <file name> : name_1 name_2
  This specifies which column to use as
  --col  <column>    :   1  or  2
  --swap             : swap one name to another
";
	return $txt;
}
# --------------------------------------------

