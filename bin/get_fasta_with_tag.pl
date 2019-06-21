#!/usr/bin/perl
# --------------------------------------------
# Alex Lomsadze
# GaTech
# --------------------------------------------

use strict;
use warnings;

use Getopt::Long;

my $Version = "1.7";

# --------------------------------------------

if ( $#ARGV == -1 ) { print PrintUsage($Version); exit 1; }

my $in = '';
my $out = '';
my $tag = '';
my $list = '';
my $swap = '';
my $rev = '';
my $v = '';

my $result = GetOptions
(
  'in=s'   => \$in,
  'out=s'  => \$out,
  'tag=s'  => \$tag,
  'list=s' => \$list,
  'rev'    => \$rev,
  'swap'   => \$swap,
  'verbose'=> \$v,
);

if ( !$result ) { print "error on cmd\n"; print PrintUsage($Version); exit 1; }
if ( !$in or !$out ) { die "error on cmd: out or in is missing\n"; }
if ( $in eq $out ) { die "error on cmd: out equals in\n"; }
if ( (!$tag and !$list) ) { die "error on cmd: tag or list is missing \n"; }
if ( $tag and $swap ) { print "warning, swap option is ignored with tag option\n"; }
if ( $rev and $swap ) { print "warning, swap option is ignored with rev option\n"; }
if ( $tag and $tag =~ /\s/ ) { die "error, white space is not allowed in tag: $tag\n"; }

# --------------------------------------------

my %h;
LoadList( $list, \%h ) if( $list );

# --------------------------------------------

my $count_in = 0;
my $count_out = 0;

$/ = ">";

open( my $IN, $in ) or die "Can't open $in: $!\n";
open( OUT, ">", $out ) or die "Can't open $out: $!\n";

while ( my $record = <$IN> )
{
	chomp $record;
	next if ( $record =~ /^\s*$/ );

	++$count_in;

	if ( $list )
	{
		if ( $record =~/^\s*(\S+)\s+/ )
		{
			my $id = $1;
			
			if ( exists( $h{$id} ) and !$rev )
			{
				print OUT ">";

				if ( $swap )
				{
					if ( $record =~ /^\s*\S+[ \t]\S+/ )
					{
						$record =~ s/^\s*\S+\s*.*?\n//;
					}
					else
					{
						$record =~ s/^\s*\S+\s*?\n//;
					}
					print OUT $h{$id} ."\n";
				}

				print OUT "$record";

				++$count_out;

				print "found ". $id ."\n" if $v;

				$h{$id} = '';
			}
			elsif( !exists( $h{$id} ) and $rev )
			{
				print OUT ">";
				print OUT "$record";

				++$count_out;

				print "not found ". $id ."\n" if $v;
			}
		}
		else
		{
			print "error, unexpected\n";
			exit 1;
		}
	}
	elsif ( $tag )
	{
		if ( !$rev and ($record =~ /^\s*$tag/) )
		{
			print OUT ">";
			print OUT "$record";

			++$count_out;
		}
		elsif ( $rev and ($record !~ /^\s*$tag/) ) 
		{
			print OUT ">";
			print OUT "$record";

			++$count_out;
		}
	}
	else
	{
		print "error, unexpected\n";
		exit 1;
	}
}

print OUT "\n";

close $IN;
close OUT;

CheckAllFound( \%h ) if $list;

print "$count_in in $in file\n";
print "$count_out in $out file\n";

# done

# --------------------------------------------
sub CheckAllFound
{
	my $ref = shift;

	for my $key (keys %$ref )
	{
		if ( $ref->{$key} )
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
	
	if (!$name)
	{
		print "error, file name is empty\n";
		exit 1;
	} 
	
	open( my $F, $name ) or die "Can't open $name: $!\n";
	while( my $line = <$F> )
	{
		next if ( $line =~ /^\s*$/ );
		next if ( $line =~ /^#/ );
		
		if ( $line =~ /^\s*(\S+)\s*$/ )
		{
			if ( exists $ref->{$1} )
			{
				print "error, duplicated entry found in the list: $name\n";
				exit 1;
			}
 
			$ref->{$1} = $1;
		}
		elsif ( $line =~ /^\s*(\S+)\s+(\S+)\s*/ )
		{
			if ( exists $ref->{$1} )
			{
				print "error, duplicated entry found in the list: $name\n";
				exit 1;
			}

			$ref->{$1} = $2;
		}
		else
		{
			print "warning, unexpected format was found: $line";
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
	my $label = shift;
	
	my $txt = "Usage: $0  --in <file name>  --out <file name>  --tag <string>  OR  --list <file name>

  This program takes as input FASTA formatted sequence and
  outputs FASTA records with definition lines that
     matches the tag
  or
     are in the list

  --tag  <string>, match this:
            '> string*'
  --list <file name>, match this:
            '> string '
  --swap swap fasta defline to one matching the list
  --rev  reverse selection, get records not matching the tag

version $label;
";
	return $txt;
}
# --------------------------------------------

