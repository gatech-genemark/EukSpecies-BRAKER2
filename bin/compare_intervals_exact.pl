#!/usr/bin/perl
# ==============================================================
# Alex Lomsadze, Tomas Bruna
# Georgia Institute of Technology, Atlanta, Georgia, US
# update 2019 
#
# This script compares intervals from GFF/GTF/GFF3 files and
# generates report suitable for drawing Venn diagram.
# Supported interval types:
# 	 "CDS"
# 	 "intron"
# 	 "star_codon"
# 	 "stop_codon"
# Comparison can be done between 2 or 3 files.
# ==============================================================

use strict;
use warnings;

use Getopt::Long qw( GetOptions );
use Data::Dumper;

# ------------------------------------------------
my $v = '';
my $debug = 0;

my $f1 = '';
my $f2 = '';
my $f3 = '';

my $compare_cds = 0;
my $compare_introns = 0;
my $compare_donors = 0;
my $compare_acceptors = 0;
my $compare_starts = 0;
my $compare_stops = 0;
my $compare_genes = 0;
my $compare_transcripts = 0;

my $no_phase = '';

my $out_file = '';
my $original = 0;

my $shared12 = '';
my $shared13 = '';
my $shared23 = '';

my $shared123 = '';

my $unique1 = '';
my $unique2 = '';
my $unique3 = '';
# ------------------------------------------------

Usage() if ( @ARGV < 1 );
ParseCMD();
CheckBeforeRun();

my %h1;
my %h2;
my %h3;

my %shared_12;
my %shared_13;
my %shared_23;

my %shared_123;

my %unique_1;
my %unique_2;
my %unique_3;

ParseGFF( $f1, \%h1 );
ParseGFF( $f2, \%h2 );
ParseGFF( $f3, \%h3 ) if $f3;

if( !$f3 )
{
	Compare2();
}
else
{
	Compare3();
}

exit 0;

# ================= subs =========================
sub Compare3
{
	my $c1 = scalar keys %h1;
	my $c2 = scalar keys %h2;
	my $c3 = scalar keys %h3;
	
	my $match = 0;
	
	my $c1_c2 = 0;
	my $c1_c3 = 0;
	my $c2_c3 = 0;
	
	my $c1_uniq = 0;
	my $c2_uniq = 0;
	my $c3_uniq = 0;
	

	foreach my $key ( keys %h1 )
	{
		if ( exists $h2{$key} and exists $h3{$key} )
		{
			$match += 1;
			$shared_123{$key} += 1 if $out_file;
		}
		elsif ( exists $h2{$key} and  ! exists $h3{$key} )
		{
			$c1_c2 += 1;
			$shared_12{$key} += 1 if $out_file;
		}
		elsif (! exists $h2{$key} and  exists $h3{$key} )
		{
			$c1_c3 += 1;
			$shared_13{$key} += 1 if $out_file;
		}
		else
		{
			$unique_1{$key} += 1 if $out_file;
		}
	}

	foreach my $key ( keys %h2 )
	{
		if ( ! exists $h1{$key} and exists $h3{$key} )
		{
			$c2_c3 += 1;
			$shared_23{$key} += 1 if $out_file;
		}
		
		if ($out_file)
		{
			if ( ! exists $h1{$key} and ! exists $h3{$key} )
			{
				$unique_2{$key} += 1 if $out_file;
			}
		}
	}

	if ( $out_file )
	{
		foreach my $key ( keys %h3 )
                {
			if ( ! exists $h1{$key} and ! exists $h2{$key} )
			{
				$unique_3{ $key } += 1;
			}
		}
	}

	$c1_uniq = $c1 - $match - $c1_c2 - $c1_c3;
        $c2_uniq = $c2 - $match - $c1_c2 - $c2_c3;
        $c3_uniq = $c3 - $match - $c1_c3 - $c2_c3;

	die "error, no intervals in file: $f1\n" if (!$c1);
	die "error, no intervals in file: $f2\n" if (!$c2);
	die "error, no intervals in file: $f3\n" if (!$c3);

	print "\n";
	print "#in\tin_all\tin_1\tin_2\tin_3\tunique\t\%match_all\tCDS\n" if ($v and !$compare_introns );
	print "#in\tin_all\tin_1\tin_2\tin_3\tunique\t\%match_all\tIntron\n" if ($v and $compare_introns );
	print $c1 ."\t". $match ."\t". "0"    ."\t". $c1_c2 ."\t". $c1_c3 ."\t". $c1_uniq ."\t". sprintf( "%.2f", 100.0*$match/$c1 ) ."\t". $f1 ."\n";
	print $c2 ."\t". $match ."\t". $c1_c2 ."\t". "0"    ."\t". $c2_c3 ."\t". $c2_uniq ."\t". sprintf( "%.2f", 100.0*$match/$c2 ) ."\t". $f2 ."\n";
	print $c3 ."\t". $match ."\t". $c1_c3 ."\t". $c2_c3 ."\t". "0"    ."\t". $c3_uniq ."\t". sprintf( "%.2f", 100.0*$match/$c3 ) ."\t". $f3 ."\n";
	print "\n";

	if ($out_file)
	{
		PrintKeys( \%unique_1,  "unique 1",   $out_file ) if $unique1;
		PrintKeys( \%unique_2,  "unique 2",   $out_file ) if $unique2;
		PrintKeys( \%unique_3,  "unique 3",   $out_file ) if $unique3;
		PrintKeys( \%shared_12, "shared 1-2", $out_file ) if $shared12;
		PrintKeys( \%shared_13, "shared 1-3", $out_file ) if $shared13;
		PrintKeys( \%shared_23, "shared 2-3", $out_file ) if $shared23;
		PrintKeys( \%shared_123, "shared 1-2-3", $out_file ) if $shared123;

		if ($debug)
		{
			TestCount( $match, \%shared_123 );
			TestCount( $c1_c2, \%shared_12 );
			TestCount( $c1_c3, \%shared_13 );
			TestCount( $c2_c3, \%shared_23 );
			TestCount( $c1_uniq, \%unique_1 );
			TestCount( $c2_uniq, \%unique_2 );
			TestCount( $c3_uniq, \%unique_3 );
			TestCount( $c1, \%shared_123, \%shared_12, \%shared_13, \%unique_1 );
			TestCount( $c2, \%shared_123, \%shared_12, \%shared_23, \%unique_2 );
			TestCount( $c3, \%shared_123, \%shared_13, \%shared_23, \%unique_3 );
        	}
	}
}
# ------------------------------------------------
sub PrintKeys
{
	my $ref = shift;
	my $label = shift;
	my $name = shift;

	open( my $OUT, ">", $name ) or die;

	print $OUT "# $label\n";

	foreach my $key (sort keys %{$ref})
	{
		if ( ! $original )
		{
			print $OUT "$key\n";
		}
		else
		{
			 print $OUT  PrintOriginal( $key, $original );
		}
	}

	close $OUT;
}
# ------------------------------------------------
sub PrintOriginal
{
	my $key = shift;
	my $id = shift;

	if ( $id == 1 )
	{
		return $h1{$key};
	}
	elsif ( $id == 2 )
	{
		return $h2{$key};
	}
	elsif ( $id == 3 )
	{
		return $h3{$key};
	}
	else
		{ die "error, unexpected value found in PrintOriginal\n"; }
}
# ------------------------------------------------
sub Compare2
{
	my $total_c1 = scalar keys %h1;
	my $total_c2 = scalar keys %h2;
	my $match = 0;
	my $c1_uniq = 0;
	my $c2_uniq = 0;
	
	foreach my $key ( keys %h1 )
	{
		if ( exists $h2{$key} )
		{
			$match += 1;
			$shared_12{$key} += 1 if $out_file;
		}
		else
		{
			$unique_1{ $key } += 1 if $out_file
		}
	}

	if ( $out_file )
	{
		foreach my $key ( keys %h2 )
		{
			if ( ! exists $h1{$key} )
			{
				$unique_2{ $key } += 1;
			}
		}
	}

	$c1_uniq = $total_c1 - $match;
	$c2_uniq = $total_c2 - $match;

	die "error, no intervals in file: $f1\n" if (!$total_c1);
	die "error, no intervals in file: $f2\n" if (!$total_c2);

	print "\n";
	print "#in\tmatch\tunique\t\%match\tCDS\n" if ($v and $compare_cds );
	print "#in\tmatch\tunique\t\%\tIntron\n" if ($v and $compare_introns );
	print "#in\tmatch\tunique\t\%\tStarts\n" if ($v and $compare_starts );
	print "#in\tmatch\tunique\t\%\tStops\n" if ($v and $compare_stops );
	print "#in\tmatch\tunique\t\%\tDonors\n" if ($v and $compare_donors );
	print "#in\tmatch\tunique\t\%\tAcceptors\n" if ($v and $compare_acceptors );
	print "\n";
	print  $total_c1 ."\t". $match ."\t". $c1_uniq ."\t". sprintf( "%.2f", 100.0*$match/$total_c1 ) ."\t". $f1 ."\n";
	print  $total_c2 ."\t". $match ."\t". $c2_uniq ."\t". sprintf( "%.2f", 100.0*$match/$total_c2 ) ."\t". $f2 ."\n";
	print "\n";

	if ( $out_file )
	{
		PrintKeys( \%unique_1,  "unique 1",   $out_file ) if $unique1;
		PrintKeys( \%unique_2,  "unique 2",   $out_file ) if $unique2;
		PrintKeys( \%shared_12, "shared 1-2", $out_file ) if $shared12;

		if ($debug)
		{
			TestCount( $match, \%shared_12 );
			TestCount( $c1_uniq, \%unique_1 );
			TestCount( $c2_uniq, \%unique_2 );
			TestCount( $total_c1, \%shared_12, \%unique_1 );
			TestCount( $total_c2, \%shared_12, \%unique_2 );
		}
	}
}
# ------------------------------------------------
sub TestCount
{
	my $count = shift;
	my $ref1 = shift;
	my $ref2 = shift;
	my $ref3 = shift;
	my $ref4 = shift;

	my $size = scalar (keys %{$ref1});

	$size += scalar (keys %{$ref2}) if ($ref2);
	$size += scalar (keys %{$ref3}) if ($ref3);
	$size += scalar (keys %{$ref4}) if ($ref4);

	die "error in counting: $count $size\n" if ( $size != $count );
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

		next if ( $line !~ /\t/ );
		
		if( $line =~ /^(\S+)\t\S+\t(\S+)\t(\d+)\t(\d+)\t\S+\t([-+.])\t(\S+)\s*/ )
		{
			my $id     = $1;
			my $type   = $2;
			my $start  = $3;
			my $end    = $4;
			my $strand = $5;
			my $ph     = $6;

			if ( $compare_cds and ( $type eq "CDS") )
			{
				;
			}
			elsif ( $compare_introns and ( $type =~ /^[Ii]ntron/) )
			{
				;
			}
			elsif ( $compare_starts and ( $type =~ /^[Ss]tart_codon/) )
			{
				;
			}
			elsif ( $compare_stops and ( $type =~ /^[Ss]top_codon/) )
			{
				;
			}
			elsif ( $compare_donors and ( $type =~ /^[Ii]ntron/) )
			{
				;
			}
			elsif ( $compare_acceptors and ( $type =~ /^[Ii]ntron/) )
			{
				;
			}
			else
				{ next; }

			if ( $v )
			{
				print "warning, strand is not defined: $line" if ($strand eq ".");
			}

			my $key = $id ."_". $start ."_". $end ."_". $strand;
			
			if ( $compare_donors )
			{
				if ( $strand eq "+" )
				{
					$key = $id ."_". $start ."_". $strand;
				}
				elsif ( $strand eq "-" )
				{
					$key = $id ."_". $end   ."_". $strand;
				}
			}

			if ( $compare_acceptors )
			{
				if ( $strand eq "+" )
				{
					$key = $id ."_". $end   ."_". $strand;
				}
				elsif ( $strand eq "-" )
				{
					$key = $id ."_". $start ."_". $strand;
				}
			}

			if ( ! $no_phase )
			{
				$key .= "_". $ph;
			}

			if ( ! $original )
			{
				$ref->{$key} += 1;
			}
			else
			{
				if ( exists $ref->{$key} and $v )
				{
					print "warning, more than one record with the same key was detected: $line";
				}

				$ref->{$key} .= $line;
			}
		}
		else
		{
			print "warning, unxpected line format found in gff: $line\n";
		}
	}	
	close $IN;
	
	if ($v)
	{
		print "# CDS in file $name: ". (scalar keys %$ref) ."\n" if $compare_cds;
		print "# Introns in file $name: ". (scalar keys %$ref) ."\n" if $compare_introns;
		print "# Starts in file $name: ". (scalar keys %$ref) ."\n" if $compare_starts;
		print "# Stops in file $name: ". (scalar keys %$ref) ."\n" if $compare_stops;
		print "# Donors in file $name: ". (scalar keys %$ref) ."\n" if $compare_donors;
		print "# Acceptors in file $name: ". (scalar keys %$ref) ."\n" if $compare_acceptors;
	}
	
	print Dumper($ref) if $debug;
}
# ------------------------------------------------
sub CheckBeforeRun
{
	die "error, file is missing $0: option --f1\n" if( ! -e $f1 );
	die "error, file is missing $0: option --f2\n" if( ! -e $f2 );

	if ( $f3 )
	{
		die "error, file is missing $0: option --f3\n" if( ! -e $f3 );
	}

	if ( $out_file )
	{
		die "error, output file name matches input file\n" if (( $out_file eq $f1 ) or ( $out_file eq $f1 ));
		die "error, output file name matches input file\n" if ( $f3 and ( $out_file eq $f3 ));
	}

	if ($original)
	{
		if ( $original == 1 or $original == 2 )
			{;}
		elsif ( $f3 and $original == 3 )
			{;}
		else
			{ die "error, value of --original parameter is out of allowed range: $original\n"; }
	}
}
# ------------------------------------------------
sub ParseCMD
{
	my $opt_results = GetOptions
	(
		'f1=s'   => \$f1,
		'f2=s'   => \$f2,
		'f3=s'   => \$f3,
		'out=s'  => \$out_file,
		'cds'       => \$compare_cds,
		'introns'   => \$compare_introns,
		'don'       => \$compare_donors,
		'acc'       => \$compare_acceptors,
		'starts'    => \$compare_starts,
		'stops'     => \$compare_stops,
		'no_phase'  => \$no_phase,
		'verbose'   => \$v,
		'debug'     => \$debug,
		'shared123' => \$shared123,
		'shared12'  => \$shared12,
		'shared13'  => \$shared13,
		'shared23'  => \$shared23,
		'unique1'   => \$unique1,
		'unique2'   => \$unique2,
		'unique3'   => \$unique3,
		'original=i' => \$original,
	);

	die "error on command line\n" if( !$opt_results );
	die "error, unexpected argument found on command line\n" if( @ARGV > 0 );

	my $count = 0;
	$count += 1 if $compare_introns;
	$count += 1 if $compare_donors;
	$count += 1 if $compare_acceptors;
	$count += 1 if $compare_starts;
	$count += 1 if $compare_stops;
	$count += 1 if $compare_cds;

	if ($count == 0 )
	{
		$compare_cds = 1;
	}
	elsif ($count > 1)
	{
		die "erros, more thean one comparision type was specifed on command line\n";
	}

	$v = 1 if $debug;
}
# ------------------------------------------------
sub Usage
{
	print qq(
Usage:
$0  --f1 [name]  --f2 [name]

Compare intervals from input files:
 
   --f1  [name]  file in GFF or GTF format
   --f2  [name]  file in GFF or GTF format

Optional:
   --f3  [name]  file in GFF or GTF format

Default comparision is done for 'CDS' type

   --cds         compare fields of 'CDS' type
   --introns     compare fields of 'intron' type
   --starts      compare fields of 'start_codon' type
   --stops       compare fields of 'stop_codon' type
   --don         compare donors using 'intron' type
   --acc         compare acceptors using 'intron' type

   --no_phase    ignore phase of record in comparision

Output subset of records into file based on the overlap status:

   --out [name]  output file name

   --shared12    output keys for shared by 1 and 2, not 3
   --shared13    output keys for shared by 1 and 3, not 2
   --shared23    output keys for shared by 2 and 3, not 1

   --shared123   output keys for shared by 1-2-3

   --unique1     output keys unique in 1
   --unique2     output keys unique in 2
   --unique3     output keys unique in 3

   --original [number]  1 or 2 of 3; print output as in original input files
                 1 corresponds to file with --f1 name, etc.
General:
   --verbose
   --bedug
);
	exit 1;
}
# ------------------------------------------------

