#!/usr/bin/perl
# --------------------------------------------
# Alex Lomsadze
# GaTech
# update 2019
#
# Enrich GFF3 file:
#     * add introns, if missing
#     * separate intorns in UTRs from introns between CDS
#     * add start and stop codons
#     * add exon type label: sinle, inititial, internal or terminal
#     * add splice site dinucleotides:  gt_ag, gc_ag, etc to intons
#     * add start and stop triplets to start/stop codon records 
#     * add complete or partial gene label to gene
# --------------------------------------------

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $VERSION = "v1_2019";

# --------------------------------------------
my $in = '';
my $out = '';
my $v = '';
my $debug = '';
my $warnings = '';
my $cds_gene_only = '';
my $min_intron = 10;
# --------------------------------------------
Usage() if ( @ARGV < 1 );
ParseCMD();
CheckBeforeRun();
# --------------------------------------------

my @record = ();

my $count_CDS_genes = 0;

open( my $IN, $in ) or die "error om open file $in: $!\n";
open( my $OUT, ">", $out ) or die "error om open file $out: $!\n";

while(my $line = <$IN>)
{
	next if( $line =~ /^\s*$/ );

	if ( $line =~ /^#[^#]/ )
	{
		print $OUT $line;
		next;
	}

	if ( $line =~ /^##[^#]/ )
	{
		print $OUT $line;
		next;
	}

	if ( $line =~ /^###\s*$/ )
	{
		if ( IsCDSgene(\@record) )
		{
			my @new_record = EnrichRecord(\@record);
			$count_CDS_genes += 1;

			PrintRecord(\@new_record);
		}
		else
		{
			PrintRecord(\@record) if ( !$cds_gene_only );
		}

		@record = ();

		print $OUT $line;
		next;
	}

	AddLineToRecord( $line, \@record );
}

if ( @record )
{
	if ( IsCDSgene(\@record) )
	{
		my @new_record = EnrichRecord(\@record);
		$count_CDS_genes += 1;

		PrintRecord(\@new_record);
	}
	else
	{
		PrintRecord(\@record) if ( !$cds_gene_only );
	}

	@record = ();

	print $OUT, "###\n";
}

close $OUT;
close $IN;

if($v)
{
	print "# CDS genes in: $count_CDS_genes\n";
}

print "done\n";

# --------------------------------------------
sub EnrichRecord
{
	my $ref = shift;

	my @new_r = ();

	my $gene_id = GetGeneID($ref, \@new_r);
	my %mrna = GetMrnaIDs($ref);
	
	SplitByMRNA( $ref, \%mrna );

	foreach my $key ( keys %mrna )
	{
		push @new_r, [ @{$mrna{$key}[0]} ];

		my @cds = SelectCDSsorted( $mrna{$key} );
		my @introns = CreateCdsIntrons( \@cds );
		my @start_codon = CreateStartCodon( \@cds );
		my @stop_codon = CreateStopCodon( \@cds );

		foreach my $entry (@cds)
		{
			$entry->[8] = "Parent=". $key;
		}

		foreach my $entry (@introns)
		{
			$entry->[8] = "Parent=". $key;
		}

		foreach my $entry (@start_codon)
		{
			$entry->[8] = "Parent=". $key;
		}

		foreach my $entry (@stop_codon)
		{
			$entry->[8] = "Parent=". $key;
		}

		push @new_r, @cds;
		push @new_r, @introns;
		push @new_r, @start_codon;
		push @new_r, @stop_codon;
	}

	return @new_r;
}
# --------------------------------------------
sub CreateStopCodon
{
	my $ref = shift;

	my @arr = ();
	my @current = @{$ref->[0]};

	if ( $current[6] eq '+' )
	{
		my $idx = (scalar @{$ref}) - 1;

		@current = @{$ref->[$idx]};

		if ( $current[4] - $current[3] + 1 >= 3 )
		{
			$current[2] = "stop_codon";
			$current[3] = $current[4] - 3 + 1;
			$current[7] = 0;

			push @arr, [ @current ];
		}
		else
		{
			if ( $current[4] - $current[3] + 1 == 2 )
			{
				$current[2] = "stop_codon";
				$current[3] = $current[4] - 1;
				$current[7] = 0;

				push @arr, [ @current ];

				@current = @{$ref->[$idx - 1]};

				$current[2] = "stop_codon";
				$current[3] = $current[4];
				$current[7] = 1;

				push @arr, [ @current ];
                        }
			elsif ( $current[4] - $current[3] + 1 == 1 )
			{
				$current[2] = "stop_codon";
				$current[3] = $current[4];
				$current[7] = 0;

				push @arr, [ @current ];

				@current = @{$ref->[$idx - 1]};

				$current[2] = "stop_codon";
				$current[3] = $current[4] - 1;
				$current[7] = 2;

				push @arr, [ @current ];
			}

			if ($warnings)
			{
				print "warning, split stop codon detected:\n";
				print Dumper(\@arr);
			}
		}
	}
	elsif ( $current[6] eq '-' )
	{
		if ( $current[4] - $current[3] + 1 >= 3 )
		{
			$current[2] = "stop_codon";
			$current[4] = $current[3] + 3 - 1;
			$current[7] = 0;

			push @arr, [ @current ];
		}
		else
		{
			if ( $current[4] - $current[3] + 1 == 2 )
			{
				$current[2] = "stop_codon";
				$current[4] = $current[3] + 1;
				$current[7] = 2;

				push @arr, [ @current ];

				@current = @{$ref->[1]};

				$current[2] = "stop_codon";
				$current[4] = $current[3];
				$current[7] = 0;

				push @arr, [ @current ];
			}
			elsif ( $current[4] - $current[3] + 1 == 1 )
			{
				$current[2] = "stop_codon";
				$current[4] = $current[3];
				$current[7] = 1;

				push @arr, [ @current ];

				@current = @{$ref->[1]};

				$current[2] = "stop_codon";
				$current[4] = $current[3] + 1;
				$current[7] = 0;

				push @arr, [ @current ];
			}

			if ($warnings)
			{
				print "warning, split stop codon detected:\n";
				print Dumper(\@arr);
			}
		}
	}

	return @arr;
}
# --------------------------------------------
sub CreateStartCodon
{
	my $ref = shift;

	my @arr = ();
	my @current = @{$ref->[0]};

	if ( $current[6] eq '+' )
	{
		if ( $current[4] - $current[3] + 1 >= 3 )
		{
			$current[2] = "start_codon";
			$current[4] = $current[3] + 3 - 1;
			$current[7] = 0;

			push @arr, [ @current ];
		}
		else
		{
			if ( $current[4] - $current[3] + 1 == 2 )
			{
				$current[2] = "start_codon";
				$current[4] = $current[3] + 1;
				$current[7] = 0;

				push @arr, [ @current ];

				@current = @{$ref->[1]};

				$current[2] = "start_codon";
				$current[4] = $current[3];
				$current[7] = 1;

				push @arr, [ @current ];
			}
			elsif ( $current[4] - $current[3] + 1 == 1 )
			{
				$current[2] = "start_codon";
				$current[4] = $current[3];
				$current[7] = 0;

				push @arr, [ @current ];

				@current = @{$ref->[1]};

				$current[2] = "start_codon";
				$current[4] = $current[3] + 1;
				$current[7] = 2;

				push @arr, [ @current ];
			}

			if ($warnings)
			{
				print "warning, split start codon detected:\n";
				print Dumper(\@arr);
			}
		}
	}
	elsif ( $current[6] eq '-' )
	{
		my $idx = (scalar @{$ref}) - 1;

		@current = @{$ref->[$idx]};

		if ( $current[4] - $current[3] + 1 >= 3 )
		{
			$current[2] = "start_codon";
			$current[3] = $current[4] - 3 + 1;
			$current[7] = 0;

			push @arr, [ @current ];
		}
		else
		{
			if ( $current[4] - $current[3] + 1 == 2 )
			{
				$current[2] = "start_codon";
				$current[3] = $current[4] - 1;
				$current[7] = 0;

				push @arr, [ @current ];

				@current = @{$ref->[$idx - 1]};

				$current[2] = "start_codon";
				$current[3] = $current[4];
				$current[7] = 1;

				push @arr, [ @current ];
                        }
                        elsif ( $current[4] - $current[3] + 1 == 1 )
                        {
                                $current[2] = "start_codon";
                                $current[3] = $current[4];
                                $current[7] = 0;

                                push @arr, [ @current ];

                                @current = @{$ref->[$idx - 1]};

                                $current[2] = "start_codon";
                                $current[3] = $current[4] - 1;
                                $current[7] = 2;

                                push @arr, [ @current ];
                        }

			if ($warnings)
			{
				print "warning, split start codon detected:\n";
				print Dumper(\@arr);
			}
		}
	} 

	return @arr;
}
# --------------------------------------------
sub PrintRecord
{
	my $ref = shift;

	foreach my $entry ( @{$ref} )
	{
		print $OUT   $entry->[0] ."\t". $entry->[1] ."\t". $entry->[2] ."\t". $entry->[3] ."\t". $entry->[4] ."\t". $entry->[5] ."\t". $entry->[6] ."\t". $entry->[7];
		if ( defined $entry->[8] )
		{
			print $OUT "\t". $entry->[8];
		}
		print $OUT  "\n";
	}
}
# --------------------------------------------
sub CreateCdsIntrons
{
	my $ref = shift;

	my $size = scalar @{$ref};

	my @arr = ();

	return @arr if ($size == 1);

	my $i = 0;
	my $j = 1;

	while( $j < $size)
	{
		if ( $ref->[$i][6] eq $ref->[$j][6] )
		{
			my @current = @{$ref->[$i]};
			$current[2] = "intron";
			$current[3] = $ref->[$i][4] + 1;
			$current[4] = $ref->[$j][3] - 1;

			if ( $ref->[$i][6] eq "+" )
			{
				$current[7] = 3 - $ref->[$i][7];
				$current[7] = 0 if ( $current[7] == 0 );
			}
			elsif ( $ref->[$i][6] eq "-" )
			{
				$current[7] = $ref->[$i][7];
			}

			if ( $current[4] - $current[3] + 1 < $min_intron )
			{
				$current[2] = "gap";

				if ($warnings)
				{
					print "warning, distance between CDS-CDS is below $min_intron: gap was introduced\n";
					print Dumper(\@current);
				}
			}

			push @arr, [@current];
		}

		$i += 1;
		$j += 1;
	}

	return @arr;
}
# --------------------------------------------
sub SelectCDSsorted
{
	my $ref = shift;

	my @arr = ();

	foreach my $entry (@{$ref})
	{
		if ( $entry->[2] eq "CDS" )
		{
			if ( $entry->[6] !~ /^[+-]$/ )
				{ die "error, CDS strand value is missing: $entry->[6]\n"; } 

			push @arr, [ @{$entry} ];
		}
	}

	@arr = sort{ $a->[0] cmp $b->[0] || $a->[3] <=> $b->[3] } @arr;

	if ( $warnings )
	{
		my $i = 0;
		my $j = 1;
		my $size = scalar @arr;

		while( $j < $size )
		{
			if ( $arr[$i][4] >= $arr[$j][3] )
			{
				print "warning, two CDS from the same mRNA overlap: $arr[$i][3] .. $arr[$i][4]  $arr[$j][3] .. $arr[$j][4]\n"; 
			}

			$i += 1;
			$j += 1;
		}

		my $strand = $arr[0][6];
		my $seq = $arr[0][0];

		foreach my $current (@arr)
		{
			if ( $strand ne $current->[6] )
			{
				print "warning, two strands detected in one mRNA: $strand  $current->[6]\n";
				last;
			}	
		}

		foreach my $current (@arr)
		{
			if ( $seq ne $current->[0] )
			{
				print "warning, two seqid detected in one mRNA: $strand  $current->[0]\n";
				last;
			}
		}
	}

	return @arr;
}
# --------------------------------------------
sub SplitByMRNA
{
	my $ref = shift;
	my $h_mrna = shift;

	my %other_transcripts = ();

	foreach my $entry (@{$ref})
	{
		next if ( $entry->[2] eq "mRNA" );
		next if ( $entry->[2] eq "gene" );

		if ( $entry->[2] eq "nc_primary_transcript")
		{
			if ( $entry->[8] =~ /ID=(\S+?);/ )
			{
				$other_transcripts{$1} = 1;
				next;
			}
		}

		if ( $entry->[8] =~ /Parent=(\S+?);/ or $entry->[8] =~ /Parent=(\S+)/)
		{
			my @pset = split( ',', $1 );

			foreach my $value (@pset)
			{
				if ( exists $h_mrna->{$value} )
				{
					push @{$h_mrna->{$value}}, [ @{$entry} ];
				}
				else
				{
					if ( exists $other_transcripts{$value} )
					{
						next;
					}
					else
						{ print Dumper($ref); die "error, feature not in mRNA:\n$value\n"; }
				}
			}
		}
		else
			{ die "error, line without the Parent in:\n$entry->[8]\n"; }
	}
}
# --------------------------------------------
sub GetMrnaIDs
{
	my $ref = shift;

	# one or many per record

	my %mRNA = ();

	foreach my $entry (@{$ref})
	{
		if ( $entry->[2] eq "mRNA" )
		{
			if ( $entry->[8] =~ /ID=(\S+?);/ )
			{
				my $ID = $1;

				if ( ! exists $mRNA{$ID} )
				{
					push @{$mRNA{$ID}}, [ @{$entry} ];
				}
				else
					{ die "error, mRNA ID duplication found: $ID\n"; }
			}
			else
				{ die "error, mRNA ID field not found in: $entry->[8]\n"; }
		}
	}

	if ( scalar (keys %mRNA ) == 0 )
	{
		print Dumper($ref);
		die "error, no mRNA in record\n";
	}

	return %mRNA;
}
# --------------------------------------------
sub GetGeneID
{
	my $ref = shift;
	my $arr_r = shift;

	# one "gene" per record

	my $gene = '';

	foreach my $entry (@{$ref})
	{
		if ( $entry->[2] eq "gene" )
		{
			if ( !$gene )
			{
				if ( $entry->[8] =~ /ID=(\S+?);/ ) 
				{
					$gene = $1;

					push @{$arr_r}, [@{$entry}];
				}
				else
					{ die "error, gene ID field not found in: $entry->[8]\n"; }
			}
			else
				{ die "error, gene entry duplication was detected in record: $gene\n"; }
		}
	}

	return $gene;
}
# --------------------------------------------
sub IsCDSgene
{
	my $ref = shift;

	my $gene_count = 0;
	my $cds_count = 0;

	foreach my $entry (@{$ref})
	{
		if ( $entry->[2] eq "CDS" )
		{
			$cds_count += 1;
		}
		elsif ( $entry->[2] eq "gene" )
		{
			$gene_count += 1;
		}
	}

	if ( $gene_count > 1 )
	{
		print Dumper($ref);
		die "error, gene label duplication was detected:\n";
	}

	if ( $gene_count == 1 and $cds_count > 0 )
	{
		return 1;
	}
	else
	{
		return 0;
	}
}
# --------------------------------------------
sub CheckForValidGFF
{
	my $ref = shift;
	my $str = shift;

	# seqid 1
	if ( $ref->[0] !~ /^\w+$/ )
		{ die "error, unexpect seqid found:\n$str\n"; }

	# start 4
	if ( $ref->[3] !~ /^\d+$/ )
		{ die "error, unexpect start found:\n$str\n"; }

	# end 5
	if ( $ref->[4] !~ /^\d+$/ )
		{ die "error, unexpect end found:\n$str\n"; }

	# start <= end
	if ( $ref->[3] > $ref->[4] )
		{ die "error, start is more than end:\n$str\n"; }

	if ( $ref->[6] !~ /^[+-.]$/ )
		{ die "error, wrong strand value:\n$str\n"; }
}
# --------------------------------------------
sub AddLineToRecord
{
	my $str = shift;
	my $ref = shift;

	chomp $str;

	my @arr = split( '\t', $str );

	my $size = @arr;

	if ( $size != 8 and $size != 9 )
		{ die "error, unexpected number of TABs found:\n$str\n"; }

	CheckForValidGFF( \@arr, $str ) if $debug;

	push @{$ref}, [@arr];
}
# ------------------------------------------------
sub CheckBeforeRun
{
	die "error, file not found: option --in $in\n" if( ! -e $in );
	die "error, output file name matches input file: $in $out\n" if ( $out eq $in );
}
# ------------------------------------------------
sub ParseCMD
{
	my $opt_results = GetOptions
	(
		'in=s'    => \$in,
		'out=s'   => \$out,
		'verbose'  => \$v,
		'debug'    => \$debug,
		'warnings' => \$warnings,
		'cds_gene_only' => \$cds_gene_only,
		'min_intron=i'  => \$min_intron,
        );

	die "error on command line\n" if( !$opt_results );
	die "error, unexpected argument found on command line\n" if( @ARGV > 0 );

	$v = 1 if $debug;
	$warnings = 1 if $debug;
}
# ------------------------------------------------
sub Usage
{
        print qq(
Usage:
$0  --in [name]  --out [name]

  This program takes as input 'nice' GFF3 formatted genome annotation file
  and enriches annotation by adding ...

  --cds_gene_only    output only protein coding genes
  --min_intron [$min_intron]  minimum length of intron to calculate from CDS-CDS

General:
  --verbose
  --debug
  --warnings

Version: $VERSION

);
	exit 1;
}
# ------------------------------------------------

