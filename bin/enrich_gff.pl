#!/usr/bin/perl
# --------------------------------------------
# Alex Lomsadze
# GaTech
# Last update 2019
#
# Enrich GFF3 file:
#     * add CDS introns
#     * add start and stop codons
#     * add exon type label: single, inititial, internal or terminal
#     * add splice site dinucleotides: gt_ag, gc_ag, etc to introns
#     * add start and stop triplets to start/stop codons
#     * compare calculated features with existing, if any
# To do:
#     * add label complete or partial to gene   -!
#     * separate intorns in UTRs and CDS   -!
#
# Input file must be in "nice" GFF3 format
# for example, file generated by "gt" code with --nice option
# --------------------------------------------

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;

my $VERSION = "v3_2019";

# --------------------------------------------
my $in = '';
my $out = '';
my $v = '';
my $debug = '';
my $warnings = '';
my $cds_gene_only = '';
my $min_intron = 10;
my $seq = '';
# --------------------------------------------
Usage() if ( @ARGV < 1 );
ParseCMD();
CheckBeforeRun();
# --------------------------------------------

# Download genome - optional
my %genome = ();
LoadGenome( $seq, \%genome ) if $seq;

# put all GFF line for one gene into this array
# genes are proceed one a time and array is reused
my @record = ();

my $count_CDS_genes = 0;

# helpers for warnings
my %mess = ();
$mess{"all_annot_introns"} = 0;
$mess{"match_annot"} = 0;
$mess{"mismatch_annot"} = 0;

print "# starting GFF3 parsing\n" if $v;

open( my $IN, $in ) or die "error on open file $in: $!\n";
open( my $OUT, ">", $out ) or die "error on open file $out: $!\n";
while(my $line = <$IN>)
{
	# skip empty lines
	next if( $line =~ /^\s*$/ );

	# print GFF3 support lines
	if ( $line =~ /^#[^#]/ or $line =~ /^##[^#]/ )
	{
		print $OUT $line;
		next;
	}

	# one full record was placed into array
	if ( $line =~ /^#{3}\s*$/ )
	{
		if ( IsCDSgene(\@record) )
		{
			my @new_record = EnrichRecord(\@record, $seq, \%mess);
			PrintRecord(\@new_record);

			print $OUT $line;
			$count_CDS_genes += 1;

			if ( $v and  $count_CDS_genes % 1000 == 0 )
			{
				print "# cds genes processed: $count_CDS_genes\n";
			}
		}
		elsif ( !$cds_gene_only )
		{
			PrintRecord(\@record);
			print $OUT $line;
		}

		@record = ();
	}
	else
	{
		AddLineToRecord( $line, \@record );
	}
}
close $OUT;
close $IN;

if ( @record )
	{ die "error, file didnot end with end-of GFF record line: \#\#\#\n"; }

if ($warnings)
{
	print "# all_annot_introns match_annot mismatch_annot\n";
	print "# $mess{'all_annot_introns'} $mess{'match_annot'} $mess{'mismatch_annot'}\n";
}

print "# CDS genes in: $count_CDS_genes\n" if $v;
print "# done\n" if $v;

# --------------------------------------------
sub LoadGenome
{
	my $name = shift;
	my $ref = shift;

	my $seqid = "";

	open( my $IN, $name ) or die "error on open file $name: $!\n";
	while(my $line = <$IN>)
	{
		next if ( $line =~ /^#/ );
		next if ( $line =~ /^\s*$/ );

		if ( $line =~ /^\s*>/ )
		{
			if ( $line =~ /^>\s*(\S+)\s+/ )
			{
				$seqid = $1;
				print "# seqid: $seqid\n" if $v;
			}
			else
				{ die "error, unexpected defline format found: $line\n"; }
		}
		else
		{
			if ( ! $seqid )
				{ die "error, seqid is missing\n"; }

			$line =~ s/\s|[0-9]//g;
			$line = uc($line);

			$ref->{$seqid} .= $line;
		}
	}
	close $IN;

	if ($v)
	{
		print "# genome sequence from file: $name\n";
		print "# number of sequences: ". (scalar keys %{$ref}) ."\n";
	}
}
# --------------------------------------------
sub EnrichRecord
{
	# ref on array of arrays of GFF3 values from one gene
	my $ref = shift;
	# parse sequence or not
	my $use_seq = shift;
	# for passing messages - warning information
	my $m = shift;

	my @new_rec = ();
	my @gene = ();

	# add gene line
	my $gene_id = GetGeneID($ref, \@gene);

	# prepare mrna lines
	my %mrna = GetMrnaIDs($ref);

	# sort other lines by mrna
	SplitByMRNA( $ref, \%mrna );

	AddCountToAttrInHashMRNA( \%mrna );

	push @new_rec, @gene;

	foreach my $key ( keys %mrna )
	{
		my @cds = SelectCDSsorted( $mrna{$key} );

		my @introns = ();
		my @start_codon = ();
		my @stop_codon = ();

		if ( @cds > 0 )
		{
			@introns = CreateCdsIntrons( \@cds );
			@start_codon = CreateStartCodon( \@cds );
			@stop_codon = CreateStopCodon( \@cds );

			foreach my $entry (@cds)          { $entry->[8] = "Parent=". $key .";"; }
			foreach my $entry (@introns)      { $entry->[8] = "Parent=". $key .";"; }
			foreach my $entry (@start_codon)  { $entry->[8] = "Parent=". $key .";"; }
			foreach my $entry (@stop_codon)   { $entry->[8] = "Parent=". $key .";"; }

			AddLabelsToCDS(\@cds);
			AddCountToAttr(\@cds);
			AddCountToAttr(\@introns);
			AddCountToAttrSemiReverse(\@start_codon, 1);
			AddCountToAttrSemiReverse(\@stop_codon, 0);

			if ( $use_seq )
			{
				AddSpliceSites( \@introns, \%genome );
				AddStartStopSites( \@start_codon, \%genome );
				AddStartStopSites( \@stop_codon, \%genome );
			}

			push @new_rec, [ @{$mrna{$key}[0]} ];
			push @new_rec, @cds;
			push @new_rec, @introns;
			push @new_rec, @start_codon;
			push @new_rec, @stop_codon;

			if ($warnings)
			{
				CompareIntrons( $ref, \@introns, $m );
			}
		}
	}

	return @new_rec;
}
# --------------------------------------------
sub CompareIntrons
{
	my $annot = shift;
	my $calc = shift;
	my $m = shift;

	my %h = ();

	# collect all annotated introns
	#
	foreach my $entry ( @{$annot} )
	{
		if ( $entry->[2] =~ /[Ii]ntron/ )
		{
			my $key = $entry->[0] ."_". $entry->[3] ."_". $entry->[4] ."_". $entry->[6];
			$h{$key} = 0;
		}
	}

	$m->{"all_annot_introns"} += scalar (keys %h);

	if ( $m->{"all_annot_introns"} > 0 )
	{
		foreach my $entry ( @{$calc} )
		{
			my $key = $entry->[0] ."_". $entry->[3] ."_". $entry->[4] ."_". $entry->[6];

			if ( exists $h{$key} )
			{
				$m->{"match_annot"} += 1;
			}
			else
			{
				$m->{"mismatch_annot"} += 1;
			}
		}
	}
}
# --------------------------------------------
sub AddStartStopSites
{
	my $ref = shift;
	my $h_genome = shift;

	foreach my $entry ( @{$ref} )
	{
		if ($warnings)
		{
			if ( $entry->[2] !~ /[Ss]tart_codon/ and $entry->[2] !~ /[Ss]top_codon/ )
				{ die "error, start/stop type is expected: $entry->[2]\n"; }
		}

		my $CODON = '';

		if ( $entry->[8] =~ /1_1/ )
		{
			$CODON = substr( $h_genome->{$entry->[0]}, $entry->[3] -1, 3 );
		}
		else
		{
			if ( $entry->[4] - $entry->[3] +1 == 1 )
			{
				$CODON = substr( $h_genome->{$entry->[0]}, $entry->[3] -1, 1 );
			}
			elsif ( $entry->[4] - $entry->[3] +1 == 2 )
			{
				$CODON = substr( $h_genome->{$entry->[0]}, $entry->[3] -1, 2 );
			}
			else
				{die;}
		}

		if ( $entry->[6] eq '+' )
		{
			;
		}
		elsif ( $entry->[6] eq '-' )
		{
			$CODON = RevComp($CODON);
		}
		else
			{ die "error, strand is required for sequence extraction:"; }

		$entry->[8] .= ("site_seq=". $CODON .";");
	}
}
# --------------------------------------------
sub AddSpliceSites
{
	my $ref = shift;
	my $h_genome = shift;

	foreach my $entry ( @{$ref} )
	{
		if ($warnings)
		{
			if ( $entry->[2] !~ /[Ii]ntron/ and $entry->[2] !~ /gap/ )
				{ die "error, intron type is expected: $entry->[2]\n"; }
		}

		my $DON = '';
		my $ACC = '';

		if ( $entry->[6] eq '+' )
		{
			$DON = substr( $h_genome->{$entry->[0]}, $entry->[3] -1, 2 );
			$ACC = substr( $h_genome->{$entry->[0]}, $entry->[4] -2, 2 );
		}
		elsif ( $entry->[6] eq '-' )
		{
			$ACC = substr( $h_genome->{$entry->[0]}, $entry->[3] -1, 2 );
			$DON = substr( $h_genome->{$entry->[0]}, $entry->[4] -2, 2 );

			$DON = RevComp($DON);
			$ACC = RevComp($ACC);
		}
		else
			{ die "error, strand is required for sequence extraction:"; }

		$entry->[8] .= ("site_seq=". $DON ."_". $ACC .";");
	}
}
# --------------------------------------------
sub RevComp
{
	my $s = shift;

	$s = reverse($s);

	$s =~ s/A/1/g;
	$s =~ s/T/A/g;
	$s =~ s/1/T/g;

	$s =~ s/C/2/g;
	$s =~ s/G/C/g;
	$s =~ s/2/G/g;

	return $s;
}
# --------------------------------------------
sub AddCountToAttr
{
	my $ref = shift;

	my $size = scalar @{$ref};

	return if ( $size == 0 );

	if ( $ref->[0][6] eq "+" )
	{
		for( my $i = 0; $i < $size; $i += 1 )
		{
			$ref->[$i][8] .= ";" if ( $ref->[$i][8] !~ m/;$/ );
			$ref->[$i][8] .= "count=". ($i+1) ."_". $size .";";
		}
	}
	elsif ( $ref->[0][6] eq "-" )
	{
		for( my $i = $size -1; $i >= 0; $i -= 1 )
		{
			$ref->[$i][8] .= ";" if ( $ref->[$i][8] !~ m/;$/ );
			$ref->[$i][8] .= "count=". ($i+1) ."_". $size .";";
		}
	}
}
# --------------------------------------------
sub AddCountToAttrSemiReverse
{
	my $ref = shift;
	my $is_start = shift;

	my $size = scalar @{$ref};

	return if ( $size == 0 );

	if ( $size == 1 )
	{
		$ref->[0][8] .= ";" if ( $ref->[0][8] !~ /;$/ );
		$ref->[0][8] .= "count=1_1;";
	}
	elsif ( $size == 2 )
	{
		$ref->[0][8] .= ";" if ( $ref->[0][8] !~ /;$/ );
		$ref->[1][8] .= ";" if ( $ref->[1][8] !~ /;$/ );

		if (( $is_start and ($ref->[0][6] eq "+" )) or ( !$is_start and ($ref->[0][6] eq "-" )))
		{
			$ref->[0][8] .= "count=1_2;";
			$ref->[1][8] .= "count=2_2;";
		}
		elsif (( !$is_start and ($ref->[0][6] eq "+" )) or ( $is_start and ($ref->[0][6] eq "-" )))
		{
			$ref->[0][8] .= "count=2_2;";
			$ref->[1][8] .= "count=1_2;";
		}
	}
}
# --------------------------------------------
sub AddCountToAttrInHashMRNA
{
	my $ref = shift;

	my $size = scalar ( keys %{$ref} );
	my $i = 0;

	foreach my $key ( keys %{$ref} )
	{
		$ref->{$key}[0][8] .= ";" if ( $ref->{$key}[0][8] !~ m/;$/ );
		$ref->{$key}[0][8] .= "count=". ($i+1) ."_". $size .";";

		$i += 1;
	}
}
# --------------------------------------------
sub AddLabelsToCDS
{
	my $ref = shift;

	my $size = scalar @{$ref};

	if ( $size == 1 )
	{
		$ref->[0][8] .= ";" if ( $ref->[0][8] !~ /;$/ );

		$ref->[0][8] .= "cds_type=Single;";
	}
	elsif ( $size > 1 )
	{
		if ( $ref->[0][6] eq "+" )
		{
			$ref->[0][8] .= ";" if ( $ref->[0][8] !~ /;$/ );
			$ref->[0][8] .= "cds_type=Initial;";

			$ref->[$size -1][8] .= ";" if ( $ref->[$size -1][8] !~ /;$/ );
			$ref->[$size -1][8] .= "cds_type=Terminal;";
		}
		elsif ( $ref->[0][6] eq "-" )
		{
			$ref->[$size -1][8] .= ";" if ( $ref->[$size -1][8] !~ /;$/ );
			$ref->[$size -1][8] .= "cds_type=Initial;";

			$ref->[0][8] .= ";" if ( $ref->[0][8] !~ /;$/ );
			$ref->[0][8]        .= "cds_type=Terminal;";
		}

		for( my $i = 1; $i < $size -1; $i += 1 )
		{
			$ref->[$i][8] .= ";" if ( $ref->[$i][8] !~ /;$/ );
			$ref->[$i][8] .= "cds_type=Internal;";
		}
	}
}
# --------------------------------------------
sub CreateStopCodon
{
	# ref on array of arrays
	my $ref = shift;

	# output
	my @arr = ();

	# first CDS from the left
	my @current = @{$ref->[0]};

	if ( $current[6] eq '+' )
	{
		my $idx = (scalar @{$ref}) - 1;

		# last CDS from the left
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
			# last CDS from the left

			if ( $current[4] - $current[3] + 1 == 2 )
			{
				$current[2] = "stop_codon";
				$current[3] = $current[4] - 1;
				$current[7] = 2;

				push @arr, [ @current ];

				# before last CDS from left
				@current = @{$ref->[$idx - 1]};

				$current[2] = "stop_codon";
				$current[3] = $current[4];
				$current[7] = 0;

				push @arr, [ @current ];
                        }
			elsif ( $current[4] - $current[3] + 1 == 1 )
			{
				$current[2] = "stop_codon";
				$current[3] = $current[4];
				$current[7] = 1;

				push @arr, [ @current ];

				@current = @{$ref->[$idx - 1]};

				$current[2] = "stop_codon";
				$current[3] = $current[4] - 1;
				$current[7] = 0;

				push @arr, [ @current ];
			}

			if ($warnings)
			{
				print "warning, split stop codon detected:\n";
				print Dumper(\@arr) if $debug;
			}
		}
	}
	elsif ( $current[6] eq '-' )
	{
		# first CDS from the left

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
				print Dumper(\@arr) if $debug;
			}
		}
	}

	@arr = sort{ $a->[3] <=> $b->[3] } @arr;

	return @arr;
}
# --------------------------------------------
sub CreateStartCodon
{
	# ref on array of array - values of CDS from one mrna - sorted
	my $ref = shift;

	# output array of arrays
	my @arr = ();

	# this is first CDS on the left side
	my @current = @{$ref->[0]};

	if ( $current[6] eq '+' )
	{
		# complete start codon
		if ( $current[4] - $current[3] + 1 >= 3 )
		{
			$current[2] = "start_codon";
			$current[4] = $current[3] + 3 - 1;
			$current[7] = 0;

			push @arr, [ @current ];
		}
		else
		{
			if ( @{$ref} < 2 )
				{ die "error, not enough data to position start codon\n"; }

			# to do : check strand

			if ( $current[4] - $current[3] + 1 == 2 )
			{
				$current[2] = "start_codon";
				$current[4] = $current[3] + 1;
				$current[7] = 0;

				push @arr, [ @current ];

				# mode to second CDS from left
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

				# mode to second CDS from left
				@current = @{$ref->[1]};

				$current[2] = "start_codon";
				$current[4] = $current[3] + 1;
				$current[7] = 2;

				push @arr, [ @current ];
			}

			if ($warnings)
			{
				print "warning, split start codon detected:\n";
				print Dumper(\@arr) if $debug;
			}
		}
	}
	elsif ( $current[6] eq '-' )
	{
		my $idx = (scalar @{$ref}) - 1;

		# this is last CDS from the left
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
			if ( @{$ref} < 2 )
				{ die "error, not enough data to position start codon\n"; }

			if ( $current[4] - $current[3] + 1 == 2 )
			{
				$current[2] = "start_codon";
				$current[3] = $current[4] - 1;
				$current[7] = 0;

				push @arr, [ @current ];

				# before last CDS from left
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

				# before last CDS from left
                                @current = @{$ref->[$idx - 1]};

                                $current[2] = "start_codon";
                                $current[3] = $current[4] - 1;
                                $current[7] = 2;

                                push @arr, [ @current ];
                        }

			if ($warnings)
			{
				print "warning, split start codon detected:\n";
				print Dumper(\@arr) if $debug;
			}
		}
	} 

	@arr = sort{ $a->[3] <=> $b->[3] } @arr;

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
	# ref on array of arrys - CDS values from GFF file - one mrna - sorted
	my $ref = shift;

	my $size = scalar @{$ref};

	# ref on arry of arrays with introns - output
	my @arr = ();

	# two CDS minimum for intron deriviation
	return @arr if ($size < 2);

	my $i = 0;
	my $j = 1;

	while( $j < $size)
	{
		# two exons must be on the same strand

		if ( $ref->[$i][6] eq $ref->[$j][6] )
		{
			my @current = @{$ref->[$i]};

			$current[2] = "intron";
			$current[3] = $ref->[$i][4] + 1;
			$current[4] = $ref->[$j][3] - 1;

			if ( $ref->[$i][6] eq "+" )
			{
				$current[7] = 3 - $ref->[$j][7];
				$current[7] = 0 if ( $current[7] == 3 );
			}
			elsif ( $ref->[$i][6] eq "-" )
			{
				$current[7] = 3 - $ref->[$i][7];
				$current[7] = 0 if ( $current[7] == 3 );
			}

			if ( $current[4] - $current[3] + 1 < $min_intron )
			{
				$current[2] = "gap";

				if ($warnings)
				{
					print "warning, distance between CDS-CDS is below $min_intron: intron was replased by gap label\n";
					print Dumper(\@current) if $debug;
				}
			}

			push @arr, [@current];
		}
		else
		{
			print "warning, oposite strand CDS were detected: intron is not assigned in such cases\n" if $warnings;
		}

		$i += 1;
		$j += 1;
	}

	return @arr;
}
# --------------------------------------------
sub SelectCDSsorted
{
	# ref on array of arrays - GFF values one gene
	my $ref = shift;

	# put here only CDS lines
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

	# is it possible to have trans-splicing from different chromosomes?
	
	@arr = sort{ $a->[0] cmp $b->[0] || $a->[3] <=> $b->[3] } @arr;

	if ( $warnings and (@arr > 0))
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
		my $seqid = $arr[0][0];

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
			if ( $seqid ne $current->[0] )
			{
				print "warning, two seqid detected in one mRNA: $seqid  $current->[0]\n";
				last;
			}
		}
	}

	return @arr;
}
# --------------------------------------------
sub SplitByMRNA
{
	# ref array of arrays with GGF values of one gene
	my $ref = shift;
	# ref on hash of arrays - with GFF values separated by mrna ID
	my $h_mrna = shift;

	# some mRNA may be non-coding; keep names here:
	my %other_transcripts = ();

	foreach my $entry (@{$ref})
	{
		next if ( $entry->[2] eq "mRNA" );
		next if ( $entry->[2] eq "gene" );
		next if ( $entry->[2] eq "nc_primary_transcript" );
		next if ( $entry->[2] eq "lnc_RNA" );
		next if ( $entry->[2] eq "transcript" );
		next if ( $entry->[2] eq "unconfirmed_transcript" );
		next if ( $entry->[2] =~ /^\S_gene_segment$/ );

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
					{
						print Dumper($ref) if $debug;
						die "error, feature not in mRNA:\n$value\n";
					}
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
	# ref on array of arrays - GFF values one gene
	my $ref = shift;

	# one or many mRNA per record 
	my %mRNA = ();

	foreach my $entry (@{$ref})
	{
		if (($entry->[2] eq "mRNA") or ($entry->[2] eq "nc_primary_transcript") or ($entry->[2] eq "lnc_RNA") or ($entry->[2] eq "transcript") or ($entry->[2] eq "unconfirmed_transcript") or ( $entry->[2] =~ /^\S_gene_segment$/))
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
		print Dumper($ref) if $debug;
		die "error, no mRNA in record\n";
	}

	return %mRNA;
}
# --------------------------------------------
sub GetGeneID
{	
	# ref on array of arrays of GFF3 values from one gene
	my $ref_in = shift;
	# ref on array of arrays - output
	my $ref_out = shift;

	# gene ID - one gene line per record
	my $gene_id = '';

	foreach my $entry (@{$ref_in})
	{
		if ( $entry->[2] eq "gene" )
		{
			if ( !$gene_id )
			{
				if ( $entry->[8] =~ /ID=(\S+?);/ ) 
				{
					$gene_id = $1;

					push @{$ref_out}, [@{$entry}];
				}
				else
					{ die "error, gene ID field not found in: $entry->[8]\n"; }
			}
			else
				{ die "error, gene entry duplication was detected in record: $gene_id\n"; }
		}
	}

	if ( ! $gene_id )
		{ die "error, gene id is missing:\n"; }

	return $gene_id;
}
# --------------------------------------------
sub IsCDSgene
{
	# ref on array of arrays of GFF3 values from one gene
	my $ref = shift;

	my $gene_count = 0;
	my $cds_count = 0;

	my $gene_name = '';

	foreach my $entry (@{$ref})
	{
		if ( $entry->[2] eq "CDS" )
		{
			$cds_count += 1;
		}
		elsif ( $entry->[2] eq "gene" )
		{
			$gene_count += 1;
			$gene_name = $entry->[8];
		}
	}

	if ( $gene_count > 1 )
	{
		print Dumper($ref) if $debug;
		die "error, gene label duplication was detected: $gene_name\n";
	}

	if ( $gene_count == 1 and $cds_count > 0 )
	{
		return 1;
	}
	elsif ( $cds_count == 0 )
	{
		return 0;
	}
	else
	{
		die "error, unexpected combination of gene and cds count was detected: $gene_name $gene_count $cds_count\n";
	}
}
# --------------------------------------------
sub CheckForValidGFF
{
	my $str = shift;
	my $ref = shift;

	# seqid 1
	if ( $ref->[0] !~ /^\w+$/ )
		{ die "error, unexpected seqid format found:\n$str\n"; }

	# start 4
	if ( $ref->[3] !~ /^\d+$/ )
		{ die "error, unexpected start format found:\n$str\n"; }

	# end 5
	if ( $ref->[4] !~ /^\d+$/ )
		{ die "error, unexpected end format found:\n$str\n"; }

	# start <= end
	if ( $ref->[3] > $ref->[4] )
		{ die "error, start is more than end:\n$str\n"; }

	# strand 6
	if ( $ref->[6] !~ /^[+-.]$/ )
		{ die "error, wrong strand value:\n$str\n"; }

	# phase 7
	if ( $ref->[7] !~ /^[.012]$/ )
		 { die "error, wrong phase value:\n$str\n"; }
}
# --------------------------------------------
sub AddLineToRecord
{
	# str - line from GFF3 file
	# ref on array of arrays - put all split lines from one gene here
	my $str = shift;
	my $ref = shift;

	chomp $str;

	my @arr = split( '\t', $str );

	my $size = @arr;

	if ( $size != 8 and $size != 9 )
		{ die "error, unexpected number of TABs found:\n$str\n"; }

	CheckForValidGFF( $str, \@arr ) if $debug;

	push @{$ref}, [@arr];
}
# ------------------------------------------------
sub CheckBeforeRun
{
	die "error, file not found: option --in $in\n" if( ! -e $in );
	die "error, output file name matches input file: $in $out\n" if ( $out eq $in );

	if ( $seq )
	{
		die "error, file not found: option --seq $seq\n" if( ! -e $seq );
		die "error, output file name matches input file: $seq $out\n" if ( $out eq $seq );
	}
}
# ------------------------------------------------
sub ParseCMD
{
	my $opt_results = GetOptions
	(
		'in=s'     => \$in,
		'out=s'    => \$out,
		'verbose'  => \$v,
		'debug'    => \$debug,
		'warnings' => \$warnings,
		'cds_gene_only' => \$cds_gene_only,
		'min_intron=i'  => \$min_intron,
		'seq=s'    => \$seq,
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
  and enriches annotation by adding introns, stop + start codons, CDS types, etc

Optional 
  --cds_gene_only    output only protein coding genes
  --min_intron [$min_intron]  minimum length of intron to calculate from CDS-CDS
  --seq [name]       file with sequence in FASTA format

General:
  --verbose
  --debug
  --warnings

Version: $VERSION

);
	exit 1;
}
# ------------------------------------------------

