#!/usr/bin/perl
# ==============================================================
# Alex Lomsadze
# 2019
# Georgia Institute of Technology, Atlanta, Georgia, US
#
# This script selects one transcript with longest CDS per gene from GTF file
# ==============================================================

use strict;
use warnings;
use Getopt::Long qw( GetOptions );
use Data::Dumper;

# ------------------------------------------------

my $in = '';
my $out = '';

my $no_phase = '';
my $stat = '';
my $v ='';

# ------------------------------------------------

Usage() if ( @ARGV < 1 );
ParseCMD();
CheckBeforeRun();

my %transcript_to_join_cds;
my %transcript_to_gene;
my %transcript_to_length;

ParseGFF( $in, \%transcript_to_join_cds, \%transcript_to_gene, \%transcript_to_length );

my %cds_to_transcripts = ReverseKeyValue( \%transcript_to_join_cds );

#print Dumper(\%);

my %transcript_to_transcripts = TrToTrs( \%cds_to_transcripts );

my %gene_to_longest_transcript = GeneToLength( \%transcript_to_length, \%transcript_to_gene, \%transcript_to_transcripts );

PrintValueKey($out, \%gene_to_longest_transcript);

if ( $stat )
{
	open( my $OUT, ">", $stat ) or die "error on open file $stat: $!\n";
	foreach my $key (keys %cds_to_transcripts)
	{
		if ( $cds_to_transcripts{$key} =~ /^(\S+)/ )
        	{
			my $tr_id = $1;

			if (!exists $transcript_to_gene{$tr_id}) {die "$tr_id";}

			my $tr_per_cds = () = $cds_to_transcripts{$key} =~ / /g;
			my $cds_per_tr = () = $key =~ / /g;

			print $OUT $transcript_to_gene{$tr_id} ."\t". $tr_id ."\t". $transcript_to_length{$tr_id} ."\t". $tr_per_cds ."\t". $cds_per_tr ."\t". $cds_to_transcripts{$key} ."\n";
		}
		else {die"error, unexpected format found: $key\n";}
	}
	close $OUT;
}


exit 0;

# ================= subs =========================
sub GeneToLength
{
	my $tr_2_length = shift;
	my $tr_2_gene = shift;
	my $tr_2_trs = shift;

	my $count_identicall = 0;
	my $count_diff_length = 0;
	my $count_match_in_length = 0;

	my %gene_2_tr;

	foreach my $key (keys %{$tr_2_gene})
	{
		if ( exists $gene_2_tr{ $tr_2_gene->{$key} } )
		{
			if ( $tr_2_trs->{$gene_2_tr{ $tr_2_gene->{$key} }} =~ $key )
			{
				$count_identicall += 1;
				next;
			}

			if ( $tr_2_length->{ $gene_2_tr{ $tr_2_gene->{$key} }  } < $tr_2_length->{ $key } )
			{
				$gene_2_tr{ $tr_2_gene->{$key} } = $key;
				$count_diff_length += 1;
			}
			elsif ( $tr_2_length->{ $gene_2_tr{ $tr_2_gene->{$key} }  } == $tr_2_length->{ $key } )
			{
				$count_match_in_length += 1;
				print "# warning match in length found: $gene_2_tr{ $tr_2_gene->{$key} } $key\n" if $v;
			}
		}
		else
		{
			$gene_2_tr{ $tr_2_gene->{$key} } = $key;
		}
	}

	if ( $v )
	{
		print "iso CDS is identicall: $count_identicall\n";
		print "iso CDS is longer: $count_diff_length\n";
		print "iso transcript match in CDS length: $count_match_in_length\n";
	}

	return %gene_2_tr;
}
# ------------------------------------------------
sub TrToTrs
{
	my $cds_to_tr = shift;

	my %h;

	foreach my $key (keys %{$cds_to_tr})
	{
		my @arr = split( ' ', $cds_to_tr->{$key} );

		foreach my $id (@arr)
		{
			$h{$id} = $cds_to_tr->{$key};
		}
	}

	return %h;
}
# ------------------------------------------------
sub PrintFirstId
{
	my $ref = shift;

	foreach my $key (keys %{$ref})
	{
		if ( $ref->{$key} =~ /^(\S+)/ )
		{
			print $1 ."\n";
		}
		else {die"error, unexpected format found: $key\n";}
	}

	exit 1;
}
# ------------------------------------------------
sub ReverseKeyValue
{
	my $ref = shift;
	my %new_h;

	foreach my $key (keys %{$ref})
	{
		$new_h{ $ref->{$key} } .= ($key ." ");
	}

	if ($v)
	{
		print "non-redundant transcripts: ". (scalar keys %new_h) ."\n";
	}

	return %new_h;
}
# ------------------------------------------------
sub PrintValueKey
{
	my $name = shift;
	my $ref = shift;

	open( my $OUT, ">", $name ) or die "error on open file $name: $!\n";
	foreach my $key (sort keys %{$ref})
	{
		print $OUT "$ref->{$key}\t$key\n";
	}
	close $OUT;
}
# ------------------------------------------------
sub PrintKeys
{
	my $ref = shift;
	my $label = shift;
	my $name = shift;

	open( my $OUT, ">", $name ) or die"error on open file";

	print $OUT "# $label\n";

	foreach my $key (sort keys %{$ref})
	{
		print $OUT "$key\n";
	}

	close $OUT;
}
# ------------------------------------------------
sub ParseGFF
{
	my ($name, $tr_2_cds, $tr_2_gene, $tr_2_length) = @_;

	open( my $IN, $name ) or die "error on open file $name: $!\n";
	while( my $line = <$IN> )
	{
		next if( $line =~ /^\s*#/ );
		next if( $line =~ /^\s*$/ );

#		if( $line =~ /^(\S+)\t\S+\t(\S+)\t(\d+)\t(\d+)\t\S+\t([-+.])\t(\S+)\ttranscript_id \"(\S+)\"; gene_id \"(\S+)\";/ )		
		if( $line =~ /^(\S+)\t\S+\t(\S+)\t(\d+)\t(\d+)\t\S+\t([-+.])\t(\S+)\s+gene_id \"(\S+)\"; transcript_id \"(\S+)\";/ )
		{
			my $id      = $1;
			my $type    = $2;
			my $start   = $3;
			my $end     = $4;
			my $strand  = $5;
			my $ph      = $6;
			my $gene_id = $7;
			my $tr_id   = $8;

			next if ( $type ne "CDS" );

			my $key = $id ."_". $start ."_". $end ."_". $strand;

			if ( ! $no_phase )
			{
				$key .= "_". $ph;

				print "warning, strand is not defined: $strand\n" if ($strand eq ".");
			}


			$key .= " ";

			$tr_2_cds->{$tr_id} .= $key;
			$tr_2_gene->{$tr_id} = $gene_id;
			$tr_2_length->{$tr_id} += ($end - $start + 1);
		}
		else
		{
			print "warning, unxpected format found in gff: $line\n";
		}
	}	
	close $IN;
	
	if ($v)
	{
		print "transcripts in file $name: ". (scalar keys %$tr_2_cds) ."\n";
		print "transcripts in file $name: ". (scalar keys %$tr_2_gene) ."\n";
		print "transcripts in file $name: ". (scalar keys %$tr_2_length) ."\n";
	}
	
#	print Dumper($tr_2_cds) if $v;
}
# ------------------------------------------------
sub CheckBeforeRun
{
	die "error, file name is missing $0: option --in\n" if( ! $out );
	die "error, file is missing $0: option --in\n" if( ! -e $in );
	die "error, file name is missing $0: option --out\n" if( ! $out );
	die "error, output file name matches input file\n" if ( $out eq $in ) ;
};
# ------------------------------------------------
sub ParseCMD
{
	my $opt_results = GetOptions
	(
		'in=s'     => \$in,
		'out=s'    => \$out,
		'no_phase' => \$no_phase,
		'stat=s'   => \$stat,
		'verbose'  => \$v,
	);

	die "error on command line\n" if( !$opt_results );
	die "error, unexpected argument found on command line\n" if( @ARGV > 0 );
};
# ------------------------------------------------
sub Usage
{
	print qq(
Usage:
$0  --in [name] --out [name]

Find longest CDS and non-redundunt set of CDS genes in GTF file.
Input file must be sorted by columns 1,n4,n5
 
   --in  [name]  file in GTF format
   --out [name]  output file name

Optional:
   --no_phase    ignore phase in CDS comparision
   --stat [name] output stat on input
   --verbose
);
	exit 1;
}
# ------------------------------------------------

