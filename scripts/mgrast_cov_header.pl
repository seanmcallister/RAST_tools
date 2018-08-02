#!/usr/bin/perl
use strict;
use Getopt::Std;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Convert fasta header to include _[cov=#] for mgrast.

#INFO on Preparing depth.txt:
##Metabat runs a script called jgi_summarize_bam_contig_depths. The output of this (*.depth.txt) is the input for this script.
##You can't just use the one you used for binning, because that depth.txt file is a bit too conservative. We want to annotate the coverage of the unbinned contigs, not throw them out due to size/coverage.
##SO, you need to set the --minContigLength to 1, and --minContigDepth to 1.
##Example: 'jgi_summarize_bam_contig_depths --outputDepth S6_qccontigs_simpname.fasta.depth.txt --pairedContigs S6_qccontigs_simpname.fasta.paired.txt --minContigLength 1 --minContigDepth 1  S6qcreads_v_S6qcspadescontigs_sorted.bam'


# - - - - - U S E R    V A R I A B L E S - - - - - - - -


# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %options=();
getopts("f:c:h", \%options);


if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-f = assembly fasta file\n";
        print "-c = coverage stats from metabat's depth.txt output. Less this script for more info in the header\n";
	print "-h = This help message\n\n";
	die;
    }

my %Sequences;
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
open(IN, "<$options{c}") or die "\n\nNADA $options{c} you FOOL!!!\n\n";
my @cov_data = <IN>; close(IN); shift(@cov_data);
foreach my $i (@cov_data)
	{	chomp($i);
		my @data = split('\t', $i);
		$Sequences{$data[0]}{'cov'} = $data[2];
	}

$/=">";
my $infile = $options{f};
open(IN2, "<$infile") or die "\n\nNADA $infile you FOOL!!!\n\n";
my @DATA = <IN2>; close(IN2); shift(@DATA); $/="\n";	
foreach my $entry (@DATA)
	{	my @data = split('\n', $entry);
		my $seq = '';
		foreach my $i (1..$#data)
		{	$seq .= $data[$i];  }
		$seq =~ s/>//;
		my $cleanhead = $data[0];
		my @shorthead = split(' ', $cleanhead);
		my $ok = $shorthead[0];
		$Sequences{$ok}{'gappy-ntseq'} = uc($seq);
	}



foreach my $i (sort keys %Sequences)
	{	print STDOUT ">".$i."_[cov=".$Sequences{$i}{'cov'}."]\n".$Sequences{$i}{'gappy-ntseq'}."\n";	
	}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
