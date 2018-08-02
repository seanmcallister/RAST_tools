#!/usr/bin/perl -w
use strict;
use Getopt::Std;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Modify the RAST tab-delimited gene files for use in import to ANVIO as annotations.
#RAST calls this the "Spreadsheet (tab-separated text format)" under downloads.

# - - - - - O P T I O N S  - - - - - -
my %options=();
getopts("i:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-i = infile\n";
	print "-h = This help message\n\n";
	die;
    }

my %GeneDict;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
open(IN, "<$options{i}") or die "\n\nFile $options{i} does not exist or was not given. Try -h for the help file.\n\n";
my @data = <IN>; close(IN); shift(@data);
my $unid = 1000001;
foreach my $line (@data)
	{   	chomp($line);
                my @geneinfo = split('\t', $line);
                my $original_contigid = $geneinfo[0];
                $original_contigid =~ s/\./_/g;
                $original_contigid =~ s/\|/_/g;
                $original_contigid =~ s/-/_/g;
                $original_contigid =~ s/_$//;
                $GeneDict{$unid}{'contig_id'} = uc($original_contigid);
                $GeneDict{$unid}{'feature_id'} = $geneinfo[1];
                $GeneDict{$unid}{'type'} = $geneinfo[2];
                my $start_original = $geneinfo[4];
                my $start_zeroindex = $start_original - 1;
                $GeneDict{$unid}{'start'} = $start_zeroindex;
                my $stop_original = $geneinfo[5];
                my $stop_zeroindex = $stop_original - 1;
                $GeneDict{$unid}{'stop'} = $stop_zeroindex;
                $GeneDict{$unid}{'strand'} = $geneinfo[6];
                $GeneDict{$unid}{'annotation'} = $geneinfo[7];
                $unid += 1;
	}

print "gene_callers_id\tsource\taccession\tfunction\te_value\n";

foreach my $i (sort keys %GeneDict)
    {if ($GeneDict{$i}{'type'} eq "peg") {
        my $geneids = $i - 1000000;
        print "$geneids\tRAST\t\t$GeneDict{$i}{'annotation'}\t0\n";
      }
    }

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
